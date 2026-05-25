#!/usr/bin/env python3
"""
CoSAG Cluster Optimizer (CoSAG-nf adapted)

Adapted for the multi-cluster round JSON produced by
prepare_cluster_json.py / prepare_round_json.py:

  {
    "version": "1.0", "round": 1, "member_type": "sag",
    "n_clusters": N, "n_sags": M,
    "clusters": [
      {
        "cluster_id": "R1_C1", "cluster_size": 5,
        "leaf_sags": [...],
        "members": [{"sag_id", "read1", "read2", "individual_contigs", "checkm2": {...}}, ...],
        "merged_reads": {"read1", "read2"},
        "coassembly_contigs": "...",
        "checkm2_coassembly": {...},
        "status": "completed"
      }, ...
    ]
  }

Algorithm (unchanged): TNF-based hierarchical clustering -> outlier removal ->
greedy iterative single-SAG removal -> SPAdes + CheckM2 evaluation.

Key adaptations vs. the original script:
  * Top-level JSON is a list of clusters (can process all, or a single --cluster-id).
  * Pre-computed coassembly_contigs/checkm2_coassembly is reused as iteration-0
    baseline, avoiding a redundant SPAdes run on the full member set.
  * fastq.gz handled in binary mode (gzip streams are concatenable).
  * current_result carries forward between iterations (no re-assembly of the
    set we just chose).
  * Iteration terminates cleanly at min_sags (no silent no-op loops).
  * CheckM2 TSV parsed by header name (not fragile column indices).
  * Canonical 4-mer TNF (136 keys) by default.
  * genome_size / contig_n50 are stored in the final result.

Runtime: expected to run inside the CoSAG-nf container (SPAdes + CheckM2 +
numpy/scipy/biopython). Call as:

  cluster_optimizer.py \\
      --json round1.json \\
      --checkm2-db uniref100.KO.1.dmnd \\
      --outdir results \\
      [--cluster-id R1_C1] [--threads N]

For ad-hoc use outside the container, set SPADES_BIN / CHECKM2_BIN env vars
to point at the executables you want to use.
"""

import argparse
import json
import logging
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime
from itertools import product

import numpy as np
from Bio import SeqIO
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist


# Tool executables. The script assumes it runs *inside* the CoSAG-nf
# container where `spades.py` and `checkm2` are both on PATH. For ad-hoc
# testing outside the container, override via environment variables.
SPADES_BIN  = os.environ.get('SPADES_BIN',  'spades.py')
CHECKM2_BIN = os.environ.get('CHECKM2_BIN', 'checkm2')


def canonical_tetranucleotides():
    """Return 136 canonical 4-mers (each 4-mer and its RC collapsed to lex-min)."""
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    seen, canon = set(), []
    for p in product('ATGC', repeat=4):
        k = ''.join(p)
        rc = ''.join(comp[b] for b in reversed(k))
        key = min(k, rc)
        if key not in seen:
            seen.add(key)
            canon.append(key)
    return canon


class ClusterOptimizer:
    def __init__(self, cluster_dict, output_dir, checkm2_db, round_num=1,
                 target_contamination=5.0, min_completeness=90.0,
                 max_iterations=10, min_sags=5, threads=20,
                 use_canonical_tnf=True):
        self.cluster_data = cluster_dict
        self.cluster_id = cluster_dict['cluster_id']
        self.round_num = round_num
        self.output_dir = output_dir
        self.checkm2_db = os.path.abspath(checkm2_db)
        self.target_contamination = target_contamination
        self.min_completeness = min_completeness
        self.max_iterations = max_iterations
        self.min_sags = min_sags
        self.threads = threads
        self.use_canonical = use_canonical_tnf

        os.makedirs(output_dir, exist_ok=True)
        self.all_contig_dir = os.path.join(output_dir, "all_contig")
        os.makedirs(self.all_contig_dir, exist_ok=True)

        self._setup_logging()

        self.tetras = (canonical_tetranucleotides()
                       if self.use_canonical
                       else [''.join(p) for p in product('ATGC', repeat=4)])
        self._rc_table = str.maketrans('ATGC', 'TACG')
        self.optimization_members = self._select_optimization_members()

        n_members = cluster_dict.get('cluster_size', len(cluster_dict.get('members', [])))
        msg = (f"Optimizing cluster {self.cluster_id} "
               f"(round {round_num}, {n_members} members)")
        self.logger.info(msg)
        self.logger.info(f"Target: Completeness ≥ {min_completeness}% AND "
                         f"Contamination ≤ {target_contamination}%")
        print(msg)

    def _setup_logging(self):
        log_file = os.path.join(
            self.output_dir,
            f"optimization_{self.cluster_id}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log",
        )
        self.logger = logging.getLogger(f'ClusterOptimizer.{self.cluster_id}')
        self.logger.setLevel(logging.INFO)
        for h in self.logger.handlers[:]:
            self.logger.removeHandler(h)
        fh = logging.FileHandler(log_file, encoding='utf-8')
        fh.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        self.logger.addHandler(fh)
        self.logger.info(f"=== Cluster {self.cluster_id} optimization started ===")
        self.logger.info(f"Log file: {log_file}")

    # ------------------------------------------------------------------
    #  TNF signatures
    # ------------------------------------------------------------------
    def calculate_tetranucleotide_signature(self, fasta_file):
        counts = {t: 0 for t in self.tetras}
        total = 0
        try:
            for rec in SeqIO.parse(fasta_file, "fasta"):
                seq = str(rec.seq).upper()
                for i in range(len(seq) - 3):
                    k = seq[i:i + 4]
                    if 'N' in k:
                        continue
                    if self.use_canonical:
                        rc = k.translate(self._rc_table)[::-1]
                        k = min(k, rc)
                    if k in counts:
                        counts[k] += 1
                        total += 1
        except Exception as e:
            self.logger.error(f"Error reading {fasta_file}: {e}")
            return None
        if total == 0:
            return None
        return np.array([counts[t] / total for t in self.tetras])

    def calculate_signatures_for_all_members(self):
        print("Step 1: Calculating tetranucleotide signatures...")
        signatures, valid_members = {}, []
        for m in self.optimization_members:
            mid = self._member_id(m)
            # Round-2+ compatibility: member may carry either
            # 'individual_contigs' (SAG) or 'coassembly_contigs' (prev-round cluster)
            contigs_path = m.get('individual_contigs') or m.get('coassembly_contigs')
            if not contigs_path or not os.path.exists(contigs_path):
                self.logger.warning(f"Missing contigs for {mid}: {contigs_path}")
                print(f"  Warning: missing contigs for {mid}")
                continue
            print(f"  Processing {mid}...")
            sig = self.calculate_tetranucleotide_signature(contigs_path)
            if sig is not None:
                signatures[mid] = sig
                valid_members.append(m)
            else:
                self.logger.warning(f"Skipping {mid} (empty/invalid contigs)")
                print(f"  Warning: skipping {mid}")
        print(f"  Got signatures for {len(signatures)} / {len(self.optimization_members)} members")
        return signatures, valid_members

    @staticmethod
    def _member_id(member):
        """
        Robust member ID across rounds:
          - round1 members: sag_id
          - round2+ members: cluster_id
        """
        return member.get('sag_id') or member.get('cluster_id') or member.get('id')

    def _select_optimization_members(self):
        """
        Round-2+ preference:
          If leaf_sags contains per-SAG contigs, optimize at SAG level.
          Otherwise fall back to cluster members.
        """
        leaves = self.cluster_data.get('leaf_sags', [])
        leaf_dicts = [x for x in leaves if isinstance(x, dict)]
        leaf_with_contigs = [x for x in leaf_dicts if x.get('individual_contigs')]
        if leaf_with_contigs:
            self.logger.info("Using leaf_sags as optimization members (%d)", len(leaf_with_contigs))
            return leaf_with_contigs
        return self.cluster_data.get('members', [])

    # ------------------------------------------------------------------
    #  Outlier detection & candidate ranking
    # ------------------------------------------------------------------
    def perform_hierarchical_clustering(self, signatures):
        # Hierarchical clustering with 3 splits needs enough members to be meaningful.
        if len(signatures) < 6:
            return list(signatures.keys()), []
        print("Step 2: Hierarchical clustering for outlier detection...")
        ids = list(signatures.keys())
        mat = np.array([signatures[i] for i in ids])
        Z = linkage(pdist(mat, metric='euclidean'), method='ward')
        labels = fcluster(Z, t=3, criterion='maxclust')
        groups = {}
        for i, c in enumerate(labels):
            groups.setdefault(c, []).append(ids[i])
        main = max(groups.values(), key=len)
        outliers = [m for g in groups.values() if g is not main for m in g]
        print(f"  Main: {len(main)} | Outliers: {len(outliers)} -> {outliers}")
        return main, outliers

    def identify_candidates_to_remove(self, signatures, member_ids, n=3):
        if len(member_ids) <= self.min_sags:
            return []
        ids = list(member_ids)
        mat = np.array([signatures[i] for i in ids])
        totals = {}
        for i, mid in enumerate(ids):
            # Self-distance is 0, so sum over all rows == sum over non-self.
            totals[mid] = float(np.linalg.norm(mat - mat[i], axis=1).sum())
        ranked = sorted(totals.items(), key=lambda x: x[1], reverse=True)
        max_cand = min(n, len(ids) - self.min_sags)
        return [mid for mid, _ in ranked[:max_cand]]

    # ------------------------------------------------------------------
    #  Assembly + evaluation
    # ------------------------------------------------------------------
    def create_merged_fastq_files(self, selected_ids):
        """Concatenate gzipped fastq files in binary mode (gzip streams are concatenable)."""
        tmp = tempfile.mkdtemp(prefix=f"cosag_{self.cluster_id}_")
        r1 = os.path.join(tmp, f"{self.cluster_id}_merged_R1.fastq.gz")
        r2 = os.path.join(tmp, f"{self.cluster_id}_merged_R2.fastq.gz")
        selected = set(selected_ids)

        # Build lookup from cluster-level leaf_sags if present
        # (round2+ JSON carries read paths here instead of in members).
        leaf_lookup = {}
        for leaf in self.cluster_data.get('leaf_sags', []):
            if isinstance(leaf, dict):
                sid = leaf.get('sag_id')
                if sid and leaf.get('read1') and leaf.get('read2'):
                    leaf_lookup[sid] = (leaf['read1'], leaf['read2'])

        copied_pairs = set()
        with open(r1, 'wb') as o1, open(r2, 'wb') as o2:
            # If selected IDs are SAG leaf IDs, copy them directly first.
            for sid in selected:
                pair = leaf_lookup.get(sid)
                if not pair or pair in copied_pairs:
                    continue
                with open(pair[0], 'rb') as f:
                    shutil.copyfileobj(f, o1)
                with open(pair[1], 'rb') as f:
                    shutil.copyfileobj(f, o2)
                copied_pairs.add(pair)

            for m in self.cluster_data['members']:
                mid = self._member_id(m)
                if mid not in selected:
                    continue

                # Round1: member directly has read1/read2.
                if m.get('read1') and m.get('read2'):
                    pair = (m['read1'], m['read2'])
                    if pair not in copied_pairs:
                        with open(pair[0], 'rb') as f:
                            shutil.copyfileobj(f, o1)
                        with open(pair[1], 'rb') as f:
                            shutil.copyfileobj(f, o2)
                        copied_pairs.add(pair)
                    continue

                # Round2+: member references leaf_sags; resolve to read pairs.
                for sid in m.get('leaf_sags', []):
                    pair = leaf_lookup.get(sid)
                    if not pair:
                        self.logger.warning(f"Missing leaf reads for {mid}/{sid}")
                        continue
                    if pair in copied_pairs:
                        continue
                    with open(pair[0], 'rb') as f:
                        shutil.copyfileobj(f, o1)
                    with open(pair[1], 'rb') as f:
                        shutil.copyfileobj(f, o2)
                    copied_pairs.add(pair)
        return r1, r2, tmp

    def run_spades_assembly(self, r1, r2, prefix):
        out = os.path.join(self.output_dir, f"{prefix}_spades")
        cmd = [
            SPADES_BIN, "--sc", "--careful",
            "--phred-offset", "33",
            "-1", r1, "-2", r2, "-o", out,
            "-t", str(self.threads),
        ]
        self.logger.info(f"SPAdes: {prefix}")
        self.logger.info("CMD: " + " ".join(cmd))
        try:
            subprocess.run(cmd, capture_output=True, text=True, check=True)
        except subprocess.CalledProcessError as e:
            err = (e.stderr or "")[-500:]
            self.logger.error(f"SPAdes failed ({prefix}): {err}")
            print(f"  SPAdes failed: {prefix}")
            return None
        contigs = os.path.join(out, "contigs.fasta")
        if not os.path.exists(contigs):
            self.logger.error(f"SPAdes contigs.fasta missing: {out}")
            return None
        self.logger.info(f"SPAdes ok: {contigs}")
        return os.path.abspath(contigs)

    def run_checkm2_evaluation(self, contigs_file, prefix):
        out = os.path.join(self.output_dir, f"{prefix}_checkm2")
        os.makedirs(out, exist_ok=True)
        cmd = [
            CHECKM2_BIN, "predict",
            "--input", contigs_file,
            "--output-directory", out,
            "--threads", str(self.threads),
            "--database_path", self.checkm2_db,
            "--force",
        ]
        self.logger.info(f"CheckM2: {prefix}")
        try:
            subprocess.run(cmd, capture_output=True, text=True, check=True)
        except subprocess.CalledProcessError as e:
            err = (e.stderr or "")[-500:]
            self.logger.error(f"CheckM2 failed ({prefix}): {err}")
            print(f"  CheckM2 failed: {prefix}")
            return None
        qf = os.path.join(out, "quality_report.tsv")
        if not os.path.exists(qf):
            self.logger.error(f"quality_report.tsv missing: {out}")
            return None
        return self._parse_checkm2_report(qf)

    @staticmethod
    def _parse_checkm2_report(qf):
        """Parse CheckM2 TSV by header name, robust to column-order changes."""
        with open(qf) as f:
            header = f.readline().rstrip('\n').split('\t')
            line = f.readline()
        if not line:
            return None
        row = dict(zip(header, line.rstrip('\n').split('\t')))

        def g(key, cast=float):
            v = row.get(key)
            if v in (None, '', 'NA'):
                return None
            try:
                return cast(v)
            except ValueError:
                return None

        return {
            'completeness':   g('Completeness'),
            'contamination':  g('Contamination'),
            'genome_size':    g('Genome_Size'),
            'contig_n50':     g('Contig_N50'),
            'gc_content':     g('GC_Content'),
            'coding_density': g('Coding_Density'),
        }

    def test_member_combination(self, selected_ids, iteration, tag=""):
        """Run SPAdes + CheckM2 on a subset. Returns result dict or None."""
        self.logger.info(f"[iter {iteration}] testing {tag}: {len(selected_ids)} members")
        print(f"    Testing {tag}: {len(selected_ids)} members")
        r1, r2, tmp = self.create_merged_fastq_files(selected_ids)
        try:
            prefix = f"iter{iteration}_{tag}_{len(selected_ids)}m"
            contigs = self.run_spades_assembly(r1, r2, prefix)
            if not contigs:
                return None
            metrics = self.run_checkm2_evaluation(contigs, prefix)
            if not metrics or metrics['completeness'] is None:
                return None
            result = {
                'iteration': iteration, 'tag': tag,
                'member_count': len(selected_ids),
                'members': list(selected_ids),
                'contigs_file': contigs,
                'output_prefix': prefix,
                **metrics,
            }
            self.logger.info(
                f"  -> Comp={metrics['completeness']:.2f}% "
                f"Cont={metrics['contamination']:.2f}% "
                f"N50={metrics['contig_n50']:.0f} Size={metrics['genome_size']:.0f}"
            )
            print(f"      Comp={metrics['completeness']:.2f}%  "
                  f"Cont={metrics['contamination']:.2f}%")
            return result
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

    # ------------------------------------------------------------------
    #  Baseline reuse (from JSON)
    # ------------------------------------------------------------------
    def baseline_from_json(self):
        """Use the JSON-provided full-set coassembly as iteration-0, if present."""
        cm = self.cluster_data.get('checkm2_coassembly')
        contigs = self.cluster_data.get('coassembly_contigs')
        if not cm or not contigs or not os.path.exists(contigs):
            return None
        if cm.get('completeness') is None or cm.get('contamination') is None:
            return None
        all_ids = [self._member_id(m) for m in self.optimization_members if self._member_id(m)]
        return {
            'iteration': 0, 'tag': 'baseline_all',
            'member_count': len(all_ids),
            'members': all_ids,
            'contigs_file': os.path.abspath(contigs),
            'output_prefix': f"{self.cluster_id}_baseline",
            'completeness':   cm.get('completeness'),
            'contamination':  cm.get('contamination'),
            'genome_size':    cm.get('genome_size'),
            'contig_n50':     cm.get('contig_n50'),
            'gc_content':     cm.get('gc_content'),
            'coding_density': cm.get('coding_density'),
        }

    # ------------------------------------------------------------------
    #  Main loop
    # ------------------------------------------------------------------
    def _is_optimal(self, r):
        return (r and r.get('completeness') is not None
                and r['completeness'] >= self.min_completeness
                and r['contamination'] <= self.target_contamination)

    def _is_acceptable(self, r):
        return (r and r.get('contamination') is not None
                and r['contamination'] <= 10.0)

    def optimize_cluster(self):
        print(f"\n=== Optimizing Cluster {self.cluster_id} ===")

        # Signatures
        signatures, valid_members = self.calculate_signatures_for_all_members()
        valid_ids = {self._member_id(m) for m in valid_members if self._member_id(m)}

        # Trackers
        best_result = None
        best_acceptable = None
        completeness_history = []
        has_acceptable = False
        acceptable_max = 0.0
        decline_threshold = 10.0
        completeness_floor = 80.0

        # -------- baseline (iteration 0, free from JSON) --------
        baseline = self.baseline_from_json()
        if baseline:
            self.logger.info(
                f"Baseline (all {baseline['member_count']} members): "
                f"Comp={baseline['completeness']:.2f}% "
                f"Cont={baseline['contamination']:.2f}%"
            )
            print(f"Baseline (full set, n={baseline['member_count']}): "
                  f"Comp={baseline['completeness']:.2f}%  "
                  f"Cont={baseline['contamination']:.2f}%")
            completeness_history.append(baseline['completeness'])
            if self._is_optimal(baseline):
                print("  ✓ Baseline already meets target — skipping optimization.")
                self._finalize(baseline, "OPTIMAL_BASELINE", iterations=0,
                               completeness_history=completeness_history,
                               acceptable_max=baseline['completeness'],
                               has_acceptable=True)
                return
            if self._is_acceptable(baseline):
                best_acceptable = baseline
                has_acceptable = True
                acceptable_max = baseline['completeness']

        # If too few valid members for iterative removal, bail to baseline.
        if len(signatures) < self.min_sags:
            msg = (f"Too few valid members ({len(signatures)} < {self.min_sags}); "
                   f"optimization skipped.")
            self.logger.warning(msg)
            print(f"⚠ {msg}")
            if baseline:
                self._finalize(baseline, "BASELINE_ONLY", iterations=0,
                               completeness_history=completeness_history,
                               acceptable_max=acceptable_max,
                               has_acceptable=has_acceptable)
            return

        # -------- hierarchical outlier removal --------
        sub_sigs = {i: signatures[i] for i in valid_ids}
        _, outliers = self.perform_hierarchical_clustering(sub_sigs)
        current_ids = [i for i in signatures if i in valid_ids and i not in outliers]
        if outliers:
            print(f"Removed {len(outliers)} hierarchical-clustering outliers.")

        # If no outliers were removed, the current set == full set,
        # so baseline IS the first current_result — skip re-assembly.
        current_result = baseline if (baseline and not outliers) else None

        # -------- greedy iterative removal --------
        iteration = 1
        while len(current_ids) >= self.min_sags and iteration <= self.max_iterations:
            print(f"\n--- Iter {iteration}: {len(current_ids)} members ---")

            # (1) evaluate current set if not carried forward
            if current_result is None:
                current_result = self.test_member_combination(
                    current_ids, iteration, "current"
                )

            if current_result:
                c = current_result['completeness']
                if iteration != 0 and (not completeness_history
                                       or completeness_history[-1] != c):
                    completeness_history.append(c)

                if self._is_acceptable(current_result):
                    if not has_acceptable:
                        has_acceptable = True
                        print("  ✓ First acceptable result (contamination ≤ 10%)")
                    if c > acceptable_max:
                        acceptable_max = c
                    if best_acceptable is None or c > best_acceptable['completeness']:
                        best_acceptable = current_result

                if self._is_optimal(current_result):
                    best_result = current_result
                    print("  ✓ Optimal target reached. Stopping.")
                    break

                # early stop: significant completeness decline
                stop = False
                if has_acceptable and (acceptable_max - c) > decline_threshold:
                    print(f"  ⚠ Completeness declined >{decline_threshold:.0f}% "
                          f"from best acceptable ({acceptable_max:.1f}% → {c:.1f}%). "
                          f"Stopping.")
                    stop = True
                if c < completeness_floor:
                    print(f"  ⚠ Completeness below floor "
                          f"{completeness_floor:.0f}%. Stopping.")
                    stop = True
                if stop:
                    break

            # (2) cannot remove further — we are at the minimum
            if len(current_ids) <= self.min_sags:
                print(f"  Reached minimum member count ({self.min_sags}). Stopping.")
                break

            # (3) test removing each of the most distant candidates
            candidates = self.identify_candidates_to_remove(sub_sigs, current_ids, n=3)
            if not candidates:
                break

            print(f"  Candidates to remove: {candidates}")
            test_results = []
            for cand in candidates:
                trial_ids = [i for i in current_ids if i != cand]
                r = self.test_member_combination(
                    trial_ids, iteration, f"rm_{cand[:8]}"
                )
                if r:
                    test_results.append((cand, r))

            if not test_results:
                print("  All removal trials failed — dropping first candidate heuristically.")
                current_ids.remove(candidates[0])
                current_result = None
                iteration += 1
                continue

            # (4) predictive stop — all trials drop completeness too much
            if has_acceptable and self._is_acceptable(current_result):
                max_trial = max(r['completeness'] for _, r in test_results)
                if (acceptable_max - max_trial) > decline_threshold:
                    print(f"  ⚠ All removal trials would drop completeness "
                          f">{decline_threshold:.0f}%. Preserving acceptable result.")
                    break

            # (5) pick the best trial: optimal > acceptable+max_comp > min_cont
            optimal_trials = [(c, r) for c, r in test_results if self._is_optimal(r)]
            if optimal_trials:
                cand, r = max(optimal_trials, key=lambda x: x[1]['completeness'])
                current_ids.remove(cand)
                best_result = r
                completeness_history.append(r['completeness'])
                print(f"  ✓ Optimal by removing {cand}. Stopping.")
                break

            acceptable_trials = [(c, r) for c, r in test_results if self._is_acceptable(r)]
            if acceptable_trials:
                cand, r = max(acceptable_trials, key=lambda x: x[1]['completeness'])
            else:
                cand, r = min(test_results, key=lambda x: x[1]['contamination'])

            print(f"  → Remove {cand}: Comp={r['completeness']:.1f}% "
                  f"Cont={r['contamination']:.1f}%")
            current_ids.remove(cand)
            current_result = r  # carry forward — skip re-assembly next iter
            completeness_history.append(r['completeness'])

            if self._is_acceptable(r) and (best_acceptable is None
                                           or r['completeness'] > best_acceptable['completeness']):
                best_acceptable = r
                if r['completeness'] > acceptable_max:
                    acceptable_max = r['completeness']
                    has_acceptable = True

            iteration += 1

        # -------- finalize: prefer optimal > acceptable > baseline --------
        if best_result and self._is_optimal(best_result):
            final, tag = best_result, "OPTIMAL"
        elif best_acceptable:
            final, tag = best_acceptable, "BEST_UNDER_10pct_CONTAMINATION"
        elif baseline:
            final, tag = baseline, "BASELINE_FALLBACK"
        else:
            self.logger.error("No valid result obtained.")
            print("No valid result obtained.")
            return

        self._finalize(final, tag, max(0, iteration - 1),
                       completeness_history, acceptable_max, has_acceptable)

    # ------------------------------------------------------------------
    #  Output writers
    # ------------------------------------------------------------------
    def _finalize(self, result, result_type, iterations,
                  completeness_history, acceptable_max, has_acceptable):
        print(f"\n=== Cluster {self.cluster_id} — Done ===")
        print(f"Type: {result_type}")
        print(f"Iterations: {iterations}")
        print(f"Members: {result['member_count']}")
        print(f"Completeness: {result['completeness']:.2f}%")
        print(f"Contamination: {result['contamination']:.2f}%")
        if result.get('genome_size'):
            print(f"Genome size: {result['genome_size']:.0f} bp   "
                  f"N50: {result.get('contig_n50') or 0:.0f}")
        if self._is_optimal(result):
            print("✓ EXCELLENT: meets both targets")
        elif self._is_acceptable(result):
            print("✓ GOOD: contamination ≤ 10%")
        else:
            print("⚠ ACCEPTABLE: best available result")

        result['result_type'] = result_type
        result['optimization_summary'] = {
            'total_iterations': iterations,
            'max_iterations': self.max_iterations,
            'target_completeness': self.min_completeness,
            'target_contamination': self.target_contamination,
            'meets_targets': self._is_optimal(result),
            'completeness_history': completeness_history,
            'best_acceptable_completeness': acceptable_max if has_acceptable else None,
        }

        result_file = os.path.join(
            self.output_dir, f"cluster_{self.cluster_id}_optimized.json"
        )
        with open(result_file, 'w') as f:
            json.dump(result, f, indent=2, default=str)
        self.logger.info(f"Wrote: {result_file}")
        print(f"Wrote: {result_file}")

        self._save_best_contigs(result)
        self._write_readme(result)

    def _save_best_contigs(self, result):
        src = result.get('contigs_file')
        if not src or not os.path.exists(src):
            self.logger.error(f"Best contigs not found: {src}")
            return
        dst = os.path.join(self.all_contig_dir, f"{self.cluster_id}_best.fasta")
        shutil.copy2(src, dst)
        sz = os.path.getsize(dst) / (1024 * 1024)
        self.logger.info(f"Saved best contigs: {dst} ({sz:.2f} MB)")
        print(f"✓ Best contigs: {dst} ({sz:.2f} MB)")

    def _write_readme(self, result):
        path = os.path.join(self.all_contig_dir, f"{self.cluster_id}_README.txt")
        with open(path, 'w') as f:
            f.write(f"Cluster {self.cluster_id} (round {self.round_num})\n")
            f.write("=" * 50 + "\n")
            f.write(f"Generated : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Type      : {result['result_type']}\n\n")
            f.write("Quality:\n")
            f.write(f"  Completeness : {result['completeness']:.2f}%\n")
            f.write(f"  Contamination: {result['contamination']:.2f}%\n")
            if result.get('genome_size'):
                f.write(f"  Genome size  : {result['genome_size']:.0f} bp\n")
            if result.get('contig_n50'):
                f.write(f"  Contig N50   : {result['contig_n50']:.0f}\n")
            if result.get('gc_content') is not None:
                f.write(f"  GC content   : {result['gc_content']:.3f}\n")
            f.write(f"  Member count : {result['member_count']}\n\n")
            f.write("Selected members:\n")
            for i, mid in enumerate(result['members'], 1):
                f.write(f"  {i}. {mid}\n")


# ----------------------------------------------------------------------
#  CLI
# ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="CoSAG cluster optimizer (CoSAG-nf adapted)."
    )
    ap.add_argument('--json', required=True,
                    help='Round JSON from prepare_cluster_json.py / prepare_round_json.py')
    ap.add_argument('--outdir', required=True, help='Output directory')
    ap.add_argument('--checkm2-db', required=True,
                    help='Path to CheckM2 diamond DB (uniref100.KO.1.dmnd). '
                         'In Nextflow, stage this as a process input.')
    ap.add_argument('--cluster-id', default=None,
                    help='Process only this cluster_id (default: process all)')
    ap.add_argument('--min-completeness', type=float, default=90.0)
    ap.add_argument('--target-contamination', type=float, default=5.0)
    ap.add_argument('--max-iterations', type=int, default=10)
    ap.add_argument('--min-sags', type=int, default=5)
    ap.add_argument('--threads', type=int, default=20)
    ap.add_argument('--no-canonical-tnf', action='store_true',
                    help='Use raw 256-mer TNF instead of 136 canonical 4-mers')
    args = ap.parse_args()

    # --- sanity checks up front ---
    if not os.path.exists(args.json):
        print(f"Error: JSON not found: {args.json}")
        sys.exit(1)
    if not os.path.exists(args.checkm2_db):
        print(f"Error: CheckM2 DB not found: {args.checkm2_db}")
        sys.exit(1)
    for tool in (SPADES_BIN, CHECKM2_BIN):
        if shutil.which(tool) is None:
            print(f"Error: '{tool}' not found on PATH. "
                  f"Run inside the CoSAG-nf container, or set SPADES_BIN/CHECKM2_BIN.")
            sys.exit(1)

    with open(args.json) as f:
        data = json.load(f)

    # Accept both shapes:
    #   * multi-cluster round JSON: {"clusters": [...], "round": N, ...}
    #   * single-cluster JSON (per-cluster scatter):
    #       {"cluster_id": "...", "members": [...], ...}
    if isinstance(data, dict) and 'clusters' in data:
        all_clusters = data['clusters']
        round_num    = int(data.get('round', 1))
        member_type  = data.get('member_type', 'sag')
    elif isinstance(data, dict) and 'cluster_id' in data and 'members' in data:
        all_clusters = [data]
        round_num    = int(data.get('round', 1))
        member_type  = data.get('member_type', 'sag')
    else:
        print("Error: JSON must have top-level 'clusters' list "
              "or be a single cluster dict with 'cluster_id' + 'members'.")
        sys.exit(1)

    if member_type not in ('sag', 'cluster'):
        print(f"Warning: unfamiliar member_type='{member_type}', proceeding anyway")

    clusters = all_clusters
    if args.cluster_id:
        clusters = [c for c in clusters if c['cluster_id'] == args.cluster_id]
        if not clusters:
            print(f"Error: cluster_id '{args.cluster_id}' not found in JSON")
            sys.exit(1)

    print("=== CoSAG Optimizer ===")
    print(f"Input   : {args.json}")
    print(f"DB      : {args.checkm2_db}")
    print(f"Round   : {round_num}  |  member_type: {member_type}")
    print(f"Clusters: {len(clusters)}")
    print(f"Target  : Completeness ≥ {args.min_completeness}% AND "
          f"Contamination ≤ {args.target_contamination}%")
    print()

    for i, cluster in enumerate(clusters, 1):
        size = cluster.get('cluster_size', len(cluster.get('members', [])))
        if size < args.min_sags:
            print(f"[{i}/{len(clusters)}] Skip {cluster['cluster_id']} "
                  f"(size {size} < min_sags {args.min_sags})")
            continue
        print(f"\n[{i}/{len(clusters)}] >>> {cluster['cluster_id']}")
        try:
            opt = ClusterOptimizer(
                cluster_dict=cluster,
                output_dir=os.path.join(args.outdir, cluster['cluster_id']),
                checkm2_db=args.checkm2_db,
                round_num=round_num,
                target_contamination=args.target_contamination,
                min_completeness=args.min_completeness,
                max_iterations=args.max_iterations,
                min_sags=args.min_sags,
                threads=args.threads,
                use_canonical_tnf=not args.no_canonical_tnf,
            )
            opt.optimize_cluster()
        except Exception as e:
            print(f"Error on {cluster['cluster_id']}: {e}")
            import traceback
            traceback.print_exc()
            continue

    print("\n=== All clusters done ===")


if __name__ == "__main__":
    main()