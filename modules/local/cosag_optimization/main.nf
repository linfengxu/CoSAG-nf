// modules/local/cosag_optimization.nf
//
// Per-cluster SAG refinement (TNF-based iterative co-assembly optimisation).
// Expects the upstream workflow to scatter a round JSON into per-cluster
// JSON files; one task per cluster.
//
// Input tuple
//   * meta          – [id, cluster_id, round, ...]
//   * cluster_json  – single-cluster JSON with embedded read paths in
//                     members[].read1 / members[].read2
//
// Outputs
//   * optimized_json  – per-cluster optimisation result
//   * best_contigs    – best CoSAG assembly FASTA
//   * selected_sags   – selected SAG IDs, one per line
//   * selected_reads  – TSV of selected_sag_id → (read1, read2)
//   * log, versions

process COSAG_OPTIMIZATION {
    tag "${meta.cluster_id ?: meta.id}"
    label 'process_high'

    container 'quay.io/xulf2022/spades_checkm2:v1'

    publishDir "${params.outdir}/cosag_optimization/round${meta.round ?: 1}",
               mode: params.publish_dir_mode ?: 'copy',
               pattern: "{${meta.id}_output/**,${meta.id}_sag_opt.log}"

    input:
    tuple val(meta), path(cluster_json)

    output:
    tuple val(meta), path("${meta.id}_output/**/cluster_*_optimized.json"),   emit: optimized_json, optional: true
    tuple val(meta), path("${meta.id}_output/**/all_contig/*_best.fasta"),    emit: best_contigs,   optional: true
    tuple val(meta), path("${meta.id}_output/${meta.id}_selected_sags.txt"),  emit: selected_sags,  optional: true
    tuple val(meta), path("${meta.id}_output/${meta.id}_selected_reads.tsv"), emit: selected_reads, optional: true
    path  "${meta.id}_sag_opt.log",                                           emit: log
    path  "versions.yml",                                                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def prefix   = meta.id
    def cluster_id = meta.cluster_id ?: meta.id
    def min_comp = params.target_completeness   ?: 90
    def tgt_cont = params.target_contamination  ?: 5
    def max_it   = params.tnf_max_iterations    ?: 10
    def min_sag  = params.tnf_min_sags          ?: 5
    """
    set -euo pipefail

    mkdir -p ${prefix}_output

    # ── Optimisation ────────────────────────────────────────────────────
    cosag_optimizer.py \\
        --json                   ${cluster_json} \\
        --checkm2-db             ${params.checkm2_db} \\
        --outdir                 ${prefix}_output \\
        --threads                ${task.cpus} \\
        --min-completeness       ${min_comp} \\
        --target-contamination   ${tgt_cont} \\
        --max-iterations         ${max_it} \\
        --min-sags               ${min_sag} \\
        ${args} \\
        > ${prefix}_sag_opt.log 2>&1

    # ── Post-process: union of selected SAGs + their read paths ─────────
    # Reads are already embedded in cluster_json.members[], so no external
    # sag_db lookup is needed.
    python3 <<PYEOF
import json, pathlib, sys

out_dir   = pathlib.Path("${prefix}_output")
opt_files = sorted(out_dir.rglob("cluster_*_optimized.json"))
if not opt_files:
    print("[post] No optimised JSON produced — skipping aggregation.", file=sys.stderr)
    sys.exit(0)

# Ensure optimised JSON filename carries cluster_id for downstream traceability.
cluster_id = "${cluster_id}"
renamed_opt_files = []
for i, p in enumerate(opt_files, start=1):
    target_name = f"cluster_{cluster_id}_optimized.json" if i == 1 else f"cluster_{cluster_id}_{i}_optimized.json"
    target = p.with_name(target_name)
    if p != target:
        p.replace(target)
    renamed_opt_files.append(target)
opt_files = renamed_opt_files

# Build lookups from cluster_json.
# Accept both single-cluster and multi-cluster JSON shapes.
cj = json.loads(pathlib.Path("${cluster_json}").read_text())
clusters = cj["clusters"] if isinstance(cj, dict) and "clusters" in cj else [cj]

# 1) SAG-level reads lookup: sag_id -> (read1, read2)
reads_map = {}
for c in clusters:
    # Round1 single-cluster shape: members already are SAGs
    for m in c.get("members", []):
        sid = m.get("sag_id")
        r1, r2 = m.get("read1"), m.get("read2")
        if sid and r1 and r2:
            reads_map[sid] = (r1, r2)
    # Round2+ shape: leaf_sags carries SAG read paths
    for leaf in c.get("leaf_sags", []):
        if isinstance(leaf, dict):
            sid = leaf.get("sag_id")
            r1, r2 = leaf.get("read1"), leaf.get("read2")
            if sid and r1 and r2:
                reads_map[sid] = (r1, r2)

# 2) Cluster-member expansion: cluster_id -> [leaf sag ids]
member_to_leaf_sags = {}
for c in clusters:
    for m in c.get("members", []):
        cid = m.get("cluster_id")
        if cid and isinstance(m.get("leaf_sags"), list):
            member_to_leaf_sags[cid] = [s for s in m["leaf_sags"] if isinstance(s, str)]

seen, selected = set(), []
for p in opt_files:
    d = json.loads(p.read_text())
    if isinstance(d, dict):
        d["cluster_id"] = cluster_id
        p.write_text(json.dumps(d, indent=2, ensure_ascii=False) + "\\n")
    for m in d.get("members", d.get("sags", [])):
        # Optimizer may output SAG IDs (round1) or cluster IDs (round2+)
        expanded = member_to_leaf_sags.get(m, [m])
        for sid in expanded:
            if sid not in seen:
                seen.add(sid)
                selected.append(sid)

(out_dir / "${prefix}_selected_sags.txt").write_text(
    "\\n".join(selected) + ("\\n" if selected else "")
)
print(f"[post] selected_sags.txt: {len(selected)} SAGs", file=sys.stderr)

rows, missing = ["sag_id\\tread1\\tread2"], 0
for sid in selected:
    r = reads_map.get(sid)
    if not r:
        missing += 1
        print(f"[post] WARN: {sid} has no reads in cluster_json", file=sys.stderr)
        continue
    rows.append(f"{sid}\\t{r[0]}\\t{r[1]}")

(out_dir / "${prefix}_selected_reads.tsv").write_text("\\n".join(rows) + "\\n")
print(f"[post] selected_reads.tsv: {len(rows)-1} rows, {missing} missing",
      file=sys.stderr)
PYEOF

# ── Versions ────────────────────────────────────────────────────────
cat <<END_VERSIONS > versions.yml
"${task.process}":
    python:    \$(python3 --version 2>&1 | sed 's/Python //')
    biopython: \$(python3 -c 'import Bio; print(Bio.__version__)'    2>/dev/null || echo NA)
    numpy:     \$(python3 -c 'import numpy; print(numpy.__version__)' 2>/dev/null || echo NA)
    scipy:     \$(python3 -c 'import scipy; print(scipy.__version__)' 2>/dev/null || echo NA)
    spades:    \$(spades.py --version 2>&1 | awk '{print \$NF; exit}')
    checkm2:   \$(checkm2 --version 2>&1 | head -1)
END_VERSIONS
    """

    stub:
    def prefix = meta.id
    def cluster_id = meta.cluster_id ?: meta.id
    """
    mkdir -p ${prefix}_output/STUB_cluster/all_contig
    touch ${prefix}_output/STUB_cluster/cluster_${cluster_id}_optimized.json
    touch ${prefix}_output/STUB_cluster/all_contig/STUB_best.fasta
    touch ${prefix}_output/${prefix}_selected_sags.txt
    touch ${prefix}_output/${prefix}_selected_reads.tsv
    touch ${prefix}_sag_opt.log

cat <<END_VERSIONS > versions.yml
"${task.process}":
    python:    stub
    biopython: stub
    numpy:     stub
    scipy:     stub
    spades:    stub
    checkm2:   stub
END_VERSIONS
    """
}
