#!/usr/bin/env python3
"""
Build Round N (N>=2) CoSAG JSON from previous round updated JSON + new cluster assignments.

Round 2:
  prev_json (Round 1 updated) 里 member_type=sag，直接从 members[].read1/read2 取 reads

Round 3+:
  prev_json (Round 2 updated) 里 member_type=cluster，members 里没有 sag_id
  需要 --round1-json 溯源原始 SAG reads

Usage:
    # Round 2
    prepare_round_json.py \
        --prev-json   round1_updated.json \
        --clusters    hierarchical_clusters.tsv \
        --round       2 \
        --outdir      . \
        --out         round2_clusters.json \
        --threads     8

    # Round 3
    prepare_round_json.py \
        --prev-json   round2_updated.json \
        --clusters    hierarchical_clusters.tsv \
        --round       3 \
        --round1-json round1_updated.json \
        --outdir      . \
        --out         round3_clusters.json \
        --threads     8
"""

import argparse
import json
import logging
import subprocess
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


# ─────────────────────────────────────────────
# Parsers
# ─────────────────────────────────────────────

def load_sag_reads_from_round1(path):
    """
    从 Round 1 JSON 建立 sag_id → {read1, read2, checkm2} 查找表
    Round 3+ 时用于溯源原始 SAG reads
    """
    if not path:
        return {}
    with open(path) as f:
        data = json.load(f)
    sag_info = {}
    for cluster in data["clusters"]:
        for member in cluster.get("members", []):
            sid = member.get("sag_id")
            if sid and member.get("read1") and member.get("read2"):
                sag_info[sid] = {
                    "read1":   member["read1"],
                    "read2":   member["read2"],
                    "checkm2": member.get("checkm2"),   # ← 新增
                    "individual_contigs": member.get("individual_contigs"),
                }
    log.info("  round1_json: %d SAGs with reads", len(sag_info))
    return sag_info


def load_prev_json(path, round1_sag_reads=None):
    """
    从上一轮 updated JSON 构建两个查找表:
      cluster_lookup : cluster_id → cluster dict
      sag_info       : sag_id    → {read1, read2, checkm2}

    Round 2: prev_json 是 Round 1（member_type=sag），直接从 members 取
    Round 3+: prev_json 是 Round 2（member_type=cluster），
              reads/checkm2 从 round1_sag_reads 溯源
    """
    with open(path) as f:
        data = json.load(f)

    cluster_lookup = {}
    sag_info = {}
    prev_round  = data.get("round", 1)
    member_type = data.get("member_type", "sag")

    for cluster in data["clusters"]:
        cluster_lookup[cluster["cluster_id"]] = cluster

        if member_type == "sag":
            # Round 1 JSON：直接从 members[].read1/read2/checkm2 取
            for member in cluster.get("members", []):
                sid = member.get("sag_id")
                if sid and member.get("read1") and member.get("read2"):
                    sag_info[sid] = {
                        "read1":   member["read1"],
                        "read2":   member["read2"],
                        "checkm2": member.get("checkm2"),   # ← 新增
                        "individual_contigs": member.get("individual_contigs"),
                    }
        else:
            # Round 2+ JSON：leaf_sags 展平，信息从 round1_sag_reads 溯源
            if round1_sag_reads:
                for entry in cluster.get("leaf_sags", []):
                    # leaf_sags 可能是字符串（旧格式）或对象（新格式）
                    sag_id = entry if isinstance(entry, str) else entry.get("sag_id")
                    if sag_id and sag_id in round1_sag_reads and sag_id not in sag_info:
                        sag_info[sag_id] = round1_sag_reads[sag_id]
            else:
                log.warning(
                    "  prev_json is Round 2+ (member_type=cluster) "
                    "but --round1-json not provided, reads cannot be traced"
                )

    log.info(
        "  prev_json: round=%d, member_type=%s, %d clusters, %d SAGs with reads",
        prev_round, member_type, len(cluster_lookup), len(sag_info)
    )
    return cluster_lookup, sag_info, prev_round


def parse_clusters(path):
    """
    新一轮聚类结果 TSV:
      Member_ID  Cluster_ID  Cluster_Size
    Member_ID 是上一轮的 cluster_id（如 R1_C1 或 R2_C1）
    """
    result = defaultdict(list)
    with open(path) as f:
        f.readline()  # header
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            member_id, cluster_id = parts[0], parts[1]
            result[cluster_id].append(member_id)
    log.info("  new clusters: %d", len(result))
    return dict(result)


# ─────────────────────────────────────────────
# Read merging
# ─────────────────────────────────────────────

def merge_one_cluster(cluster_id, r1_list, r2_list, outdir, round_num):
    cluster_dir = (Path(outdir) / "cluster_{}".format(cluster_id)).resolve()
    cluster_dir.mkdir(parents=True, exist_ok=True)
    r1_out = cluster_dir / "merged_R1.fastq.gz"
    r2_out = cluster_dir / "merged_R2.fastq.gz"

    for inputs, out in [(r1_list, r1_out), (r2_list, r2_out)]:
        valid = [p for p in inputs if p]
        if not valid:
            return cluster_id, False, "no valid read paths"
        cmd = "cat {} > {}".format(
            " ".join('"{}"'.format(p) for p in valid), out
        )
        ret = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
        if ret.returncode != 0:
            return cluster_id, False, ret.stderr.decode().strip()

    return cluster_id, True, None


def merge_all(cluster_reads, outdir, round_num, threads):
    log.info("Merging reads for %d clusters (threads=%d)...", len(cluster_reads), threads)
    results = {}
    with ThreadPoolExecutor(max_workers=threads) as pool:
        futures = {
            pool.submit(merge_one_cluster, cid, d["r1"], d["r2"], outdir, round_num): cid
            for cid, d in cluster_reads.items()
        }
        for future in as_completed(futures):
            cid, ok, err = future.result()
            results[cid] = ok
            if ok:
                log.info("  [ok] cluster_%s", cid)
            else:
                log.error("  [fail] cluster_%s: %s", cid, err)

    n_ok = sum(results.values())
    log.info("Merge done: %d/%d succeeded", n_ok, len(results))
    if n_ok < len(results):
        log.error("%d clusters failed to merge", len(results) - n_ok)
        sys.exit(1)
    return results


# ─────────────────────────────────────────────
# JSON builder
# ─────────────────────────────────────────────

def build_json(new_clusters, cluster_lookup, sag_info, round_num, outdir):
    cluster_list  = []
    cluster_reads = {}

    prefix = "R{}_C".format(round_num)

    missing_members = []
    missing_sags    = []

    for cluster_id in sorted(new_clusters.keys(),
                             key=lambda x: int(x) if x.isdigit() else x):
        member_ids = new_clusters[cluster_id]

        members       = []
        leaf_sags     = []          # 现在存 dict，而非 str
        leaf_sag_ids  = set()       # 去重用
        r1_list       = []
        r2_list       = []

        for mid in member_ids:
            if mid not in cluster_lookup:
                missing_members.append(mid)
                log.warning("  member %s not found in prev_json", mid)
                continue

            prev_cluster = cluster_lookup[mid]

            # 展平 leaf_sags（兼容旧格式字符串列表 & 新格式对象列表）
            for entry in prev_cluster.get("leaf_sags", []):
                sag_id = entry if isinstance(entry, str) else entry.get("sag_id")
                if not sag_id or sag_id in leaf_sag_ids:
                    continue
                leaf_sag_ids.add(sag_id)

                # 构造 leaf_sag 对象，包含原始 read1/read2/checkm2
                if sag_id in sag_info:
                    info = sag_info[sag_id]
                    leaf_sags.append({
                        "sag_id":  sag_id,
                        "read1":   info["read1"],
                        "read2":   info["read2"],
                        "checkm2": info.get("checkm2"),
                        "individual_contigs": info.get("individual_contigs"),
                    })
                    r1_list.append(info["read1"])
                    r2_list.append(info["read2"])
                else:
                    leaf_sags.append({
                        "sag_id": sag_id,
                        "read1": None,
                        "read2": None,
                        "checkm2": None,
                        "individual_contigs": None,
                    })
                    if sag_id not in missing_sags:
                        missing_sags.append(sag_id)

            members.append({
                "cluster_id":         mid,
                "cluster_size":       prev_cluster.get("cluster_size", 0),
                "leaf_sags":          prev_cluster.get("leaf_sags", []),
                "coassembly_contigs": prev_cluster.get("coassembly_contigs"),
                "checkm2":            prev_cluster.get("checkm2_coassembly"),
                # Carry forward previous-round CoSAG taxonomy so round2 can
                # compare member-cluster taxonomy consistency directly.
                "gtdbtk_classification": (
                    (prev_cluster.get("co_assembly_results") or {}).get("gtdbtk_classification")
                ),
            })

        new_cid     = "{}{}".format(prefix, cluster_id)
        cluster_dir = (Path(outdir) / "cluster_{}".format(cluster_id)).resolve()
        cluster_reads[cluster_id] = {"r1": r1_list, "r2": r2_list}

        cluster_list.append({
            "cluster_id":         new_cid,
            "cluster_size":       len(member_ids),
            "members":            members,
            "leaf_sags":          leaf_sags,      # ← 现在是带详细信息的对象列表
            "merged_reads": {
                "read1": str(cluster_dir / "merged_R1.fastq.gz"),
                "read2": str(cluster_dir / "merged_R2.fastq.gz"),
            },
            "coassembly_contigs": None,
            "checkm2_coassembly": None,
            "status":             "pending",
        })

    if missing_members:
        log.warning("  %d members not in prev_json (first 5): %s",
                    len(missing_members), missing_members[:5])
    if missing_sags:
        log.warning("  %d SAGs missing reads (first 5): %s",
                    len(missing_sags), missing_sags[:5])

    # 兼容新格式（对象列表）统计 n_sags
    all_leaf_sag_ids = set()
    for c in cluster_list:
        for entry in c["leaf_sags"]:
            sid = entry if isinstance(entry, str) else entry.get("sag_id")
            if sid:
                all_leaf_sag_ids.add(sid)

    data = {
        "version":     "1.0",
        "round":       round_num,
        "member_type": "cluster",
        "n_clusters":  len(cluster_list),
        "n_sags":      len(all_leaf_sag_ids),
        "clusters":    cluster_list,
    }

    return data, cluster_reads


# ─────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prev-json",   required=True,  dest="prev_json")
    parser.add_argument("--clusters",    required=True)
    parser.add_argument("--round",       required=True,  type=int)
    parser.add_argument("--round1-json", required=False, default=None,
                        dest="round1_json",
                        help="Round 1 JSON，Round 3+ 时用于溯源原始 SAG reads")
    parser.add_argument("--outdir",      required=True)
    parser.add_argument("--out",         default="round_clusters.json")
    parser.add_argument("--threads",     type=int, default=8)
    args = parser.parse_args()

    log.info("Loading inputs (round %d)...", args.round)

    # Round 3+ 需要 round1_json
    if args.round >= 3 and not args.round1_json:
        log.error("--round1-json is required for round >= 3")
        sys.exit(1)

    # 先加载 Round 1 的 SAG 完整信息（Round 3+ 用）
    round1_sag_reads = load_sag_reads_from_round1(args.round1_json)

    # 加载上一轮 JSON
    cluster_lookup, sag_info, prev_round = load_prev_json(
        args.prev_json, round1_sag_reads
    )

    new_clusters = parse_clusters(args.clusters)

    log.info("Building JSON...")
    data, cluster_reads = build_json(
        new_clusters, cluster_lookup, sag_info, args.round, args.outdir
    )

    log.info("Merging reads...")
    merge_all(cluster_reads, args.outdir, args.round, args.threads)

    out_path = Path(args.out)
    with open(str(out_path), "w") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)

    log.info("Done -> %s  (clusters: %d, SAGs: %d)",
             out_path, data["n_clusters"], data["n_sags"])


if __name__ == "__main__":
    main()