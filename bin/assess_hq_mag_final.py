#!/usr/bin/env python3
"""
Final HQ MAG assessment after GTDB merge: attach per-cluster hq_mag_assessment
using strict QC (completeness > threshold, contamination < threshold) plus
optional barrnap (16S / 23S / 5S rRNA) and optional tRNAscan-SE (tRNA count).

Does not change upstream extraction or taxonomy logic — metadata only.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple


def _to_float(v: Any) -> Optional[float]:
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def choose_contig_metrics(
    cluster: Dict[str, Any],
) -> Tuple[Optional[float], Optional[float], Optional[str]]:
    """Same source priority as extract_cosag_contigs.py."""
    coasm = cluster.get("co_assembly_results", {}) or {}
    tnf = coasm.get("tnf_optimized", {}) or {}
    opt = tnf.get("checkm2_optimized", {}) or {}

    opt_comp = _to_float(opt.get("completeness"))
    opt_contam = _to_float(opt.get("contamination"))

    if opt_comp is not None and opt_contam is not None:
        return opt_comp, opt_contam, "optimized"

    base = cluster.get("checkm2_coassembly", {}) or {}
    base_comp = _to_float(base.get("completeness"))
    base_contam = _to_float(base.get("contamination"))

    if base_comp is not None and base_contam is not None:
        return base_comp, base_contam, "coassembly"

    return None, None, None


def parse_barrnap_rrna(gff_path: Path) -> Tuple[bool, bool, bool]:
    """Return (has_16s, has_23s, has_5s) from barrnap GFF attributes."""
    text = gff_path.read_text(encoding="utf-8", errors="replace")
    h16 = h23 = h5 = False
    for line in text.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 9:
            continue
        attr = parts[8]
        al = attr.lower()
        if re.search(r"16s[_ ]?rrna", al) or "16s ribosomal" in al:
            h16 = True
        if re.search(r"23s[_ ]?rrna", al) or "23s ribosomal" in al:
            h23 = True
        if re.search(r"(?:^|[;=])5s[_ ]?rrna", al) or "5s ribosomal" in al:
            h5 = True
    return h16, h23, h5


def count_trna_trnascan(out_path: Path) -> int:
    """Count predicted tRNAs from tRNAscan-SE tabular output (best-effort)."""
    lines = out_path.read_text(encoding="utf-8", errors="replace").splitlines()
    n = 0
    for line in lines:
        ls = line.strip()
        if not ls or ls.startswith("#"):
            continue
        low = ls.lower()
        if "sequence" in low and "trna" in low and "bound" in low:
            continue
        parts = ls.split("\t")
        if len(parts) >= 4 and parts[0] and not parts[0].startswith("-"):
            n += 1
    return n


def load_annotation_maps(
    barrnap_dir: Path,
    trna_dir: Optional[Path],
    valid_ids: Set[str],
) -> Tuple[Dict[str, Path], Dict[str, Path]]:
    barr: Dict[str, Path] = {}
    trna: Dict[str, Path] = {}
    pat_b = re.compile(r"^cluster_(.+)_barrnap\.gff$")
    pat_t = re.compile(r"^cluster_(.+)_trnascan\.out$")
    if barrnap_dir.is_dir():
        for p in sorted(barrnap_dir.glob("*_barrnap.gff")):
            m = pat_b.match(p.name)
            if m and m.group(1) in valid_ids:
                barr[m.group(1)] = p
    if trna_dir is not None and trna_dir.is_dir():
        for p in sorted(trna_dir.glob("*_trnascan.out")):
            m = pat_t.match(p.name)
            if m and m.group(1) in valid_ids:
                trna[m.group(1)] = p
    return barr, trna


def run(
    json_in: Path,
    barrnap_dir: Path,
    trna_dir: Optional[Path],
    json_out: Path,
    completeness_gt: float,
    contamination_lt: float,
    min_trna: int,
    require_trna: bool,
) -> None:
    data = json.loads(json_in.read_text(encoding="utf-8"))
    clusters: List[Dict[str, Any]] = data.get("clusters", [])
    valid_ids = {str(c.get("cluster_id")) for c in clusters if c.get("cluster_id") is not None}

    barr_map, trna_map = load_annotation_maps(barrnap_dir, trna_dir, valid_ids)

    data["hq_mag_standard_definition"] = {
        "description_zh": (
            "高质量 SAG/MAG：completeness 严格大于阈值，contamination 严格小于阈值；"
            "具备 16S、23S、5S rRNA（barrnap）。"
            + (
                f" tRNAscan 预测 tRNA ≥ {min_trna}。"
                if require_trna
                else " 当前流程未运行 tRNAscan，不考核 tRNA 数量。"
            )
        ),
        "completeness_gt_percent": completeness_gt,
        "contamination_lt_percent": contamination_lt,
        "require_rrna_16s_23s_5s": True,
        "require_trna_evidence": require_trna,
        "min_trna_count": min_trna if require_trna else None,
        "criteria_version": "1.1",
    }

    for c in clusters:
        cid = str(c.get("cluster_id", ""))
        comp, contam, src = choose_contig_metrics(c)

        assessment: Dict[str, Any] = {
            "qc_metrics_source": src,
            "completeness_percent": comp,
            "contamination_percent": contam,
        }

        qc_comp = comp is not None and comp > completeness_gt
        qc_cont = contam is not None and contam < contamination_lt
        assessment["completeness_gt_threshold"] = qc_comp
        assessment["contamination_lt_threshold"] = qc_cont

        bf = barr_map.get(cid)
        tf = trna_map.get(cid) if require_trna else None
        assessment["barrnap_gff_available"] = bf is not None
        assessment["trnascan_out_available"] = (tf is not None) if require_trna else False

        if bf is not None:
            h16, h23, h5 = parse_barrnap_rrna(bf)
            assessment["rrna_16s_detected"] = h16
            assessment["rrna_23s_detected"] = h23
            assessment["rrna_5s_detected"] = h5
            assessment["rrna_full_set_16s_23s_5s"] = bool(h16 and h23 and h5)
        else:
            assessment["rrna_16s_detected"] = False
            assessment["rrna_23s_detected"] = False
            assessment["rrna_5s_detected"] = False
            assessment["rrna_full_set_16s_23s_5s"] = False

        if require_trna and tf is not None:
            tc = count_trna_trnascan(tf)
            assessment["trna_count"] = tc
            assessment["trna_min_threshold_met"] = tc >= min_trna
        else:
            assessment["trna_count"] = None
            assessment["trna_min_threshold_met"] = None if not require_trna else False

        if require_trna:
            full = bool(
                qc_comp
                and qc_cont
                and assessment["rrna_full_set_16s_23s_5s"]
                and assessment["trna_min_threshold_met"]
            )
        else:
            full = bool(qc_comp and qc_cont and assessment["rrna_full_set_16s_23s_5s"])
        assessment["passes_full_hq_mag_standard"] = full
        c["hq_mag_assessment"] = assessment

    json_out.parent.mkdir(parents=True, exist_ok=True)
    json_out.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--json", required=True, type=Path)
    p.add_argument("--barrnap_dir", required=True, type=Path)
    p.add_argument("--trna_dir", type=Path, default=None, help="Optional; only used if --require_trna")
    p.add_argument("--out_json", required=True, type=Path)
    p.add_argument("--completeness_gt", type=float, default=90.0)
    p.add_argument("--contamination_lt", type=float, default=5.0)
    p.add_argument("--min_trna", type=int, default=18)
    p.add_argument(
        "--require_trna",
        type=lambda x: str(x).lower() in ("1", "true", "yes"),
        default=False,
        help="If true, require tRNAscan outputs and min_trna for full HQ pass.",
    )
    return p.parse_args()


def main() -> int:
    args = parse_args()
    run(
        args.json,
        args.barrnap_dir,
        args.trna_dir,
        args.out_json,
        args.completeness_gt,
        args.contamination_lt,
        args.min_trna,
        args.require_trna,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
