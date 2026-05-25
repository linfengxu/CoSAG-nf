#!/usr/bin/env python3
"""
update_gtdbtk_results.py
-------------------------
Update cluster JSON with GTDB-Tk taxonomic classification results.

Supports both individual SAG results (SAG_ASSEMBLY stage) and
CoSAG co-assembly results (QUALITY_TAXONOMY stage).

The script determines which section to update based on --stage:
    sag        → cluster.members[*].gtdbtk_classification
    coassembly → cluster.co_assembly_results.gtdbtk_classification
"""

import argparse
import csv
import json
import logging
import re
import sys
from datetime import datetime
from pathlib import Path
from typing import Optional


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def load_json(path: Path) -> dict:
    with path.open(encoding="utf-8") as fh:
        return json.load(fh)


def find_cluster(data: dict, cluster_id: str) -> Optional[dict]:
    for c in data.get("clusters", []):
        if str(c.get("cluster_id")) == str(cluster_id):
            return c
    return None


def parse_taxonomy(classification_string: str) -> dict:
    """
    Parse GTDB-Tk classification string into structured dict.
    e.g. 'd__Bacteria;p__Firmicutes;c__Bacilli;...'
    """
    levels = ["domain", "phylum", "class", "order", "family", "genus", "species"]
    result = {l: "" for l in levels}
    if not classification_string:
        return result
    for i, part in enumerate(classification_string.split(";")):
        if i < len(levels):
            result[levels[i]] = re.sub(r"^[a-z]__", "", part.strip())
    return result


def load_gtdbtk_summary(path: Path) -> dict:
    """
    Load gtdbtk.bac120.summary.tsv or gtdbtk.ar53.summary.tsv.
    Returns { user_genome: { classification, taxonomy, metadata } }.
    """
    results = {}
    with path.open(encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            genome = row.get("user_genome", "").strip()
            if not genome:
                continue
            classification = row.get("classification", "")
            results[genome] = {
                "user_genome":      genome,
                "classification":   classification,
                "taxonomy":         parse_taxonomy(classification),
                "metadata": {
                    "fastani_reference":          row.get("fastani_reference", ""),
                    "fastani_ani":                row.get("fastani_ani", ""),
                    "fastani_af":                 row.get("fastani_af", ""),
                    "closest_placement_reference":row.get("closest_placement_reference", ""),
                    "closest_placement_ani":      row.get("closest_placement_ani", ""),
                    "pplacer_taxonomy":            row.get("pplacer_taxonomy", ""),
                    "classification_method":       row.get("classification_method", ""),
                    "msa_percent":                row.get("msa_percent", ""),
                    "red_value":                  row.get("red_value", ""),
                    "warnings":                   row.get("warnings", ""),
                    "note":                       row.get("note", ""),
                },
                "updated_at": datetime.now().isoformat(),
            }
    log.info("Loaded %d GTDB-Tk classifications from %s", len(results), path.name)
    return results


def extract_cluster_id(genome_name: str) -> Optional[str]:
    """Extract cluster ID from genome name like 'cluster_144_CoSAG'."""
    m = re.search(r"cluster_(\w+)", genome_name)
    return m.group(1) if m else None


def extract_sag_id(genome_name: str) -> str:
    """
    Extract SAG ID from genome name.
    Convention: <sag_id>_contigs  (e.g. AACAACACGAAACCACCACGH_contigs)
    """
    return genome_name.replace("_contigs", "").strip()


# ---------------------------------------------------------------------------
# Update functions
# ---------------------------------------------------------------------------
def update_coassembly(data: dict, gtdbtk: dict) -> list:
    """Add GTDB-Tk results to cluster.co_assembly_results.gtdbtk_classification."""
    updated = []
    for genome, info in gtdbtk.items():
        cid = extract_cluster_id(genome)
        if cid is None:
            log.warning("Cannot extract cluster ID from: %s", genome)
            continue
        cluster = find_cluster(data, cid)
        if cluster is None:
            log.warning("Cluster %s not found", cid)
            continue
        cluster.setdefault("co_assembly_results", {})
        cluster["co_assembly_results"]["gtdbtk_classification"] = info
        updated.append(cid)
        log.info(
            "Cluster %s: %s → %s %s",
            cid,
            info["taxonomy"]["phylum"],
            info["taxonomy"]["genus"],
            info["taxonomy"]["species"],
        )
    return updated


def update_sag(data: dict, gtdbtk: dict) -> list:
    """
    Add GTDB-Tk SAG results.

    Compatible with:
      - round1: cluster.members[*].sag_id
      - round2+: cluster.leaf_sags[*].sag_id (and members[*].leaf_sags lists)
    """
    # Build lookup: sag_id → classification
    sag_lookup = {}
    for genome, info in gtdbtk.items():
        sag_id = extract_sag_id(genome)
        sag_lookup[sag_id] = info

    updated_sags = []
    for cluster in data.get("clusters", []):
        cluster_updated = []
        # round2+ primary location: cluster.leaf_sags
        for leaf in cluster.get("leaf_sags", []):
            if isinstance(leaf, dict):
                sag_id = leaf.get("sag_id", "")
                if sag_id in sag_lookup:
                    leaf["gtdbtk_classification"] = sag_lookup[sag_id]
                    updated_sags.append(sag_id)
                    cluster_updated.append(sag_id)

        for member in cluster.get("members", []):
            # round1 members
            sag_id = member.get("sag_id", "")
            if sag_id in sag_lookup:
                member["gtdbtk_classification"] = sag_lookup[sag_id]
                updated_sags.append(sag_id)
                cluster_updated.append(sag_id)

            # round2+ nested member leaf list (usually list[str])
            # Keep a single source of truth at cluster.leaf_sags[*].gtdbtk_classification.
            # Do not mirror these annotations into member.leaf_gtdbtk.

        if cluster_updated:
            sid = cluster_updated[-1]
            info = sag_lookup.get(sid)
            if info:
                log.info("Cluster %s SAG annotations updated (%d)", cluster.get("cluster_id", "NA"), len(set(cluster_updated)))
    return updated_sags


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def run(
    main_json:   Path,
    bac_summary: Optional[Path],
    ar_summary:  Optional[Path],
    out_json:    Path,
    stage:       str,
) -> None:

    data = load_json(main_json)
    gtdbtk: dict = {}

    for path in [bac_summary, ar_summary]:
        if path and path.exists():
            gtdbtk.update(load_gtdbtk_summary(path))

    if not gtdbtk:
        # Round-2 standalone may not provide SAG summaries.
        # In that case, keep JSON unchanged and continue.
        log.warning("No GTDB-Tk summary files provided or found; writing unchanged JSON")
        out_json.write_text(json.dumps(data, indent=2, ensure_ascii=False))
        log.info("Written unchanged JSON: %s", out_json)
        return

    if stage == "coassembly":
        updated = update_coassembly(data, gtdbtk)
        log.info("Updated %d clusters (co-assembly stage)", len(updated))
    else:
        updated = update_sag(data, gtdbtk)
        log.info("Updated %d SAG members (sag stage)", len(updated))

    out_json.write_text(json.dumps(data, indent=2, ensure_ascii=False))
    log.info("Written: %s", out_json)

    if not updated:
        log.warning("No entries updated — check genome name conventions")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--main_json",    required=True, type=Path)
    p.add_argument("--bac_summary",  type=Path, default=None,
                   help="gtdbtk.bac120.summary.tsv")
    p.add_argument("--ar_summary",   type=Path, default=None,
                   help="gtdbtk.ar53.summary.tsv")
    p.add_argument("--out_json",     required=True, type=Path)
    p.add_argument("--stage",        required=True,
                   choices=["sag", "coassembly"],
                   help="sag = individual SAG members; "
                        "coassembly = co-assembly cluster results")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run(
        main_json   = args.main_json,
        bac_summary = args.bac_summary,
        ar_summary  = args.ar_summary,
        out_json    = args.out_json,
        stage       = args.stage,
    )
