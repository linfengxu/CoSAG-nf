#!/usr/bin/env python3
"""
filter_sags.py
--------------
Stage 1: Pre-clustering filter for individual single-cell amplified genomes (SAGs).

Strategy:
    Individual SAGs suffer from MDA amplification bias → low completeness
    is expected and is NOT used as a filter criterion by default.
    By default only contamination is checked. All other filters are disabled
    by default but can be enabled via CLI arguments.

Filter criteria (default active):
    - contamination  <= --max_contamination   (default 10%)

Filter criteria (default disabled):
    - completeness   >= --min_completeness    (default 0.0,  disabled)
    - coding_density >= --min_coding_density  (default 0.0,  disabled)
    - genome_size    >= --min_genome_size     (default 0,    disabled)
    - total_contigs  <= --max_contigs         (default 99999, disabled)
"""

import argparse
import csv
import json
import logging
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Defaults — only contamination is active by default
# ---------------------------------------------------------------------------
DEFAULTS = {
    "max_contamination":  10.0,
    "min_completeness":    0.0,   # 0.0   = disabled (MDA bias: do not filter)
    "min_coding_density":  0.0,   # 0.0   = disabled
    "min_genome_size":        0,  # 0     = disabled
    "max_contigs":        99999,  # 99999 = disabled
}

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------
def classify(row: dict, cfg: dict) -> str:
    try:
        comp   = float(row["Completeness"])
        cont   = float(row["Contamination"])
        coding = float(row.get("Coding_Density", 1.0))
        gsize  = int(row.get("Genome_Size", 0))
        nctg   = int(row.get("Total_Contigs", 0))
    except (ValueError, KeyError) as e:
        return f"discard:parse_error({e})"

    # Default active
    if cont > cfg["max_contamination"]:
        return f"discard:high_contamination({cont:.2f}%)"

    # Optional filters — only active when user overrides defaults
    if cfg["min_completeness"] > 0 and comp < cfg["min_completeness"]:
        return f"discard:low_completeness({comp:.2f}%)"
    if cfg["min_coding_density"] > 0 and coding < cfg["min_coding_density"]:
        return f"discard:low_coding_density({coding:.3f})"
    if cfg["min_genome_size"] > 0 and gsize < cfg["min_genome_size"]:
        return f"discard:genome_too_small({gsize:,}bp)"
    if cfg["max_contigs"] < 99999 and nctg > cfg["max_contigs"]:
        return f"discard:too_fragmented({nctg}_contigs)"

    return "pass"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def run(
    quality_report: Path,
    output_pass:    Path,
    output_fail:    Path,
    stats_json:     Path,
    cfg:            dict,
) -> None:

    if not quality_report.exists():
        log.error("Input file not found: %s", quality_report)
        sys.exit(1)

    pass_rows, fail_rows = [], []
    discard_reasons: dict = {}

    with quality_report.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            log.error("Empty or malformed input TSV.")
            sys.exit(1)

        for row in reader:
            result = classify(row, cfg)
            row["filter_result"] = result
            if result == "pass":
                pass_rows.append(row)
            else:
                fail_rows.append(row)
                reason_key = result.split(":")[1].split("(")[0] \
                    if ":" in result else result
                discard_reasons[reason_key] = \
                    discard_reasons.get(reason_key, 0) + 1

    if pass_rows:
        _write_tsv(pass_rows, output_pass)
    else:
        output_pass.write_text("")
        log.warning("No SAGs passed filters.")

    if fail_rows:
        _write_tsv(fail_rows, output_fail)

    total = len(pass_rows) + len(fail_rows)

    # Report which filters were actually active
    active_filters   = {"max_contamination": cfg["max_contamination"]}
    disabled_filters = {}

    for key, label in [
        ("min_completeness",   "min_completeness"),
        ("min_coding_density", "min_coding_density"),
        ("min_genome_size",    "min_genome_size_bp"),
        ("max_contigs",        "max_contigs"),
    ]:
        val = cfg[key]
        is_disabled = (val == 0.0 or val == 0 or val == 99999)
        if is_disabled:
            disabled_filters[label] = "disabled (default)"
        else:
            active_filters[label] = val

    stats = {
        "total":            total,
        "passed":           len(pass_rows),
        "discarded":        len(fail_rows),
        "pass_rate_pct":    round(len(pass_rows) / total * 100, 2) if total else 0,
        "discard_reasons":  discard_reasons,
        "active_filters":   active_filters,
        "disabled_filters": disabled_filters,
    }
    stats_json.write_text(json.dumps(stats, indent=2))

    log.info("─" * 50)
    log.info("Total SAGs : %d", total)
    log.info("  Passed   : %d (%.1f%%)", len(pass_rows), stats["pass_rate_pct"])
    log.info("  Discarded: %d", len(fail_rows))
    for reason, count in sorted(discard_reasons.items(), key=lambda x: -x[1]):
        log.info("    %-35s %d", reason, count)
    log.info("Active filters  : %s", list(active_filters.keys()))
    log.info("Disabled filters: %s", list(disabled_filters.keys()))
    log.info("─" * 50)
    log.info("NOTE: completeness = 0.0 (disabled) — MDA amplification bias expected.")

    if len(pass_rows) == 0:
        log.error("No SAGs passed. Check thresholds or input quality.")
        sys.exit(1)


def _write_tsv(rows: list, path: Path) -> None:
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=list(rows[0].keys()),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--input",       required=True, type=Path,
                   help="CheckM2 quality_report.tsv")
    p.add_argument("--output_pass", required=True, type=Path,
                   help="SAGs passing filter → proceed to clustering")
    p.add_argument("--output_fail", required=True, type=Path,
                   help="Discarded SAGs with discard reason")
    p.add_argument("--stats_json",  required=True, type=Path,
                   help="Summary statistics JSON")

    # Default active
    p.add_argument(
        "--max_contamination", type=float,
        default=DEFAULTS["max_contamination"],
        help=f"Max contamination %% (default: {DEFAULTS['max_contamination']})"
    )

    # Default disabled
    p.add_argument(
        "--min_completeness", type=float,
        default=DEFAULTS["min_completeness"],
        help="Min completeness %%, 0=disabled (default: 0.0, "
             "disabled due to MDA amplification bias)"
    )
    p.add_argument(
        "--min_coding_density", type=float,
        default=DEFAULTS["min_coding_density"],
        help="Min coding density, 0=disabled (default: 0.0)"
    )
    p.add_argument(
        "--min_genome_size", type=int,
        default=DEFAULTS["min_genome_size"],
        help="Min genome size in bp, 0=disabled (default: 0)"
    )
    p.add_argument(
        "--max_contigs", type=int,
        default=DEFAULTS["max_contigs"],
        help="Max total contigs, 99999=disabled (default: 99999)"
    )

    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    cfg  = {k: getattr(args, k) for k in DEFAULTS}
    run(
        quality_report = args.input,
        output_pass    = args.output_pass,
        output_fail    = args.output_fail,
        stats_json     = args.stats_json,
        cfg            = cfg,
    )
