#!/usr/bin/env python3
"""
parse_nextclade_summary.py
--------------------------

Parse nextclade_result.tsv into a concise nextclade_short.tsv summary.

Input:  nextclade_result.tsv
        Columns: sample_name  segment  dataset  clade  coverage

Output: nextclade_short.tsv
        Columns: sample_name  coverage  Clade

Rules:
  - For single-segment pathogens (RSV, SARS-CoV-2, etc.):
        Clade = <clade>           e.g. "A.D.3" or "25B"
  - For multi-segment pathogens (influenza):
        Clade = SEG(clade) SEG(clade) ...
        e.g. "HA(J.2.2) NA(B.4.2.1)"
        Only segments with a non-empty clade are listed.
  - Coverage comes from the TOTAL row.

Usage:
    python parse_nextclade_summary.py nextclade_result.tsv
    python parse_nextclade_summary.py nextclade_result.tsv -o custom_output.tsv
"""

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path


# Known influenza segments (used to detect multi-segment pathogens)
FLU_SEGMENTS = {"HA", "NA", "MP", "NS", "NP", "PA", "PB1", "PB2"}

# If a sample has more than this many unique segments (excluding TOTAL),
# it is considered multi-segment.
MULTI_SEGMENT_THRESHOLD = 2


def parse_nextclade_result(tsv_path: Path) -> list[dict]:
    """
    Parse nextclade_result.tsv and return a list of summary rows.

    Each row: {"sample_name": str, "coverage": str, "Clade": str}
    """
    # Group rows by sample_name
    sample_data: dict[str, dict] = defaultdict(
        lambda: {"segments": [], "total_coverage": ""}
    )

    with tsv_path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            sample = row.get("sample_name", "").strip()
            segment = row.get("segment", "").strip()
            clade = row.get("clade", "").strip()
            coverage = row.get("coverage", "").strip()

            if not sample:
                continue

            if segment == "TOTAL":
                sample_data[sample]["total_coverage"] = coverage
            else:
                sample_data[sample]["segments"].append(
                    {"segment": segment, "clade": clade}
                )

    results: list[dict] = []
    for sample, data in sample_data.items():
        segments = data["segments"]
        total_coverage = data["total_coverage"]

        # Determine if this is a multi-segment pathogen
        unique_segments = {s["segment"] for s in segments}
        is_multi_segment = len(unique_segments) > MULTI_SEGMENT_THRESHOLD

        if is_multi_segment:
            # Multi-segment (e.g. influenza): "HA(J.2.2) NA(B.4.2.1)"
            # Only include segments that have a clade assignment
            clade_parts = []
            for seg_info in segments:
                seg = seg_info["segment"]
                clade = seg_info["clade"]
                if clade:
                    clade_parts.append(f"{seg}({clade})")
            clade_str = " ".join(clade_parts) if clade_parts else "-"
        else:
            # Single-segment: just use the clade value directly
            clade_str = "-"
            for seg_info in segments:
                if seg_info["clade"]:
                    clade_str = seg_info["clade"]
                    break

        results.append(
            {
                "sample_name": sample,
                "coverage": total_coverage,
                "Clade": clade_str,
            }
        )

    # Sort by sample name for deterministic output
    results.sort(key=lambda r: r["sample_name"])
    return results


def write_short_tsv(rows: list[dict], out_path: Path) -> None:
    """Write summary rows to a TSV file."""
    columns = ["sample_name", "coverage", "Clade"]
    with out_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    print(f"Output written to: {out_path} ({len(rows)} samples)")


def main():
    parser = argparse.ArgumentParser(
        description="Summarize nextclade_result.tsv into nextclade_short.tsv"
    )
    parser.add_argument(
        "input",
        type=Path,
        help="Input nextclade_result.tsv file",
    )
    parser.add_argument(
        "-o", "--output",
        type=Path,
        default=None,
        help="Output TSV file (default: nextclade_short.tsv in same directory as input)",
    )
    args = parser.parse_args()

    input_path: Path = args.input.resolve()
    if not input_path.is_file():
        print(f"Error: input file not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    output_path: Path = (
        args.output.resolve() if args.output else input_path.parent / "nextclade_short.tsv"
    )

    rows = parse_nextclade_result(input_path)
    if not rows:
        print("Warning: no data rows found in input.", file=sys.stderr)
        sys.exit(1)

    write_short_tsv(rows, output_path)


if __name__ == "__main__":
    main()
