#!/usr/bin/env python3
"""
irma_nextclade.py - IRMA assembly + Nextclade typing automation pipeline.

Usage:
    python irma_nextclade.py --fq_list fq.list --database nextclade_data --outdir output_dir

Input TSV format (tab-separated, with header):
    sample    MODULE    R1    R2

Supported IRMA modules: FLU (others reserved for future).
"""

import argparse
import csv
import json
import logging
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Optional

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# IRMA module -> {segment_number: segment_name}
# Empty dict means single-segment genome, segment name = module name
MODULE_SEGMENT_MAP: dict[str, dict[int, str]] = {
    "FLU": {
        1: "PB2",
        2: "PB1",
        3: "PA",
        4: "HA",
        5: "NP",
        6: "NA",
        7: "MP",
        8: "NS",
    },
    "RSV": {},  # single segment, uses module name "RSV"
    # Reserved for future support
    "CoV": {},
    "EBOLA": {},
    "FLU_AD": {},
}

SUPPORTED_MODULES = {"FLU", "RSV", "CoV", "HMPV","HPIV","CHIKV"}

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

logger = logging.getLogger("irma_nextclade")


def setup_logging(outdir: Path) -> None:
    log_file = outdir / "irma_nextclade.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(str(log_file), encoding="utf-8"),
            logging.StreamHandler(sys.stdout),
        ],
    )


# ---------------------------------------------------------------------------
# Helper: run external command
# ---------------------------------------------------------------------------


def run_cmd(cmd: str, description: str = "") -> subprocess.CompletedProcess:
    """Run a shell command, log it, and raise on failure."""
    logger.info(f"CMD: {cmd}")
    try:
        result = subprocess.run(
            cmd, shell=True, capture_output=True, text=True, check=True
        )
        if result.stdout.strip():
            logger.debug(result.stdout.strip())
        return result
    except subprocess.CalledProcessError as e:
        logger.error(f"FAILED ({description}): {cmd}")
        logger.error(f"  stderr: {e.stderr.strip()}")
        raise


# ---------------------------------------------------------------------------
# Step 1: fastp QC
# ---------------------------------------------------------------------------


def parse_fastp_json(json_path: Path) -> dict:
    """Extract QC metrics from fastp JSON output."""
    with json_path.open() as fh:
        data = json.load(fh)

    summary = data.get("summary", {})
    before = summary.get("before_filtering", {})
    after = summary.get("after_filtering", {})
    dup = data.get("duplication", {})
    insert_size = data.get("insert_size", {})

    # Library size: try common key names
    lib_size = insert_size.get("peak", 0)
    if not lib_size:
        lib_size = data.get("adapter_cutting", {}).get("insert_size_peak", 0)

    return {
        "Reads_Num": before.get("total_reads", 0),
        "Base_Num": before.get("total_bases", 0),
        "Q20": round(before.get("q20_rate", 0) * 100, 2),
        "Q30": round(before.get("q30_rate", 0) * 100, 2),
        "Clean_Reads_Num": after.get("total_reads", 0),
        "Clean_Base_Num": after.get("total_bases", 0),
        "Clean_Q20": round(after.get("q20_rate", 0) * 100, 2),
        "Clean_Q30": round(after.get("q30_rate", 0) * 100, 2),
        "Duplication_Rate": round(dup.get("rate", 0) * 100, 2),
        "Library_Size": lib_size,
    }


def step_fastp(
    sample: str, r1: Path, r2: Path, fastp_dir: Path
) -> Optional[tuple[dict, Path, Path]]:
    """
    Run fastp and return (QC metrics dict, C1 path, C2 path).

    fastp outputs to fastp_dir (a temporary directory) to avoid
    polluting the IRMA output directory (IRMA requires an empty or
    non-existent target directory).

    Returns None on failure.
    """
    fastp_dir.mkdir(parents=True, exist_ok=True)
    c1 = fastp_dir / "C1.fq.gz"
    c2 = fastp_dir / "C2.fq.gz"
    json_out = fastp_dir / "fastp.json"
    html_out = fastp_dir / "fastp.html"

    cmd = f"fastp -i {r1} -I {r2} -o {c1} -O {c2} -j {json_out} -h {html_out} --detect_adapter_for_pe"
    try:
        run_cmd(cmd, f"fastp QC for {sample}")
    except subprocess.CalledProcessError:
        logger.error(f"fastp failed for sample {sample}, skipping.")
        return None

    qc_data = parse_fastp_json(json_out)
    return qc_data, c1, c2


# ---------------------------------------------------------------------------
# Step 2: IRMA assembly
# ---------------------------------------------------------------------------


def step_irma(
    sample: str, module: str, irma_dir: Path, c1: Path, c2: Path,
    config_path: Path,
) -> bool:
    """
    Run IRMA assembly. Returns True on success.

    IRMA requires its output directory to not exist or be empty (otherwise
    it creates a versioned suffix like '-V2'). So we ensure irma_dir
    does not exist before running.
    """
    # IRMA will create irma_dir itself. Remove it if it already exists
    # (from a previous run), otherwise IRMA appends '-V2' suffix.
    if irma_dir.exists():
        shutil.rmtree(irma_dir)

    cmd = f"IRMA {module} {c1} {c2} {irma_dir} -c {config_path}"
    try:
        run_cmd(cmd, f"IRMA assembly for {sample}")
    except subprocess.CalledProcessError:
        logger.error(f"IRMA failed for sample {sample}.")
        return False

    # Verify output exists
    consensus_dir = irma_dir / "amended_consensus"
    if not consensus_dir.is_dir():
        logger.error(f"IRMA output not found: {consensus_dir}")
        return False

    return True


# ---------------------------------------------------------------------------
# Step 3: FASTA rename & save
# ---------------------------------------------------------------------------


def extract_segment_number(filename: str) -> Optional[int]:
    """Extract segment number from IRMA filename like 'flu_b1_3.fa'."""
    match = re.search(r"_(\d+)\.fa$", filename)
    if match:
        return int(match.group(1))
    return None


def step_rename_fasta(
    sample: str, module: str, irma_dir: Path, fasta_dir: Path
) -> Optional[Path]:
    """
    Read IRMA amended_consensus FASTAs, rename seqIDs, write per-sample FASTA.

    For multi-segment modules (e.g. FLU), files are numbered (sample_1.fa)
    and mapped via MODULE_SEGMENT_MAP.
    For single-segment modules (e.g. RSV, empty segment map), pick
    {sample}.pad.fa and {sample}.fa

    Returns path to the renamed FASTA, or None on failure.
    """
    seg_map = MODULE_SEGMENT_MAP.get(module, {})
    is_multi_segment = bool(seg_map)

    consensus_dir = irma_dir / "amended_consensus"
    if not consensus_dir.is_dir():
        logger.error(f"amended_consensus not found: {consensus_dir}")
        return None

    renamed_records: list[SeqRecord] = []

    if is_multi_segment:
        fa_files = sorted(f for f in consensus_dir.iterdir() if f.suffix == ".fa")
        for fa_path in fa_files:
            seg_num = extract_segment_number(fa_path.name)
            if seg_num is None or seg_num not in seg_map:
                logger.warning(f"Unknown segment number in {fa_path.name}, skipping.")
                continue
            seg_name = seg_map[seg_num]
            for record in SeqIO.parse(str(fa_path), "fasta"):
                record.id = f"{sample}_{seg_name}"
                record.name = f"{sample}_{seg_name}"
                record.description = ""
                renamed_records.append(record)
    else:
        pad_fa = consensus_dir / f"{sample}.pad.fa"
        fa = consensus_dir / f"{sample}.fa"

        use_fa = pad_fa if pad_fa.is_file() else fa
        if not use_fa.is_file():
            logger.error(f"No consensus FASTA found for {sample} in {consensus_dir}")
            return None
        else:
            logger.info(f"Using consensus FASTA for {sample}: {use_fa.name}")
        for record in SeqIO.parse(str(use_fa), "fasta"):
            record.id = f"{sample}_{module}"
            record.name = f"{sample}_{module}"
            record.description = ""
            renamed_records.append(record)

    if not renamed_records:
        logger.warning(f"No consensus sequences found for {sample}.")
        return None

    fasta_dir.mkdir(parents=True, exist_ok=True)
    out_path = fasta_dir / f"{sample}.fasta"
    SeqIO.write(renamed_records, str(out_path), "fasta")
    logger.info(f"Renamed FASTA written: {out_path} ({len(renamed_records)} sequences)")
    return out_path


# ---------------------------------------------------------------------------
# Step 4: Collect coverage diagrams
# ---------------------------------------------------------------------------


def collect_coverage_plots(sample: str, irma_dir: Path, plot_dir: Path) -> None:
    """
    Copy IRMA coverage diagram PDFs to outdir/coverage_plot/.

    Searches irma_dir/figures/*coverageDiagram.pdf, copies each with
    sample name prefixed to avoid filename collisions across samples.
    """
    figures_dir = irma_dir / "figures"
    if not figures_dir.is_dir():
        logger.warning(f"figures directory not found: {figures_dir}")
        return

    pdfs = sorted(figures_dir.glob("*coverageDiagram.pdf"))
    if not pdfs:
        logger.warning(f"No coverageDiagram.pdf found in {figures_dir}")
        return

    plot_dir.mkdir(parents=True, exist_ok=True)
    for src in pdfs:
        dst = plot_dir / f"{sample}_{src.name}"
        shutil.copy2(src, dst)

    logger.info(f"Copied {len(pdfs)} coverage plots to {plot_dir}")


# ---------------------------------------------------------------------------
# Step 7: Alignment rate (minimap2 + samtools flagstat)
#
# Uses consensus sequences from IRMA output (amended_consensus/*.fa).
# ---------------------------------------------------------------------------


def step_alignment_rate(
    sample: str,
    run_dir: Path,
    c1: Path,
    c2: Path,
    irma_dir: Path,
) -> Optional[dict]:
    """
    Collect all FASTA files from IRMA output directory,
    merge into one reference FASTA sorted by sequence length,
    align clean reads with minimap2,
    and compute alignment rate using samtools flagstat.

    Args:
        sample: Sample name.
        run_dir: Per-sample working directory (outdir/run_{sample}/).
        c1, c2: Clean reads from fastp.
        irma_dir: IRMA output directory (outdir/run_{sample}/{sample}/).

    Returns dict with alignment stats, or None on failure.
    """
    align_dir = run_dir / "alignment"
    align_dir.mkdir(parents=True, exist_ok=True)

    # Collect all .fasta files directly under irma_dir (not recursive)
    ref_fa = align_dir / "reference.fasta"
    records: list[SeqRecord] = []

    for fa_path in irma_dir.glob("*.fasta"):
        for rec in SeqIO.parse(str(fa_path), "fasta"):
            records.append(rec)

    if not records:
        logger.error(f"No reference sequences found in {irma_dir} for {sample}")
        return None

    # Sort by sequence length (descending)
    records.sort(key=lambda r: len(str(r.seq)), reverse=True)

    SeqIO.write(records, str(ref_fa), "fasta")
    logger.info(f"Reference from IRMA: {ref_fa} ({len(records)} sequences)")

    # minimap2 index
    ref_mmi = align_dir / "reference.mmi"
    cmd = f"minimap2 -d {ref_mmi} {ref_fa}"
    try:
        run_cmd(cmd, f"minimap2 index for {sample}")
    except subprocess.CalledProcessError:
        logger.error(f"minimap2 index failed for {sample}")
        return None

    # minimap2 align + sort BAM
    bam_out = align_dir / "aligned.bam"
    cmd = (
        f"minimap2 -ax sr {ref_mmi} {c1} {c2} | "
        f"samtools sort -o {bam_out}"
    )
    try:
        run_cmd(cmd, f"minimap2 align for {sample}")
    except subprocess.CalledProcessError:
        logger.error(f"minimap2 align failed for {sample}")
        return None

    # samtools index
    cmd = f"samtools index {bam_out}"
    try:
        run_cmd(cmd, f"samtools index for {sample}")
    except subprocess.CalledProcessError:
        logger.warning(f"samtools index failed for {sample}, continuing.")

    # samtools flagstat --output-fmt json
    flagstat_json = align_dir / "flagstat.json"
    cmd = f"samtools flagstat --output-fmt json {bam_out} > {flagstat_json}"
    try:
        run_cmd(cmd, f"samtools flagstat for {sample}")
    except subprocess.CalledProcessError:
        logger.error(f"samtools flagstat failed for {sample}")
        return None

    # samtools coverage -Q 20
    coverage_tsv = align_dir / "coverage.tsv"
    cmd = f"samtools coverage -Q 20 {bam_out} > {coverage_tsv}"
    try:
        run_cmd(cmd, f"samtools coverage for {sample}")
    except subprocess.CalledProcessError:
        logger.warning(f"samtools coverage failed for {sample}")
        coverage_tsv = None

    # Parse flagstat JSON
    with flagstat_json.open() as fh:
        stats = json.load(fh)

    qc_passed = stats.get("QC-passed reads", {})
    alignment_rate = qc_passed.get("mapped %", 0) or 0

    # Parse samtools coverage -> overall coverage = sum(covbases) / sum(endpos)
    overall_coverage = 0.0
    if coverage_tsv is not None and coverage_tsv.is_file():
        total_covbases = 0
        total_endpos = 0
        with coverage_tsv.open() as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                total_covbases += int(row.get("covbases", 0))
                total_endpos += int(row.get("endpos", 0))
        if total_endpos > 0:
            overall_coverage = round(total_covbases / total_endpos * 100, 2)

    result = {
        "Alignment_Rate": alignment_rate,
        "Coverage": overall_coverage,
    }
    logger.info(f"Alignment rate for {sample}: {alignment_rate}%, Coverage: {overall_coverage}%")
    return result


# ---------------------------------------------------------------------------
# Step 5: nextclade sort
# ---------------------------------------------------------------------------


def parse_sort_tsv(tsv_path: Path) -> dict[str, str]:
    """Parse nextclade sort TSV -> {seqName: dataset}."""
    result = {}
    with tsv_path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            seq_name = row.get("seqName", "")
            dataset = row.get("dataset", "")
            if seq_name and dataset:
                result[seq_name] = dataset
    return result


def ensure_minimizer_index(database_path: Path) -> Optional[Path]:
    """
    Ensure the global nextclade minimizer_index.json exists locally.

    Downloads from Nextstrain data server if missing.
    Returns the path to the minimizer index, or None on failure.
    """
    minimizer_index_path = database_path / "minimizer_index.json"
    if minimizer_index_path.is_file():
        return minimizer_index_path

    logger.warning(f"Downloading minimizer_index.json to {database_path} to speed up sort")
    database_path.mkdir(parents=True, exist_ok=True)
    cmd = (
        f"curl -fsSLo {minimizer_index_path} "
        f"https://data.clades.nextstrain.org/v3/minimizer_index.json"
    )
    try:
        run_cmd(cmd, "download minimizer index")
    except subprocess.CalledProcessError:
        logger.error("Failed to download minimizer index")
        return None

    if not minimizer_index_path.is_file():
        logger.error("minimizer_index.json not found after download")
        return None

    return minimizer_index_path


def step_nextclade_sort(
    sample: str, run_dir: Path, input_fasta: Path, database_path: Path
) -> Optional[dict[str, str]]:
    """
    Run nextclade sort on per-sample FASTA.

    Downloads the global minimizer_index.json if absent and passes it
    via ``-m`` to accelerate dataset assignment.

    Returns {seqName: dataset}, or None on failure.
    """
    sort_tsv = run_dir / "sort.tsv"

    minimizer_index = ensure_minimizer_index(database_path)
    if minimizer_index is not None:
        cmd = f"nextclade sort -m {minimizer_index} -r {sort_tsv} {input_fasta}"
    else:
        logger.warning("Falling back to nextclade sort without minimizer index")
        cmd = f"nextclade sort -r {sort_tsv} {input_fasta}"

    try:
        run_cmd(cmd, f"nextclade sort for {sample}")
    except subprocess.CalledProcessError:
        logger.error(f"nextclade sort failed for {sample}.")
        return None

    if not sort_tsv.is_file():
        logger.error(f"sort TSV not found: {sort_tsv}")
        return None

    seq_to_dataset = parse_sort_tsv(sort_tsv)
    logger.info(f"nextclade sort: {len(seq_to_dataset)} sequences assigned to datasets")
    return seq_to_dataset


# ---------------------------------------------------------------------------
# Step 6: nextclade dataset management
# ---------------------------------------------------------------------------


def safe_dataset_name(dataset_name: str) -> str:
    """Convert dataset path to filesystem-safe name for filenames."""
    return dataset_name.replace("/", "_")


def ensure_dataset(database_path: Path, dataset_name: str) -> Optional[Path]:
    """
    Ensure a nextclade dataset exists locally.

    Checks database_path/dataset_name/reference.fasta.
    Downloads if missing.

    Returns the dataset directory path, or None on failure.
    """
    dataset_dir = database_path / dataset_name
    ref_fasta = dataset_dir / "reference.fasta"

    if ref_fasta.is_file():
        logger.debug(f"Dataset exists: {dataset_dir}")
        return dataset_dir

    logger.info(f"Dataset not found, downloading: {dataset_name}")
    dataset_dir.mkdir(parents=True, exist_ok=True)

    cmd = (
        f"nextclade dataset get "
        f"--name {dataset_name} "
        f"--output-dir {dataset_dir}"
    )
    try:
        run_cmd(cmd, f"download dataset {dataset_name}")
    except subprocess.CalledProcessError:
        logger.error(f"Failed to download dataset: {dataset_name}")
        return None

    if not ref_fasta.is_file():
        logger.error(f"reference.fasta still missing after download: {dataset_dir}")
        return None

    return dataset_dir


# ---------------------------------------------------------------------------
# Step 6: nextclade run
# ---------------------------------------------------------------------------


def step_nextclade_run(
    sample: str,
    run_dir: Path,
    input_fasta: Path,
    seq_to_dataset: dict[str, str],
    database_path: Path,
) -> list[dict]:
    """
    Run nextclade run for each dataset group.

    Returns list of result dicts with keys:
        seqName, segment, dataset, clade, coverage, lenAligned
    """
    # Group sequences by dataset
    dataset_seqs: dict[str, list[str]] = defaultdict(list)
    for seq_name, dataset_name in seq_to_dataset.items():
        dataset_seqs[dataset_name].append(seq_name)

    nextclade_dir = run_dir / "nextclade"
    nextclade_dir.mkdir(parents=True, exist_ok=True)

    all_results: list[dict] = []

    # Read all sequences from the input FASTA into memory for subsetting
    seq_records: dict[str, SeqRecord] = {}
    for record in SeqIO.parse(str(input_fasta), "fasta"):
        seq_records[record.id] = record

    for dataset_name, seq_names in dataset_seqs.items():
        # Ensure dataset is available
        dataset_dir = ensure_dataset(database_path, dataset_name)
        if dataset_dir is None:
            logger.warning(
                f"Skipping dataset {dataset_name} for sample {sample}: download failed"
            )
            continue

        # Collect sequences belonging to this dataset
        records_for_dataset: list[SeqRecord] = []
        for sn in seq_names:
            if sn in seq_records:
                records_for_dataset.append(seq_records[sn])
            else:
                logger.warning(f"seqName {sn} not found in FASTA, skipping.")

        if not records_for_dataset:
            continue

        # Write temporary FASTA for this dataset group
        safe_name = safe_dataset_name(dataset_name)
        temp_fasta = nextclade_dir / f"{safe_name}_input.fasta"
        SeqIO.write(records_for_dataset, str(temp_fasta), "fasta")

        tsv_out = nextclade_dir / f"{safe_name}.tsv"
        json_out = nextclade_dir / f"{safe_name}.json"

        cmd = (
            f"nextclade run "
            f"-D {dataset_dir} "
            f"-t {tsv_out} "
            f"-J {json_out} "
            f"{temp_fasta}"
        )
        try:
            run_cmd(cmd, f"nextclade run {dataset_name} for {sample}")
        except subprocess.CalledProcessError:
            logger.error(f"nextclade run failed: {dataset_name} for {sample}")
            continue

        # Parse JSON results
        if not json_out.is_file():
            logger.error(f"nextclade JSON output not found: {json_out}")
            continue

        with json_out.open() as fh:
            json_data = json.load(fh)

        for r in json_data.get("results", []):
            seq_name = r.get("seqName", "")
            clade = r.get("clade") or ""
            coverage = r.get("coverage", 0)
            len_aligned = r.get("lenAligned", 0)

            # Extract segment name from seqName (e.g., "flu_b1_PB2" -> "PB2")
            segment = seq_name.split("_")[-1] if "_" in seq_name else seq_name

            all_results.append(
                {
                    "seqName": seq_name,
                    "segment": segment,
                    "dataset": dataset_name,
                    "clade": clade,
                    "coverage": round(coverage, 4),
                    "lenAligned": len_aligned,
                }
            )

    return all_results


# ---------------------------------------------------------------------------
# Step 8: Aggregate results
# ---------------------------------------------------------------------------


def compute_total_coverage(results: list[dict]) -> float:
    """
    Compute total coverage across all segments.

    total_coverage = sum(lenAligned_i) / sum(lenAligned_i / coverage_i)

    This is equivalent to a reference-length-weighted average of per-segment coverage.
    """
    total_aligned = 0
    total_ref = 0.0
    for r in results:
        aligned = r["lenAligned"]
        cov = r["coverage"]
        if aligned > 0 and cov > 0:
            ref_len = aligned / cov
            total_aligned += aligned
            total_ref += ref_len
        elif aligned > 0:
            total_aligned += aligned
            total_ref += aligned

    if total_ref == 0:
        return 0.0
    return round(total_aligned / total_ref, 4)


def write_qc_tsv(all_qc: list[dict], out_path: Path) -> None:
    """Write aggregated QC TSV."""
    columns = [
        "sample_name",
        "Reads_Num",
        "Base_Num",
        "Q20",
        "Q30",
        "Clean_Reads_Num",
        "Clean_Base_Num",
        "Clean_Q20",
        "Clean_Q30",
        "Duplication_Rate",
        "Library_Size",
        "Alignment_Rate",
        "Coverage",
    ]
    with out_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for row in all_qc:
            writer.writerow(row)
    logger.info(f"QC results written: {out_path}")


def write_nextclade_tsv(
    sample: str, results: list[dict], out_path: Path
) -> None:
    """Write nextclade results for one sample to its own TSV file."""
    columns = ["sample_name", "segment", "dataset", "clade", "coverage"]

    with out_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t")
        writer.writeheader()

        # Per-segment rows
        for r in results:
            writer.writerow(
                {
                    "sample_name": sample,
                    "segment": r["segment"],
                    "dataset": r["dataset"],
                    "clade": r["clade"],
                    "coverage": r["coverage"],
                }
            )

        # Total coverage row
        total_cov = compute_total_coverage(results)
        writer.writerow(
            {
                "sample_name": sample,
                "segment": "TOTAL",
                "dataset": "-",
                "clade": "-",
                "coverage": total_cov,
            }
        )


def merge_nextclade_tsv(sample_tsvs: list[Path], merged_path: Path) -> None:
    """Merge per-sample nextclade TSVs into one final TSV."""
    with merged_path.open("w") as out_fh:
        for i, tsv_path in enumerate(sample_tsvs):
            if not tsv_path.is_file():
                continue
            with tsv_path.open() as in_fh:
                header = in_fh.readline()  # always skip header
                if i == 0:
                    out_fh.write(header)
                for line in in_fh:
                    out_fh.write(line)
    logger.info(f"Merged nextclade results: {merged_path}")


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------


def read_fq_list(fq_list_path: Path) -> list[dict]:
    """Read input sample list TSV."""
    samples = []
    with fq_list_path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            # Normalize column names (strip whitespace, case-insensitive)
            norm_row = {k.strip().upper(): v.strip() for k, v in row.items()}
            sample = norm_row.get("SAMPLE", "")
            module = norm_row.get("MODULE", "")
            r1 = norm_row.get("R1", "")
            r2 = norm_row.get("R2", "")
            if not all([sample, module, r1, r2]):
                logger.warning(f"Skipping incomplete row: {row}")
                continue
            samples.append(
                {"sample": sample, "module": module, "R1": r1, "R2": r2}
            )
    return samples


def process_sample(
    sample_info: dict,
    outdir: Path,
    database_path: Path,
    irma_config: Path,
) -> Optional[tuple[dict, list[dict]]]:
    """
    Process a single sample through the full pipeline.

    Returns (qc_data, nextclade_results) or None on critical failure.
    """
    sample = sample_info["sample"]
    module = sample_info["module"]
    r1 = Path(sample_info["R1"])
    r2 = Path(sample_info["R2"])

    # Directory layout:
    #   run_dir    = outdir/run_{sample}/           (working dir)
    #   irma_dir   = outdir/run_{sample}/{sample}/  (IRMA output only)
    #   fastp_dir  = outdir/run_{sample}/fastp/     (fastp output)
    #   sort.tsv   = outdir/run_{sample}/sort.tsv
    #   nextclade/ = outdir/run_{sample}/nextclade/  (nextclade run results)
    run_dir = outdir / f"run_{sample}"
    irma_dir = run_dir / sample
    fastp_dir = run_dir / "fastp"
    fasta_dir = outdir / "fasta"
    run_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"{'=' * 60}")
    logger.info(f"Processing sample: {sample} (module={module})")

    # Validate module
    if module not in SUPPORTED_MODULES:
        logger.error(f"Unsupported IRMA module '{module}' for sample {sample}. Supported: {SUPPORTED_MODULES}")
        return None

    # Step 1: fastp (output to temp dir to keep irma_dir clean for IRMA)
    logger.info(f"[Step 1/7] fastp QC for {sample}")
    fastp_result = step_fastp(sample, r1, r2, fastp_dir)
    if fastp_result is None:
        logger.error(f"Skipping sample {sample} due to fastp failure.")
        if fastp_dir.is_dir():
            shutil.rmtree(fastp_dir)
        return None
    qc_data, c1, c2 = fastp_result
    qc_data["sample_name"] = sample

    # Step 2: IRMA (irma_dir = run_dir/sample, must not exist)
    logger.info(f"[Step 2/7] IRMA assembly for {sample}")
    if not step_irma(sample, module, irma_dir, c1, c2, irma_config):
        logger.error(f"Skipping sample {sample} due to IRMA failure.")
        return None

    # Step 3: Rename FASTA
    logger.info(f"[Step 3/7] Renaming consensus FASTA for {sample}")
    renamed_fasta = step_rename_fasta(sample, module, irma_dir, fasta_dir)
    if renamed_fasta is None:
        logger.error(f"No consensus sequences for {sample}, skipping nextclade.")
        return qc_data, []

    # Step 4: Collect coverage diagrams
    logger.info(f"[Step 4/7] Collecting coverage plots for {sample}")
    collect_coverage_plots(sample, irma_dir, outdir / "coverage_plot")

    # Step 5: nextclade sort (output to run_dir)
    logger.info(f"[Step 5/7] nextclade sort for {sample}")
    seq_to_dataset = step_nextclade_sort(sample, run_dir, renamed_fasta, database_path)
    if seq_to_dataset is None:
        logger.warning(f"nextclade sort failed for {sample}, skipping nextclade run + alignment.")
        qc_data["Alignment_Rate"] = 0
        qc_data["Coverage"] = 0
        return qc_data, []

    # Step 6: nextclade run (output to run_dir/nextclade)
    logger.info(f"[Step 6/7] nextclade run for {sample}")
    nextclade_results = step_nextclade_run(
        sample, run_dir, renamed_fasta, seq_to_dataset, database_path
    )

    # Step 7: Alignment rate — uses consensus sequences from IRMA results
    logger.info(f"[Step 7/7] Alignment rate for {sample}")
    align_result = step_alignment_rate(
        sample, run_dir, c1, c2, irma_dir
    )
    if align_result is not None:
        qc_data["Alignment_Rate"] = align_result["Alignment_Rate"]
        qc_data["Coverage"] = align_result["Coverage"]
    else:
        qc_data["Alignment_Rate"] = 0
        qc_data["Coverage"] = 0

    return qc_data, nextclade_results


def main():
    parser = argparse.ArgumentParser(
        description="IRMA assembly + Nextclade typing automation pipeline"
    )
    parser.add_argument(
        "--fq_list", required=True, help="Input sample list TSV file"
    )
    parser.add_argument(
        "--database",
        default=None,
        help="Path to nextclade dataset directory (default: nextclade_data in script dir)",
    )
    parser.add_argument(
        "--outdir", required=True, help="Output directory"
    )
    args = parser.parse_args()

    # Resolve default database path
    if args.database is None:
        script_dir = Path(__file__).resolve().parent
        args.database = script_dir / "nextclade_data"
    database_path = Path(args.database).resolve()

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    # Create shared TMPDIR and IRMA config (redirect TMP away from /tmp)
    tmp_dir = outdir / "TMPDIR"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    irma_config = outdir / "irma_config.cfg"
    irma_config.write_text(f"TMP={tmp_dir}\n")
    logger.info(f"IRMA config: {irma_config} (TMP={tmp_dir})")

    setup_logging(outdir)
    logger.info(f"Pipeline started")
    logger.info(f"  fq_list  : {args.fq_list}")
    logger.info(f"  database : {database_path}")
    logger.info(f"  outdir   : {outdir}")

    # Read input
    samples = read_fq_list(Path(args.fq_list))
    if not samples:
        logger.error("No valid samples found in fq_list.")
        sys.exit(1)
    logger.info(f"Found {len(samples)} sample(s) to process")

    # Output file paths
    qc_tsv_path = outdir / "QC.tsv"
    nextclade_tsv_path = outdir / "nextclade_result.tsv"

    all_qc: list[dict] = []
    sample_nextclade_tsvs: list[Path] = []

    for sample_info in samples:
        result = process_sample(sample_info, outdir, database_path, irma_config)

        if result is None:
            logger.warning(f"Sample {sample_info['sample']} failed completely.")
            continue

        qc_data, nextclade_results = result

        # Collect QC
        all_qc.append(qc_data)

        # Write per-sample nextclade TSV
        sample = sample_info["sample"]
        sample_nc_tsv = outdir / f"run_{sample}" / "nextclade.tsv"
        write_nextclade_tsv(sample, nextclade_results, sample_nc_tsv)
        sample_nextclade_tsvs.append(sample_nc_tsv)

    # Write aggregated QC
    write_qc_tsv(all_qc, qc_tsv_path)

    # Merge all per-sample nextclade TSVs
    if sample_nextclade_tsvs:
        merge_nextclade_tsv(sample_nextclade_tsvs, nextclade_tsv_path)

    logger.info(f"{'=' * 60}")
    logger.info(f"Pipeline finished. {len(all_qc)}/{len(samples)} samples processed successfully.")
    logger.info(f"  QC results      : {qc_tsv_path}")
    logger.info(f"  Nextclade results: {nextclade_tsv_path}")


if __name__ == "__main__":
    main()
