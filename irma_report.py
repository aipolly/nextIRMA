#!/usr/bin/env python3
"""
irma_report.py - Generate HTML report (EN + CN) from IRMA+Nextclade pipeline output.

Usage:
    python irma_report.py --outdir output_dir

Outputs:
    outdir/report.html     (English)
    outdir/report_ch.html  (Chinese)
"""

import argparse
import csv
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Language texts: key -> (en, ch)
# ---------------------------------------------------------------------------

TEXT = {
    "title": (
        "IRMA + Nextclade Pipeline Report",
        "IRMA + Nextclade 分析流程报告",
    ),
    "samples_processed": ("Samples Processed", "处理样本数"),
    "segments_analyzed": ("Segments Analyzed", "分析片段数"),
    "pipeline_overview_title": ("Pipeline Overview", "分析流程概述"),
    "pipeline_overview": (
        "This report summarizes the results of the IRMA + Nextclade automated analysis pipeline for viral pathogen genomic data. The pipeline consists of the following steps:",
        "本报告汇总了 IRMA + Nextclade 自动化分析流程的病毒病原体基因组数据结果。流程包含以下步骤：",
    ),
    "step_fastp": (
        "Read quality control, adapter trimming, and quality filtering of raw paired-end FASTQ data.",
        "对原始双端 FASTQ 数据进行质控、接头切除和质量过滤。",
    ),
    "step_irma": (
        "Reference-guided iterative assembly of viral genomes, producing per-segment consensus sequences.",
        "基于参考序列的迭代组装，生成各基因组片段的一致性序列。",
    ),
    "step_align": (
        "Alignment of clean reads to the reference sequences from Nextclade datasets to calculate alignment rate and genome coverage.",
        "将过滤后的 clean reads 比对到 Nextclade 数据集的标准参考序列，计算比对率和基因组覆盖度。",
    ),
    "step_sort": (
        "Automatic dataset assignment for each segment based on sequence similarity.",
        "基于序列相似性自动为每个片段分配 Nextclade 数据集。",
    ),
    "step_run": (
        "Clade assignment, mutation calling, and quality checks using the assigned dataset.",
        "使用分配的数据集进行分支鉴定、突变检测和质量检查。",
    ),
    "qc_summary": ("QC Summary", "质量控制汇总"),
    "qc_explain_title": ("QC Metrics Explanation", "QC 指标说明"),
    "qc_reads_num": (
        "Total read count and base count of raw paired-end data (R1 + R2 combined).",
        "原始双端数据（R1 + R2 合计）的总 read 数和碱基数。",
    ),
    "qc_q20q30": (
        "Percentage of bases with Phred quality score &ge;20 / &ge;30 in raw data.",
        "原始数据中 Phred 质量值 &ge;20 / &ge;30 的碱基百分比。",
    ),
    "qc_clean": (
        "Read and base counts after fastp quality filtering and adapter trimming.",
        "经 fastp 质量过滤和接头切除后的 read 数和碱基数。",
    ),
    "qc_clean_q": (
        "Q20/Q30 of the filtered clean data.",
        "过滤后 clean 数据的 Q20/Q30。",
    ),
    "qc_dup": (
        "Estimated PCR duplication rate reported by fastp. High duplication (&gt;30%) may indicate over-amplification.",
        "fastp 报告的 PCR 重复率估算值。重复率过高（&gt;30%）可能表明扩增过度。",
    ),
    "qc_libsize": (
        "Insert size peak estimated by fastp from paired-end read overlap analysis.",
        "fastp 通过双端 reads 重叠分析估算的插入片段峰值。",
    ),
    "qc_align_rate": (
        "Percentage of clean reads successfully aligned to the reference sequences from Nextclade datasets (minimap2 + samtools flagstat, QC-passed reads mapped %). Reflects how well the data matches the standard reference genomes.",
        "成功比对到 Nextclade 数据集标准参考序列的 clean reads 百分比（minimap2 + samtools flagstat，QC-passed reads mapped %）。反映数据与标准参考基因组的匹配程度。",
    ),
    "qc_coverage": (
        "Genome coverage calculated as <code>sum(covbases) / sum(endpos)</code> from <code>samtools coverage -Q 20</code>. Represents the proportion of the reference genome covered by at least one base with MAPQ &ge;20.",
        "通过 <code>samtools coverage -Q 20</code> 计算的基因组覆盖度，公式为 <code>sum(covbases) / sum(endpos)</code>。表示参考基因组中被至少一个 MAPQ &ge;20 的碱基覆盖的比例。",
    ),
    "good": ("Good", "良好"),
    "acceptable": ("Acceptable", "可接受"),
    "warning_low": ("Warning", "偏低"),
    "seg_coverage_title": ("Per-Segment Coverage", "各片段覆盖度"),
    "nc_results": ("Nextclade Results", "Nextclade 分型结果"),
    "nc_explain_title": ("Nextclade Results Explanation", "Nextclade 结果说明"),
    "nc_segment": (
        "Genome segment name (e.g., PB2, HA, NA for influenza; RSV for respiratory syncytial virus). <code>TOTAL</code> row shows the weighted average coverage across all segments.",
        "基因组片段名称（如流感的 PB2、HA、NA；RSV 为全基因组）。<code>TOTAL</code> 行显示所有片段的加权平均覆盖度。",
    ),
    "nc_dataset": (
        "Nextclade dataset automatically assigned by <code>nextclade sort</code> based on sequence similarity. Determines the reference and clade nomenclature used for analysis.",
        "由 <code>nextclade sort</code> 基于序列相似性自动分配的数据集，决定分析所用的参考序列和分支命名。",
    ),
    "nc_clade": (
        "Phylogenetic clade assigned by Nextclade. Empty if the dataset does not define clades for that segment (common for internal genes). For HA/NA segments, clade names follow WHO nomenclature.",
        "Nextclade 分配的系统发育分支。若该片段的数据集未定义分支（常见于内部基因），则为空。HA/NA 片段的分支名遵循 WHO 命名规范。",
    ),
    "nc_coverage": (
        "Fraction of the reference genome covered by the consensus sequence, as calculated by Nextclade.",
        "一致性序列覆盖参考基因组的比例，由 Nextclade 计算。",
    ),
    "nc_total_coverage": (
        "Weighted average: <code>sum(lenAligned) / sum(lenAligned / coverage)</code>, equivalent to reference-length-weighted mean of per-segment coverage.",
        "加权平均值：<code>sum(lenAligned) / sum(lenAligned / coverage)</code>，即按参考长度加权的各片段覆盖度均值。",
    ),
    "no_data": ("No QC data available.", "暂无 QC 数据。"),
    "no_nc": ("No Nextclade data available.", "暂无 Nextclade 数据。"),
}


def t(key: str, lang: str) -> str:
    """Get text by key and language. lang='en' or 'ch'."""
    idx = 0 if lang == "en" else 1
    return TEXT[key][idx]


def badge(level: str, lang: str) -> str:
    cls = {"good": "badge-green", "acceptable": "badge-orange", "warning": "badge-red"}[level]
    thresh = {"good": "&ge;95%", "acceptable": "&ge;80%", "warning": "&lt;80%"}[level]
    label = t(level if level != "warning" else "warning_low", lang)
    return f'<span class="badge {cls}">{thresh}: {label}</span>'


# ---------------------------------------------------------------------------
# Data reading
# ---------------------------------------------------------------------------


def read_tsv(path: Path) -> tuple[list[str], list[dict]]:
    """Read a TSV file, return (columns, rows)."""
    if not path.is_file():
        return [], []
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        columns = list(reader.fieldnames or [])
        rows = list(reader)
    return columns, rows


# ---------------------------------------------------------------------------
# Formatting helpers
# ---------------------------------------------------------------------------


def fmt_num(val: str) -> str:
    try:
        n = int(val)
        return f"{n:,}"
    except (ValueError, TypeError):
        try:
            f = float(val)
            if f == int(f):
                return f"{int(f):,}"
            return f"{f:,.2f}"
        except (ValueError, TypeError):
            return str(val)


def fmt_pct(val: str) -> str:
    try:
        f = float(val)
    except (ValueError, TypeError):
        return f"<td>{val}</td>"
    color = "#2e7d32" if f >= 95 else "#f57f17" if f >= 80 else "#c62828"
    return f'<td style="color:{color};font-weight:bold">{f}%</td>'


def fmt_coverage(val: str) -> str:
    try:
        f = float(val)
    except (ValueError, TypeError):
        return f"<td>{val}</td>"
    color = "#2e7d32" if f >= 0.95 else "#f57f17" if f >= 0.80 else "#c62828"
    pct = f * 100
    return f'<td style="color:{color};font-weight:bold">{pct:.2f}%</td>'


# ---------------------------------------------------------------------------
# HTML generators
# ---------------------------------------------------------------------------

CSS = """
body {
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto,
                 "Helvetica Neue", Arial, sans-serif;
    margin: 0; padding: 20px; background: #f5f5f5; color: #333;
}
h1 { text-align: center; color: #1565c0; border-bottom: 3px solid #1565c0; padding-bottom: 10px; }
h2 { color: #1976d2; border-left: 4px solid #1976d2; padding-left: 10px; margin-top: 30px; }
.section { background: #fff; border-radius: 8px; padding: 20px; margin: 15px 0;
           box-shadow: 0 1px 3px rgba(0,0,0,0.12); }
.table { width: 100%; border-collapse: collapse; font-size: 14px; }
.table th { background: #1976d2; color: #fff; padding: 10px 8px; text-align: left;
            position: sticky; top: 0; }
.table td { padding: 8px; border-bottom: 1px solid #e0e0e0; }
.table tbody tr:hover { background: #e3f2fd; }
.coverage-grid { display: flex; flex-wrap: wrap; gap: 20px; }
.sample-block { flex: 1; min-width: 300px; max-width: 500px; }
.sample-block h3 { margin: 0 0 10px 0; color: #333; font-size: 15px; }
.bar-container { width: 100%; }
.bar-row { display: flex; align-items: center; margin: 3px 0; }
.bar-label { width: 50px; font-size: 12px; font-weight: bold; color: #555; flex-shrink: 0; }
.bar-track { flex: 1; height: 18px; background: #e0e0e0; border-radius: 3px;
             margin: 0 8px; overflow: hidden; }
.bar-fill { height: 100%; border-radius: 3px; transition: width 0.3s ease; }
.bar-value { width: 55px; font-size: 12px; text-align: right; color: #555; flex-shrink: 0; }
.summary { display: flex; gap: 15px; flex-wrap: wrap; margin-bottom: 15px; }
.summary-card { flex: 1; min-width: 150px; background: #e3f2fd; border-radius: 6px;
                padding: 12px; text-align: center; }
.summary-card .num { font-size: 28px; font-weight: bold; color: #1565c0; }
.summary-card .label { font-size: 12px; color: #666; margin-top: 4px; }
.method { background: #fafafa; border-left: 4px solid #90caf9; padding: 12px 16px;
          margin-top: 12px; font-size: 13px; color: #555; line-height: 1.7; border-radius: 4px; }
.method strong { color: #333; }
.method ul { margin: 4px 0; padding-left: 20px; }
.method li { margin: 2px 0; }
.badge { display: inline-block; padding: 1px 6px; border-radius: 3px;
         font-size: 11px; font-weight: bold; margin: 0 2px; }
.badge-green { background: #e8f5e9; color: #2e7d32; }
.badge-orange { background: #fff3e0; color: #f57f17; }
.badge-red { background: #ffebee; color: #c62828; }
"""


def generate_qc_table(columns: list[str], rows: list[dict], lang: str) -> str:
    if not rows:
        return f"<p>{t('no_data', lang)}</p>"
    pct_cols = {"Q20", "Q30", "Clean_Q20", "Clean_Q30", "Duplication_Rate", "Alignment_Rate", "Coverage"}
    num_cols = {"Reads_Num", "Base_Num", "Clean_Reads_Num", "Clean_Base_Num", "Library_Size"}
    header = "<tr>" + "".join(f"<th>{c.replace('_', ' ')}</th>" for c in columns) + "</tr>"
    body = []
    for row in rows:
        cells = []
        for c in columns:
            val = row.get(c, "")
            if c in pct_cols:
                cells.append(fmt_pct(val))
            elif c in num_cols:
                cells.append(f"<td>{fmt_num(val)}</td>")
            else:
                cells.append(f"<td>{val}</td>")
        body.append("<tr>" + "".join(cells) + "</tr>")
    return f'<table class="table"><thead>{header}</thead><tbody>{"".join(body)}</tbody></table>'


def generate_nextclade_table(rows: list[dict], lang: str) -> str:
    if not rows:
        return f"<p>{t('no_nc', lang)}</p>"
    columns = ["sample_name", "segment", "dataset", "clade", "coverage"]
    header = "<tr>" + "".join(f"<th>{c.replace('_', ' ')}</th>" for c in columns) + "</tr>"
    body = []
    for row in rows:
        is_total = row.get("segment") == "TOTAL"
        style = ' style="background-color:#e3f2fd;font-weight:bold"' if is_total else ""
        cells = []
        for c in columns:
            val = row.get(c, "")
            if c == "coverage" and val:
                if is_total:
                    try:
                        f = float(val)
                        cells.append(f'<td style="font-weight:bold">{f*100:.2f}%</td>')
                    except (ValueError, TypeError):
                        cells.append(f"<td>{val}</td>")
                else:
                    cells.append(fmt_coverage(val))
            elif c == "dataset" and val and val != "-":
                short = val.replace("nextstrain/", "")
                cells.append(f'<td title="{val}">{short}</td>')
            else:
                cells.append(f"<td>{val}</td>")
        body.append(f'<tr{style}>{"".join(cells)}</tr>')
    return f'<table class="table"><thead>{header}</thead><tbody>{"".join(body)}</tbody></table>'


def generate_coverage_bars(rows: list[dict], lang: str) -> str:
    sample_segments: dict[str, list[tuple[str, float]]] = {}
    for row in rows:
        seg = row.get("segment", "")
        if seg == "TOTAL":
            continue
        sample = row.get("sample_name", "")
        try:
            cov = float(row.get("coverage", 0))
        except (ValueError, TypeError):
            cov = 0.0
        sample_segments.setdefault(sample, []).append((seg, cov))
    if not sample_segments:
        return ""
    bars = []
    for sample, segments in sample_segments.items():
        seg_bars = []
        for seg, cov in segments:
            pct = cov * 100
            bar_color = "#4caf50" if pct >= 95 else "#ff9800" if pct >= 80 else "#f44336"
            seg_bars.append(
                f'<div class="bar-row">'
                f'<span class="bar-label">{seg}</span>'
                f'<div class="bar-track"><div class="bar-fill" style="width:{pct:.1f}%;background:{bar_color}"></div></div>'
                f'<span class="bar-value">{pct:.2f}%</span></div>'
            )
        bars.append(
            f'<div class="sample-block"><h3>{sample}</h3>'
            f'<div class="bar-container">{"".join(seg_bars)}</div></div>'
        )
    return (
        f'<div class="section">'
        f'<h2>{t("seg_coverage_title", lang)}</h2>'
        f'<div class="coverage-grid">{"".join(bars)}</div></div>'
    )


def build_html(qc_columns: list[str], qc_rows: list[dict],
               nc_rows: list[dict], lang: str) -> str:
    """Build HTML report. lang='en' or 'ch'."""
    L = lang
    coverage_section = generate_coverage_bars(nc_rows, L)

    return f"""<!DOCTYPE html>
<html lang="{'en' if L == 'en' else 'zh-CN'}">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>{t('title', L)}</title>
<style>{CSS}</style>
</head>
<body>
<h1>{t('title', L)}</h1>

<div class="section">
<div class="summary">
  <div class="summary-card">
    <div class="num">{len(qc_rows)}</div>
    <div class="label">{t('samples_processed', L)}</div>
  </div>
  <div class="summary-card">
    <div class="num">{sum(1 for r in nc_rows if r.get('segment') != 'TOTAL')}</div>
    <div class="label">{t('segments_analyzed', L)}</div>
  </div>
</div>
<div class="method">
  <strong>{t('pipeline_overview_title', L)}:</strong> {t('pipeline_overview', L)}
  <ul>
    <li><strong>fastp</strong> &mdash; {t('step_fastp', L)}</li>
    <li><strong>IRMA</strong> &mdash; {t('step_irma', L)}</li>
    <li><strong>minimap2 + samtools</strong> &mdash; {t('step_align', L)}</li>
    <li><strong>Nextclade sort</strong> &mdash; {t('step_sort', L)}</li>
    <li><strong>Nextclade run</strong> &mdash; {t('step_run', L)}</li>
  </ul>
</div>
</div>

<div class="section">
<h2>{t('qc_summary', L)}</h2>
{generate_qc_table(qc_columns, qc_rows, L)}
<div class="method">
  <strong>{t('qc_explain_title', L)}:</strong>
  <ul>
    <li><strong>Reads_Num / Base_Num:</strong> {t('qc_reads_num', L)}</li>
    <li><strong>Q20 / Q30:</strong> {t('qc_q20q30', L)}
        {badge('good', L)} {badge('acceptable', L)} {badge('warning', L)}
    </li>
    <li><strong>Clean_Reads_Num / Clean_Base_Num:</strong> {t('qc_clean', L)}</li>
    <li><strong>Clean_Q20 / Clean_Q30:</strong> {t('qc_clean_q', L)}</li>
    <li><strong>Duplication_Rate:</strong> {t('qc_dup', L)}</li>
    <li><strong>Library_Size:</strong> {t('qc_libsize', L)}</li>
    <li><strong>Alignment_Rate:</strong> {t('qc_align_rate', L)}</li>
    <li><strong>Coverage:</strong> {t('qc_coverage', L)}
        {badge('good', L)} {badge('acceptable', L)} {badge('warning', L)}
    </li>
  </ul>
</div>
</div>

{coverage_section}

<div class="section">
<h2>{t('nc_results', L)}</h2>
{generate_nextclade_table(nc_rows, L)}
<div class="method">
  <strong>{t('nc_explain_title', L)}:</strong>
  <ul>
    <li><strong>segment:</strong> {t('nc_segment', L)}</li>
    <li><strong>dataset:</strong> {t('nc_dataset', L)}</li>
    <li><strong>clade:</strong> {t('nc_clade', L)}</li>
    <li><strong>coverage:</strong> {t('nc_coverage', L)}
        {badge('good', L)} {badge('acceptable', L)} {badge('warning', L)}
    </li>
    <li><strong>TOTAL coverage:</strong> {t('nc_total_coverage', L)}</li>
  </ul>
</div>
</div>

</body>
</html>"""


def main():
    parser = argparse.ArgumentParser(
        description="Generate HTML reports (EN + CN) from IRMA+Nextclade pipeline output"
    )
    parser.add_argument("--outdir", required=True, help="Pipeline output directory")
    args = parser.parse_args()

    outdir = Path(args.outdir).resolve()
    qc_columns, qc_rows = read_tsv(outdir / "QC.tsv")
    _, nc_rows = read_tsv(outdir / "nextclade_result.tsv")

    if not qc_rows and not nc_rows:
        print(f"No data found in {outdir}", file=sys.stderr)
        sys.exit(1)

    for lang, filename in [("en", "report.html"), ("ch", "report_ch.html")]:
        path = outdir / filename
        path.write_text(build_html(qc_columns, qc_rows, nc_rows, lang), encoding="utf-8")
        print(f"Report generated: {path}")


if __name__ == "__main__":
    main()
