# IRMA + Nextclade 自动化流程设计文档

## 1. 概述

对病毒 NGS 数据（流感、RSV、CoV 等）进行自动化分析：

```
fastp QC → IRMA 组装 → FASTA 重命名 → 覆盖度图收集 → 比对率计算 → Nextclade sort → Nextclade run → HTML 报告生成
```

## 2. CLI 接口

### 主流程

```bash
python irma_nextclade.py \
    --fq_list fq.list \
    --database database_path \
    --outdir output_dir
```

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--fq_list` | 输入样本列表 TSV | 必填 |
| `--database` | nextclade dataset 存储路径 | 脚本目录下 `nextclade_data` |
| `--outdir` | 输出目录 | 必填 |

### 报告生成

```bash
python irma_report.py --outdir output_dir
```

| 参数 | 说明 |
|------|------|
| `--outdir` | 流程输出目录（同 irma_nextclade.py 的 outdir） |

输出 `report.html`（英文）和 `report_ch.html`（中文）。

## 3. 输入格式

TSV 文件，含表头，列名：

```
sample	MODULE	R1	R2
flu_b1	FLU	/path/to/flu_b1_R1.fq.gz	/path/to/flu_b1_R2.fq.gz
h3n2_b2	FLU	/path/to/h3n2_b2_R1.fq.gz	/path/to/h3n2_b2_R2.fq.gz
rsv_b4	RSV	/path/to/rsv_b4.R1.fq.gz	/path/to/rsv_b4.R2.fq.gz
cov_b5	CoV	/path/to/cov_b5_R1.fq.gz	/path/to/cov_b5_R2.fq.gz
```

- `MODULE`: IRMA 模块名，目前支持 `FLU`、`RSV`、`CoV`，可选 `EBOLA`、`FLU_AD`（预留）

## 4. 输出文件

### 4.1 QC.tsv — 所有样本 QC 汇总

```
sample_name	Reads_Num	Base_Num	Q20	Q30	Clean_Reads_Num	Clean_Base_Num	Clean_Q20	Clean_Q30	Duplication_Rate	Library_Size	Alignment_Rate	Coverage
flu_b1	240798	36360498	95.55	94.35	237340	31733696	97.17	96.45	20.42	138	20.36	100.0
```

| 字段 | 说明 |
|------|------|
| Reads_Num / Base_Num | 原始数据 R1+R2 合计 read 数和碱基数 |
| Q20 / Q30 | 原始数据 Phred &ge;20/&ge;30 碱基百分比 |
| Clean_Reads_Num / Clean_Base_Num | fastp 过滤后的 read 数和碱基数 |
| Clean_Q20 / Clean_Q30 | 过滤后数据 Q20/Q30 |
| Duplication_Rate | PCR 重复率（%） |
| Library_Size | fastp 估算的 insert size peak |
| Alignment_Rate | clean reads 比对到 Nextclade 数据集标准参考序列的比率（%），来自 `samtools flagstat` |
| Coverage | 基因组覆盖度（%），来自 `samtools coverage -Q 20`，公式 `sum(covbases) / sum(endpos)` |

### 4.2 fasta/{sample}.fasta — 每样本重命名后的 consensus FASTA

- 每样本一个文件，存放于 `outdir/fasta/` 目录
- seqID 格式：`{sample_name}_{SegmentName}`
- 示例：`flu_b1_PB2`、`flu_b1_HA`、`rsv_b4_RSV`

### 4.3 nextclade_result.tsv — Nextclade 分型结果汇总

```
sample_name	segment	dataset	clade	coverage
flu_b1	PB2	nextstrain/flu/b/pb2		0.981
flu_b1	PB1	nextstrain/flu/b/pb1		0.950
flu_b1	HA	nextstrain/flu/b/ha/KX058884	C.5.7	0.933
...
flu_b1	TOTAL	-	-	0.965
```

- 每个样本每个 segment 一行，最后一行为该样本的 TOTAL coverage
- `coverage` 为小数（0~1）
- `clade` 可能为空（取决于 dataset 是否定义了 clade，常见于内部基因）
- `total_coverage = sum(lenAligned_i) / sum(lenAligned_i / coverage_i)`（参考长度加权平均）

### 4.4 report.html / report_ch.html — HTML 报告

- 英文版 `report.html`，中文版 `report_ch.html`
- 包含：流程概述、QC 汇总表、各片段覆盖度柱状图、Nextclade 分型表、指标说明
- 颜色编码：&ge;95% 绿色，&ge;80% 橙色，&lt;80% 红色

## 5. 输出目录结构

```
outdir/
├── irma_config.cfg                      # 共享 IRMA 配置（TMP 重定向）
├── TMPDIR/                              # IRMA 临时目录
├── irma_nextclade.log                   # 运行日志
├── QC.tsv                               # 所有样本 QC 汇总
├── nextclade_result.tsv                 # 所有样本 nextclade 汇总
├── report.html                          # 英文 HTML 报告
├── report_ch.html                       # 中文 HTML 报告
├── fasta/                               # 每样本 renamed consensus
│   ├── flu_b1.fasta
│   └── h3n2_b2.fasta
├── coverage_plot/                       # 覆盖度图（sample 前缀）
│   ├── flu_b1_A_HA_H3-coverageDiagram.pdf
│   ├── flu_b1_A_MP-coverageDiagram.pdf
│   └── ...
├── run_flu_b1/                          # 每样本工作目录
│   ├── fastp/                           # fastp 输出（独立子目录）
│   │   ├── C1.fq.gz                     # clean R1
│   │   ├── C2.fq.gz                     # clean R2
│   │   ├── fastp.json
│   │   └── fastp.html
│   ├── flu_b1/                          # IRMA 纯输出目录
│   │   ├── amended_consensus/           # 编号命名的 consensus
│   │   │   ├── flu_b1_1.fa
│   │   │   ├── flu_b1_2.fa
│   │   │   └── ...
│   │   ├── figures/                     # 覆盖度图 PDF
│   │   ├── B_HA.fasta                   # IRMA 按 segment 名的 FASTA
│   │   └── ...                          # IRMA 其他输出（bam/vcf/logs 等）
│   ├── alignment/                       # 比对率计算中间文件
│   │   ├── reference.fasta              # Nextclade dataset 参考序列合并
│   │   ├── reference.mmi               # minimap2 索引
│   │   ├── aligned.bam                  # 比对结果
│   │   ├── aligned.bam.bai
│   │   ├── flagstat.json               # samtools flagstat JSON
│   │   └── coverage.tsv                # samtools coverage TSV
│   ├── sort.tsv                         # nextclade sort 结果
│   ├── nextclade.tsv                    # 该样本 nextclade 汇总（用于合并）
│   └── nextclade/                       # nextclade run 结果（按 dataset）
│       ├── nextstrain_flu_b_ha_KX058884_input.fasta
│       ├── nextstrain_flu_b_ha_KX058884.tsv
│       ├── nextstrain_flu_b_ha_KX058884.json
│       └── ...
└── run_h3n2_b2/
    └── ...
```

**关键目录设计说明：**

- `run_{sample}/` = 每样本工作目录（含 fastp、sort、alignment 等中间文件）
- `run_{sample}/{sample}/` = IRMA 纯输出（IRMA 要求目标目录不存在，否则追加 `-V2` 后缀）
- `run_{sample}/fastp/` = fastp 输出独立子目录，避免污染 IRMA 目标目录

## 6. 处理流程

### Step 1: fastp QC

```bash
fastp -i {R1} -I {R2} \
      -o {run_dir}/fastp/C1.fq.gz -O {run_dir}/fastp/C2.fq.gz \
      -j {run_dir}/fastp/fastp.json -h {run_dir}/fastp/fastp.html \
      --detect_adapter_for_pe
```

- 输出到 `run_{sample}/fastp/` 独立子目录（不直接放在 IRMA 目标目录下）
- `--detect_adapter_for_pe` 自动检测双端接头

**解析 fastp.json 提取字段：**

| 目标字段 | JSON 路径 |
|----------|----------|
| Reads_Num | `summary.before_filtering.total_reads` |
| Base_Num | `summary.before_filtering.total_bases` |
| Q20 | `summary.before_filtering.q20_rate × 100` |
| Q30 | `summary.before_filtering.q30_rate × 100` |
| Clean_Reads_Num | `summary.after_filtering.total_reads` |
| Clean_Base_Num | `summary.after_filtering.total_bases` |
| Clean_Q20 | `summary.after_filtering.q20_rate × 100` |
| Clean_Q30 | `summary.after_filtering.q30_rate × 100` |
| Duplication_Rate | `duplication.rate × 100` |
| Library_Size | `insert_size.peak` |

### Step 2: IRMA 组装

```bash
IRMA {MODULE} {run_dir}/fastp/C1.fq.gz {run_dir}/fastp/C2.fq.gz {irma_dir} -c {config_path}
```

- `irma_dir` = `run_{sample}/{sample}/`（IRMA 目标目录）
- 共享配置文件 `outdir/irma_config.cfg` 内容为 `TMP={outdir}/TMPDIR`，将 IRMA 临时文件重定向
- IRMA 要求 `irma_dir` 不存在，否则自动追加 `-V2` 后缀；流程会先清除已有目录
- `amended_consensus/` 下为编号命名的 consensus FASTA

### Step 3: FASTA 重命名与保存

1. 读取 `{irma_dir}/amended_consensus/` 下的 `.fa` 文件
2. **多片段模块（FLU）**：根据文件名末尾数字编号映射为 segment name
3. **单片段模块（RSV/CoV）**：选取 `{sample}.pad.fa`
4. 重命名 seqID 为 `{sample}_{SegmentName}`（单片段为 `{sample}_{MODULE}`）
5. 写入 `outdir/fasta/{sample}.fasta`

**FLU segment 编号→名称映射（实测确认）：**

```python
FLU_SEGMENT_MAP = {
    1: "PB1", 2: "PB2", 3: "PA", 4: "HA",
    5: "NP",  6: "NA",  7: "MP", 8: "NS"
}
```

### Step 4: 收集覆盖度图

从 `{irma_dir}/figures/*coverageDiagram.pdf` 复制到 `outdir/coverage_plot/`，添加样本名前缀避免跨样本文件名冲突：

```
{irma_dir}/figures/A_HA_H3-coverageDiagram.pdf
→ outdir/coverage_plot/flu_b1_A_HA_H3-coverageDiagram.pdf
```

### Step 5: 比对率计算（minimap2 + samtools）

1. 收集 nextclade run 使用的各 dataset 中的 `reference.fasta`，合并为一个参考 FASTA
2. minimap2 建立 MMI 索引
3. minimap2 比对 clean reads → samtools sort → BAM
4. `samtools flagstat --output-fmt json` → 提取 `QC-passed reads` → `mapped %` 作为 Alignment_Rate
5. `samtools coverage -Q 20` → 计算 Coverage = `sum(covbases) / sum(endpos)`

**注意：**
- `samtools flagstat --output-fmt json` 输出到 stdout，需重定向 `> flagstat.json`
- flagstat JSON 结构：`{"QC-passed reads": {"mapped %": 35.15, ...}}`

### Step 6: Nextclade sort

```bash
nextclade sort -r {run_dir}/sort.tsv {outdir}/fasta/{sample}.fasta
```

- 输入为 Step 3 产出的每样本 renamed FASTA
- 输出 TSV 格式：`index\tseqName\tdataset\tscore\tnumHits`
- 从中提取每条序列对应的 dataset 名

### Step 7: Nextclade run

**Dataset 下载逻辑：**

```
sort 返回 dataset 名: "nextstrain/flu/b/pb2"
本地路径: {database_path}/nextstrain/flu/b/pb2/
```

检查 `database_path/{dataset_name}/reference.fasta` 是否存在：
- 存在 → 直接使用
- 不存在 → 下载：

```bash
nextclade dataset get --name {dataset_name} --output-dir {database_path}/{dataset_name}
```

**按 dataset 分组运行：**

同一样本内，按 dataset 分组的序列合并为临时 FASTA，一次运行：

```bash
nextclade run -D {database_path}/{dataset_name} \
    -t {run_dir}/nextclade/{safe_dataset_name}.tsv \
    -J {run_dir}/nextclade/{safe_dataset_name}.json \
    {temp_input.fasta}
```

其中 `safe_dataset_name` 将 `/` 替换为 `_`（如 `nextstrain_flu_b_ha_KX058884`）。

### Step 8: 解析 nextclade JSON 汇总结果

从 JSON 中提取：

```python
for result in json_data["results"]:
    seq_name = result["seqName"]          # flu_b1_PB2
    clade = result.get("clade", "")       # C.5.7 或空字符串
    coverage = result["coverage"]         # 0.9326 (0~1)
    len_aligned = result["lenAligned"]    # 1885
```

**Total coverage 计算：**

```
total_coverage = sum(lenAligned_i) / sum(lenAligned_i / coverage_i)
```

等价于按参考长度加权的 coverage 平均值。

## 7. Dataset 路径设计

`database_path` 下的目录结构与 nextclade sort 返回的 dataset 名完全一致：

```
nextclade_data/
└── nextstrain/
    └── flu/
        ├── b/
        │   ├── pb2/
        │   │   ├── reference.fasta
        │   │   ├── tree.json
        │   │   ├── sequences.fasta
        │   │   ├── genome_annotation.gff3
        │   │   └── pathogen.json
        │   ├── pb1/
        │   ├── pa/
        │   ├── ha/
        │   │   └── KX058884/
        │   │       └── ... (同上)
        │   ├── np/
        │   ├── na/
        │   │   └── CY073894/
        │   ├── mp/
        │   └── ns/
        ├── h1n1pdm/
        │   └── ...
        └── h3n2/
            └── ...
```

- 按需下载，首次使用某种 dataset 时自动下载
- 后续运行复用已有 dataset

## 8. 报告生成（irma_report.py）

`irma_report.py` 读取 `outdir/QC.tsv` 和 `outdir/nextclade_result.tsv`，生成双语 HTML 报告。

- 通过 `lang` 参数控制语言（`"en"` / `"ch"`）
- 所有文案定义在 `TEXT` 字典中，key → `(en_text, ch_text)`
- 报告包含：
  - 流程概述（5 步说明）
  - QC 汇总表（带颜色编码）
  - 各片段覆盖度柱状图（CSS 横向条形图）
  - Nextclade 分型结果表（TOTAL 行高亮）
  - QC 指标说明（含阈值标签：&ge;95% Good / &ge;80% Acceptable / &lt;80% Warning）

## 9. 模块支持

```python
MODULE_SEGMENT_MAP = {
    "FLU": {
        1: "PB1", 2: "PB2", 3: "PA", 4: "HA",
        5: "NP",  6: "NA",  7: "MP", 8: "NS"
    },
    "RSV": {},   # 单片段，使用 {sample}.pad.fa，seqID = {sample}_RSV
    "CoV": {},   # 单片段，使用 {sample}.pad.fa，seqID = {sample}_CoV
    "EBOLA": {}, # 预留
    "FLU_AD": {},# 预留
}

SUPPORTED_MODULES = {"FLU", "RSV", "CoV"}
```

**多片段 vs 单片段逻辑：**
- `seg_map` 非空（FLU）：遍历 `amended_consensus/*.fa`，按编号映射 segment name
- `seg_map` 为空（RSV/CoV）：只取 `{sample}.pad.fa`

## 10. 错误处理

| 场景 | 处理 |
|------|------|
| fastp 失败 | 记录错误，跳过该样本 |
| IRMA 失败 | 记录错误，跳过该样本 |
| 无 consensus 序列 | 跳过 nextclade，保留已有 QC 数据 |
| nextclade sort 失败 | 记录警告，跳过 nextclade run |
| dataset 下载失败 | 记录错误，跳过相关序列 |
| nextclade run 失败（单个 dataset） | 记录错误，继续其他 dataset |
| 部分 segment 的 clade 为空 | 输出空字符串 |
| 不支持的 MODULE | 报错并跳过该样本 |

## 11. 依赖

- Python >= 3.10
- BioPython
- fastp
- IRMA
- minimap2
- samtools
- nextclade >= 3.21
