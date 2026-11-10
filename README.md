# 🧬 PCHiC-Forge — From reads to robust promoter–enhancer maps
[![snakemake-check](https://github.com/yousra291987/PCHiC_workflow/actions/workflows/snakemake-check.yaml/badge.svg)](https://github.com/yousra291987/PCHiC_workflow/actions/workflows/snakemake-check.yaml)

_End-to-end Snakemake pipeline for PCHi-C and Hi-C: QC, significant interaction detection, .hic generation, insulation scores, and TAD calling._

---

## 📘 Overview

This repository provides a **modular and reproducible Snakemake pipeline** for analyzing **Promoter Capture-HiC (PCHi-C)** and **Hi-C** data.  
It performs all major processing and analysis steps — from read mapping and QC, through significant interaction detection (Chicago / ChiCMaxima), to 3D genome structure generation and insulation score computation.

---

## 🧩 Workflow Summary

| Step | Description | Main Tool |
|------|--------------|------------|
| **1. HiCUP QC** | Maps paired reads, filters invalid di-tags | [`HiCUP`](https://www.bioinformatics.babraham.ac.uk/projects/hicup/) |
| **2. Capture Efficiency** | Computes P0/P1/P2 statistics from `.mat` | R / data.table |
| **3. Chicago** | Detects significant promoter interactions | [`Chicago`](https://bioconductor.org/packages/release/bioc/html/Chicago.html) |
| **4. ChiCMaxima** | Generates interaction BED (IBED) for visualization | [`ChiCMaxima`](https://github.com/yousra291987/ChiCMaxima) |
| **5. Hi-C generation** | Converts matrices to `.hic` for Juicebox / Arrowhead | [`juicer_tools`](https://github.com/aidenlab/juicer/wiki/Juicer-Tools-Quick-Start) |
| **6. Insulation Score** | Calculates domain insulation using FAN-C | [`FAN-C`](https://github.com/vaquerizaslab/fanc) |

Each stage is implemented as a **Snakemake rule** and wrapped in **R scripts** for reproducibility and transparency.

---

## 📁 Repository Structure
```
├── Snakefile
├── config/
│ ├── config.yaml
│ └── envs/
│ ├── r.yaml
│ ├── fanc.yaml
│ ├── chicago.yaml
│ └── juicer.yaml
├── scripts/
│ ├── 00_CaptureHiC_library_qualityControl.R
│ ├── 00_CaptureHiC_CaptureEfficiency_Control.R
│ ├── 01_Prepare_Chicago.R
│ ├── 02_ChiCMaxima.R
│ ├── 03_Generate_hic_and_TADs.R
│ ├── 04_InsulationScore.R
│ ├── BamToFragments.pl
│ ├── Bam2Juicer.pl
│ └── align2ibed.pl
├── design/
│ ├── mm10_HindIII/
│ │ ├── MM10_HindIII.baitmap
│ │ ├── MM10_HindIII.rmap
│ │ └── MM10_HindIII_ProbeInfo.txt
└── README.md
```

---

## ⚙️ Dependencies

The workflow uses **Snakemake** and manages environments via **conda**.

### Core requirements

| Tool | Version (tested) | Purpose |
|------|------------------|----------|
| Snakemake | ≥ 7.8 | Workflow engine |
| R | ≥ 4.2 | Statistical scripting |
| Python | ≥ 3.9 | FAN-C & auxiliary scripts |
| Perl | ≥ 5.30 | Bam2Juicer / BamToFragments |
| Java | ≥ 8 | Juicer tools |

### Conda environments

All environments are provided in `config/envs/`:

- `r.yaml` → core R dependencies (optparse, data.table, tidyverse, etc.)
- `chicago.yaml` → Bioconductor Chicago
- `fanc.yaml` → FAN-C
- `juicer.yaml` → openjdk + juicer_tools

You can create them manually or let Snakemake handle them automatically.

---

##  Installation

1. **Clone the repository**
   ```bash
   git clone https://github.com/yousra291987/PCHiC_workflow.git
   cd PCHiC_workflow
   
2. **Install conda (if not available)**
   ```bash 
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   bash Miniconda3-latest-Linux-x86_64.sh
   
3. **Install Snakemake**
   ```bash 
   conda install -c bioconda -c conda-forge snakemake
 
4. **Set up environments**
 ```bash 
   Set up environments (optional)
```

##  Configuration

All runtime parameters are defined in config/config.yaml.
```yaml
datadir: "results"
basedir: "/path/to/project"
SampleName: "SampleA"

design:
dir: "design/mm10_HindIII"
rmap: "design/mm10_HindIII/MM10_HindIII.rmap"
baitmap: "design/mm10_HindIII/MM10_HindIII.baitmap"
probeinfo: "design/mm10_HindIII/MM10_HindIII_ProbeInfo.txt"

tools:
rscript: "Rscript"
perl: "perl"
python: "python"
java: "java"
juicer_tools: "juicer_tools.jar"
hicup: "/path/to/hicup"
bam_to_fragments: "scripts/BamToFragments.pl"
bam2juicer: "scripts/Bam2Juicer.pl"
align2ibed: "scripts/align2ibed.pl"

params:
juicer_mem_gb: 10
juicer_q: 10
juicer_genome_id: "mm10"
fanc_resolution: "10kb"
fanc_balancing: "KR"
fanc_windows: ["0.5mb","1mb","1.5mb","2mb","2.5mb"]
```
## 🧬 Adjusting the configuration for your data
Before running the workflow, open the file config/config.yaml and update it with your own paths:
Under the samples: section, specify the absolute or relative paths to your paired-end FASTQ files (r1 and r2) for each sample.
Update the hicup: section with the correct paths to your Bowtie2 binary, Bowtie2 index, and genome digest file corresponding to your restriction enzyme and genome assembly.
You can add multiple samples; the pipeline will automatically process each in parallel.
```yaml
Example:
samples:
  MySample01:
    r1: "data/MySample01_R1.fastq.gz"
    r2: "data/MySample01_R2.fastq.gz"
  MySample02:
    r1: "data/MySample02_R1.fastq.gz"
    r2: "data/MySample02_R2.fastq.gz"

hicup:
  bowtie2_bin: "/usr/bin/bowtie2"
  index: "/refs/bowtie2/mm10"
  digest: "design/mm10_HindIII/mm10_HindIII_digest.txt"
```

### 🔁 Technical replicates (multiple FASTQ pairs per sample)

A **technical replicate** is the *same* library sequenced more than once (e.g. across lanes or
sequencing runs). You list all FASTQ pairs for the sample in the config, and the pipeline produces a
single `mapped.bam` so everything downstream (QC, Chicago, `.hic`, insulation) treats the sample as one.

> **How the merge actually happens.** HiCUP itself does **not** combine pairs: per the
> [HiCUP docs](https://www.bioinformatics.babraham.ac.uk/projects/hicup/read_the_docs/html/index.html),
> "a file pair generates one output file", so *N* pairs give *N* `*.hicup.bam` files. The `hicup_qc`
> step then **merges those per-replicate BAMs** into one `mapped.bam` with `samtools merge`.
> Note: HiCUP de-duplicates each pair separately, so duplicates spanning different pairs are not
> removed (normally negligible for technical replicates of one library). If you need de-duplication
> across everything, concatenate the FASTQs into a single pair instead and run HiCUP once on that.

A sample may be declared in any of the following ways:

**1. Single pair (legacy form) — still works:**
```yaml
samples:
  MySample01:
    r1: "data/MySample01_R1.fastq.gz"
    r2: "data/MySample01_R2.fastq.gz"
```

**2. Several technical replicates — preferred form (`fastqs:`):**
```yaml
samples:
  MySample02:
    fastqs:
      - r1: "data/MySample02_lane1_R1.fastq.gz"
        r2: "data/MySample02_lane1_R2.fastq.gz"
      - r1: "data/MySample02_lane2_R1.fastq.gz"
        r2: "data/MySample02_lane2_R2.fastq.gz"
      # add one r1/r2 entry per sequencing run — the pipeline merges them into one mapped.bam
```

**3. Equal-length `r1` / `r2` lists — matched element-wise:**
```yaml
samples:
  MySample03:
    r1: ["data/MySample03_L1_R1.fastq.gz", "data/MySample03_L2_R1.fastq.gz"]
    r2: ["data/MySample03_L1_R2.fastq.gz", "data/MySample03_L2_R2.fastq.gz"]
```

You can mix forms across samples. A few rules the pipeline enforces:

- Every pair needs **both** `r1` and `r2`.
- In form 3, the `r1` and `r2` lists must have the **same length** — each R1 needs a matching R2.
- R1 and R2 are matched **by position**, so keep the entries in the same order.
- Paths may be relative or absolute (they are resolved to absolute paths internally).
- **File names are arbitrary.** Use whatever paths your FASTQs actually have
  (`testis_R1.fastq.gz`, `miti.chicken.testis.R1.fastq.gz`, …); there is no
  required naming convention — the values you write here are passed to HiCUP verbatim.

## 🔧 Downstream tool setup (CHiCAGO & juicer_tools)

QC and capture-efficiency run with the bundled conda environments, but two
downstream steps need external tools that are **not** shipped with this repo.
Configure them in `config/config.yaml` before running those steps.

**CHiCAGO interaction calling.** The `Chicago` R package is now installed from
the conda env (`envs/chicago.yaml`), but CHiCAGO also needs the **chicagoTools**
helper scripts (`bam2chicago.sh`, `makeDesignFiles_py3.py`), which are
distributed separately from the R package. Point the config at them:

```yaml
tools:
  chicago_bam2chicago: "/path/to/chicagoTools/bam2chicago.sh"
  chicago_make_design: "/path/to/chicagoTools/makeDesignFiles_py3.py"
  python: "python"

chicago:
  # true on the FIRST run (or whenever rmap/baitmap change) to build the
  # CHiCAGO design files (.npb/.nbpb/.poe) into design.dir. Needs chicago_make_design.
  force_make_design: true
```

If `chicago_bam2chicago` is left empty, the CHiCAGO step stops early with a clear
message instead of running on missing inputs. The design directory must contain
the CHiCAGO design files (`.npb`, `.nbpb`, `.poe`); set `force_make_design: true`
once to generate them, then you can set it back to `false`.

**`.hic` generation (juicer_tools).** Download `juicer_tools` and set an
**absolute** path (the default `"juicer_tools.jar"` is only a placeholder):

```yaml
tools:
  juicer_tools: "/path/to/juicer_tools.1.22.01.jar"
```

ChiCMaxima and capture-efficiency do **not** need either of these, and (as of the
latest version) ChiCMaxima no longer waits on CHiCAGO, so it runs even if CHiCAGO
is not configured.

### ⏯️ Turning steps on/off

If a tool isn't available to you (e.g. you can't install CHiCAGO, or you don't
have `juicer_tools`), switch the corresponding step **off** in `config.yaml`
instead of leaving tool paths empty — an empty path does not skip the step, it
just makes it fail. Each step defaults to `true`:

```yaml
steps:
  capture_efficiency: true
  chicago: false          # skip CHiCAGO entirely (its rule is never requested)
  chicmaxima: true
  hic_and_tads: true      # set false if you don't have juicer_tools
  insulation: true        # only runs when hic_and_tads is also true
```

A step set to `false` is simply not requested, so Snakemake never enters its
rule. When the workflow starts it prints the resolved state of every step, e.g.

```
[steps] capture_efficiency=on, chicago=OFF, chicmaxima=on, hic_and_tads=OFF, insulation=on
```

Check that line to confirm your toggle took effect. Two common reasons a
`false` is ignored: the key is not indented **under** `steps:`, or it is
**quoted** (`hic_and_tads: "false"` is the string `"false"`, which is truthy —
write it unquoted: `hic_and_tads: false`).

To skip a step for a single run **without editing the config** (this always
works, regardless of the toggles), use `--omit-from <rule>` or `--until <rule>`.
You can repeat `--omit-from`:

```bash
# Run QC + capture efficiency + ChiCMaxima only; skip CHiCAGO and .hic/TADs:
snakemake -c 4 --use-conda --omit-from chicago --omit-from hic_and_tads
```

> **Non-standard genomes & `.hic`:** `juicer_tools pre` needs its last argument to
> be a known genome ID (`hg19`, `mm10`, …) **or a chromosome-sizes file**. For a
> custom assembly (e.g. `ggal7`), set `params.juicer_genome_id` to the path of a
> `chrom.sizes` file for that genome, otherwise the `hic_and_tads` step will size
> chromosomes incorrectly or fail.

## 🚀 Running the Pipeline

1. **Dry-run (check commands)**

```bash
snakemake --dry-run
```
2. **Execute with 4 cores**

```bash
snakemake -c 4 --use-conda
```
3. **Create a detailed HTML report**

```bash
snakemake --report pipeline_report.html
```
4. **Restart after interruption**

```bash
snakemake -c 4 --rerun-incomplete
```

## 📤 Outputs
All results are written under results/<SampleName>/.

| Output                     | Description                                |
| -------------------------- | ------------------------------------------ |
| `*.hicup.bam`              | Filtered BAM from HiCUP                    |
| `*.mat`                    | Contact matrix for Chicago / QC            |
| `*_CaptureEfficiency.stat` | Capture efficiency statistics              |
| `*_Chicago_output/`        | Chicago outputs + significant interactions |
| `*_Cis.ibed`               | ChiCMaxima IBED file                       |
| `*.hic`                    | Final Juicebox-compatible contact map      |
| `*.insulation`             | Domain insulation profile (FAN-C)          |


## 🧠 Notes
You may adapt the design files for any restriction enzyme and genome build.
The pipeline can run on Linux, macOS, or HPC clusters (Snakemake supports SLURM/LSF/etc.).
Ensure the paths to HiCUP, juicer_tools, and auxiliary Perl scripts are correct in your config.

## 👩‍💻 Author
Yousra Ben Zouari

## 💬 Contact & Support
For questions, feature requests, or bug reports, please open an issue on GitHub or contact
📧 yousra.benzouari@gmail.com
