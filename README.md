# 🧬 Capture-HiC Analysis Pipeline
[![snakemake-check](https://github.com/yousra291987/PCHiC_workflow/actions/workflows/snakemake-check.yaml/badge.svg)](https://github.com/yousra291987/PCHiC_workflow/actions/workflows/snakemake-check.yaml)

_End-to-end Snakemake pipeline for PCHi-C and Hi-C: QC, significant interaction detection, .hic generation, insulation scores, and TAD calling._

---

## 📘 Overview

This repository provides a **modular and reproducible Snakemake pipeline** for analyzing **Promoter Capture-HiC (PCHi-C)** and **Hi-C** data.  
It performs all major processing and analysis steps — from read mapping and QC, through significant interaction detection (Chicago / ChiCMaxima), to 3D genome structure generation and insulation score computation.

The pipeline was developed and maintained by **Yousra Ben Zouari**  
during her work on chromatin architecture and enhancer regulation.

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
