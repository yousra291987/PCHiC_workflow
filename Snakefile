# capture-hic-pipeline/Snakefile
# -------------------------------------------------------------------
# Capture-HiC workflow: QC, interaction calling (Chicago/ChiCMaxima),
# .hic generation and insulation (FAN-C)
# Portable (conda), multi-sample, species-agnostic via config.yaml
# -------------------------------------------------------------------

import os
from pathlib import Path

# Default for users; CI overrides with: --configfile config/config-ci.yaml
configfile: "config/config.yaml"

# Accept samples as dict {sample: {...}} or list [sample, ...]
def _resolve_samples(cfg):
    s = cfg.get("samples", {})
    if isinstance(s, dict):
        return list(s.keys())
    if isinstance(s, list):
        return s
    return []

# Optional: per-sample metadata (only present if samples is a dict)
def _sample_meta(sample_id):
    _s = config.get("samples", {})
    if isinstance(_s, dict):
        return _s.get(sample_id, {})
    return {}

SAMPLES = _resolve_samples(config)
OUTDIR  = config.get("outdir", "results")
THREADS = config.get("threads", {}).get("default", 4)

DESIGN = config.get("design", {})        # rmap, baitmap, probeinfo, dir
PARAMS = config.get("params", {})
TOOLS  = config.get("tools", {"rscript": "Rscript"})  # tool paths with safe default

# -------------------------------------------------------------------
# Default targets (precomputed; no wildcards to infer)
# -------------------------------------------------------------------
TARGETS = []
for s in SAMPLES:
    TARGETS += [
        f"{OUTDIR}/{s}/hicup_{s}.config",
        f"{OUTDIR}/{s}/{s}.hic",
        f"{OUTDIR}/{s}/{s}.insulation",
        f"logs/hicup/{s}.Rout",
        f"logs/capture_eff/{s}.Rout",
        f"logs/chicago/{s}.Rout",
        f"logs/chicmaxima/{s}.Rout",
    ]

rule all:
    input: TARGETS

# -------------------------------------------------------------------
# 1) HICUP config File
# -------------------------------------------------------------------
rule hicup_config:
    """
    Render per-sample HiCUP config from `config` + `_sample_meta`.
    Produces: {outdir}/{sample}/hicup_{sample}.config
    """
    output:
        cfg = OUTDIR + "/{sample}/hicup_{sample}.config"
    params:
        outdir   = OUTDIR,
        # HiCUP global settings (with sane defaults)
        quiet    = lambda w: config.get('hicup', {}).get('quiet', 1),
        keep     = lambda w: config.get('hicup', {}).get('keep', 0),
        zip      = lambda w: config.get('hicup', {}).get('zip', 0),
        bowtie2  = lambda w: config.get('hicup', {}).get('bowtie2_bin', ''),
        index    = lambda w: config.get('hicup', {}).get('index', ''),     # bowtie2 index prefix (no .bt2)
        digest   = lambda w: config.get('hicup', {}).get('digest', ''),    # hicup_digester output
        fmt      = lambda w: config.get('hicup', {}).get('format', ''),    # "", "Sanger", etc.
        longest  = lambda w: config.get('hicup', {}).get('longest', 2000),
        shortest = lambda w: config.get('hicup', {}).get('shortest', 100),
        # Per-sample FASTQs from your helper (fallbacks are CI-friendly)
        r1       = lambda w: os.path.abspath(_sample_meta(w.sample).get('r1', 'R1.fastq.gz')),
        r2       = lambda w: os.path.abspath(_sample_meta(w.sample).get('r2', 'R2.fastq.gz')),
    threads: 1
    message:
        "Writing HiCUP config for {wildcards.sample} -> {output.cfg}"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir}/{wildcards.sample}

        # Quick guardrails: ensure r1/r2 look non-empty
        if [ -z "{params.r1}" ] || [ -z "{params.r2}" ]; then
          echo "[hicup_config] r1/r2 missing for {wildcards.sample}. Check config['samples']." >&2
          exit 2
        fi

        cat > {output.cfg} <<CFG
# Auto-generated HiCUP config for {wildcards.sample}
Outdir:{params.outdir}/{wildcards.sample}
Threads:{threads}
Quiet:{params.quiet}
Keep:{params.keep}
Zip:{params.zip}
Bowtie2:{params.bowtie2}
Index:{params.index}
Digest:{params.digest}
Format:{params.fmt}
Longest:{params.longest}
Shortest:{params.shortest}
{params.r1}
{params.r2}
CFG
"""

# -------------------------------------------------------------------
# 1) HICUP QC (pure shell so conda is allowed)
# -------------------------------------------------------------------
rule hicup_qc:
    input:
        script    = "scripts/00_CaptureHiC_library_qualityControl.R",
        probeinfo = DESIGN.get("probeinfo", ""),
        rmap      = DESIGN.get("rmap", ""),
        baitmap   = DESIGN.get("baitmap", ""),
        hicup_cfg = OUTDIR + "/{sample}/hicup_{sample}.config",   # from hicup_config
    output:
        rout      = "logs/hicup/{sample}.Rout",
        bam       = OUTDIR + "/{sample}/mapped.bam",
        bed       = temp(OUTDIR + "/{sample}/{sample}.hicup.bed"),
        mat       = OUTDIR + "/{sample}/{sample}.mat",
    params:
        outdir     = OUTDIR,
        hicup_bin  = TOOLS.get("hicup", "hicup"),
        bedtools   = TOOLS.get("bedtools", "bedtools"),
        perl       = TOOLS.get("perl", "perl"),
        bam2frags  = TOOLS.get("bam_to_fragments", "scripts/BamToFragments.pl"),
        rscript    = TOOLS.get("rscript", "Rscript"),
    threads: config.get("hicup", {}).get("threads", THREADS)
    conda:
        "envs/r.yaml"
    log:
        "logs/hicup/{sample}.log",
    shell:
        r"""
        set -euo pipefail
        mkdir -p logs/hicup {params.outdir}/{wildcards.sample}

        {params.rscript} "{input.script}" \
          --output "{params.outdir}/{wildcards.sample}" \
          --sample "{wildcards.sample}" \
          --probeinfo "{input.probeinfo}" \
          --rmap "{input.rmap}" \
          --baitmap "{input.baitmap}" \
          --hicup_config "{input.hicup_cfg}" \
          --hicup_bin "{params.hicup_bin}" \
          --bedtools_bin "{params.bedtools}" \
          --perl_bin "{params.perl}" \
          --bam_to_fragments "{params.bam2frags}" \
          --mapped_bam "{output.bam}" \
          --bed_out "{output.bed}" \
          --mat_out "{output.mat}" \
          > {log} 2>&1

        echo "OK" > {output.rout}
        """

# -------------------------------------------------------------------
# 2) Capture efficiency control
# -------------------------------------------------------------------
rule capture_efficiency:
    input:
        script = "scripts/00_CaptureHiC_CaptureEfficiency_Control.R",
        mat    = OUTDIR + "/{sample}/{sample}.mat",
        prev   = "logs/hicup/{sample}.Rout",
    output:
        rout  = "logs/capture_eff/{sample}.Rout",
        stats = OUTDIR + "/{sample}/{sample}_CaptureEfficiency.stat",
    params:
        outdir  = OUTDIR,
        rscript = TOOLS.get("rscript","Rscript"),
    threads: THREADS
    conda:
        "envs/r.yaml"
    log:
        "logs/capture_eff/{sample}.log",
    shell:
        r"""
        mkdir -p logs/capture_eff
        {params.rscript} "{input.script}" \
          --dir "{params.outdir}/{wildcards.sample}" \
          --sample "{wildcards.sample}" \
          --mat "{input.mat}" \
          --stats "{output.stats}" \
          > {log} 2>&1
        echo "OK" > {output.rout}
        """

# -------------------------------------------------------------------
# 3) CHICAGO interaction calling
# -------------------------------------------------------------------
rule chicago:
    input:
        script  = "scripts/01_Prepare_Chicago.R",
        rmap    = DESIGN.get("rmap", ""),
        baitmap = DESIGN.get("baitmap", ""),
        bam     = OUTDIR + "/{sample}/mapped.bam",
        prev    = "logs/capture_eff/{sample}.Rout",
    params:
        outdir            = OUTDIR,
        design_dir        = DESIGN.get("dir", ""),
        bam2chicago       = TOOLS.get("chicago_bam2chicago", ""),
        make_design       = TOOLS.get("chicago_make_design", ""),
        python_bin        = TOOLS.get("python", "python"),
        force_make_design = config.get("chicago", {}).get("force_make_design", False),
        force_flag        = "TRUE" if config.get("chicago", {}).get("force_make_design", False) else "FALSE",
        rscript           = TOOLS.get("rscript","Rscript"),
    output:
        rout = "logs/chicago/{sample}.Rout",
        out  = OUTDIR + "/{sample}/chicago_interactions.tsv",
    threads: THREADS
    conda:
        "envs/chicago.yaml"
    log:
        "logs/chicago/{sample}.log",
    shell:
        r"""
        set -euo pipefail
        mkdir -p logs/chicago

        # Ensure BiocManager and CHiCAGO are present (works on macOS ARM too)
        {params.rscript} -e 'if (!requireNamespace("BiocManager", quietly=TRUE)) {{
                               install.packages("BiocManager", repos="https://cloud.r-project.org")
                             }}
                             if (!requireNamespace("Chicago", quietly=TRUE)) {{
                               BiocManager::install("Chicago", ask=FALSE, update=FALSE)
                             }}' >> {log} 2>&1

        # Run your pipeline step
        {params.rscript} "{input.script}" \
          --dir "{params.outdir}/{wildcards.sample}" \
          --sample "{wildcards.sample}" \
          --rmap "{input.rmap}" \
          --baitmap "{input.baitmap}" \
          --bam "{input.bam}" \
          --design_dir "{params.design_dir}" \
          --out "{output.out}" \
          --bam2chicago "{params.bam2chicago}" \
          --make_design "{params.make_design}" \
          --python "{params.python_bin}" \
          --force_make_design {params.force_flag} \
          >> {log} 2>&1

        echo "OK" > {output.rout}
        """
# -------------------------------------------------------------------
# 4) ChiCMaxima interaction calling
# -------------------------------------------------------------------
rule chicmaxima:
    input:
        script  = "scripts/02_ChiCMaxima.R",
        rmap    = DESIGN.get("rmap", ""),
        baitmap = DESIGN.get("baitmap", ""),
        prev    = "logs/chicago/{sample}.Rout",
    params:
        outdir     = OUTDIR,
        rscript    = TOOLS.get("rscript","Rscript"),
        perl       = TOOLS.get("perl","perl"),
        align2ibed = TOOLS.get("align2ibed","scripts/align2ibed.pl"),
    output:
        rout = "logs/chicmaxima/{sample}.Rout",
        ibed = OUTDIR + "/{sample}/{sample}_Cis.ibed",
    threads: THREADS
    conda:
        "envs/r.yaml"
    log:
        "logs/chicmaxima/{sample}.log",
    shell:
        r"""
        mkdir -p logs/chicmaxima
        {params.rscript} "{input.script}" \
          --dir "{params.outdir}/{wildcards.sample}" \
          --sample "{wildcards.sample}" \
          --rmap "{input.rmap}" \
          --baitmap "{input.baitmap}" \
          --perl_bin "{params.perl}" \
          --align2ibed "{params.align2ibed}" \
          --ibed_out "{output.ibed}" \
          > {log} 2>&1
        echo "OK" > {output.rout}
        """

# -------------------------------------------------------------------
# 5) Generate .hic and optional TADs; then FAN-C insulation
# -------------------------------------------------------------------
rule hic_and_tads:
    input:
        script = "scripts/03_Generate_hic_and_TADs.R",
        rmap   = DESIGN.get("rmap", ""),
        mat    = OUTDIR + "/{sample}/{sample}.mat",
        bam    = OUTDIR + "/{sample}/mapped.bam",
    params:
        outdir     = OUTDIR,
        java_bin   = TOOLS.get("java", "java"),
        juicer_jar = TOOLS.get("juicer_tools", ""),
        perl_bin   = TOOLS.get("perl", "perl"),
        bam2juicer = TOOLS.get("bam2juicer", ""),
        genome_id  = PARAMS.get("juicer_genome_id", ""),
        mem_gb     = PARAMS.get("juicer_mem_gb", 8),
        q_threshold= PARAMS.get("juicer_q", 1),
        rscript    = TOOLS.get("rscript","Rscript"),
    output:
        rout = "logs/hic/{sample}.Rout",
        hic  = OUTDIR + "/{sample}/{sample}.hic",
    threads: THREADS
    conda:
        "envs/r.yaml"
    log:
        "logs/hic/{sample}.log",
    shell:
        r"""
        mkdir -p logs/hic {params.outdir}/{wildcards.sample}
        {params.rscript} "{input.script}" \
          --dir "{params.outdir}/{wildcards.sample}" \
          --sample "{wildcards.sample}" \
          --mat "{input.mat}" \
          --rmap "{input.rmap}" \
          --pre "{params.outdir}/{wildcards.sample}/{wildcards.sample}_pre.txt" \
          --hic "{output.hic}" \
          --java "{params.java_bin}" \
          --juicer_tools "{params.juicer_jar}" \
          --perl "{params.perl_bin}" \
          --bam2juicer "{params.bam2juicer}" \
          --genome_id "{params.genome_id}" \
          --mem_gb {params.mem_gb} \
          --q {params.q_threshold} \
          > {log} 2>&1
        echo "OK" > {output.rout}
        """

rule insulation:
    input:
        script = "scripts/04_InsulationScore.R",
        hic    = OUTDIR + "/{sample}/{sample}.hic",
        prev   = "logs/hic/{sample}.Rout",
    params:
        outdir      = OUTDIR,
        resolution  = PARAMS.get("fanc_resolution", "25kb"),
        balancing   = PARAMS.get("fanc_balancing", "ICE"),
        windows     = PARAMS.get("fanc_windows", []),
        windows_str = " ".join(PARAMS.get("fanc_windows", [])) if isinstance(PARAMS.get("fanc_windows", []), (list, tuple)) else str(PARAMS.get("fanc_windows", "")),
        rscript     = TOOLS.get("rscript","Rscript"),
    output:
        rout = "logs/insulation/{sample}.Rout",
        ins  = OUTDIR + "/{sample}/{sample}.insulation",
    threads: THREADS
    conda:
        "envs/fanc.yaml"
    log:
        "logs/insulation/{sample}.log",
    shell:
        r"""
        mkdir -p logs/insulation
        {params.rscript} "{input.script}" \
          --dir "{params.outdir}/{wildcards.sample}" \
          --sample "{wildcards.sample}" \
          --hic "{input.hic}" \
          --out "{output.ins}" \
          --resolution "{params.resolution}" \
          --balancing "{params.balancing}" \
          --windows "{params.windows_str}" \
          > {log} 2>&1
        echo "OK" > {output.rout}
        """
