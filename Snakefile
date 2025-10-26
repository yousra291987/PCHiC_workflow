# capture-hic-pipeline/Snakefile
# -------------------------------------------------------------------
# Capture-HiC workflow: QC, interaction calling (Chicago/ChiCMaxima),
# .hic generation and insulation (FAN-C)
# Portable (conda), multi-sample, species-agnostic via config.yaml
# -------------------------------------------------------------------

import os
from pathlib import Path

configfile: "config/config.yaml"

SAMPLES = config["samples"]
OUTDIR  = config.get("outdir", "results")
THREADS = config.get("threads", {}).get("default", 4)

DESIGN = config["design"]  # rmap, baitmap, probeinfo
PARAMS = config.get("params", {})
TOOLS  = config.get("tools", {"rscript": "Rscript"})

# -------------------------------------------------------------------
# Default targets
# -------------------------------------------------------------------
rule all:
    input:
        expand(f"{OUTDIR}/{{sample}}/{{sample}}.hic", sample=SAMPLES),
        expand(f"{OUTDIR}/{{sample}}/{{sample}}.insulation", sample=SAMPLES),
        expand("logs/hicup/{sample}.Rout", sample=SAMPLES),
        expand("logs/capture_eff/{sample}.Rout", sample=SAMPLES),
        expand("logs/chicago/{sample}.Rout", sample=SAMPLES),
        expand("logs/chicmaxima/{sample}.Rout", sample=SAMPLES)

# -------------------------------------------------------------------
# Utilities
# -------------------------------------------------------------------
def ensure_dirs(*paths):
    for p in paths:
        Path(p).parent.mkdir(parents=True, exist_ok=True)

# -------------------------------------------------------------------
# 1) HICUP QC (placeholder wrapper around your R QC + mapping step)
#    Expects your script to write a BAM (or symlink) per sample.
# -------------------------------------------------------------------

rule hicup_qc:
    input:
        script     = "scripts/00_CaptureHiC_library_qualityControl.R",
        probeinfo  = DESIGN["probeinfo"],
        rmap       = DESIGN["rmap"],
        baitmap    = DESIGN["baitmap"],
        hicup_cfg  = config["hicup"]["config"],
        digest     = config["hicup"]["digest"]
    output:
        rout = "logs/hicup/{sample}.Rout",
        bam  = f"{OUTDIR}/{{sample}}/mapped.bam",
        bed  = temp(f"{OUTDIR}/{{sample}}/{{sample}}.hicup.bed"),
        mat  = f"{OUTDIR}/{{sample}}/{{sample}}.mat"
    threads: THREADS
    conda: "envs/r.yaml"
    log: "logs/hicup/{sample}.log"
    shell:
        r"""
        mkdir -p logs/hicup {OUTDIR}/{wildcards.sample}
        {TOOLS[rscript]} "{input.script}" \
          --output "{OUTDIR}/{wildcards.sample}" \
          --sample "{wildcards.sample}" \
          --probeinfo "{input.probeinfo}" \
          --rmap "{input.rmap}" \
          --baitmap "{input.baitmap}" \
          --hicup_config "{input.hicup_cfg}" \
          --digest "{input.digest}" \
          --hicup_bin "{TOOLS[hicup]}" \
          --bedtools_bin "{TOOLS[bedtools]}" \
          --perl_bin "{TOOLS[perl]}" \
          --bam_to_fragments "{TOOLS[bam_to_fragments]}" \
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
        mat    = f"{OUTDIR}/{{sample}}/{{sample}}.mat",
        prev   = "logs/hicup/{sample}.Rout"
    output:
        rout  = "logs/capture_eff/{sample}.Rout",
        stats = f"{OUTDIR}/{{sample}}/{{sample}}_CaptureEfficiency.stat"
    threads: THREADS
    conda: "envs/r.yaml"
    log: "logs/capture_eff/{sample}.log"
    shell:
        r"""
        mkdir -p logs/capture_eff
        {TOOLS[rscript]} "{input.script}" \
          --dir "{OUTDIR}/{wildcards.sample}" \
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
        script   = "scripts/01_Prepare_Chicago.R",
        rmap     = DESIGN["rmap"],
        baitmap  = DESIGN["baitmap"],
        bam      = f"{OUTDIR}/{{sample}}/mapped.bam",
        prev     = "logs/capture_eff/{sample}.Rout"
    params:
        design_dir = config["design"]["dir"],
        bam2chicago = config["tools"].get("chicago_bam2chicago", ""),
        make_design = config["tools"].get("chicago_make_design", ""),
        python_bin  = config["tools"].get("python", "python"),
        force_make_design = config.get("chicago", {}).get("force_make_design", False)
    output:
        rout = "logs/chicago/{sample}.Rout",
        out  = f"{OUTDIR}/{{sample}}/chicago_interactions.tsv"
    threads: THREADS
    conda: "envs/chicago.yaml"
    log: "logs/chicago/{sample}.log"
    shell:
        r"""
        mkdir -p logs/chicago
        {TOOLS[rscript]} "{input.script}" \
          --dir "{OUTDIR}/{wildcards.sample}" \
          --sample "{wildcards.sample}" \
          --rmap "{input.rmap}" \
          --baitmap "{input.baitmap}" \
          --bam "{input.bam}" \
          --design_dir "{params.design_dir}" \
          --out "{output.out}" \
          --bam2chicago "{params.bam2chicago}" \
          --make_design "{params.make_design}" \
          --python "{params.python_bin}" \
          --force_make_design {"TRUE" if params.force_make_design else "FALSE"} \
          > {log} 2>&1
        echo "OK" > {output.rout}
        """

# -------------------------------------------------------------------
# 4) ChiCMaxima interaction calling
# -------------------------------------------------------------------
rule chicmaxima:
    input:
        script  = "scripts/02_ChiCMaxima.R",
        rmap    = DESIGN["rmap"],
        baitmap = DESIGN["baitmap"],
        # We’ll glob the BED inside the script, but bam/Rout ensure order:
        prev    = "logs/chicago/{sample}.Rout"
    output:
        rout = "logs/chicmaxima/{sample}.Rout",
        ibed = f"{OUTDIR}/{{sample}}/{{sample}}_Cis.ibed"
    threads: THREADS
    conda: "envs/r.yaml"
    log: "logs/chicmaxima/{sample}.log"
    shell:
        r"""
        mkdir -p logs/chicmaxima
        {TOOLS[rscript]} "{input.script}" \
          --dir "{OUTDIR}/{wildcards.sample}" \
          --sample "{wildcards.sample}" \
          --rmap "{input.rmap}" \
          --baitmap "{input.baitmap}" \
          --perl_bin "{TOOLS[perl]}" \
          --align2ibed "{TOOLS[align2ibed]}" \
          --ibed_out "{output.ibed}" \
          > {log} 2>&1
        echo "OK" > {output.rout}
        """

# -------------------------------------------------------------------
# 5) Generate .hic and optional TADs; then FAN-C insulation
# -------------------------------------------------------------------
rule hic_and_tads:
    input:
        script  = "scripts/03_Generate_hic_and_TADs.R",
        rmap    = DESIGN["rmap"],
        mat     = f"{OUTDIR}/{{sample}}/{{sample}}.mat",
        bam     = f"{OUTDIR}/{{sample}}/mapped.bam"   # not used by the script but keeps ordering
    params:
        java_bin     = config["tools"].get("java", "java"),
        juicer_jar   = config["tools"]["juicer_tools"],
        perl_bin     = config["tools"].get("perl", "perl"),
        bam2juicer   = config["tools"]["bam2juicer"],
        genome_id    = config["params"]["juicer_genome_id"],
        mem_gb       = config["params"]["juicer_mem_gb"],
        q_threshold  = config["params"]["juicer_q"]
    output:
        rout = "logs/hic/{sample}.Rout",
        hic  = f"{OUTDIR}/{{sample}}/{{sample}}.hic"
    threads: THREADS
    # Use a small env with R + Java + Perl; or split into a dedicated env if you prefer.
    conda: "envs/r.yaml"
    log: "logs/hic/{sample}.log"
    shell:
        r"""
        mkdir -p logs/hic {OUTDIR}/{wildcards.sample}
        {TOOLS[rscript]} "{input.script}" \
          --dir "{OUTDIR}/{wildcards.sample}" \
          --sample "{wildcards.sample}" \
          --mat "{input.mat}" \
          --rmap "{input.rmap}" \
          --pre "{OUTDIR}/{wildcards.sample}/{wildcards.sample}_pre.txt" \
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
        hic    = f"{OUTDIR}/{{sample}}/{{sample}}.hic",
        prev   = "logs/hic/{sample}.Rout"
    params:
        resolution = config["params"]["fanc_resolution"],
        balancing  = config["params"]["fanc_balancing"],
        windows    = config["params"]["fanc_windows"]   # list in YAML
    output:
        rout = "logs/insulation/{sample}.Rout",
        ins  = f"{OUTDIR}/{{sample}}/{{sample}}.insulation"
    threads: THREADS
    conda: "envs/fanc.yaml"
    log: "logs/insulation/{sample}.log"
    shell:
        r"""
        mkdir -p logs/insulation
        {TOOLS[rscript]} "{input.script}" \
          --dir "{OUTDIR}/{wildcards.sample}" \
          --sample "{wildcards.sample}" \
          --hic "{input.hic}" \
          --out "{output.ins}" \
          --resolution "{params.resolution}" \
          --balancing "{params.balancing}" \
          --windows "{' '.join(params.windows)}" \
          > {log} 2>&1
        echo "OK" > {output.rout}
        """