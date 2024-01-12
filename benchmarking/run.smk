import os
import json
import yaml
import shutil

from pathlib import Path
from snakemake.utils import validate
from io import StringIO
import pandas as  pd

manif_base = [f"{rep}_depth10000000_spike{perc}" for rep in range(1,6) for perc in ["0.00000", "0.00001", "0.00010", "0.00100", "0.01000", "0.10000", "1.00000"]]

pathbase="/lila/data/brinkvd/users/watersn/spikequant/benchmarking/data/"

manif = "R1s,R2s,bindir,assembly\n" + "\n".join([f"{pathbase}{x}_R1.fastq.gz,{pathbase}{x}_R2.fastq.gz,{pathbase}metawrap/refined_binning_{x}/metawrap_70_10_bins/,{pathbase}spades_{x}.assembly.fasta" for x in manif_base])

df = pd.read_csv(StringIO(manif), sep=",")
df["sample"] = manif_base
df = df.set_index("sample")

onstart:
    print(df)

configfile: os.path.join(str(workflow.current_basedir), "../config/config.yaml")


envvars:
    "TMPDIR",
    "SNAKEMAKE_PROFILE"

profile = os.getenv("SNAKEMAKE_PROFILE")


rule all:
    input:
        expand("{sample}/{sample}_all_coverage.tsv", sample=manif_base),
        "merged_benchmarking_coverage.tsv",

rule run_benchmarking_pipeline:
    input:
        R1s=lambda wildcards: df.loc[wildcards.sample, "R1s"].split(","),
        R2s=lambda wildcards: df.loc[wildcards.sample, "R2s"].split(","),
        bindir=lambda wildcards: df.loc[wildcards.sample, "bindir"].split(","),
        assembly=lambda wildcards: df.loc[wildcards.sample, "assembly"].split(","),
    output:
        "{sample}/{sample}_all_coverage.tsv"
    params:
        profile=profile,
        workflow_dir=workflow.current_basedir,
    resources:
        walltime=5*60
    shell:"""
    export SNAKEMAKE_PROFILE={params.profile}
    # downloading references can fail.
    if [ -d "{wildcards.sample}/D6331.refseq" ]
    then
    rm -r {wildcards.sample}/D6331.refseq*
    fi
    snakemake --snakefile {params.workflow_dir}/Snakefile --directory {wildcards.sample}/ --config R1=[{input.R1s}] \
    R2=[{input.R2s}] bindir={input.bindir} assembly={input.assembly} sample={wildcards.sample} -f {wildcards.sample}_all_coverage.tsv \
    --rerun-incomplete
    """

rule merge_tabulated_tables:
    input:
        expand("{sample}/{sample}_all_coverage.tsv", sample=manif_base),
    output:
        f"merged_benchmarking_coverage.tsv"
    shell: """
    head -n 1 {input[0]} > {output}
    for i in {input}
    do
    cat $i | grep -v "Relative_Abundance_Perc" >> {output}
    done
    """
