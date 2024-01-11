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

manif = "R1s,R2s,bindir,assembly\n" + "\n".join([f"{pathbase}{x}_R1.fastq.gz,{pathbase}{x}_R2.fastq.gz,{pathbase}rawbinning_{x}/concoct/concoct_bins/,{pathbase}spades_{x}.assembly.fasta" for x in manif_base])

df = pd.read_csv(StringIO(manif), sep=",")
df["sample"] = manif_base
df = df.set_index("sample")
print(df)
onstart:
    print(df)

configfile: os.path.join(str(workflow.current_basedir), "../config/config.yaml")




rule all:
    input:
        expand("{sample}_all_coverage.tsv", sample=manif_base)

rule:
    input:
        R1s=lambda wildcards: df.loc[wildcards.sample, "R1s"].split(","),
        R2s=lambda wildcards: df.loc[wildcards.sample, "R2s"].split(","),
        bindir=lambda wildcards: df.loc[wildcards.sample, "bindir"].split(","),
        assembly=lambda wildcards: df.loc[wildcards.sample, "assembly"].split(","),
    output:
        "{sample}_all_coverage.tsv"
    shell: """
    snakemake --snakefile benchmarking/Snakefile --directory $PWD/benchmarking_results/ --config R1={input.R1} \
    R2={input.R2} bindir={input.bindir} assembly={input.assembly}
    """
