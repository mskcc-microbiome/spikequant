import os
import json
import yaml
import shutil

from pathlib import Path
from snakemake.utils import validate
from io import StringIO
import pandas

manif_base = [f"{rep}_depth10000000_spike{perc}" for rep in range(1,6) for perc in ["0", "0.00001", "0.00010", "0.00100", "0.01000", "0.10000", "1.00000"]]

bindir: "../testdata/data/metawrap/rawbinning_1_depth100000_spike0.10000/concoct/concoct_bins/"
assembly: "../testdata/data/spades_1_depth100000_spike0.10000.assembly.fasta"
pathbase="/lila/data/brinkvd/users/watersn/spikequant/benchmarking/data/"

manif = "R1s,R2s,bindir,assembly\n" + "\n".join([f"{pathbase}{x}_R1.fastq.gz,{pathbase}{x}_R2.fastq.gz,{pathbase}rawbinning_{x}/concoct/concoct_bins/,{pathbase}{x}.assembly.fasta" for x in manif_base])

df = pd.read_csv(StringIO(manif), sep=",")
df["sample"] = manif_base


configfile: os.path.join(str(workflow.current_basedir), "../config/config.yaml")



covermreports = expand('coverm/{sample}_and_{dbtype}_via_{mapper}_bins.coverage_mqc.tsv',
                       sample=config["sample"],
                       mapper=mappers,
                       dbtype=dbtypes)


rule all:
    input:
        despiked_reads,
        ani_results,
        covermreports,
        "all_coverage.tsv"
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
