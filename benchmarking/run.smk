import os
import json
import yaml
import shutil

from pathlib import Path
from snakemake.utils import validate
from io import StringIO
import pandas as  pd

mappers = ["minimap2-sr", "bwa-mem"]
zymo = config["zymo_genomes"]

# make this a dict for easier lookup
spike_manifests = {os.path.splitext(os.path.basename(x))[0]: x for x in config["spike_manifests"]}

#    "raw": "/data/brinkvd/users/watersn/spikequant/benchmarking/manifests/raw.csv",
#    "rrnamasked": "/data/brinkvd/users/watersn/spikequant/benchmarking/manifests/rrnamasked.csv",
#    "split": "/data/brinkvd/users/watersn/spikequant/benchmarking/manifests/split.csv"
#}
#dbtypes = [f"{x}-{y}" for x in ["spike", "assembly", "bins", "mock"] for y in quanttype]

# read in the config actually used to create the test data
with open(config["testdata_config"], "r") as inf:
    testdata_config = yaml.safe_load(inf)
    testdata_config["total_spike_fractions"] =  [f"{frac:.5f}" for frac in testdata_config["total_spike_fractions"]]


#manif_base = [f"{rep}_depth10000000_spike{perc}" for rep in range(1,6) for perc in testdata_config["total_spike_fractions]["0.00000", "0.00001", "0.00010", "0.00100", "0.01000", "0.10000", "1.00000"]]
manif_base = expand("{rep}_depth10000000_spike{perc}",
                    rep=testdata_config["reps"],
                    depth=testdata_config["total_depths"],
                    perc=testdata_config["total_spike_fractions"])
samples = expand("{manif_base}_offtarget{offtarget}_mapper{mapper}_quanttype{quanttype}",
                 manif_base = manif_base,
                 offtarget = ["none", "zymo", "bins", "assembly"],
                 mapper=mappers,
                 quanttype=spike_manifests.keys())

pathbase=config["benchmarkdata"]
manif = "R1s,R2s,bindir,assembly\n" + "\n".join([f"{pathbase}{x}_R1.fastq.gz,{pathbase}{x}_R2.fastq.gz,{pathbase}metawrap/refined_binning_{x}/metawrap_70_10_bins/,{pathbase}spades_{x}.assembly.fasta" for x in manif_base])

df = pd.read_csv(StringIO(manif), sep=",")
df["sample"] = manif_base
df = df.set_index("sample")

# sanity check the existence of these files
for i in df.R1s:
    assert os.path.exists(i), f"{i} doesn't exist"
for i in df.R2s:
    assert os.path.exists(i), f"{i} doesn't exist"
for i in df.bindir:
    assert os.path.isdir(i), f"{i} doesn't exist"
for i in df.assembly:
    assert os.path.exists(i), f"{i} doesn't exist"

onstart:
    print(df)

configfile: os.path.join(str(workflow.current_basedir), "../config/config.yaml")


envvars:
    "TMPDIR",
    "SNAKEMAKE_PROFILE"

profile = os.getenv("SNAKEMAKE_PROFILE")

# "raw|split|rrnamasked",
wildcard_constraints:
    quanttype="|".join(spike_manifests.keys()),
    offtarget="none|zymo|bins|assembly",
    mapper="minimap2-sr|bwa-mem",


rule all:
    input:
        expand("{sample}/coverm/{sample}.coverage_mqc.tsv", sample=samples),
        f"all_coverage.tsv"

def get_sample_base(sample):
    """ sample_base is the part of the sample name relevant to data generation,
    as opposed to the benchmarking parameters like the mapper of the offtarget.
    We need this to look up against the file manifest "df"
    """
    return(sample.split("_offtarget")[0])

def get_offtarget_from_sample(wildcards):
    offtarget = wildcards.sample.split("offtarget")[1].split("_")[0]
    if offtarget == "none":
        return "none"
    elif offtarget == "bins":
        return df.loc[get_sample_base(wildcards.sample), "bindir"].split(",")
    elif offtarget == "assembly":
        return df.loc[get_sample_base(wildcards.sample), "assembly"].split(",")
    elif offtarget == "zymo":
        return zymo
    else:
        raise ValueError("offtarget not found in sample")

def get_spiketable_from_sample(wildcards):
    qtype=wildcards.sample.split("quanttype")[1].split("_")[0]
    return(spike_manifests[qtype])

rule run_benchmarking_pipeline:
    # same as in main pipeline, the offtarget is given as a param so it can be
    # optional when specified as "none"
    # the target is specified relative to --directory
    input:
        R1s=lambda wildcards: df.loc[get_sample_base(wildcards.sample), "R1s"].split(","),
        R2s=lambda wildcards: df.loc[get_sample_base(wildcards.sample), "R2s"].split(","),
        spiketable=get_spiketable_from_sample,
    output:
        coverm="{sample}/coverm/{sample}.coverage_mqc.tsv"
    params:
#        offtarget=lambda wc: get_offtarget_from_sample(wc),
        target = lambda wc, output: os.path.join("coverm", os.path.basename(output.coverm)),
        offtarget=get_offtarget_from_sample,
        profile=profile,
        workflow_dir=workflow.current_basedir,
    resources:
        walltime=5*60
    shell:"""
    export SNAKEMAKE_PROFILE={params.profile}
    snakemake --snakefile {params.workflow_dir}/../workflow/Snakefile --directory {wildcards.sample}/ --config R1=[{input.R1s}] \
    R2=[{input.R2s}] offtarget={params.offtarget} spiketable={input.spiketable} -f {params.target} \
    --rerun-incomplete
    """


rule tabulate:
    input:
        expand("{sample}/coverm/{sample}.coverage_mqc.tsv", sample=samples)
    output:
        f"all_coverage.tsv"
    shell: """
    echo -e "sample\tGenome\tMean\tRelative_Abundance_Perc\tTrimmed_Mean\tCovered_Bases\tVariance\tLength\tRead_Count\tReads_per_base\tRPKM\tTPM" > {output[0]}
    for rep in {input}
    do
    base=$( basename $rep  | sed "s|.coverage_mqc.tsv||g")
    tail -n+4 $rep | sed "s|^|${{base}}\t|g" >> {output[0]}
    done
    """


# rule merge_tabulated_tables:
#     input:
#         expand("{sample}/{sample}_all_coverage.tsv", sample=manif_base),
#     output:
#         f"merged_benchmarking_coverage.tsv"
#     shell: """
#     head -n 1 {input[0]} > {output}
#     for i in {input}
#     do
#     cat $i | grep -v "Relative_Abundance_Perc" >> {output}
#     done
#     """
