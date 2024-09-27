import os
import json
import yaml
import shutil
import glob
from pathlib import Path
from snakemake.utils import validate
from io import StringIO
import pandas as  pd

#mappers = ["minimap2-sr", "bwa-mem", "bowtie2"]
#mappers = ["minimap2-sr", "bwa-mem"]
zymo = config["zymo_genomes"]

offtargets = ["none", "zymo", "assembly"]
offtargets = config["offtargets"]
mappers = config["mappers"]



spike_manifests = {os.path.splitext(os.path.basename(x))[0]: x for x in config["spike_manifests"]}


# read in the config actually used to create the test data
with open(config["testdata_config"], "r") as inf:
    testdata_config = yaml.safe_load(inf)
    testdata_config["total_spike_fractions"] =  [f"{frac:.5f}" for frac in testdata_config["total_spike_fractions"]]
    testdata_config["even_coverages"] =  [f"{frac:.5f}" for frac in testdata_config["even_coverages"]]
    #testdata_config["total_spike_fractions"] =  [f"{frac:.5f}" for frac in [testdata_config["total_spike_fractions"][0]]]

covtypes = expand("depth{x}_spike{y}",
                  x=testdata_config["total_depths"],
                  y=testdata_config["total_spike_fractions"])
covtypes.extend([f"evencoverage{x}"  for x in testdata_config["even_coverages"]])

#manif_base = [f"{rep}_depth10000000_spike{perc}" for rep in range(1,6) for perc in testdata_config["total_spike_fractions]["0.00000", "0.00001", "0.00010", "0.00100", "0.01000", "0.10000", "1.00000"]]
manif_base = expand("{rep}_{cov}",
                    rep=testdata_config["reps"],
                    cov=covtypes
                    )


samples = expand("{manif_base}_offtarget{offtarget}_mapper{mapper}_quanttype{quanttype}_len{lengththresh}_pctid{pctidthresh}",
                 manif_base = manif_base,
                 offtarget = offtargets,
                 mapper=mappers,
                 lengththresh = config["length_thresholds"],
                 pctidthresh = config["pctid_thresholds"],
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
if "bins" in offtargets:
    for i in df.bindir:
        assert os.path.isdir(i), f"{i} doesn't exist"
if "assembly" in offtargets:
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
    mapper="minimap2-sr|bwa-mem|bowtie2",


rule all:
    input:
        f"all_coverage.tsv",
        f"all_spike_counts.tsv",
        f"all_nospike_counts.tsv",

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

def get_mapper_from_sample(wildcards):
     return wildcards.sample.split("mapper")[1].split("_")[0]

def get_length_threshold_from_sample(wildcards):
     return wildcards.sample.split("len")[1].split("_")[0]
def get_pctid_threshold_from_sample(wildcards):
     return wildcards.sample.split("pctid")[1].split("_")[0]

rule run_benchmarking_pipeline:
    # same as in main pipeline, the offtarget is given as a param so it can be
    # optional when specified as "none"
    # the target is specified relative to --directory
    input:
        R1s=lambda wildcards: df.loc[get_sample_base(wildcards.sample), "R1s"].split(","),
        R2s=lambda wildcards: df.loc[get_sample_base(wildcards.sample), "R2s"].split(","),
        spiketable=get_spiketable_from_sample,
    output:
        counts="{sample}/{sample}.spike_reads.tsv",
        r1="{sample}/{sample}_nospike_R1.fastq.gz",
        r2="{sample}/{sample}_nospike_R2.fastq.gz",
        r1spike="{sample}/{sample}_spike_R1.fastq.gz",
        r2spike="{sample}/{sample}_spike_R2.fastq.gz",
    params:
        offtarget=get_offtarget_from_sample,
        profile=profile,
        mapper=get_mapper_from_sample,
        lenthresh = get_length_threshold_from_sample,
        pctidthresh = get_pctid_threshold_from_sample,
        workflow_dir=workflow.current_basedir,
    resources:
        walltime=5*60
    shell:"""
    export SNAKEMAKE_PROFILE={params.profile}
    snakemake --snakefile {params.workflow_dir}/../workflow/Snakefile --directory {wildcards.sample}/ --config R1=[{input.R1s}] \
    R2=[{input.R2s}] offtarget={params.offtarget} mapper={params.mapper} spiketable={input.spiketable} sample={wildcards.sample} \
    pctid_threshold={params.pctidthresh} length_threshold={params.lenthresh} \
    --rerun-incomplete
    """

rule tabulate:
    input:
        counts=expand("{sample}/{sample}.spike_reads.tsv", sample=samples),
    output:
        f"all_coverage.tsv"
    shell: """
    echo -e "sample\tGenome\ttotal_length\tcounts" > {output[0]}
    for rep in {input}
    do
    base=$( basename $rep  | sed "s|.spike_reads.tsv||g")
    cat $rep | sed "s|^|${{base}}\t|g" >> {output[0]}
    done
    """

rule count_missed_spikes:
    input:
        r1="{sample}/{sample}_nospike_R1.fastq.gz",
        r2="{sample}/{sample}_nospike_R2.fastq.gz",
    output:
        counts = "{sample}/{sample}_nospike_counts"
    params:
        tag="nospike"
    shell: """
    # try to catch empty gz fastq
    if [ -n "$(gunzip < {input.r1} | head -c 1 | tr '\\0\\n' __)" ]; then
        zcat  {input.r1} | grep "@" | cut -f1 -d- | uniq -c | sed 's|^|{wildcards.sample}\tR1\t{params.tag}\t|g' > {output.counts}
        zcat  {input.r2} | grep "@" | cut -f1 -d- | uniq -c | sed 's|^|{wildcards.sample}\tR2\t{params.tag}\t|g' >> {output.counts}
    else
        echo -e "{wildcards.sample}\tR1\tnospike\t0\tempty\n{wildcards.sample}\tR2\tnospike\t0\tempty" > {output.counts}
    fi
    """

use rule count_missed_spikes as count_false_spikes with:
    input:
        r1="{sample}/{sample}_spike_R1.fastq.gz",
        r2="{sample}/{sample}_spike_R2.fastq.gz",
    output:
        counts = "{sample}/{sample}_spike_counts"
    params:
        tag="spike"

rule aggregate_missed_spikes:
    input:
        counts=expand("{sample}/{sample}_nospike_counts", sample=samples),
    output:
        f"all_nospike_counts.tsv"
    shell:"""
    cat {input} > {output}
    """

rule aggregate_false_counts:
    input:
        counts=expand("{sample}/{sample}_spike_counts", sample=samples),
    output:
        f"all_spike_counts.tsv"
    shell:"""
    cat {input} > {output}
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
