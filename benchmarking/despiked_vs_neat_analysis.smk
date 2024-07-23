import os
import json
import yaml
import shutil
import glob
from pathlib import Path
from snakemake.utils import validate
from io import StringIO
import pandas as  pd



#manifest = pd.read_csv(os.path.join(str(workflow.current_basedir), "../benchmarking/v5/manifest_R1"), header=None, names=["R1"])
manifest = pd.read_csv(config["manifest"], header=None, names=["R1"])
manifest["R2"] = manifest["R1"].str.replace("_R1", "_R2")
# remove samples we cant evaluate (eg all-spike samples should have empty depiked fastqs)
minsize = 10000
manifest["sizes"] = [os.path.getsize(x) for x in manifest["R1"]]
print(manifest.shape)
manifest = manifest[manifest["sizes"] > minsize]
print(manifest.shape)
#############
manifest["present"] = [os.path.exists(x) and os.path.exists(y) for x,y in zip(manifest["R1"], manifest["R2"])]
manifest = manifest[manifest["present"]]
print(manifest.shape)
#############
manifest["sample"] = manifest["R1"].apply(lambda x: Path(x).stem.split("_R1")[0])
manifest["tmp"] = manifest["sample"] # lazy
manifest = manifest.set_index("tmp")
#manifest = manifest.tail(20)
print(manifest.tail(10))


metaphlan_results = [f"metaphlan/{x}_metaphlan3_profile.txt" for x in manifest["sample"]]
humann_gf_results = [f"humann/{x}_humann3_genefamilies.tsv" for x in manifest["sample"]]

rule all:
    input:
        metaphlan_results,
        humann_gf_results,




def get_R1_R2(wildcards):
    return {
        "R1": manifest.loc[wildcards.sample, "R1"],
        "R2": manifest.loc[wildcards.sample, "R2"],
    }

# If your reads lack /1 and /2 then humann will miscount
# the cat paired end reads and metaphlan and humann3 part
rule cat_pair_and_reformat:
    input:
        unpack(get_R1_R2)
    output:
        joined=temp("tmp/{sample}_cat.fastq.gz"),
    container: "docker://staphb/bbtools:38.97"
    shell:
        """
        reformat.sh in={input.R1} in2={input.R2}  out={output.joined}  addslash int
        """


rule humann3_run_uniref90:
    input:
        fastq="tmp/{sample}_cat.fastq.gz",
        metaphlan_profile="metaphlan/{sample}_metaphlan3_profile.txt",
        choco_db=config["choco_db"],
        uniref90_db=config["uniref90_db"],
    output:
        ab="humann/{sample}_humann3_pathabundance.tsv",
        genefam="humann/{sample}_humann3_genefamilies.tsv",
        stats="humann/{sample}_humann3_stats.txt",
    params:
        out_prefix="{sample}_humann3",
        out_dir=lambda w, output: os.path.dirname(output[0]),
    container:
        config["docker_biobakery"]
    resources:
        mem_mb=lambda wildcards, attempt, input: attempt
        * 1024
        * max(input.fastq.size // 1000000000, 1)
        * 12,
        runtime=12 * 60,
    threads: 48
    # we have an extra log in case there is an error with humann.  Cause
    # we skip the built in logging because
    # they dont actually log errors to their --o-log :(

    log:
        e="logs/humann_{sample}.e",
        o="logs/humann_{sample}.o",
    shell:
        """
        # see https://forum.biobakery.org/t/metaphlan-v4-0-2-and-huma-3-6-metaphlan-taxonomic-profile-provided-was-not-generated-with-the-expected-database/4296/8
        cat {input.metaphlan_profile} | cut -f 1-4 > {wildcards.sample}_tmp_metaphlan.tsv
        humann \
            --input {input.fastq} \
            --output {params.out_dir} \
            --output-basename {params.out_prefix} \
            --search-mode uniref90 \
            --remove-column-description-output \
            --protein-database {input.uniref90_db} \
            --nucleotide-database {input.choco_db} \
            --taxonomic-profile {wildcards.sample}_tmp_metaphlan.tsv  \
            --remove-temp-output \
            --output-max-decimals 5 \
            --threads {threads} \
            > {log.o} 2> {log.e}

        cat {log.o} | \
            grep "reads" | \
            grep "\: [0-9]" \
            > {output.stats}
        """





# Metaphlan3 and Strainphlan
rule metaphlan_run:
    input:
        fastq="tmp/{sample}_cat.fastq.gz",
        db=config["metaphlan_db"],
    output:
        outfile="metaphlan/{sample}_metaphlan3_profile.txt",
        sam="metaphlan/{sample}.sam.bz2",
    container:
        config["docker_biobakery"]
    resources:
        mem_mb=16 * 1024,
        runtime=20 * 60,
    params:
        dbbase = os.path.basename(os.path.dirname(config["metaphlan_db"]))
    threads: 48
    log:
        e="logs/metaphlan_{sample}.e",
    shell:
        """
        # the presense of this file causes an error from metaphlan
        # which makes rerunning irritating
        if [ -f "{input.fastq}.bowtie2out.txt" ]
        then
            rm {input.fastq}.bowtie2out.txt
        fi
        export METAPHLAN_BOWTIE2_DB={input.db}
        metaphlan {input.fastq} \
            --bowtie2db {input.db} \
            --index {params.dbbase} \
            --input_type fastq \
            --sample_id  {wildcards.sample} \
            -s {output.sam} \
            --add_viruses \
            --unclassified_estimation \
            --nproc {threads} \
            -t rel_ab_w_read_stats \
            -o {output.outfile}
        """















# mappers = ["minimap2-sr", "bwa-mem", "bowtie2"]
# zymo = config["zymo_genomes"]
# # make this a dict for easier lookup
# spike_manifests = {os.path.splitext(os.path.basename(x))[0]: x for x in config["spike_manifests"]}


# # read in the config actually used to create the test data
# with open(config["testdata_config"], "r") as inf:
#     testdata_config = yaml.safe_load(inf)
#     testdata_config["total_spike_fractions"] =  [f"{frac:.5f}" for frac in testdata_config["total_spike_fractions"]]
#     #testdata_config["total_spike_fractions"] =  [f"{frac:.5f}" for frac in [testdata_config["total_spike_fractions"][0]]]


# #manif_base = [f"{rep}_depth10000000_spike{perc}" for rep in range(1,6) for perc in testdata_config["total_spike_fractions]["0.00000", "0.00001", "0.00010", "0.00100", "0.01000", "0.10000", "1.00000"]]
# manif_base = expand("{rep}_depth10000000_spike{perc}",
#                     rep=testdata_config["reps"],
#                     depth=testdata_config["total_depths"],
#                     perc=testdata_config["total_spike_fractions"])

# sample_bases = expand("{manif_base}_offtarget{offtarget}_mapper{mapper}",
#                  manif_base = manif_base,
#                  offtarget = ["none", "zymo", "bins", "assembly"],
#                  mapper=mappers)

# samples = expand("{manif_base}_offtarget{offtarget}_mapper{mapper}_quanttype{quanttype}",
#                  manif_base = manif_base,
#                  offtarget = ["none", "zymo", "bins", "assembly"],
#                  mapper=mappers,
#                  quanttype=spike_manifests.keys())

# pathbase=config["benchmarkdata"]
# manif = "R1s,R2s,bindir,assembly\n" + "\n".join([f"{pathbase}{x}_R1.fastq.gz,{pathbase}{x}_R2.fastq.gz,{pathbase}metawrap/refined_binning_{x}/metawrap_70_10_bins/,{pathbase}spades_{x}.assembly.fasta" for x in manif_base])

# df = pd.read_csv(StringIO(manif), sep=",")
# df["sample"] = manif_base
# df = df.set_index("sample")

# # sanity check the existence of these files
# for i in df.R1s:
#     assert os.path.exists(i), f"{i} doesn't exist"
# for i in df.R2s:
#     assert os.path.exists(i), f"{i} doesn't exist"
# for i in df.bindir:
#     assert os.path.isdir(i), f"{i} doesn't exist"
# for i in df.assembly:
#     assert os.path.exists(i), f"{i} doesn't exist"

# onstart:
#     print(df)

# configfile: os.path.join(str(workflow.current_basedir), "../config/config.yaml")


# envvars:
#     "TMPDIR",
#     "SNAKEMAKE_PROFILE"

# profile = os.getenv("SNAKEMAKE_PROFILE")

# # "raw|split|rrnamasked",
# wildcard_constraints:
#     quanttype="|".join(spike_manifests.keys()),
#     offtarget="none|zymo|bins|assembly",
#     mapper="minimap2-sr|bwa-mem|bowtie2",


# rule all:
#     input:
#         f"all_coverage.tsv",

# def get_sample_base(sample):
#     """ sample_base is the part of the sample name relevant to data generation,
#     as opposed to the benchmarking parameters like the mapper of the offtarget.
#     We need this to look up against the file manifest "df"
#     """
#     return(sample.split("_offtarget")[0])

# def get_offtarget_from_sample(wildcards):
#     offtarget = wildcards.sample.split("offtarget")[1].split("_")[0]
#     if offtarget == "none":
#         return "none"
#     elif offtarget == "bins":
#         return df.loc[get_sample_base(wildcards.sample), "bindir"].split(",")
#     elif offtarget == "assembly":
#         return df.loc[get_sample_base(wildcards.sample), "assembly"].split(",")
#     elif offtarget == "zymo":
#         return zymo
#     else:
#         raise ValueError("offtarget not found in sample")

# def get_spiketable_from_sample(wildcards):
#     qtype=wildcards.sample.split("quanttype")[1].split("_")[0]
#     return(spike_manifests[qtype])

# def get_mapper_from_sample(wildcards):
#      return wildcards.sample.split("mapper")[1].split("_")[0]


# rule run_benchmarking_pipeline:
#     # same as in main pipeline, the offtarget is given as a param so it can be
#     # optional when specified as "none"
#     # the target is specified relative to --directory
#     input:
#         R1s=lambda wildcards: df.loc[get_sample_base(wildcards.sample), "R1s"].split(","),
#         R2s=lambda wildcards: df.loc[get_sample_base(wildcards.sample), "R2s"].split(","),
#         spiketable=get_spiketable_from_sample,
#     output:
#         counts="{sample}/{sample}.spike_reads.tsv",
#     params:
#         offtarget=get_offtarget_from_sample,
#         profile=profile,
#         mapper=get_mapper_from_sample,
#         workflow_dir=workflow.current_basedir,
#     resources:
#         walltime=5*60
#     shell:"""
#     export SNAKEMAKE_PROFILE={params.profile}
#     snakemake --snakefile {params.workflow_dir}/../workflow/Snakefile --directory {wildcards.sample}/ --config R1=[{input.R1s}] \
#     R2=[{input.R2s}] offtarget={params.offtarget} mapper={params.mapper} spiketable={input.spiketable} sample={wildcards.sample} \
#     --rerun-incomplete
#     """

# rule tabulate:
#     input:
#         counts=expand("{sample}/{sample}.spike_reads.tsv", sample=samples),
#     output:
#         f"all_coverage.tsv"
#     shell: """
#     echo -e "sample\tGenome\ttotal_length\tcounts" > {output[0]}
#     for rep in {input}
#     do
#     base=$( basename $rep  | sed "s|.spike_reads.tsv||g")
#     cat $rep | sed "s|^|${{base}}\t|g" >> {output[0]}
#     done
#     """
