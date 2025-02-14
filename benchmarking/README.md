# Generating Full Benchmark simulated dataset
A yaml config file is used to define the number of replicates the test data workflow generates, the total depth, and the proportion of the spikes in the sample. For benchmarking downstream taxonomic profiling, we generated test data with and without the T. reesi spike.

```
snakemake --snakefile testdata/Snakefile --directory $PWD/benchmarking/data/ --configfile benchmarking/benchmarking_sim_data_config.yaml
snakemake --snakefile testdata/Snakefile --directory $PWD/benchmarking/data_no_treesei/ --configfile benchmarking/benchmarking_sim_data_no_treesei_config.yaml

head -n 1 benchmarking/data/1_depth10000000_spike0.10000.statsfastq |  tr -s " " "\t" > benchmarking/data_summary.tsv && \
  cat benchmarking/data/*_depth*_spike*.statsfastq | grep -v num_seqs |  tr -s " " "\t" >> benchmarking/data_summary.tsv &&  \
  cat benchmarking/data/*evenc*.statsfastq | grep -v num_seqs | tr -s " " "\t" >> benchmarking/data_summary.tsv

head -n 1 benchmarking/data/1_depth10000000_spike0.10000.statsfastq |  tr -s " " "\t" > benchmarking/data_summary_no_treesei.tsv && \
  cat benchmarking/data_no_treesei/*_depth*_spike*.statsfastq | grep -v num_seqs |  tr -s " " "\t" >> benchmarking/data_summary_no_treesei.tsv &&  \
  cat benchmarking/data_no_treesei/*evenc*.statsfastq | grep -v num_seqs | tr -s " " "\t" >> benchmarking/data_summary_no_treesei.tsv

```

## Setting up benchmarking conditions
There are 3 main things to vary: what spikes are used, what alignment method, and what set of off-target sequences should be considered (if any).

For the off targets, we considered 4 options:
- no offtarget sequences
- the sample's assembly
- the sample's binned contigs
- the GutZymo mock

The Gutzymo mock can be downloaded as follows:


```sh
wget -N https://s3.amazonaws.com/zymo-files/BioPool/D6331.refseq.zip
unzip D6331.refseq.zip
```

### Raw Spikes
```
bash benchmarking/scripts/get_spikes.sh benchmarkingspikes/raw/
```

### BED files for primary, no-rRNA, and marker quantification
```
bash benchmarking/scripts/write_bed_files.sh benchmarking/beds/
```


<!-- ### rRNA Masked Spikes -->
<!-- Bedtools is needed for this script; It masks the rRNA from the references in benchmarking/spikes/split/. -->

<!-- ``` -->
<!-- bash benchmarking/scripts/get_spikes_and_mask_rRNAs.sh benchmarking/spikes/rrnamasked/ -->

<!-- ``` -->


### Metaphlan Marker Masked Spikes

<TBD>

### Creaing spike Manifest

The manifest file used by spikequant takes 3 columns: spike name , whether it should be used for quantification or just for depletion, and a path to the fasta.

At the end of this, you should have the following files:
```sh
benchmarking/spikes/
├── raw
│   ├── Haloarcula_hispanica.fa
│   ├── Haloarcula_hispanica.fa.fai
│   ├── Salinibacter_ruber.fa
│   └── Trichoderma_reesei.fa
├── rrnamasked
│   ├── Haloarcula_hispanica.bed
│   ├── Haloarcula_hispanica_other_masked.fa
│   ├── Haloarcula_hispanica_primary_masked.fa
│   ├── Salinibacter_ruber.bed
│   ├── Salinibacter_ruber_other_masked.fa
│   ├── Salinibacter_ruber_primary_masked.fa
│   ├── Trichoderma_reesei_other.fa
│   └── Trichoderma_reesei_primary.fa
└── split
    ├── Haloarcula_hispanica_other.fa
    ├── Haloarcula_hispanica_primary.fa
    ├── Salinibacter_ruber_other.fa
    ├── Salinibacter_ruber_primary.fa
    ├── Trichoderma_reesei_other.fa
    └── Trichoderma_reesei_primary.fa
```

A manifest for the `split` spikes would look like the following, replacing $PWD so that the paths are absolute:
```
Trichoderma_reesei,no,$PWD/benchmarking/spikes/split/Trichoderma_reesei_other.fa
Trichoderma_reesei,yes,$PWD/benchmarking/spikes/split/Trichoderma_reesei_primary.fa
Salinibacter_ruber,no,$PWD/benchmarking/spikes/split/Salinibacter_ruber_other.fa
Salinibacter_ruber,yes,$PWD/benchmarking/spikes/split/Salinibacter_ruber_primary.fa
Haloarcula_hispanica,no,$PWD/benchmarking/spikes/split/Haloarcula_hispanica_other.fa
Haloarcula_hispanica,yes,$PWD/benchmarking/spikes/split/Haloarcula_hispanica_primary.fa
```


## Executing the Benchmarking dataset

```
snakemake  --snakefile benchmarking/run.smk --directory $PWD/benchmarking/v5/ --config benchmarkdata=$PWD/benchmarking/data/ spike_manifests=[$PWD/benchmarking/manifests/raw.csv,$PWD/benchmarking/manifests/rrnamasked.csv,$PWD/benchmarking/manifests/split.csv,$PWD/benchmarking/manifests/markers.csv] testdata_config=$PWD/benchmarking/benchmarking_sim_data_config.yaml zymo_genomes=$PWD/D6331.refseq/genomes/ length_thresholds=[.97] pctid_thresholds=[.97] --keep-going  --rerun-triggers mtime -f all_coverage.tsv  --rerun-incomplete

snakemake  --snakefile benchmarking/run.smk --directory $PWD/benchmarking/v5_no_treesei/ --config benchmarkdata=$PWD/benchmarking/data_no_treesei/ spike_manifests=[$PWD/benchmarking/manifests/raw.csv,$PWD/benchmarking/manifests/rrnamasked.csv,$PWD/benchmarking/manifests/split.csv,$PWD/benchmarking/manifests/markers.csv] testdata_config=$PWD/benchmarking/benchmarking_sim_data_no_treesei_config.yaml zymo_genomes=$PWD/D6331.refseq/genomes/ length_thresholds=[.97] pctid_thresholds=[.97] --keep-going  --rerun-triggers mtime -f all_coverage.tsv  --rerun-incomplete
```



# Benchmarking performance against undepleted samples
```
find  $PWD/benchmarking/v5/ -name "*_R1.fastq.gz" | grep offtargetnone_mapperminimap2 | grep nospike | grep evencoverage > benchmarking/v5/manifest_R1 && \
find $PWD/benchmarking/data/ -name "*_R1.fastq.gz"  >> benchmarking/v5/manifest_R1

find  $PWD/benchmarking/v5_no_treesei/ -name "*_R1.fastq.gz" | grep offtargetnone_mapperminimap2 | grep nospike | grep evencoverage > benchmarking/v5/manifest_R1_no_treesei
find $PWD/benchmarking/data_no_treesei/ -name "*_R1.fastq.gz"  >> benchmarking/v5/manifest_R1_no_treesei

snakemake --snakefile benchmarking/despiked_vs_neat_analysis.smk --directory benchmarking/despike_v_neat/  --config metaphlan_db=/data/brinkvd/resources/dbs/metaphlan/mpa_vJun23_CHOCOPhlAnSGB_202403/ docker_biobakery=docker://ghcr.io/vdblab/biobakery-profiler:20240513a choco_db=/data/brinkvd/resources/dbs/chocophlan/v201901_v31/ uniref90_db=/data/brinkvd/resources/dbs/uniref90_diamond/v201901b/ manifest=$PWD/benchmarking/v5/manifest_R1

snakemake --snakefile benchmarking/despiked_vs_neat_analysis.smk --directory benchmarking/despike_v_neat_no_treesei/  \
    --config metaphlan_db=/data/brinkvd/resources/dbs/metaphlan/mpa_vJun23_CHOCOPhlAnSGB_202403/ docker_biobakery=docker://ghcr.io/vdblab/biobakery-profiler:20240513a choco_db=/data/brinkvd/resources/dbs/chocophlan/v201901_v31/ uniref90_db=/data/brinkvd/resources/dbs/uniref90_diamond/v201901b/ manifest=$PWD/benchmarking/v5/manifest_R1_no_treesei
```


# Further parameter sweeps:
```
snakemake  --snakefile benchmarking/run.smk --directory $PWD/benchmarking/parameter_sweep/ --config benchmarkdata=$PWD/benchmarking/data/ spike_manifests=[$PWD/benchmarking/manifests/raw.csv,$PWD/benchmarking/manifests/rrnamasked.csv,$PWD/benchmarking/manifests/split.csv,$PWD/benchmarking/manifests/markers.csv] testdata_config=$PWD/benchmarking/benchmarking_sim_data_config.yaml zymo_genomes=$PWD/D6331.refseq/genomes/ mappers=["minimap2-sr"] offtargets=["none"] length_thresholds=[0.5,0.85,0.9,0.95,0.97,0.99] pctid_thresholds=[0.5,0.85,0.95,0.97,.99] --keep-going  --rerun-triggers mtime -f all_coverage.tsv  --rerun-incomplete

```
