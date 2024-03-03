# Generating Full Benchmark simulated dataset
A yaml config file is used to define the number of replicates the test data workflow generates, the total depth, and the proportion of the spikes in the sample.

```
snakemake --snakefile testdata/Snakefile --directory $PWD/benchmarking/data/ --config benchmarking/benchmarking_sim_data_config.yaml

head -n 1 benchmarking/data/3_depth10000000_spike0.10000.statsfastq |  tr -s " " "\t" > benchmarking/data_summary.tsv && cat benchmarking/data/*_depth*_spike*.statsfastq | grep -v num_seqs |  tr -s " " "\t" >> benchmarking/data_summary.tsv

```

## Setting up benchmarking conditions
There are 3 main things to vary: what spikes are used, what alignment method, and what set of off-target sequences should be considered (if any).

For the off targets, we considered 4 options:
- no offtarget sequnces
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

### Primary/Other Split Spikes
Given that some spikes might contain sub-par contigs or may contain multi-copy plasmids, we also split the spikes into two files: the primary sequence, and all other sequences
```
for f in benchmarking/spikes/raw/*.fa
do
	python benchmarking/scripts/split_primary_other.py $f benchmarking/spikes/split/
done

```

### rRNA Masked Spikes
Bedtools is needed for this script; It masks the rRNA from the references in benchmarking/spikes/split/.

```
bash benchmarking/scripts/get_spikes_and_mask_rRNAs.sh benchmarking/spikes/rrnamasked/

```


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
snakemake --snakefile benchmarking/run.smk --directory $PWD/benchmarking/v2/ --config benchmarkdata=$PWD/benchmarking/data/ spike_manifests=[$PWD/benchmarking/manifests/raw.csv,$PWD/benchmarking/manifests/rrnamasked.csv,$PWD/benchmarking/manifests/split.csv] testdata_config=$PWD/benchmarking/benchmarking_sim_data_config.yaml zymo_genomes=$PWD/D6331.refseq/genomes/
```
