# Generating Full Benchmark simulated dataset
```
snakemake --snakefile testdata/Snakefile --directory $PWD/benchmarking/data/ --config benchmarking/benchmarking_sim_data_config.yaml

head -n 1 benchmarking/data/3_depth10000000_spike0.10000.statsfastq |  tr -s " " "\t" > benchmarking/data_summary.tsv && cat benchmarking/data/*_depth*_spike*.statsfastq | grep -v num_seqs |  tr -s " " "\t" >> benchmarking/data_summary.tsv

```

## Setting up benchmarking conditions
There are several core things to vary: what spikes are used, what alignment method, and what set of off-target sequences should be considered (if any).



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

A manifest for the `split` spikes would look like the following:
```
Trichoderma_reesei,no,benchmarking/spikes/split/Trichoderma_reesei_other.fa
Trichoderma_reesei,yes,benchmarking/spikes/split/Trichoderma_reesei_primary.fa
Salinibacter_ruber,no,benchmarking/spikes/split/Salinibacter_ruber_other.fa
Salinibacter_ruber,yes,benchmarking/spikes/split/Salinibacter_ruber_primary.fa
Haloarcula_hispanica,no,benchmarking/spikes/split/Haloarcula_hispanica_other.fa
Haloarcula_hispanica,yes,benchmarking/spikes/split/Haloarcula_hispanica_primary.fa
```

## Executing the Benchmarking dataset

```
snakemake --snakefile benchmarking/run.smk --config benchmarkdata=$PWD/benchmarking/data/
```
