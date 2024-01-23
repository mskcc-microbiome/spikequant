# Generating Full Benchmark simulated dataset
```
snakemake --snakefile testdata/Snakefile --directory $PWD/benchmarking/data/ --config benchmarking/benchmarking_sim_data_config.yaml

head -n 1 benchmarking/data/3_depth10000000_spike0.10000.statsfastq |  tr -s " " "\t" > benchmarking/data_summary.tsv && cat benchmarking/data/*_depth*_spike*.statsfastq | grep -v num_seqs |  tr -s " " "\t" >> benchmarking/data_summary.tsv

```

## Executing the Benchmarking dataset

```
snakemake --snakefile benchmarking/run.smk
```
