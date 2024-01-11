# Generating Full Benchmark dataset
```
snakemake --snakefile testdata/Snakefile --directory $PWD/benchmarking/data/ --config benchmarking/benchmarking_sim_data_config.yaml
```

## Executing the Benchmarking dataset

```
snakemake --snakefile testdata/Snakefile --directory $PWD/benchmarking/data/ --configfile benchmarking/benchmarking_sim_data_config.yaml
snakemake --snakefile benchmarking/run.smk
```
