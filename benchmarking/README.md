# Generating Full Benchmark simulated dataset
```
snakemake --snakefile testdata/Snakefile --directory $PWD/benchmarking/data/ --config benchmarking/benchmarking_sim_data_config.yaml
```

## Executing the Benchmarking dataset

```
snakemake --snakefile benchmarking/run.smk
```
