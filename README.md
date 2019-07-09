# snakemake-bismark-methyldackel

[![Build Status](https://travis-ci.org/dohlee/snakemake-bismark-methyldackel.svg?branch=master)](https://travis-ci.org/dohlee/snakemake-bismark-methyldackel)
![Reproduction Status](https://img.shields.io/endpoint.svg?url=https://gist.githubusercontent.com/dohlee/39755d8246a88cea530fa72706397478/raw/bismark-methyldackel.json)

Bismark-MethylDackel pipeline in snakemake.

## Reproducible pipeline

This pipeline guarantees reproducible results, which means that it guarantees same output (DNA methylation levels) with same input (Bisulfite sequencing reads). The reproducibility is continuously being tested and the result is shown as a badge above.

## Quickstart

1. Clone the repo.

```
$ git clone https://github.com/dohlee/snakemake-bismark-methyldackel.git
$ cd snakemake-bismark-methyldackel
```

2. Modify the configurations manifest file as you want.

3. Run the pipeline.

If you already have snakemake installed, it might be useful to just use `--use-conda` option. Tweak `-j` parameter according to the number of available cores on your system.

```
$ snakemake --use-conda -p -j 32
```

Or you can create separate environment for this pipeline and run it.

```
$ conda env create -f environment.yaml
$ snakemake -p -j 32
```
