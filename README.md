# next-gwas

The main objective of this project is to develop a parallel distributed computational pipeline that uses maxT permutation testing to determine which genomic variants are implicated in a disease.

---
## Contents
- Pre-requisites
- Installation
- Usage
---
## Pre-requisites

To get started you needed the following tools on your system:

- [nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation)
- [plink](http://zzz.bwh.harvard.edu/plink/download.shtml ) 

We particularly recommend installing them in a directory accesible by your by your `$PATH` variable i.e, `$HOME/bin`

## Installation

To install the next-gwas project, clone this git repository

```bash
git clone https://github.com/tino-sibanda/next-gwas.git
cd next-gwas/
```

## Usage

To run the pipeline on the example datasets , use the command  below:
```bash
nextflow run pipeline.nf
```








