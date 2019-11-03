# next-gwas

The main objective of this project is to develop a parallel distributed computational pipeline that uses maxT permutation testing to determine which genomic variants are implicated in a disease.

---
## Contents
- Pre-requisites
- Installation
- Usage
---
## Pre-requisites

Ensure that you have both [nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation) and [plink](http://zzz.bwh.harvard.edu/plink/download.shtml ) installed on your machine by following the instructions provided in the links. We highly recommend that you install these two tools in a directory accesible by your by your `$PATH` variable.

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








