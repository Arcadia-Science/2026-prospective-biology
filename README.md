# 2026-prospective-biology

[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/projects/miniconda/en/latest/)

## Purpose

Code associated with the pub "Biology needs to become prospective".

## Installation and Setup

This repository uses conda to manage software environments and installations. You can find operating system-specific instructions for installing miniconda [here](https://docs.conda.io/projects/miniconda/en/latest/). After installing conda and [mamba](https://mamba.readthedocs.io/en/latest/), run the following command to create the pipeline run environment.

```{bash}
TODO: Replace <NAME> with the name of your environment
mamba env create -n <NAME> --file envs/dev.yml
conda activate <NAME>
```

## Data

All associated data are available on [Zenodo](https://zenodo.org/records/18474710). This includes AFDB clusters and metadata from "Clustering predicted structures at the scale of the known protein universe" (Hernandez et al. 2023), a time-calibrated phylogeny from TimeTree (Kumar et al. 2022), taxonomic hierarchies accessed from the NCBI Taxonomy database, and genome/proteome statistics from the NCBI Genome resource. 

## Overview

### Description of the folder structure
```
├── code
│   ├── analysis
│   │   ├── figure_1.R
│   │   ├── figures_2_3.R
│   │   └── figures_4_5.R
│   └── utils
│       ├── process_data.R
│       └── utils.R
```
### Methods

1. Download data from Zenodo and preprocess for analyzes using `process_data.R`.
2. Run analyses and generate plots associated with Figure 1 using `Figure_1.R`.
3. Run analyses and generate plots associated with Figures 2 and 3 using `Figures_2_3.R`.
4. Run analyses and generate plots associated with Figures 4 and 5 using `Figures_4_5.R`.

### Compute Specifications

All analyses were done on an Apple MacBook Pro running macOS Montery with 32GB RAM, 10 cores, and 1TB of storage.

## Contributing

See how we recognize [feedback and contributions to our code](https://github.com/Arcadia-Science/arcadia-software-handbook/blob/main/guides-and-standards/guide-credit-for-contributions.md).

---
## For Developers

This section contains information for developers who are working off of this template. Please adjust or edit this section as appropriate when you're ready to share your repo.

### GitHub templates
This template uses GitHub templates to provide checklists when making new pull requests. These templates are stored in the [.github/](./.github/) directory.

### `.gitignore`
This template uses a `.gitignore` file to prevent certain files from being committed to the repository.

### Linting
This template automates linting using GitHub Actions and the [`lintr` linter](https://cran.r-project.org/web/packages/lintr/vignettes/lintr.html). When you push changes to your repository, GitHub will automatically run the linter and report any errors, blocking merges until they are resolved. 
