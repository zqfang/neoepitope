# Snakemake Workflow for Label-free proteomics
This project contains a snakemake workflow that uses the OpenMS (ver 2.3.0) toolset to obtain relativer protein intensities and spectral counts for MS/MS data.

https://bitbucket.org/snakemake/snakemake

## Authors

* Michael Freitas

## Requirements
The workflow requires preinstallation of snakemake and singularity.  The workflow uses a prebuild singularity container with that has OpenMS 2.3.0 and thirdparty applicaitons.  Consult the documentation for these tools on their installation.

## Usage

### Step 1: Install workflow
Download the workflow. If you use this workflow in a paper, please cite the URL of this repository.

### Step 2: Download the example mzML and fasta files
A script is provided to download example fasta and  mxML files.

### Step 3: Configure job parameters
Configure the workflow via the `config.yaml`

### Step 3: Execute workflow
Test workflow by performing a dry-run

    snakemake -np --snakefile Snakefile.py

Execute the workflow.  Specify the number of cores by replacing the $N with the desired number of processors to use.

    snakemake --cores $N --snakefile Snakefile.py

Consult the official [Snakemake documentation](https://snakemake.readthedocs.io) for more information.
