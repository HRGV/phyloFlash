---
layout: home
---

<img src="phyloFlash_logo.png" alt="phyloFlash logo" style="width: 200px; display: block; margin-left: auto; margin-right:auto;"/>

**phyloFlash** is a pipeline to rapidly reconstruct the SSU rRNAs and explore
phylogenetic composition of an Illumina (meta)genomic or transcriptomic
dataset.

***NOTE*** Version 3 changes some input options and also how mapping-based taxa
(NTUs) are handled. Please download the last release of v2.0 ([tar.gz
archive](https://github.com/HRGV/phyloFlash/archive/v2.0-beta6.tar.gz)) for the
old implementation. No changes have been made to the database setup, so
databases prepared for v2.0 can still be used for v3.0.

This manual explains how to install and use phyloFlash. Navigate from the menu
bar above or the table of contents below.

## What does phyloFlash do?

 - Summarize taxonomic diversity of a metagenome/transcriptome library from SSU
   rRNA read affiliations
 - Assemble/reconstruct full-length SSU rRNA sequences suitable for
   phylogenetic analysis
 - Quick comparison of multiple samples by their taxonomic composition using a
   heatmap

You may read more about the pipeline design and application in our
[paper](https://doi.org/10.1128/mSystems.00920-20).

## Quick-start

### Download via Conda

[Conda](https://conda.io/docs/) is a package manager that will also install
dependencies that are required if you don't have them already. phyloFlash is
distributed through the [Bioconda](http://bioconda.github.io/) channel on
Conda.

According to the [Conda
documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html),
it is recommended to install all packages at the same time to avoid dependency
conflicts, and to create new environments instead of installing to the base
environment.

We also recommend using [Mamba](https://mamba.readthedocs.io/en/latest/) as a
drop-in substitute for Conda. It implements a more effective dependency solver
and is also the default Conda frontend for the pipeline managers Snakemake.
Conda sometimes fails to solve the environment, and in these cases Mamba
usually works.

```bash
# If you haven't set up Bioconda already
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
# Try the following step if "solving environment" does not terminate
conda config --set channel_priority strict
# Create new environment named "pf" with phyloflash
# sortmerna is an optional dependency
conda create -n pf phyloflash sortmerna=2.1b
# If Conda is unable to solve the environment; requires mamba in base env
mamba create -n pf phyloflash sortmerna=2.1b
```

### Download from GitHub

If you prefer not to use Conda, or are interested in a specific version that is
not distributed there, you can download releases from the
[releases](https://github.com/HRGV/phyloFlash/releases) page on GitHub.

If you clone the repository directly off GitHub you might end up with a version
that is still under development.

```bash
# Download latest release
wget https://github.com/HRGV/phyloFlash/archive/pf3.4.tar.gz
tar -xzf pf3.4.tar.gz

# Check for dependencies and install them if necessary
cd phyloFlash-pf3.4
./phyloFlash.pl -check_env
```

### Set up database and run

This assumes that the phyloFlash scripts are already in your path.

```bash
# Install reference database (takes some time)
phyloFlash_makedb.pl --remote

# Run with test data and 16 processors (default is to use all processors available)
phyloFlash.pl -lib TEST -CPUs 16 -read1 test_files/test_F.fq.gz -read2 test_files/test_R.fq.gz

# Run with interleaved reads
phyloFlash.pl -lib LIB -read1 reads_FR.fq.gz -interleaved

# Additionally run EMIRGE for 16S rRNA sequence reconstruction
phyloFlash.pl -lib LIB -emirge -read1 reads_F.fq.gz -read2 reads_R.fq.gz

# Compress output into tar.gz archive and write a log file
phyloFlash.pl -lib LIB -zip -log -read1 reads_F.fq.gz -read2 reads_R.fq.gz

# Run both SPAdes and EMIRGE and produce all optional outputs
phyloFlash.pl -lib LIB -everything -read1 reads_F.fq.gz -read2 reads_R.fq.gz

# Run SPAdes (skip EMIRGE) and produce all optional outputs (recommended)
phyloFlash.pl -lib LIB -almosteverything -read1 reads_F.fq.gz -read2 reads_R.fq.gz

# Supply trusted contigs containing SSU rRNA sequences to screen vs reads
phyloFlash.pl -lib LIB -read1 reads_F.fq.gz -read2 reads_R.fq.gz -trusted contigs.fasta

# Use SortMeRNA instead of BBmap for initial mapping (slower, but more sensitive)
phyloFlash.pl -lib LIB -read1 reads_F.fq.gz -read2 reads_R.fq.gz -sortmerna
```

Use the `-help` option to display a brief help and the `-man` option to display
the full help message.

Use the `-sc` switch for MDA datasets (single cell) or other hard to assemble
read sets.

Use the `-zip` switch to compress output files into tar.gz archive, and `-log`
to save run messages to a log file

Example phyloFlash report from the provided test data can be viewed
[here](test.phyloFlash.html).

## Contents

 - [Installation](install.html)
 - [Usage](usage.html)
 - [Utilities](utilities.html)
 - [Output report](output.html)
 - [FAQ](FAQ.html)

## About

phyloFlash is written by Harald Gruber-Vodicka ([Google
Scholar](https://scholar.google.de/citations?user=imYEnqMAAAAJ&hl=en&oi=ao),
[GitHub](https://github.com/HRGV)), Elmar A. Pruesse ([Google
Scholar](https://scholar.google.de/citations?user=F-yGwRIAAAAJ&hl=en&oi=ao),
[GitHub](https://github.com/epruesse)), and Brandon Seah ([Google
Scholar](https://scholar.google.de/citations?user=3l8G5BwAAAAJ&hl=en&oi=ao),
[GitHub](https://github.com/kbseah))

You can find the source code for phyloFlash at GitHub:
[HRGV/phyloFlash](https://github.com/HRGV/phyloFlash/)

[Max Planck Institute for Marine Microbiology](http://www.mpi-bremen.de/)

## Citation

If you use phyloFlash for a publication, please cite our paper in _mSystems_:

Harald R Gruber-Vodicka, Brandon KB Seah, Elmar Pruesse. phyloFlash: Rapid SSU
rRNA profiling and targeted assembly from metagenomes. [ *mSystems* 5 :
e00920-20](https://doi.org/10.1128/mSystems.00920-20);
doi:10.1128/mSystems.00920-20

and also remember to cite the dependencies used, which are listed in each
phyloFlash report file.
