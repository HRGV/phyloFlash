---
layout: home
---

<img src="phyloFlash_logo.png" alt="phyloFlash logo" style="width: 200px; display: block; margin-left: auto; margin-right:auto;"/>

**phyloFlash** is a pipeline to rapidly reconstruct the SSU rRNAs and explore
phylogenetic composition of an Illumina (meta)genomic or transcriptomic
dataset.

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

We recommend installing phyloFlash and its dependencies using Conda or Mamba.
[Conda](https://conda.io/docs/) is a package manager that will also install
dependencies that are required if you don't have them already.

phyloFlash is distributed through the [Bioconda](http://bioconda.github.io/)
channel on Conda.

According to the [Conda documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html),
avoid installing new packages to your base environment but create new
environments for them as required. Also, specify all desired packages at the
same time when creating a new environment, instead of adding them sequentially,
to avoid dependency conflicts.

We also suggest using [Mamba](https://mamba.readthedocs.io/en/latest/) as a
drop-in substitute for Conda. It implements a more effective dependency solver
and is also the default Conda frontend for the pipeline manager Snakemake.
Simply replace `conda` with `mamba` in the commands below. Note that the
`defaults` channel should be enabled.

```bash
# If you haven't set up Bioconda already
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
# Try the following step if "solving environment" does not terminate
conda config --set channel_priority strict
# Create new environment named "pf" with phyloflash
conda create -n pf phyloflash
# Activate environment
conda activate pf
# Check that dependencies all installed properly
phyloFlash.pl -check_env
```


### Download pre-formatted database

Pre-formatted databases derived from SILVA releases 138 onwards are available
from the following Zenodo archives:

 * [SILVA 138.1](https://doi.org/10.5281/zenodo.7892521) (latest)
 * [SILVA 138](https://doi.org/10.5281/zenodo.7890453)

Download, checksum, and unpack:

```bash
wget https://zenodo.org/record/7892522/files/138.1.tar.gz # 5.5 GB download
tar -xzf 138.1.tar.gz # unpacks folder 138.1/ in the current location
```

Specify path to the database folder with the option `-dbhome` when running
phyloFlash (see below).

Older versions of the SILVA database have a more restrictive license, so we are
unable to distribute pre-formatted versions. You will have to download the
original SILVA files and run the `phyloFlash_makedb.pl` script yourself (see
Manual).


### Test phyloFlash with test dataset

Test data are included with phyloFlash. The following assumes that you
installed phyloFlash to a Conda environment called `pf`, and that the database
files have been unpacked to a folder `/path/to/138.1`. By default, phyloFlash
will look for the database folder in the folder where it is installed. If it is
located somewhere else, specify this to the `-dbhome` option.

```bash
conda activate pf # If Conda environment not already activated
phyloFlash.pl -dbhome /path/to/138.1 -lib TEST -CPUs 16 \
 -read1 ${CONDA_PREFIX}/lib/phyloFlash/test_files/test_F.fq.gz \
 -read2 ${CONDA_PREFIX}/lib/phyloFlash/test_files/test_R.fq.gz \
 -almosteverything
```


### Example phyloFlash commands

```bash
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
