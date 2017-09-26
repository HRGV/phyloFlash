---
layout: home
---

<img src="phyloFlash_logo.png" alt="phyloFlash logo" style="width: 200px; display: block; margin-left: auto; margin-right:auto;"/>

**phyloFlash** is a pipeline to rapidly reconstruct the SSU rRNAs and explore phylogenetic composition of an Illumina (meta)genomic or transcriptomic dataset.

***NOTE*** Version 3.0 changes some input options and also how mapping-based taxa (NTUs) are handled. Please download the last release of v2.0 ([tar.gz archive](https://github.com/HRGV/phyloFlash/archive/v2.0-beta6.tar.gz)) for the old implementation. No changes have been made to the database setup, so databases prepared for v2.0 can still be used for v3.0.

This manual explains how to install and use phyloFlash. Navigate from the menu bar above or the table of contents below.

## What does phyloFlash do?

 - Summarize taxonomic diversity of a metagenome/transcriptome library from SSU rRNA read affiliations
 - Assemble/reconstruct full-length SSU rRNA sequences suitable for phylogenetic analysis
 - Quick comparison of multiple samples by their taxonomic composition using a heatmap

## Contents

 - [Installation](install.md)
 - [Usage](usage.md)
 - [Comparing multiple samples](multiple-samples.md)
 - [Output report](output.md)
 - [FAQ](FAQ.md)

## Quick-start

```bash
# Download phyloFlash
wget https://github.com/HRGV/phyloFlash/archive/vXXX.tar.gz  
tar -xzf vXXX.tar.gz

# Check for dependencies
cd phyloFlash-XXX
./phyloFlash.pl -check_env

# Install reference database
./phyloFlash_makedb.pl --remote

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
```

Use the `-help` option to display a brief help and the `-man` option to display the full help message.

Use the `-sc` switch for MDA datasets (single cell) or other hard to assemble read sets.

Use the `-zip` switch to compress output files into tar.gz archive, and `-log` to save run messages to a log file

Example phyloFlash report from the provided test data can be viewed [here](test.phyloFlash.html).

## About

phyloFlash is written by Harald Gruber-Vodicka ([Google Scholar](https://scholar.google.de/citations?user=imYEnqMAAAAJ&hl=en&oi=ao), [GitHub](https://github.com/HRGV)), Elmar A. Pruesse ([Google Scholar](https://scholar.google.de/citations?user=F-yGwRIAAAAJ&hl=en&oi=ao), [GitHub](https://github.com/epruesse)), and Brandon Seah ([Google Scholar](https://scholar.google.de/citations?user=3l8G5BwAAAAJ&hl=en&oi=ao), [GitHub](https://github.com/kbseah))

You can find the source code for phyloFlash at GitHub:
[HRGV/phyloFlash](https://github.com/HRGV/phyloFlash/)

[Max Planck Institute for Marine Microbiology](http://www.mpi-bremen.de/)

## Citation

If you use phyloFlash for a publication, please cite as

Gruber-Vodicka HR, Pruesse E, Seah B. 2017. phyloFlash. Online:https://github.com/HRGV/phyloFlash

and also remember to cite the dependencies used.
