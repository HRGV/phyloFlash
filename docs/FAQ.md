---
layout: page
title: FAQ
permalink: /faq/
order: 4
---

# What is phyloFlash for?

Quick overview of taxonomic composition and library quality in metagenomic and metatranscriptomic samples, including multiple sample comparison.

# What input data can be analyzed?

Metagenomic and metatranscriptomic short-read sequencing, e.g. Illumina. Raw reads are okay, without quality trimming or filtering.

Transcriptomes typically have a much higher proportion of rRNA reads than genomes. Limit the number of mapped reads to < 1 million with `-readlimit 1000000`.

**NOT** suitable for:
 - Sanger or PacBio libraries
 - Amplicon libraries - use a dedicated tool like mothur

# What output will I get?

The default **phyloFlash** pipeline provides the following output:

An overview of the taxonomic composition in three ways:
 - Taxonomic affiliation of reference DB sequences which have hits
 - Full-length SSU rRNA sequences assembled by SPAdes, and their closest reference DB hits
 - Full-length SSU rRNA sequences reconstructed by EMIRGE, and their closest reference DB hits

You'll also get output parameters that will help you evaluate quality of your sequences (see below)

# How does phyloFlash show library quality?

 - Insert size histogram - should be more or less unimodal
 - Identify possible contamination, e.g. *Propionibacterium* is a common contaminant from human skin

# How can I use the results?

You can browse the results as a text file or HTML formatted report (use `-html` flag), with graphics, to look at the results from one sample.

To compare multiple samples you can [generate a heatmap](multiple-samples.md) of taxa vs. samples, which is especially useful with large numbers of samples.

Assembled and/or reconstructed full-length rRNA seqeuences can be used for downstream analyses, for phylogenetic analysis.

Reads which map to SSU rRNA are extracted into Fastq formatted read files by phyloFlash. These can also be used for further downstream analyses, e.g. assembly to full-length sequences by other software, or mapping to other reference databases.

# Should I run SPAdes and EMIRGE, or skip them?

phyloFlash incorporates two different tools for getting full-length rRNA sequences from unassembled metagenomic/transcriptomic reads - one is the SPAdes assembler, the other is EMIRGE which uses an expectation maximization method to reconstruct rRNA sequences.

If you don't need full-length sequences, you may want to turn them off `-skip_spades` and `-skip_emirge` so that phyloFlash runs faster. However full-length sequences are useful to have because NTUs alone may give a misleadingly high impression of taxonomic diversity.

EMIRGE in particular can be tricky to install because of its own software dependencies, and you when you first try the pipeline you may want to skip EMIRGE if you don't already have it installed on your system.

# How do I get help if something doesn't work?

First check if you have all the dependencies installed properly with `phyloFlash.pl -check_env`.

To get in touch, please use the [phyloFlash Google group](https://groups.google.com/forum/#!forum/phyloflash), or the [issues tracker](https://github.com/HRGV/phyloFlash/issues) on GitHub to submit bug reports, and supply us with information that will help us diagnose the problem. We'll try to help if we can but please understand that phyloFlash is for academic, non-commercial use and we may not always be able to respond promptly.
