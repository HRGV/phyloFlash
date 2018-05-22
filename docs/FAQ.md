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

The **phyloFlash** pipeline provides the following output:

An overview of the taxonomic composition in three ways:
 - Taxonomic affiliation of reference database sequences which have hits
 - Full-length SSU rRNA sequences assembled by SPAdes, and their closest reference database hits
 - Full-length SSU rRNA sequences reconstructed by EMIRGE, and their closest reference database hits (if `-emirge` option is used)

The extracted SSU rRNA reads will be mapped back to the full-length sequences to determine the proportion that was successfully assembled or reconstructed. A low ratio of assembled:unassembled could be caused by high taxonomic diversity in a sample, or other problems with assembly (e.g. several closely-related taxa that make assembly difficult).

You'll also get output parameters that will help you evaluate quality of your sequences (see below)

# How does phyloFlash show library quality?

 - Insert size histogram - should be more or less unimodal. Multiple peaks may be caused by contamination from other sequencing libraries.
 - Identify possible contamination, e.g. *Propionibacterium* is a common contaminant from human skin
 - Coverage on SSU rRNA HMM models (18S and 16S rRNA) - should be relatively even across the gene. High coverage or gaps in coverage may be a sign of contamination from amplicon libraries.
 - Mapping identity to database sequences - identities below 90% suggest the presence of novel or divergent taxa that are not well-represented in the database.

# Why use the SSU rRNA gene, and not single-copy markers / whole genomes / kmers?

Many other tools use genomic reference data to summarize the taxonomic composition of microbial metagenomes. phyloFlash uses only the SSU rRNA marker because it is the most well-represented and well-curated phylogenetic marker gene. Tools using genome-based reference data are limited by the underrepresentation of uncultivated, environmental taxa. These taxa are also underrepresented in marker gene databases, but the latter have had several decade's head start in accumulating comparative data. Using the SSU rRNA gene also allows us to take advantage of the curated taxonomy from the [SILVA](https://www.arb-silva.de/) project, whereas the taxonomic annotation of genome sequences deposited in public databases is typically user-reported and not always consistent.

# How can I use the results?

You can browse the results as a text file or HTML formatted report (with embedded graphics) to look at the results from one sample.

To compare multiple samples you can [generate a heatmap](multiple-samples.md) of taxa vs. samples, which is especially useful with large numbers of samples.

Assembled and/or reconstructed full-length rRNA seqeuences can be used for downstream analyses, for phylogenetic analysis.

Reads which map to SSU rRNA are extracted into Fastq formatted read files by phyloFlash. These can also be used for further downstream analyses, e.g. assembly to full-length sequences by other software, or mapping to other reference databases.

# Should I run SPAdes and EMIRGE, or skip them?

phyloFlash incorporates two different tools for getting full-length rRNA sequences from unassembled metagenomic/transcriptomic reads - one is the SPAdes assembler (on by default), the other is EMIRGE (off by default) which uses an expectation maximization method to reconstruct rRNA sequences.

If you don't need full-length sequences, you may want to turn off SPAdes (`-skip_spades`) so that phyloFlash runs faster. In that case, only the initial mapping to SILVA and the taxonomic summary will be produced. However full-length sequences are useful to have because NTUs alone may give a misleadingly high impression of taxonomic diversity.

EMIRGE in particular can be tricky to install because of its own software dependencies, and you when you first try the pipeline you may want to skip EMIRGE if you don't already have it installed on your system.

# How is the taxonomic summary produced?

The taxonomic affiliation of the SSU rRNA reads is taken from the reference sequence to which it maps in the SILVA database. From v3.2b1 onwards, the taxonomy reported is the last-common-ancestor consensus of the top hits, i.e. if a read has more than one best-scoring hit, it will report the lowest taxonomic level which they have in common. In previous versions, only the single best hit was retained, but this would result in overly-preecise taxonomic assignments especially for divergent sequences. To replicate the old behavior, use the `tophit` option.

# How do I get help if something doesn't work?

First check if you have all the dependencies installed properly with `phyloFlash.pl -check_env`.

See the built-in help messages and manuals for any script in the phyloFlash folder by using either the `--help` or `--man` options.

To get in touch, please use the [phyloFlash Google group](https://groups.google.com/forum/#!forum/phyloflash), or the [issues tracker](https://github.com/HRGV/phyloFlash/issues) on GitHub to submit bug reports, and supply us with information that will help us diagnose the problem. We'll try to help if we can but please understand that phyloFlash is for academic, non-commercial use and we may not always be able to respond promptly.
