---
layout: page
title: Utilities
order: 4
---

In addition to the main `phyloFlash.pl` tool, the phyloFlash folder contains a number of other tools for working with phyloFlash and analyzing the results.

To see usage instructions for any phyloFlash tool, simply run the script without arguments (short synopsis), or with the options `--help` (short help) or `--man` (full detailed manual page). 

## 1. Database setup

```bash
# download databases via FTP
phyloFlash_makedb.pl --remote

# use local copies
phyloFlash makedb.pl --silva_file path/to/silva_db --univec_file path/to/univec_db
```

Set up phyloFlash databases. See [Installation](install.html).

## 2. Compare samples

```bash
phyloFlash_compare.pl --zip LIB1.phyloFlash.tar.gz,LIB2.phyloFlash.tar.gz --task barplot

phyloFlash_compare.pl --zip LIB1.phyloFlash.tar.gz,LIB2.phyloFlash.tar.gz --task heatmap

# Do both heatmap and barplot, and output in PNG instead of PDF format
phyloFlash_compare.pl --zip LIB1.phyloFlash.tar.gz,LIB2.phyloFlash.tar.gz --task heatmap,barplot --outfmt png

# Find all *.phyloFlash.tar.gz archives in the current folder
phyloFlash_compare.pl --allzip --task barplot,heatmap
```

Compare the taxonomic composition multiple metagenomic/transcriptomic libraries using the phyloFlash NTU abundance results. Three types of comparison are available: heatmap (taxa vs samples), barplot (relative taxon abundance by sample), or distance matrix (Unifrac-like abundance-weighted taxonomic distances).

The phyloFlash.pl pipeline rapidly screens metagenomic/transcriptomic libraries for SSU rRNA reads by mapping against the SILVA SSU Ref NR database. The top reference hits per read are used to report an approximate taxonomic affiliation. Read counts per NTU (nearest taxonomic units) are reported by the pipeline to give an overview of the taxonomic diversity in the sample.

These NTU abundances are the basis for the comparison tools in phyloFlash_compare.pl. Users should be aware that taxonomic affiliations are only approximate and are probably inaccurate at lower taxonomic levels, and that taxonomic groups are not necessarily monophyletic. For more accurate (but slower) taxonomic or phylogenetic classifications, one should reanalyze the extracted reads with a dedicated method, e.g. using a phylogenetic placement algorithm on a curated reference tree.

This script is a convenient wrapper for the R scripts `phyloFlash_heatmap.R` ([more info](multiple-samples.html)) and phyloFlash_barplot.R. More options are available when running those scripts separately (see help messages by running the commands without arguments), but it is recommended to use `phyloFlash_compare.pl` if you are not familiar with the internal workings.

### 2.1 R package dependencies

The following R packages are required: "optparse", "ggplot2", "reshape2", "ggdendro", "gtable", "plyr", "RColorBrewer", "methods", "grid".

### 2.2 Input files

This script uses the taxonomic summaries from the phyloFlash runs to perform the comparison. This can be either in the form of the summary table (LIBNAME.phyloFlash.NTUabundance.csv or LIBNAME.phyloFlash.NTUfull_abundance.csv), or can be recalculated from the original SAM mapping files in a phyloFlash.tar.gz archive. The latter options is given because the way that the taxonomic summary is made was changed between phyloFlash versions, so older results can still be used.

To use the CSV files, pass the filenames as a comma-separated list to the option `--csv`, like so: `--csv LIB1.phyloFlash.NTUabundance.csv,LIB2.phyloFlash.NTUabundance.csv`.

To use the tar.gz archives, pass the filenames as a comma-separated list to the option `--zip`, like so: `--zip LIB1.phyloFlash.tar.gz,LIB2.phyloFlash.tar.gz`.

To compare all phyloFlash.tar.gz archives in the same folder, use the option `--allzip`. Use `use_SAM` to recalculate the taxonomic summary from the SAM files in the archives.

### 2.3 Analysis options

Three types of analyses can be performed, and are chosen with the option `--task`.

The taxonomic composition for each sample can be summarized as a `barplot`, displayed side-by-side in the order given.

The option `heatmap` produces a plot where both samples and taxa are clustered by similarity, and abundances are represented by a color scale.

To produce only a plain-text, tab-separated pairwise distance matrix of samples, use option `matrix`. The distances used are abundance-weighted taxonomic Unifrac-like distances (see below). More than one option can be supplied at the same time, separated with commas, e.g. `--task barplot,heatmap`.

Further options can be viewed with the `--help` or `--man` options.

### 2.3 Abundance-weighted taxonomic Unifrac-like distance

The default options for `heatmap` cluster samples by their taxonomic composition, but treat each taxon independently. However, two genera from the same family more similar than two genera from different families. The abundance-weighted taxonomic Unifrac-like distance is used to take this into account.

The original [weighted Unifrac](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1828774/) is a distance measure designed for comparing communities using a phylogenetic tree of their members, so that the phylogenetic distance is also taken into account, while weighting the result by their abundance in the different samples. The taxonomic Unifrac-like distance simply uses a taxonomy tree (based on the SILVA taxonomy) instead of a phylogenetic tree.

This distance measure can be used in two ways in the phyloFlash_compare.pl script. Either report only the raw distance matrix using the option `--task matrix`, or plot a heatmap with the "custom" distance metric for clustering samples: `--task heatmap --cluster-samples custom`

This feature is currently experimental and should be used with caution!

## 3. Download public data

```bash
# Download run given accession number and run phyloFlash
ENA_phyloFlash.pl --acc ACCESSION --phyloFlash

# Specify HTTP proxy, delete read files when done
ENA_phyloFlash.pl --acc ACCESSION --phyloFlash --http_proxy="http://proxy.server/" --cleanup
```

Given the run accession number, download read data from the [European Nucleotide Archive](http://www.ebi.ac.uk/ena/) and run phyloFlash on it. By default it only downloads the read Fastq files to the current folder. Run phyloFlash afterwards with the `--phyloFlash` option, using the `--almosteverything` option by default. Alternatively you can pass other options to phyloFlash like so: `--phyloFlash="--option1 --option2"`.

Proxy server users should specify the address to the `--http_proxy` option.

Note that only RUN accession numbers are valid. For an explanation of the different types of accession numbers used by ENA, see their [metadata documentation](http://ena-docs.readthedocs.io/en/latest/meta_01.html#metadata-model). Run accessions typically have the prefix "ERR", "SRR", or "DRR". 

## 4. Genome binning

```bash
# Megahit assembly
phyloFlash_fastgFishing.pl --fasta [Fasta] --fastg [Fastg] --out [PREFIX]

# SPAdes assembly
phyloFlash_fastgFishing.pl --fasta [Fasta] --fastg [Fastg] --paths [Paths] --out [PREFIX]
```

From (meta)genome assembly graph in Fastg format, predict SSU rRNA sequences, extract contig clusters with total length > cutoff (default 100 kbp), and match them to phyloFlash SSU rRNA-targeted assembly results.

This is intended to aid binning of microbial genomes from metagenomes. Each cluster in a Fastg graph is likely to originate from a single genome, and so represents a putative genome bin. 

If only a Fasta and Fastg file are specified, predict SSU rRNA sequences in the Fasta contigs and report all contig clusters containing an SSU rRNA with total cluster length above the cutoff.

If `--clusteronly` option is used, then the SSU rRNA are not predicted, and the script simply returns all contig clusters above the length cutoff.

If a Fasta file containing SSU rRNA sequences (e.g. from phyloFlash) is specified with option `--compare-ssu`, these will be searched against the Fasta contigs, and a table of which bin each sequence occurs in will also be reported. If you have a phyloFlash (v3.1+) output tar.gz archive, you can specify it with the option `--compare-zip` and the assembled SSU file will be automatically extracted for you.

### 4.1 Background

The conventional Fasta sequence format only allows sequences to be represented linearly. The [Fastg format](http://fastg.sourceforge.net/) allows additional information, such as uncertainties and non-linear assembly structures, to be represented. This can be useful for genome binning, as an assembler that resolves a genome into multiple linear contigs may nonetheless have information on how those contigs are connected to each other.

This is the basis for the Fastg-'fishing' implemented here. Input metagenomes assemblies (Fasta) are screened for SSU rRNA sequences, which are then used as 'bait' to fish out all connected contigs, using information from the corresponding Fastg file.

Unfortunately the implementation of the Fastg format is not well-standardized, and different assemblers (and versions thereof) may report it in different ways. Therefore only the output from assemblers MegaHIT and SPAdes are supported at the moment. The SPAdes output requires an additional file, called the "paths" file, which translates between 'EDGE' and 'NODE' names.

## 5. Other files

The following other files are used internally by phyloFlash or its associated utilities:

`PhyloFlash.pm` - Perl local module of shared subroutines

`phyloFlash_plotscript_svg.pl` - Generate plots for HTML report

`phyloFlash_report_template.html` - Template to build HTML report

`phyloFlash_heatmap.R` and `phyloFlash_barplot.R` - R scripts to generate the comparison plots; called by phyloFlash_compare.pl
