<img align="right" src="docs/phyloFlash_logo.png" width="200" alt="phyloFlash logo"/>

# phyloFlash v3.3b3

[![GitHub (pre-)release](https://img.shields.io/github/release/HRGV/phyloflash/all.svg?label=Latest%20Version)]()
[![Bioconda](https://img.shields.io/conda/vn/Bioconda/phyloFlash.svg)](https://bioconda.github.io/recipes/phyloflash/README.html)
[![](https://img.shields.io/conda/dn/Bioconda/phyloflash.svg)](https://anaconda.org/bioconda/phyloflash/files)

by Harald Gruber-Vodicka, Elmar A. Pruesse, and Brandon Seah.

*phyloFlash* is a pipeline to rapidly reconstruct the SSU rRNAs and explore phylogenetic composition of an Illumina (meta)genomic or transcriptomic dataset. **[Manual](https://hrgv.github.io/phyloFlash)**

***NOTE*** Version 3 changed some input options and also how mapping-based taxa (NTUs) are handled. Please download the last release of v2.0 ([tar.gz archive](https://github.com/HRGV/phyloFlash/archive/v2.0-beta6.tar.gz)) for the old implementation. No changes have been made to the database setup, so databases prepared for v2.0 can still be used for v3.0.

Read [our preprint](https://doi.org/10.1101/521922) on phyloFlash.

## Quick-start

### Download via Conda

[Conda](https://conda.io/docs/) is a package manager that will also install dependencies that are required if you don't have them already.

phyloFlash is distributed through the [Bioconda](http://bioconda.github.io/) channel on Conda.

```bash
# If you haven't set up Bioconda already
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
# Try the following step if "solving environment" does not terminate
conda config --set channel_priority strict
# Install packages to current environment
conda install sortmerna=2.1b # Optional - if you want to use SortMeRNA option
conda install phyloflash
```

### Download from GitHub

If you prefer not to use Conda, or are interested in a specific version that is not distributed there, you can download releases from the [releases](https://github.com/HRGV/phyloFlash/releases) page on GitHub.

If you clone the repository directly off GitHub you might end up with a version that is still under development.

```bash
# Download latest release
wget https://github.com/HRGV/phyloFlash/archive/pf3.3b3.tar.gz
tar -xzf pf3.3b3.tar.gz

# Check for dependencies and install them if necessary
cd phyloFlash-pf3.3b3
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

Use the `-help` option to display a brief help and the `-man` option to display the full help message.

Use the `-sc` switch for MDA datasets (single cell) or other hard to assemble read sets.

Use the `-zip` switch to compress output files into tar.gz archive, and `-log` to save run messages to a log file

## Output

phyloFlash screens metagenomic or metatranscriptomic reads for SSU rRNA sequences by mapping against the SILVA SSU Ref database.

Extracted reads are assembled and/or reconstructed into full-length sequences, and the proportion of reads assembled is estimated by re-mapping them to the full-length seqeunces.

phyloFlash reports a taxonomic summary of the reads from the initial mapping, the best database matches of the full-length sequences, and a second taxonomic summary of the unassembled reads left over.

Plain text and HTML-formatted reports are produced, reporting summary statistics from each run. The HTML report includes an interactive graphical summary.

## Going further

The phyloFlash suite also includes other tools for SSU rRNA-centric metagenome analyses. Run the commands without arguments to see help messages.

 * `ENA_phyloFlash.pl` - Automatically download read files from [ENA](https://www.ebi.ac.uk/ena) given a read accession number, and run phyloFlash on them
 * `phyloFlash_compare.pl` - Compare the taxonomic composition of multiple samples from their phyloFlash results. This produces a barplot, heatmap, or distance matrix based on the NTU abundances in two or more samples.
 * `phyloFlash_fastgFishing.pl` - Given a metagenomic assembly graph in [Fastg](http://fastg.sourceforge.net/) format, identify SSU rRNA sequences and extract contigs connected to them. Optionally compare to phyloFlash results from the same library.

## Manual

For further information **please refer to the [Manual](https://hrgv.github.io/phyloFlash)**.

## Versions and changes

* v3.3 beta 2
  * New options to graphical comparison scripts, and other small bug fixes
  * Fix bug due to change in SILVA project file naming convention with SILVA 138 onwards
* v3.3 beta 1
  * Add support for using SortMeRNA instead of BBmap for initial mapping step
  * Changes to how mapping data are hashed; process SAM file of initial mapping to fix known bugs with bitflag and read name reporting in BBmap and SortMeRNA
* v3.2 beta 1
  * Report ambiguous hits during mapping step, use consensus of top hits to assign taxonomy instead of single best hit
  * Add utility `phyloFlash_compare.pl` to compare taxonomic composition of multiple libraries from phyloFlash output
  * Add utility `phyloFlash_fastgFishing.pl` to extract genome bins from Fastg files
* v3.1 beta 2
  * Fix bug in Fasta headers with changed output from Bedtools v2.26+
  * Make bbmap and reformat.sh overwrite existing output files of same name
  * Rearrange elements in HTML report file
* v3.1 beta 1
  * Allow user to supply "trusted contigs" of sequence assemblies containing SSU rRNA which will also be screened against the read libraries
  * Fix bugs in tree plotting
* v3.0 beta 1
  * Re-map extracted SSU reads onto assembled sequences to check proportion assembled
  * Revamp of HTML report output. Embed interactive graphical summary, use SVG-formatted graphics, remove dependency on R packages for report graphics.
  * Changes to how mapping-based NTUs are calculated. Now count all reads (not only unambiguously-mapped) and count segments of read pairs separately.
  * No change to heatmap script for comparing multiple samples
* v2.0 complete rewrite

## Contact

Please report any problems to the [phyloFlash Google group](https://groups.google.com/forum/#!forum/phyloflash) or with the GitHub issue tracker.

(Pull requests with suggested fixes are of course also always welcome.)

We also welcome any feedback on the software and its documentation, especially suggestions for improvement!

## Acknowledgements

We thank colleagues and phyloFlash users who have contributed to phyloFlash development by testing the software, reporting bugs, and suggesting new features.

## Citation

If you use phyloFlash for a publication, please cite our preprint on BioRxiv:

Harald R Gruber-Vodicka, Brandon KB Seah, Elmar Pruesse. phyloFlash — Rapid SSU rRNA profiling and targeted assembly from metagenomes. [bioRxiv 521922](https://doi.org/10.1101/521922); doi: https://doi.org/10.1101/521922

and also remember to cite the dependencies used, which are listed in each phyloFlash report file.
