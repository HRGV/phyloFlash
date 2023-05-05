---
layout: page
title: Installation
order: 1
---

## 1. System requirements

To use **phyloFlash** you will need a GNU/Linux system with Perl, R and Python
installed. (OS X is for the brave, we have not tested this!)


## 2. Download package

### 2.1 Download via Conda

We recommend installing phyloFlash and its dependencies using Conda or Mamba.
[Conda](https://conda.io/docs/) is a package manager that will also install
dependencies that are required if you don't have them already.

phyloFlash is distributed through the [Bioconda](http://bioconda.github.io/)
channel on Conda.

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

 * Avoid installing new packages to your base environment. Instead, create new
   environments with required packages as you need them.
 * Install packages to a new environment simultaneously, instead of adding them
   sequentially. This will prevent dependency conflicts.
 * In some cases, `conda install` can hang on the "Solving environment" step.
   This appears to be because of ambiguities in dependency specifications in
   packages on different channels (see this
   [issue](https://github.com/conda/conda/issues/8197) on GitHub). Setting the
   `channel_priority` to `strict` asks Conda to always pick the higher-priority
   channel first when installing packages. This requires conda version to be
   4.6 and above.
 * We also suggest using [Mamba](https://mamba.readthedocs.io/en/latest/) as a
   drop-in substitute for Conda. It implements a more effective dependency
   solver and is also the default Conda frontend for the pipeline manager
   Snakemake. Simply replace `conda` with `mamba` in the commands. Note that
   the `defaults` channel should be enabled.
 * If you wish to use Sortmerna (optional) for extracting rRNA reads, specify
   version 2.1b: `conda create -n pf_sortmerna phyloflash sortmerna=2.1b`


### 2.2 Download from GitHub

If you wish to modify the source code, you can clone the repository from GitHub

```bash
git clone https://github.com/HRGV/phyloFlash.git
cd phyloFlash
git status
```


## 3. Check and install dependencies

Check that dependencies are available:

```bash
phyloFlash.pl -check_env
```

If you downloaded via Conda they should already be installed, otherwise you
will need to do it yourself.

phyloFlash relies on the following software:

 - [Perl >= 5.13.2](http://www.perl.org/get.html)
 - [EMIRGE](https://github.com/csmiller/EMIRGE) and its dependencies
 - [BBmap](http://sourceforge.net/projects/bbmap/)
 - [Vsearch >=2.5.0](https://github.com/torognes/vsearch)
 - [SPAdes](http://bioinf.spbau.ru/spades)
 - [Bedtools](https://github.com/arq5x/bedtools2)
 - [Mafft](http://mafft.cbrc.jp/alignment/software/)
 - [Barrnap](https://github.com/tseemann/barrnap) (customized version is provided with phyloFlash)
 - Optional: [SortMeRNA](https://github.com/biocore/sortmerna) v2.1b, if you want to use it as alternative to BBmap

These tools need to be in your `$PATH` environment variable, so that phyloFlash
can find them.

In addition, you will need [R](https://www.r-project.org/) and the following R
packages for plotting if you use the `phyloFlash_compare.pl` script for
comparing multiple samples:

 - ggdendro
 - gtable
 - reshape2
 - ggplot2
 - optparse

Within R, run the command

```R
install.packages(c("ggdendro","gtable","reshape2","ggplot2","optparse"))
```

## 4. Set up the reference database

phyloFlash uses modified versions of the SILVA SSU database of small-subunit
ribosomal RNA sequences that is maintained by the [ARB SILVA
project](www.arb-silva.de).


### 4.1. Download pre-formatted database

Pre-formatted databases derived from SILVA releases 138 onwards are available
from the following Zenodo archives:

 * [SILVA 138.1](https://doi.org/10.5281/zenodo.7892521) (latest)
 * [SILVA 138](https://doi.org/10.5281/zenodo.7890453)

NOTE: Prebuilt databases are not provided for SILVA versions before 138,
because these are released under different license(s) that prohibit usage of
the SILVA databases or parts of them within a non-academic/commercial
environment beyond a 48 h test period. SILVA version 138 onwards is released
under a more permissive Creative Commons Attribution 4.0 license.

Download, checksum, and unpack (example for release 138.1):

```bash
wget https://zenodo.org/record/7892522/files/138.1.tar.gz # 5.5 GB download
tar -xzf 138.1.tar.gz # unpacks folder 138.1/ in the current location
```

Specify path to the database folder with the option `-dbhome` when running
phyloFlash (see below).


### 4.2. Format database locally

If you wish to use earlier versions of the SILVA database, or a custom database
file, you will have to format and index them. This is done with the script
`phyloFlash_makedb.pl`. Known contamination sequences from cloning vectors are
removed, repeat regions which can have an adverse effect on sequence
reconstruction are masked, the database is clustered at 99% and 96% identity to
speed up mapping/searching, and finally indexed for the read mapper.

A full description of options for the database setup can be seen with

```bash
phyloFlash_makedb.pl --help
```

Download the desired version of the SILVA SSURef NR99 database from [the SILVA
website](https://www.arb-silva.de/download/archive/) (in Fastsa format) under the `Exports` subfolder of the respective release. The filename should be `SILVA_XXX_SSURef_Nr99_tax_silva_trunc.fasta.gz` where
`XXX` is the version number. Links to the last five releases:
 * [138.1](https://www.arb-silva.de/fileadmin/silva_databases/release_138.1/Exports/SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta.gz)
 * [138](https://www.arb-silva.de/fileadmin/silva_databases/release_138/Exports/SILVA_138_SSURef_NR99_tax_silva_trunc.fasta.gz)
 * [132](https://www.arb-silva.de/fileadmin/silva_databases/release_132/Exports/SILVA_132_SSURef_Nr99_tax_silva_trunc.fasta.gz)
 * [128](https://www.arb-silva.de/fileadmin/silva_databases/release_128/Exports/SILVA_128_SSURef_Nr99_tax_silva_trunc.fasta.gz)
 * [123.1](https://www.arb-silva.de/fileadmin/silva_databases/release_123.1/Exports/SILVA_123.1_SSURef_Nr99_tax_silva_trunc.fasta.gz)

Also download the UniVec database [from NCBI](https://www.ncbi.nlm.nih.gov/tools/vecscreen/univec/).

Specify the paths to the SILVA and UniVec files wtih the `--silva_file` and `--univec_file` options respectively to build the database locally, example below.

```bash
phyloFlash_makedb.pl --univec_file /path/to/Univec --silva_file /path/to/SILVA_128_SSURef_Nr99_tax_silva_trunc.fasta.gz
# Creates a new folder ./128
```

 * A new folder containing the database files will be created. The folder name
   will correspond to the SILVA release number and is parsed from the input
   file name (which should follow the SILVA file naming convention exactly).
 * The `--remote` option is no longer supported.
 * If you wish to use SortMeRNA in addition to or instead of BBmap for
   filtering rRNA reads, pass the option `--sortmerena` to
   `phyloFlash_makedb.pl`. This requires `sortmerna` and `indexdb_rna` to be in
   your path. At the moment only SortMeRNA v2.1b is supported.
 * When you run the main `phyloFlash.pl` script, it will by default look in the
   folder where it is installed for the subfolder with the highest SILVA
   version number. You can change this by specifying the path with the
   `-dbhome` option in `phyloFlash.pl`.


### 4.3. Set up a custom database with your own sequences

Users can supply their own databases of SSU rRNA sequences, or even other
genes, in place of the SILVA SSU database, as long as they are formatted in the
following way:

 - Sequences should be in Fasta format
 - Fasta headers should have the format `{IDENTIFIER}.{INTEGER}.{INTEGER}
   {TAXONOMY-STRING}` where:
   - `IDENTIFIER` is a unique sequence identifier which does not have spaces or
     periods
   - The difference between the two `INTEGER`s should be the length of the
     sequence, e.g. 1.1700 for a 1700 bp sequence
   - `TAXONOMY-STRING` is in SILVA or NCBI format, delimited by semicolons with
     no spaces (but spaces in taxon names allowed)
   - There is a single space before the `TAXONOMY-STRING`
 - The name of the Fasta file should begin with `SILVA_{DBNAME}_` where
   `DBNAME` is the name of the database (e.g. `CustomDB`), and will also be the
   name of the output folder containing the formatted database files. This is
   because the Fasta filename is parsed by the script.

The database setup script automatically trims cloning vectors and other
potential contaminants, and discards sequences shorter than 800 bp. If your
custom database contains a gene of interest that is a different average length,
you can change the minimum sequence length with the `--ref_minlength`
parameter.
