---
layout: page
title: Installation
order: 1
---

## Quick-start

```bash
# Install via Conda
conda install sortmerna=2.1b # Only if you want to use Sortmerna (optional dependency)
conda install phyloflash
# Check for dependencies
phyloFlash.pl -check_env
# Download and set up database in current folder (takes some time)
phyloFlash_makedb.pl --remote
```

## 1. System requirements

To use **phyloFlash** you will need a GNU/Linux system with Perl, R and Python
installed. (OS X is for the brave, we have not tested this!)

## 2. Download package

### 2.1 Download via Conda

[Conda](https://conda.io/docs/) is a package manager that will also install
dependencies that are required if you don't have them already. phyloFlash is
distributed through the [Bioconda](http://bioconda.github.io/) channel on
Conda.

According to the [Conda
documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html),
it is recommended to install all packages at the same time to avoid dependency
conflicts, and to create new environments instead of installing to the base
environment.

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
```

In some cases, `conda install` can hang on the "Solving environment" step. This
appears to be because of ambiguities in dependency specifications in packages
on different channels (see this
[issue](https://github.com/conda/conda/issues/8197) on GitHub). Setting the
`channel_priority` to `strict` asks Conda to always pick the higher-priority
channel first when installing packages. This requires conda version to be 4.6
and above.

### 2.2 Download from GitHub

If you prefer not to use Conda, or are interested in a specific version that is
not distributed there, you can download releases from the
[releases](https://github.com/HRGV/phyloFlash/releases) page on GitHub.

If you clone the repository directly off GitHub you might end up with a version
that is still under development.

```bash
# Download latest release
wget https://github.com/HRGV/phyloFlash/archive/pf3.4.tar.gz
tar -xzf pf3.4.tar.gz
```

Alternatively clone the latest development version with Git:

```bash
git clone https://github.com/HRGV/phyloFlash.git
ls phyloFlash
```

## 3. Check and install prerequisites

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
 - [Vsearch](https://github.com/torognes/vsearch)
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

## 4. Setting up the reference database

phyloFlash uses modified versions of the SILVA SSU database of small-subunit
ribosomal RNA sequences that is maintained by the [ARB SILVA
project](www.arb-silva.de).

*NOTE: The [SILVA
license](http://www.arb-silva.de/fileadmin/silva_databases/current/LICENSE.txt)
prohibits usage of the SILVA databases or parts of them within a
non-academic/commercial environment beyond a 48h test period. If you want to
use the SILVA databases with phyloFlash in a non-academic/commercial
environment please contact them at contact(at)arb-silva.de.*

The database has to be reformatted for use by phyloFlash. This is done with the
script `phyloFlash_makedb.pl`. Known contamination sequences from cloning
vectors are removed, repeat regions which can have an adverse effect on
sequence reconstruction are masked, the database is clustered at 99% and 96%
identity to speed up mapping/searching, and finally indexed for the read
mapper.

*NOTE: A .udb indexed database will be created with Vsearch if version v2.5.0+
is detected. However, the file will only be readable by the user running the
database setup script. If you wish to make it available for other users, please
change the file permissions for the .udb file accordingly.*

The final disk space required for the default SILVA SSU database is about 5 Gb.
An additional 5 Gb is required for the `.udb` indexed database for Vsearch
v2.5.0+. An additional 2.5 Gb is required for the SortMeRNA indexed database if
requested.

If you wish to use SortMeRNA in addition to or instead of BBmap for filtering
rRNA reads, pass the option `--sortmerena` to `phyloFlash_makedb.pl`. This
requires `sortmerna` and `indexdb_rna` to be in your path. At the moment only
SortMeRNA v2.1b is supported.

A full description of options for the database setup can be seen with

```bash
phyloFlash_makedb.pl --help
```

### 4.1. Downloading database automatically

To create a suitable database, just run

```bash
phyloFlash_makedb.pl --remote
```

in the directory where you unpacked phyloFlash. The script will download the
most current source databases and prepare the files required by
`phyloFlash.pl`.

*NOTE: This currently only works if you are not behind a proxy*

If you are behind a proxy and cannot download the database via the script, you
can download the current version of the SILVA database from [the SILVA
website](https://www.arb-silva.de/no_cache/download/archive/current/Exports/).
The filename should be `SILVA_XXX_SSURef_Nr99_tax_silva_trunc.fasta.gz` where
`XXX` is the current version number. You should also download the UniVec
database [from NCBI](https://www.ncbi.nlm.nih.gov/tools/vecscreen/univec/).
Then proceed with the instructions in section 4.2 below.

### 4.2. Set up database from local copy of SILVA SSU NR99

If you already have a local copy of the SILVA SSU NR99 database (in Fasta
format), and the NCBI Univec database, you can supply the paths:

```bash
phyloFlash_makedb.pl --univec_file /path/to/Univec --silva_file /path/to/SILVA_128_SSURef_Nr99_tax_silva_trunc.fasta.gz
```
By default, `phyloFlash.pl` will look in the folder where it is installed for
the subfolder with the highest SILVA version number. You can change this by
passing the `-dbhome <dir>` switch to phyloFlash.pl or by modifying the
`DBHOME` variable in `phyloFlash.pl`.

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
