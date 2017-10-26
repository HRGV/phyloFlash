---
layout: page
title: Installation
permalink: /install/
order: 1
---

## Quick-start

```bash
# Download package
wget https://github.com/HRGV/phyloFlash/archive/vXXX.tar.gz  
tar -xzf vXXX.tar.gz
# cd into folder
cd phyloFlash
# Check for dependencies
./phyloFlash.pl -check_env
# Download and set up database
./phyloFlash_makedb.pl --remote
```

## 1. System requirements

To use **phyloFlash** you will need a GNU/Linux system with Perl, R and
Python installed. (OSX is for the brave, we have not tested this!)

## 2. Download package

Download the most recent release [from Github](https://github.com/HRGV/phyloFlash/releases), either from the website or at the command line with `wget`.

```bash
wget https://github.com/HRGV/phyloFlash/archive/vXXX.tar.gz  
tar -xzf vXXX.tar.gz
```

Alternatively clone the latest development version with Git:

```bash
git clone https://github.com/HRGV/phyloFlash.git
ls phyloFlash
```

## 3. Check and install prerequisites

Install the tools phyloFlash uses:

 - [Perl >= 5.13.2](http://www.perl.org/get.html)
 - [EMIRGE](https://github.com/csmiller/EMIRGE) and its dependencies
 - [BBmap](http://sourceforge.net/projects/bbmap/)
 - [Vsearch](https://github.com/torognes/vsearch)
 - [SPAdes](http://bioinf.spbau.ru/spades)
 - [Bedtools](https://github.com/arq5x/bedtools2)
 - [Mafft](http://mafft.cbrc.jp/alignment/software/)
 - [Barrnap](https://github.com/tseemann/barrnap) (currently provided with phyloFlash)

These tools need to be in your `$PATH` environment variable, so that phyloFlash can find them. To see whether all required tools are available, just use the `-check_env` flag:

```bash
cd phyloFlash-3.0
./phyloFlash.pl -check_env
```

A quick way to install dependencies is with [Bioconda](https://bioconda.github.io/), which uses the Conda package management system.

```bash
conda config --add channels bioconda

conda install emirge
conda install bbmap
conda install vsearch
conda install spades
conda install mafft
conda install bedtools
```

## 4. Setting up the reference database

phyloFlash uses modified versions of the SILVA SSU database of small-subunit ribosomal RNA sequences that is maintained by the [ARB SILVA project](www.arb-silva.de).

*NB:* The [SILVA license](http://www.arb-silva.de/fileadmin/silva_databases/current/LICENSE.txt) prohibits usage of the SILVA databases or parts of them within a non-academic/commercial environment beyond a 48h test period. If you want to use the SILVA databases with phyloFlash in a non-academic/commercial environment please contact them at contact(at)arb-silva.de.

The database has to be reformatted for use by phyloFlash. This is done with the script `phyloFlash_makedb.pl`. Known contamination sequences from cloning vectors are removed, repeat regions which can have an adverse effect on sequence reconstruction are masked, the database is clustered at 99% and 96% identity to speed up mapping/searching, and finally indexed for the read mapper.

The final disk space required for the default SILVA SSU database is about 5 Gb.

### 4.1. Downloading database automatically

To create a suitable database, just run

```bash
./phyloFlash_makedb.pl --remote
```

in the directory where you unpacked phyloFlash. The script will download the most current source databases and prepare the files required by `phyloFlash.pl`.

*NOTE: This currently only works if you are not behind a proxy*

If you are behind a proxy and cannot download the database via the script, you can download the current version of the SILVA database from [the SILVA website](https://www.arb-silva.de/no_cache/download/archive/current/Exports/). The filename should be `SILVA_XXX_SSURef_Nr99_tax_silva_trunc.fasta.gz` where `XXX` is the current version number. You should also download the UniVec database [from NCBI](https://www.ncbi.nlm.nih.gov/tools/vecscreen/univec/). Then proceed with the instructions in section 4.2 below.

### 4.2. Set up database from local copy of SILVA SSU NR99

If you already have a local copy of the SILVA SSU NR99 database (in Fasta format), and the NCBI Univec database, you can supply the paths:

```bash
./phyloFlash_makedb.pl --univec_file /path/to/Univec --silva_file /path/to/SILVA_128_SSURef_Nr99_tax_silva_trunc.fasta.gz
```
By default, `phyloFlash.pl` will look in the folder where it is installed for the subfolder with the highest SILVA version number. You can change this by passing the `-dbhome <dir>` switch to phyloFlash.pl or by modifying the `DBHOME` variable in `phyloFlash.pl`.

### 4.3. Set up a custom database with your own sequences

Users can supply their own databases of SSU rRNA sequences, or even other genes, in place of the SILVA SSU database, as long as they are formatted in the following way:

 - Sequences should be in Fasta format
 - Fasta headers should have the format `{IDENTIFIER}.{INTEGER}.{INTEGER} {TAXONOMY-STRING}` where:
   - IDENTIFIER is a unique sequence identifier which does not have spaces or periods
   - The difference between the two INTEGERS should be the length of the sequence, e.g. 1.1700 for a 1700 bp sequence
   - TAXONOMY-STRING is in SILVA or NCBI format, delimited by semicolons with no spaces (but spaces in taxon names allowed)
   - There is a single space before the TAXONOMY-STRING
 - The name of the Fasta file should begin with `SILVA_{DBNAME}_` where DBNAME is the name of the database (e.g. `CustomDB`), and will also be the name of the output folder containing the formatted database files. This is because the Fasta filename is parsed by the script.
