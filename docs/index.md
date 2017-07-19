# phyloFlash
## Version 2.0

![](phyloFlash_logo.png)

by Harald Gruber-Vodicka and Elmar A. Pruesse with Brandon Seah.

*phyloFlash* is a pipeline to rapidly reconstruct the SSU rRNAs and
explore phylogenetic composition of an illumina (meta)genomic dataset.

> NOTE: This software is still under development and might not work as
> expected or might not work at all.

## Updates

## What phyloFlash does


## Install phyloFlash and dependencies

To use phyloFlash you will need a GNU/Linux system with Perl, R and
Python installed. (OSX is for the brave, we have not tested this!)

1. Download and unpack phyloFlash. You can find the most recent release
   here https://github.com/HRGV/phyloFlash/releases

   ```bash
   wget https://github.com/HRGV/phyloFlash/archive/v2.0-beta5.tar.gz  
   tar -xzf v2.0-beta5.tar.gz
   ```

2. Install the tools phyloFlash uses:

 - Perl >= 5.13.2 (http://www.perl.org/get.html)
 - EMIRGE and its dependencies (https://github.com/csmiller/EMIRGE)
 - bbmap (http://sourceforge.net/projects/bbmap/)
 - vsearch (https://github.com/torognes/vsearch)
 - spades (http://bioinf.spbau.ru/spades)
 - bedtools (https://github.com/arq5x/bedtools2)
 - mafft (http://mafft.cbrc.jp/alignment/software/)
 - barrnap (currently provided with phyloFlash)

 These tools need to be "in your $PATH" so that *phyloFlash* can find
 them. To see whether all required tools are available, just run
 *phyloFlash* with the option ```check_env```:

 ```bash
 cd phyloFlash-2.0
 ./phyloFlash.pl -check_env
 ```
## Setting up the reference database

phyloFlash uses modified versions of the SILVA SSU database that is maintained
by the SILVA team - for more information visit their webpage www.arb-silva.de.

The SILVA license prohibits usage of the SILVA databases or parts of them within
a non-academic/commercial environment beyond a 48h test period. If you want to use
the SILVA databases with phyloFlash in a non-academic/commercial environment please
contact them at contact(at)arb-silva.de.

The detailed SILVA license can be found at
http://www.arb-silva.de/fileadmin/silva_databases/current/LICENSE.txt.

To create a suitable database, just run

```bash
./phyloFlash_makedb.pl --remote
```

in the directory where you unpacked *phyloFlash*. The script will download
the most current source databases and prepare the files required by ```phyloFlash.pl```.

> NOTE: This currently only works if you are not behind a proxy

If you have a local copy of the SILVA SSU NR99 database (Fasta format)
and the NCBI Univec database, you can supply the paths:

```bash
./phyloFlash_makedb.pl --univec_file /path/to/univec --silva_file /path/to/silva_db
```

By default, ```phyloFlash.pl``` will look for a directory named "119"
(the most recent SILVA release as of March 2015). You can change this
by passing the ```-dbhome <dir>``` switch to phyloFlash.pl or
by modifying the ```DBHOME``` variable in phyloFlash.pl

## Use phyloFlash

```bash
./phyloFlash.pl -lib LIBRARYNAME -read1 READFILE_F.fq(.gz) -read2 READFILE_R.fq(.gz) [options]
./phyloFlash.pl -help # Help page
./phyloFlash.pl -man # Manual page in pager
```

## phyloFlash output

## phyloFlash vs other pipelines

## FAQ

## Contact
