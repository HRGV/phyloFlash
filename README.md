phyloFlash 2.0
==============

*phyloFlash* is a pipeline to rapidly reconstruct the SSU rRNAs and
explore phylogenetic composition of an illumina (meta)genomic dataset.

> NOTE: This software is still under development and might not work as
> expected or might not work at all.

Description
-----------

*phyloFlash* accepts single read files or paired end fastq read data
provided in separate files. The files do not need to be unzipped (.gz
is ok).

Raw right off the illumina machine works fine. Usually, neither
trimming nor quality filtering are necessary.

If you are using transcriptomic data with a high proportion of rRNA,
please limit the read number reads to < 1 million reads for reasonable
runtimes and results.


Installing phyloFlash
---------------------

To use phyloFlash you will need a GNU/Linux system with Perl, R and
Python installed. (OSX is for the brave, we have not tested this!)

1. Download and unpack phyloFlash. You can find the most recent release
   here https://github.com/HRGV/phyloFlash/releases

   ```bash
   wget https://github.com/HRGV/phyloFlash/releases/phyloflash-2.0.tgz
   tar -xzf phyloFlash-2.0.tgz
   ```

2. Install the tools phyloFlash uses:

 - EMIRGE and its dependencies (https://github.com/csmiller/EMIRGE)
 - bbmap (http://sourceforge.net/projects/bbmap/)
 - vsearch (https://github.com/torognes/vsearch)
 - spades (http://bioinf.spbau.ru/spades)
 - barrnap (http://www.vicbioinformatics.com/software.barrnap.shtml)
 - fastaFromBed (https://github.com/arq5x/bedtools2/releases)

 These tools need to be "in your $PATH" so that *phyloFlash* can find
 them. To see whether all required tools are available, just run
 *phyloFlash* without arguments:

 ```bash
 cd phyloFlash-2.0
 ./phyloFlash.pl
 ```

Preparing the Reference Database
--------------------------------

To create a suitable database, just run

```bash
./phyloFlash_makedb.pl
```

in the directory where you unpacked *phyloFlash*. The script will download
the most current source databases and prepare the files required by ```phyloFlash.pl```.

By default, ```phyloFlash.pl``` will look for a directory named "119" 
(the most recent SILVA release as of March 2015). You can change this
by passing the "-dbhome <dir>" switch to phyloFlash.pl or 
by modifying the "DBHOME" variable in phyloFlash.pl 

Synopsis
--------
```bash
./phyloFlash.pl -lib LIBRARYNAME -read1 READFILE_F.fq(.gz) -read2 READFILE_R.fq(.gz) [options]
./phyloFlash.pl -help
./phyloFlash.pl -man
```

Use the ```-help``` option to display a brief help and the ```-man``` 
option to display a man-file. 

Expected Performance
--------------------

10Mio 100bp PE-reads of a metagenomic library are processed in less
than 5 minutes on a normal 2014 desktop PC with 4 CPU cores and 4 GB
of RAM.

*phyloFlash* usually detects most lifeforms on earth that have a SSU
rRNA sequence that is at least 70% identical to anything in the
databases, more exotic organisms might be problematic, but test data
is hard to come by. If you happen to have such a test case and you are
willing to share please drop me a line... If you think *phyloFlash* is
not detecting a certain organism that is very distant from the known
SSU rRNA sequences please try lowering the minimal id by using
e.g. -minid .63

Test data
---------

The test data provided contains subsampled SSU reads from SRA
ERR138446, a Caenorhabditis sample with associated bacteria.

run with

```bash
perl phyloFlash.pl -lib test -read1 test_F.fq.gz -read2 test_R.fq.gz 
```

Output files
------------

YOUR_LIBRARY_NAME.phyloFlash - *phyloFlash* report file, a simple text file with reconstructed SSU composition report, performance and SSU metrics, and rough taxon list based on read mappings 

YOUR_LIBRARY_NAME.test_F.fq.gz.SSU.1.fq - the filtered SSU reads and their paired read, forward read file
YOUR_LIBRARY_NAME.test_F.fq.gz.SSU.2.fq - the filtered SSU reads and their paired read, reverse read file
YOUR_LIBRARY_NAME.inserthistogram - a csv file of the detected insert sizes, is plotted in a diagram if the -html option is set
YOUR_LIBRARY_NAME.hitstats - the mapping statistics reporting the hits accross the SSU database, this is analyzed for the unambigous mapping hit counts
YOUR_LIBRARY_NAME.spades.final.phyloFlash - assembled OTUs results list including taxonomic affiliation
YOUR_LIBRARY_NAME.spades.final.phyloFlash.dbhits.fa - database hits of the assembled OTUs
YOUR_LIBRARY_NAME.spades.final.phyloFlash.notmatched.fa - assembled OTUs that do not have a database hit
YOUR_LIBRARY_NAME.spades.out - SPAdes log
YOUR_LIBRARY_NAME.spades_rRNAs.fasta - assembled OTUs with original headers
YOUR_LIBRARY_NAME.spades_rRNAs.final.fasta  - assembled OTUs with *phyloFlash* simplified headers
YOUR_LIBRARY_NAME.barrnap.out - barrnap log
YOUR_LIBRARY_NAME.scaffolds.arch.gff - 16S predictions for assembled OTUs based on archeal SSU rRNA hmm profile
YOUR_LIBRARY_NAME.scaffolds.bac.gff - 16S predictions for assembled OTUs based on bacterial SSU rRNA hmm profile
YOUR_LIBRARY_NAME.scaffolds.euk.gff - 18S predictions for assembled OTUs based on eukaryote SSU rRNA hmm profile
YOUR_LIBRARY_NAME.emirge.final.phyloFlash.notmatched.fa - a fasta file with the reconstructed SSU sequences with no significant hit in the provided SSU database
YOUR_LIBRARY_NAME.emirge.final.fa - a fasta file with the Emirge reconstructed SSU OTUs
YOUR_LIBRARY_NAME.emirge.out - the log of the EMIRGE run
YOUR_LIBRARY_NAME.emirge.final.phyloFlash.dbhits.fa - a fasta file with the best hits for the reconstructed SSU sequences in the provided SSU database
YOUR_LIBRARY_NAME.taxa - the list of uniqe higher level taxa (e.g. orders for bacteria) in the order of their apperance 
YOUR_LIBRARY_NAME.bbmap.out - the bbmap log
YOUR_LIBRARY_NAME.all.dbhits.NR97.fa - a fasta file with the collected reference hits clustered at 97% identity
YOUR_LIBRARY_NAME.SSU.collection.alignment.fasta - an aligned multifasta of all the predicted OTUs and the references 
YOUR_LIBRARY_NAME.SSU.collection.fasta - a multifasta of all the predicted OTUs and the references 
YOUR_LIBRARY_NAME.SSU.collection.fasta.tree - an NJ tree of the mafft alignment of all the predicted OTUs and the references 

Versions and changes
--------------------

2.00 complete rewrite 

Known Bugs
----------

None as yet. 

Please report any problems to *phyloFlash* google
group at https://groups.google.com/forum/#!forum/phyloflash

(Pull requests with suggested fixes are of course also always welcome.)
