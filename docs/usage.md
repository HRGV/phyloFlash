---
layout: page
title: Usage
permalink: /usage/
order: 2
---

## Quick-start

```bash
phyloFlash.pl -lib LIBNAME -read1 READFILE_F.fq(.gz) -read2 READFILE_R.fq(.gz) [options]
phyloFlash.pl -help # Help page
phyloFlash.pl -man # Manual page in pager
```

Running `phyloFlash.pl` without arguments will show the basic help message.

## 1. Basic usage

To screen paired-end 100 bp read files named `reads_F.fq.gz` and `reads_R.fq.gz` for SSU rRNA sequences, and have output files labeled as "run01":
```bash
phyloFlash.pl -lib run01 -read1 reads_F.fq.gz -read2 -reads_R.fq.gz
```

Interleaved reads:
```bash
phyloFlash.pl -lib run01 -read1 reads_FR.fq.gz -interleaved
```

Longer read lengths:
```bash
phyloFlash.pl -lib run01 -read1 reads_F.fq.gz -read2 -reads_R.fq.gz -readlength 150
```

Limit number of processors used to 8:
```bash
phyloFlash.pl -lib run01 -read1 reads_F.fq.gz -read2 -reads_R.fq.gz -CPUs 8
```

Format results prettily as HTML:
```
phyloFlash.pl -lib run01 -read1 reads_F.fq.gz -read2 -reads_R.fq.gz -html
```

## 2. Full description of program options

### 2.1. Standard input arguments

`-lib LIBNAME` Library name for the phyloFlash run. This must be a word consisting of only letters, numbers, underscore and/or hyphen. This is used as the filename prefix for all output files.

`-read1 FILENAME` Forward reads in FASTA or FASTQ formats, maybe compressed with Gzip (.gz extension). If interleaved reads are provided, please use `--interleaved` flag in addition

`-read2 FILENAME` Reverse reads, for paired-end reads. If this option is omitted, phyloFlash will run in *experimental* single-end mode.

`-interleaved` Use this flag if read file is in interleaved format

`-readlength N` Set expected readlength (between 50 and 500). Always use if your read length differs from 100. Default: 100.

### 2.2. Performance-related

`-CPUs N` Number of threads to use. Defaults to all available CPU cores.

`-readlimit N` Limits processing to the first N reads in each input file that map to the reference database. Use this for transcriptomes with a lot of rRNA reads, and use values <1000000. Default: unlimited.

`-amplimit N` Set the limit of SSU read pairs to switch from emirge.py to emirge_amplicon.py. This feature is not reliable as emirge_amplicon.py has been problematic to run (use values >100000). Default: 500000.

## 2.3. Customizing the run

`-id N`

`-clusterid N`

`-maxinsert N`

`-sc` Use if data are from single-cell MDA libraries, option is passed to the SPAdes assembler.


### 2.3. Localization and compatibility options

`-crlf` Use CRLF as the line terminator in CSV output, to be RFC4180 compliant

`-decimalcomma` Use decimal comma instead of decimal point to fix locale problems for some European systems

### 2.4. Configuring output report

Prettify the output with the following options:

`-html` Produce an HTML-formatted version of the report file. This helps improve readability and individual sections of the report can be collapsed.

`-treemap` Include an interactive treemap of the NTU counts in the HTML report. This uses the Google Visualization API, which requires an Internet connection and that you agree to their [terms of service](https://developers.google.com/chart/terms), and is not open-source although it is free to use. (Experimental feature)

## 3. Testing phyloFlash

You will find test data in the `test_files` folder. The test data provided contains subsampled SSU reads from SRA ERR138446, a *Caenorhabditis* sample with associated bacteria. You can test if phyloFlash is working properly with these files:

```bash
phyloFlash.pl -lib TEST -read1 test_files/test_F.fq.gz -read2 test_files/test_R.fq.gz
```

## 4. Expected performance

10 million 100 bp paired-end-reads of a metagenomic library are processed in less than 5 minutes on a normal 2014 desktop PC with 4 CPU cores and 8 GB of RAM.

*phyloFlash* usually detects most lifeforms on earth that have a SSU rRNA sequence that is at least 70% identical to anything in the databases, more exotic organisms might be problematic, but test data is hard to come by. If you happen to have such a test case and you are willing to share please drop me a line... If you think *phyloFlash* is not detecting a certain organism that is very distant from the known SSU rRNA sequences please try lowering the minimal id by using e.g. `-minid 0.63`

## 5. Detailed description of output files

### 5.1. Report files

These are the main human-readable output from phyloFlash.

 - `LIBNAME.phyloFlash` *phyloFlash* report file, a simple text file with reconstructed SSU composition  report, performance and SSU metrics, and rough taxon list based on read mappings  
 - `LIBNAME.phyloFlash.html` Report file formatted in HTML, if `-html` flag was used

### 5.2. Unassembled sequence files

Reads that map to the reference database are extracted to these files in Fastq format

 - `LIBNAME.test_F.fq.gz.SSU.1.fq` the filtered SSU reads and their paired read, forward read file  
 - `LIBNAME.test_F.fq.gz.SSU.2.fq` the filtered SSU reads and their paired read, reverse read file  

### 5.3. Assembled/reconstructed sequence files

Assembled or reconstructed full-length SSU rRNA reads are output unless the `-skip_spades` or `-skip_emirge` options are used.

 - `LIBNAME.spades_rRNAs.final.fasta ` assembled OTUs from SPAdes with *phyloFlash* simplified headers

 - `LIBNAME.emirge.final.phyloFlash.notmatched.fa` a fasta file with the reconstructed SSU sequences with no significant hit in the provided SSU database
 - `LIBNAME.emirge.final.fa` a fasta file with the Emirge reconstructed SSU OTUs
 - `LIBNAME.emirge.final.phyloFlash.dbhits.fa` a fasta file with the best hits for the reconstructed SSU sequences in the provided SSU database

 - `LIBNAME.all.final.fasta` All assembled and reconstructed sequences from SPAdes and/or EMIRGE in a single file
 - `LIBNAME.all.final.phyloFlash.dbhits.fa`
 - `LIBNAME.all.final.phyloFlash.notmatched.fa`

 - `LIBNAME.all.dbhits.NR97.fa` Reference sequences from database with hits from the supplied reads, clustered at 97% identity

### 5.4. Alignments

 - `LIBNAME.SSU.collection.alignment.fasta` an aligned multifasta of all the predicted OTUs and the references
 - `LIBNAME.SSU.collection.fasta` a multifasta of all the predicted OTUs and the references
 - `LIBNAME.SSU.collection.fasta.tree` an NJ tree of the mafft alignment of all the predicted OTUs and the references. PDF and PNG versions are created for the HTML report if the `-html` option is set

### 5.5. Other statistics

 - `LIBNAME.inserthistogram` Histogram of detected insert sizes in tab-separated format, if paired-end reads were input. PDF and PNG versions are created for the HTML report if the `-html` option is set
 - `LIBNAME.idhistogram` Histogram of the % identity of reads vs. reference database sequences, in tab-separated format. PDF and PNG versions are created for the HTML report if the `-html` option is set
 - `LIBNAME.hitstats` Mapping statistics of reads mapping to the reference database, in tab-separated format. The unambiguous mapping hit counts are used to generate the NTU report
 - `LIBNAME.phyloFlash.NTUabundance.csv` the list of uniqe higher level taxa (e.g. orders for bacteria) in the order of  their apperance

 - `LIBNAME.scaffolds.arch.gff` 16S rRNA gene predictions for assembled OTUs based on archaeal SSU rRNA hmm profile  
 - `LIBNAME.scaffolds.bac.gff` 16S rRNA gene predictions for assembled OTUs based on bacterial SSU rRNA hmm profile  
 - `LIBNAME.scaffolds.euk.gff` 18S rRNA gene predictions for assembled OTUs based on eukaryote SSU rRNA hmm profile

### 5.6. CSV files used for multiple-sample comparison

These files are used by the `phyloFlash_heatmap.R` script if you wish to compare multiple samples by their taxonomic composition.

 - `LIBNAME.phyloFlash.NTUabundance.csv`
 - `LIBNAME.phyloFlash.report.csv`

### 5.7. Log files

Log files from the various tools used by the pipeline. If individual steps fail these can help to diagnose the problem.

 - `LIBNAME.barrnap.out` barrnap log  
 - `LIBNAME.bbmap.out` the bbmap log
 - `LIBNAME.spades.out` SPAdes log  
 - `LIBNAME.emirge.out` EMIRGE log
