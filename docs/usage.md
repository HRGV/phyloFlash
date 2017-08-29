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

Compress all the output files into a tar.gz archive:
```
phyloFlash.pl -lib run01 -read1 reads_F.fq.gz -read2 -reads_R.fq.gz -zip
```

## 2. Full description of program options

You can access help on the command line with the following options:

`-check_env` Invokes checking of working environment and dependencies without data input. Use to test setup.

`-help` Print brief help message

`-man` Show manual

`-outfiles` Show detailed list of output and temporary files and exit.

### 2.1. Standard input arguments

`-lib LIBNAME` Library name to use as a filename prefix for the output files for this phyloFlash run. The name must be one word comprising only letters, numbers and `_` or `-` (no whitespace or other punctuation).

`-read1 FILENAME` Forward reads in FASTA or FASTQ formats. May be compressed with Gzip (.gz extension). If interleaved reads are provided, please use `--interleaved` flag in addition for paired-end processing.

`-read2 FILENAME` Reverse reads, for paired-end reads. If this option is omitted, phyloFlash will run in *experimental* single-end mode.

`-interleaved` Use this flag if read file is in interleaved format

`-readlength N` Set expected readlength (between 50 and 500). Always use if your read length differs from 100. Default: 100.

### 2.2. Performance-related

`-CPUs N` Number of threads to use. Defaults to all available CPU cores.

`-readlimit N` Limits processing to the first N reads in each input file that map to the reference database. Use this for transcriptomes with a lot of rRNA reads, and use values <1000000. Default: unlimited.

`-amplimit N` Set the limit of SSU read pairs to switch from emirge.py to emirge_amplicon.py. This feature is not reliable as emirge_amplicon.py has been problematic to run (use values >100000). Default: 500000.

### 2.3. Customizing the run

`-skip_spades` Do not use SPAdes to assemble full-length sequences from extracted reads

`-emirge` Use EMIRGE to reconstruct full-length sequences from extracted reads. (Default: Off)

`-poscov` Use Nhmmer to find positional coverage of reads across Barrnap's HMM model of the 16S and 18S rRNA genes from a subsample of reads, as an estimate of coverage evenness. (Default: Off)

`-id N` Minimum % identity of reads to map against reference database. Must be between 63 and 98. Set to a lower value for very divergent taxa. Default: 70.

`-clusterid N` % identity threshold for reference sequence clustering step. Must be between 50 and 100. Default: 97.

`-taxlevel N` Level in the taxonomy string to use for taxonomic units (NTUs), for the taxonomic summary and to estimate diversity. Must be an integer, and starts with 1 for the highest taxonomic level (Domain). Default: 4.

`-maxinsert N` Maximum insert size allowed for paired end read mapping. Must be between 0 and 1200. Default: 1200.

`-sc` Use if data are from single-cell MDA libraries, option is passed to the SPAdes assembler. (Default: Off)

`-dbhome DIR` Directory containing phyloFlash reference databases, prepared with `phyloFlash_makedb.pl`. (Default: Look in phyloFLash folder for highest SILVA version number)

### 2.4. Localization and compatibility options

`-crlf` Use CRLF as the line terminator in CSV output, to be RFC4180 compliant (Default: Off)

`-decimalcomma` Use decimal comma instead of decimal point to fix locale problems for some European systems (Default: Off)

### 2.5. Configuring output

`-html` Produce an HTML-formatted version of the report file. This helps improve readability and individual sections of the report can be collapsed. (Default: On, turn off with `-nohtml`)

`-treemap` Include an interactive treemap of the NTU counts in the HTML report. This uses the Google Visualization API, which requires an Internet connection and that you agree to their [terms of service](https://developers.google.com/chart/terms), and is not open-source although it is free to use. (Default: Off)

`-log` Write status messages printed to STDERR also to a log file (Default: Off)

`-zip` Compress output into a tar.gz archive file (Default: Off)

`-keeptmp` Keep temporary/intermediate files (Default: Off)

`-everything` Turn on all the optional analyses and output options. Options without defaults and any local settings must still be specified. Equivalent to `-emirge -poscov -treemap -zip -log`

`-almosteverything` Like `-everything` except without `-emirge`

## 3. Testing phyloFlash

You will find test data in the `test_files` folder. The test data provided contains subsampled SSU reads from SRA ERR138446, a *Caenorhabditis* sample with associated bacteria. You can test if phyloFlash is working properly with these files:

```bash
phyloFlash.pl -lib TEST -read1 test_files/test_F.fq.gz -read2 test_files/test_R.fq.gz
```

## 4. Expected performance

10 million 100 bp paired-end-reads of a metagenomic library are processed in less than 5 minutes on a normal 2014 desktop PC with 4 CPU cores and 8 GB of RAM.

*phyloFlash* usually detects most lifeforms on earth that have a SSU rRNA sequence that is at least 70% identical to anything in the databases, more exotic organisms might be problematic, but test data is hard to come by. If you happen to have such a test case and you are willing to share please drop me a line... If you think *phyloFlash* is not detecting a certain organism that is very distant from the known SSU rRNA sequences please try lowering the minimum sequence identity for a mapping hit by using e.g. `-id 0.63`.
