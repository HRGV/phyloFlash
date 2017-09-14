phyloFlash v3.0 beta 1
==============

<img src="docs/phyloFlash_logo.png" width="200" alt="phyloFlash logo" />

by Harald Gruber-Vodicka, Elmar A. Pruesse, and Brandon Seah.

*phyloFlash* is a pipeline to rapidly reconstruct the SSU rRNAs and explore phylogenetic composition of an Illumina (meta)genomic or transcriptomic dataset.

***NOTE*** Version 3.0 changes some input options and also how mapping-based taxa (NTUs) are handled. Please download the last release of v2.0 ([tar.gz archive](https://github.com/HRGV/phyloFlash/archive/v2.0-beta6.tar.gz)) for the old implementation. No changes have been made to the database setup, so databases prepared for v2.0 can still be used for v3.0.

Quick-start
-----------

```bash
# Download phyloFlash
wget https://github.com/HRGV/phyloFlash/archive/pf3.0b1.tar.gz
tar -xzf pf3.0b1.tar.gz

# Check for dependencies
cd phyloFlash-pf3.0b1
./phyloFlash.pl -check_env

# Get missing depencies - the easiest way is to install conda/bioconda - https://conda.io/miniconda.html
# First add bioconda to the conda channels and then grab what you need
conda config --add channels bioconda
 
conda install emirge
conda install bbmap
conda install vsearch
conda install spades
conda install mafft
conda install bedtools

# Install reference database
./phyloFlash_makedb.pl --remote

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
```

Use the `-help` option to display a brief help and the `-man` option to display the full help message.

Use the `-sc` switch for MDA datasets (single cell) or other hard to assemble read sets.

Use the `-zip` switch to compress output files into tar.gz archive, and `-log` to save run messages to a log file

Output
------

phyloFlash screens metagenomic or metatranscriptomic reads for SSU rRNA sequences by mapping against the SILVA SSU Ref database.

Extracted reads are assembled and/or reconstructed into full-length sequences, and the proportion of reads assembled is estimated by re-mapping them to the full-length seqeunces.

phyloFlash reports a taxonomic summary of the reads from the initial mapping, the best database matches of the full-length sequences, and a second taxonomic summary of the unassembled reads left over.

Plain text and HTML-formatted reports are produced, reporting summary statistics from each run. The HTML report includes an interactive graphical summary.

Manual
------

For further information **please refer to the [Manual](https://hrgv.github.io/phyloFlash)**.

Versions and changes
--------------------

* v3.0 beta 1
  * Re-map extracted SSU reads onto assembled sequences to check proportion assembled
  * Revamp of HTML report output. Embed interactive graphical summary, use SVG-formatted graphics, remove dependency on R packages for report graphics.
  * Changes to how mapping-based NTUs are calculated. Now count all reads (not only unambiguously-mapped) and count segments of read pairs separately.
  * No change to heatmap script for comparing multiple samples
* v2.0 complete rewrite

Contact
-------

Please report any problems to the [phyloFlash Google group](https://groups.google.com/forum/#!forum/phyloflash) or with the GitHub issue tracker.

(Pull requests with suggested fixes are of course also always welcome.)
