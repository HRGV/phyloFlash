phyloFlash 2.0
==============

<img src="docs/phyloFlash_logo.png" style="width:200px" alt="phyloFlash logo" />

by Harald Gruber-Vodicka and Elmar A. Pruesse with Brandon Seah.

*phyloFlash* is a pipeline to rapidly reconstruct the SSU rRNAs and explore phylogenetic composition of an illumina (meta)genomic dataset.

> NOTE: This software is still under development and might not work as expected or might not work at all.


Quick-start
-----------

```bash
# Download phyloFlash
wget https://github.com/HRGV/phyloFlash/archive/v2.0-beta6.tar.gz  
tar -xzf v2.0-beta6.tar.gz

# Check for dependencies
cd phyloFlash-2.0
./phyloFlash.pl -check_env

# Install reference database
./phyloFlash_makedb.pl --remote

# Run with test data
phyloFlash.pl -lib TEST -read1 test_files/test_F.fq.gz -read2 test_files/test_R.fq.gz

# Run with interleaved reads
phyloFlash.pl -lib run01 -read1 your_reads.fq.gz -interleaved
```

Use the `-help` option to display a brief help and the `-man` option to display a man-file.

Use the `-skip_spades` and/or `-skip_emirge` options to turn off SSU sequence reconstruction with SPAdes assembler or EMIRGE respectively.

Use the `-sc` switch for MDA datasets (single cell) or other hard to assemble read sets.

Use the `-html` switch to generate HTML-formatted report, and `-treemap` to draw an interactive treemap of taxonomic classification in the HTML report.

Manual
------

For further information **please refer to the [Manual](https://kbseah.github.io/phyloFlash)**.


Versions and changes
--------------------

2.00 complete rewrite


Contact
-------

Please report any problems to the [phyloFlash Google group](https://groups.google.com/forum/#!forum/phyloflash) or with the GitHub issue tracker.

(Pull requests with suggested fixes are of course also always welcome.)
