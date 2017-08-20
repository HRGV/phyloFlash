phyloFlash v3.0 beta 1
==============

<img src="docs/phyloFlash_logo.png" width="200" alt="phyloFlash logo" />

by Harald Gruber-Vodicka and Elmar A. Pruesse with Brandon Seah.

*phyloFlash* is a pipeline to rapidly reconstruct the SSU rRNAs and explore phylogenetic composition of an Illumina (meta)genomic or transcriptomic dataset.

***NOTE*** Version 3.0 changes some input options and also how mapping-based taxa (NTUs) are handled. Please download the last release of v2.0 ([tar.gz archive](https://github.com/HRGV/phyloFlash/archive/v2.0-beta6.tar.gz)) for the old implementation.

Quick-start
-----------

```bash
# Download phyloFlash
wget https://github.com/HRGV/phyloFlash/archive/vXXX.tar.gz  
tar -xzf vXXX.tar.gz

# Check for dependencies
cd phyloFlash-XXX
./phyloFlash.pl -check_env

# Install reference database
./phyloFlash_makedb.pl --remote

# Run with test data
phyloFlash.pl -lib TEST -read1 test_files/test_F.fq.gz -read2 test_files/test_R.fq.gz

# Run with interleaved reads
phyloFlash.pl -lib run01 -read1 your_reads.fq.gz -interleaved
```

Use the `-help` option to display a brief help and the `-man` option to display a man-file.

Use the `-skip_spades` to turn off SSU rRNA sequence assembly with SPAdes, and `-emirge` to use EMIRGE for SSU rRNA sequence reconstruction.

Use the `-sc` switch for MDA datasets (single cell) or other hard to assemble read sets.

Use the `-zip` switch to compress output files into tar.gz archive

Plain text and HTML-formatted reports are produced, reporting summary statistics from each run. Use `-treemap` to draw an interactive treemap of taxonomic classification in the HTML report.

Manual
------

For further information **please refer to the [Manual](https://hrgv.github.io/phyloFlash)**.

Versions and changes
--------------------

* v3.0 beta 1
 * Re-map extracted SSU reads onto assembled sequences to check proportion assembled
 * Revamp of HTML report output. Embed interactive graphical summary, use SVG-formatted graphics, remove dependency on R packages for report graphics.
 * Changes to how mapping-based NTUs are calculated. Now count all reads (not only unambiguously-mapped) and count segments of read pairs separately.
* v2.0 complete rewrite

Contact
-------

Please report any problems to the [phyloFlash Google group](https://groups.google.com/forum/#!forum/phyloflash) or with the GitHub issue tracker.

(Pull requests with suggested fixes are of course also always welcome.)
