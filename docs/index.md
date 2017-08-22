---
layout: home
---

<img src="phyloFlash_logo.png" alt="phyloFlash logo" style="width: 200px; display: block; margin-left: auto; margin-right:auto;"/>

**phyloFlash** is a pipeline to rapidly reconstruct the SSU rRNAs and explore phylogenetic composition of an illumina (meta)genomic dataset.

This manual explains how to install and use phyloFlash. Navigate from the menu bar above or the table of contents below.

## Quick-start

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

## Contents

 - [Installation](install.md)
 - [Usage](usage.md)
 - [Comparing multiple samples](multiple-samples.md)
 - [FAQ](FAQ.md)

## About

phyloFlash is written by Harald Gruber-Vodicka ([Google Scholar](https://scholar.google.de/citations?user=imYEnqMAAAAJ&hl=en&oi=ao), [GitHub](https://github.com/HRGV)) and Elmar A. Pruesse ([Google Scholar](https://scholar.google.de/citations?user=F-yGwRIAAAAJ&hl=en&oi=ao), [GitHub](https://github.com/epruesse)) with contributions from Brandon Seah ([Google Scholar](https://scholar.google.de/citations?user=3l8G5BwAAAAJ&hl=en&oi=ao), [GitHub](https://github.com/kbseah))

You can find the source code for phyloFlash at GitHub:
[HRGV/phyloFlash](https://github.com/HRGV/phyloFlash/)

[Max Planck Institute for Marine Microbiology](http://www.mpi-bremen.de/)

## Citation

If you use phyloFlash for a publication, please cite as

Gruber-Vodicka HR, Pruesse E, Seah B. 2017. phyloFlash. Online:https://github.com/HRGV/phyloFlash

and also remember to cite the dependencies used.
