<img align="right" src="docs/phyloFlash_logo.png" width="200" alt="phyloFlash logo"/>

# phyloFlash

[![GitHub (pre-)release](https://img.shields.io/github/release/HRGV/phyloflash/all.svg?label=Latest%20Version)]()
[![Bioconda](https://img.shields.io/conda/vn/Bioconda/phyloFlash.svg)](https://bioconda.github.io/recipes/phyloflash/README.html)
[![](https://img.shields.io/conda/dn/Bioconda/phyloflash.svg)](https://anaconda.org/bioconda/phyloflash/files)

by Harald Gruber-Vodicka, Elmar A. Pruesse, and Brandon Seah.

***NOTE: This software is only being sporadically maintained. We regret that we
are unable to respond to issues on a regular basis.***

*phyloFlash* is a pipeline to rapidly reconstruct the SSU rRNAs and explore
phylogenetic composition of an Illumina (meta)genomic or transcriptomic
dataset. **[Manual](https://hrgv.github.io/phyloFlash)**

Read [our paper](https://doi.org/10.1128/mSystems.00920-20) on phyloFlash.


## Quick-start

### Download via Conda

We recommend installing phyloFlash and its dependencies using Conda or Mamba.
[Conda](https://conda.io/docs/) is a package manager that will also install
dependencies that are required if you don't have them already.

phyloFlash is distributed through the [Bioconda](http://bioconda.github.io/)
channel on Conda.

According to the [Conda documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html),
avoid installing new packages to your base environment but create new
environments for them as required. Also, specify all desired packages at the
same time when creating a new environment, instead of adding them sequentially,
to avoid dependency conflicts.

We also suggest using [Mamba](https://mamba.readthedocs.io/en/latest/) as a
drop-in substitute for Conda. It implements a more effective dependency solver
and is also the default Conda frontend for the pipeline manager Snakemake.
Simply replace `conda` with `mamba` in the commands below. Note that the
`defaults` channel should be enabled.

```bash
# If you haven't set up Bioconda already
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
# Try the following step if "solving environment" does not terminate
conda config --set channel_priority strict
# Create new environment named "pf" with phyloflash
conda create -n pf phyloflash
# Activate environment
conda activate pf
# Check that dependencies all installed properly
phyloFlash.pl -check_env
```


### Download pre-formatted database

Pre-formatted databases derived from SILVA releases 138 onwards are available
from the following Zenodo archives:

 * [SILVA 138.1](https://doi.org/10.5281/zenodo.7892521) (latest)
 * [SILVA 138.1, taxonomy with main ranks only](https://doi.org/10.5281/zenodo.10047346) (see details in repository)
 * [SILVA 138](https://doi.org/10.5281/zenodo.7890453)

Download, checksum, and unpack:

```bash
wget https://zenodo.org/record/7892522/files/138.1.tar.gz # 5.5 GB download
tar -xzf 138.1.tar.gz # unpacks folder 138.1/ in the current location
```

Specify path to the database folder with the option `-dbhome` when running
phyloFlash (see below).

Older versions of the SILVA database have a more restrictive license, so we are
unable to distribute pre-formatted versions. You will have to download the
original SILVA files and run the `phyloFlash_makedb.pl` script yourself (see
Manual).


### Test phyloFlash with test dataset

Test data are included with phyloFlash. The following assumes that you
installed phyloFlash to a Conda environment called `pf`, and that the database
files have been unpacked to a folder `/path/to/138.1`. By default, phyloFlash
will look for the database folder in the folder where it is installed. If it is
located somewhere else, specify this to the `-dbhome` option.

```bash
conda activate pf # If Conda environment not already activated
phyloFlash.pl -dbhome /path/to/138.1 -lib TEST -CPUs 16 \
 -read1 ${CONDA_PREFIX}/lib/phyloFlash/test_files/test_F.fq.gz \
 -read2 ${CONDA_PREFIX}/lib/phyloFlash/test_files/test_R.fq.gz \
 -almosteverything
```


### Example phyloFlash commands

```bash
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

# Run SPAdes (skip EMIRGE) and produce all optional outputs (recommended)
phyloFlash.pl -lib LIB -almosteverything -read1 reads_F.fq.gz -read2 reads_R.fq.gz

# Supply trusted contigs containing SSU rRNA sequences to screen vs reads
phyloFlash.pl -lib LIB -read1 reads_F.fq.gz -read2 reads_R.fq.gz -trusted contigs.fasta

# Use SortMeRNA instead of BBmap for initial mapping (slower, but more sensitive)
phyloFlash.pl -lib LIB -read1 reads_F.fq.gz -read2 reads_R.fq.gz -sortmerna
```

Use the `-help` option to display a brief help and the `-man` option to display
the full help message.

Use the `-sc` switch for MDA datasets (single cell) or other hard to assemble
read sets.

Use the `-zip` switch to compress output files into tar.gz archive, and `-log`
to save run messages to a log file


## Output

phyloFlash screens metagenomic or metatranscriptomic reads for SSU rRNA
sequences by mapping against the SILVA SSU Ref database.

Extracted reads are assembled and/or reconstructed into full-length sequences,
and the proportion of reads assembled is estimated by re-mapping them to the
full-length seqeunces.

phyloFlash reports a taxonomic summary of the reads from the initial mapping,
the best database matches of the full-length sequences, and a second taxonomic
summary of the unassembled reads left over.

Plain text and HTML-formatted reports are produced, reporting summary
statistics from each run. The HTML report includes an interactive graphical
summary.


## Going further

The phyloFlash suite also includes other tools for SSU rRNA-centric metagenome
analyses. Run the commands without arguments to see help messages.

 * `ENA_phyloFlash.pl` - Automatically download read files from
   [ENA](https://www.ebi.ac.uk/ena) given a read accession number, and run
   phyloFlash on them
 * `phyloFlash_compare.pl` - Compare the taxonomic composition of multiple
   samples from their phyloFlash results. This produces a barplot, heatmap, or
   distance matrix based on the NTU abundances in two or more samples.
 * `phyloFlash_fastgFishing.pl` - Given a metagenomic assembly graph in
   [Fastg](http://fastg.sourceforge.net/) format, identify SSU rRNA sequences
   and extract contigs connected to them. Optionally compare to phyloFlash
   results from the same library.


## Manual

For further information please refer to the
[Manual](https://hrgv.github.io/phyloFlash) as well as the command-line help
page `phyloFlash.pl -man`.


## Versions and changes

* v3.4
  * Automatic detection of read lengths
  * More informative error messages
  * Optionally change minimum reference length when building database, e.g. for
    custom non-SSU rRNA databases
* v3.3 beta 4
  * Fix bug with double counting of reads that was introduced in v3.3b3
* v3.3 beta 3
  * Single-end reads can now be used as input
* v3.3 beta 2
  * New options to graphical comparison scripts, and other small bug fixes
  * Fix bug due to change in SILVA project file naming convention with SILVA
    138 onwards
* v3.3 beta 1
  * Add support for using SortMeRNA instead of BBmap for initial mapping step
  * Changes to how mapping data are hashed; process SAM file of initial mapping
    to fix known bugs with bitflag and read name reporting in BBmap and
    SortMeRNA
* v3.2 beta 1
  * Report ambiguous hits during mapping step, use consensus of top hits to
    assign taxonomy instead of single best hit
  * Add utility `phyloFlash_compare.pl` to compare taxonomic composition of
    multiple libraries from phyloFlash output
  * Add utility `phyloFlash_fastgFishing.pl` to extract genome bins from Fastg
    files
* v3.1 beta 2
  * Fix bug in Fasta headers with changed output from Bedtools v2.26+
  * Make bbmap and reformat.sh overwrite existing output files of same name
  * Rearrange elements in HTML report file
* v3.1 beta 1
  * Allow user to supply "trusted contigs" of sequence assemblies containing
    SSU rRNA which will also be screened against the read libraries
  * Fix bugs in tree plotting
* v3.0 beta 1
  * Re-map extracted SSU reads onto assembled sequences to check proportion
    assembled
  * Revamp of HTML report output. Embed interactive graphical summary, use
    SVG-formatted graphics, remove dependency on R packages for report
    graphics.
  * Changes to how mapping-based NTUs are calculated. Now count all reads (not
    only unambiguously-mapped) and count segments of read pairs separately.
  * No change to heatmap script for comparing multiple samples
* v2.0 complete rewrite


## Contact

Please report any problems to the [phyloFlash Google
group](https://groups.google.com/forum/#!forum/phyloflash) or with the GitHub
issue tracker.

(Pull requests with suggested fixes are of course also always welcome.)

We also welcome any feedback on the software and its documentation, especially
suggestions for improvement!


## Acknowledgements

We thank colleagues and phyloFlash users who have contributed to phyloFlash
development by testing the software, reporting bugs, and suggesting new
features.


## Citation

If you use phyloFlash for a publication, please cite our paper in _mSystems_:

Harald R Gruber-Vodicka, Brandon KB Seah, Elmar Pruesse. phyloFlash: Rapid SSU
rRNA profiling and targeted assembly from metagenomes. [ *mSystems* 5 :
e00920-20](https://doi.org/10.1128/mSystems.00920-20);
doi:10.1128/mSystems.00920-20

and also remember to cite the dependencies used, which are listed in each
phyloFlash report file.
