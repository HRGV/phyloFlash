---
layout: page
title: FAQ
order: 5
---

# What is phyloFlash for?

Quick overview of taxonomic composition and library quality in metagenomic and
metatranscriptomic samples, including multiple sample comparison.

# What input data can be analyzed?

Metagenomic and metatranscriptomic short-read paired-end sequencing, e.g.
Illumina. Raw reads are okay, without quality trimming or filtering.

Transcriptomes typically have a much higher proportion of rRNA reads than
genomes, even if libraries have been prepared with poly-A selection or rRNA
depletion. Limit the number of mapped reads to below 1 million with `-readlimit
1000000`.

**NOT** suitable for:
 - Sanger or PacBio libraries
 - Amplicon libraries - use a dedicated tool like mothur

# What output will I get?

The **phyloFlash** pipeline provides the following output:

An overview of the taxonomic composition in three ways:
 - Taxonomic affiliation of reference database sequences which have hits
 - Full-length SSU rRNA sequences assembled by SPAdes, and their closest
   reference database hits
 - Full-length SSU rRNA sequences reconstructed by EMIRGE, and their closest
   reference database hits (if `-emirge` option is used)

The extracted SSU rRNA reads will be mapped back to the full-length sequences
to determine the proportion that was successfully assembled or reconstructed. A
low ratio of assembled:unassembled could be caused by high taxonomic diversity
in a sample, or other problems with assembly (e.g. several closely-related taxa
that make assembly difficult).

You'll also get output parameters that will help you evaluate quality of your
sequences (see below)

# Conda gets stuck on "Solving environment" when I try to install phyloFlash - help!

In some cases, `conda install` or `conda create` can hang on the "Solving
environment" step. This appears to be because of ambiguities in dependency
specifications in packages on different channels (see this
[issue](https://github.com/conda/conda/issues/8197) on GitHub). Setting the
`channel_priority` to `strict` asks Conda to always pick the higher-priority
channel first when installing packages. This requires conda version to be 4.6
and above.

Run the following command before installing phyloFlash with Conda. It will
modify the `.condarc` configuration file of the current user by default
(usually `~/.condarc`).

```bash
conda config --set channel_priority strict
```

Add either the `--system` or `--env` flags to modify the system or environment
config respectively. See `conda config --help` and `conda config --describe`
for more information.

# How does phyloFlash show library quality?

 - Insert size histogram - should be more or less unimodal. Multiple peaks may
   be caused by contamination from other sequencing libraries.
 - Identify possible contamination, e.g. *Propionibacterium* is a common
   contaminant from human skin
 - Coverage on SSU rRNA HMM models (18S and 16S rRNA) - should be relatively
   even across the gene. High coverage or gaps in coverage may be a sign of
   contamination from amplicon libraries.
 - Mapping identity to database sequences - identities below 90% suggest the
   presence of novel or divergent taxa that are not well-represented in the
   database.

# Why use the SSU rRNA gene, and not single-copy markers / whole genomes / kmers?

Many other tools use genomic reference data to summarize the taxonomic
composition of microbial metagenomes. phyloFlash uses only the SSU rRNA marker
because it is the most well-represented and well-curated phylogenetic marker
gene. Tools using genome-based reference data are limited by the
underrepresentation of uncultivated, environmental taxa. These taxa are also
underrepresented in marker gene databases, but the latter have had several
decade's head start in accumulating comparative data. Using the SSU rRNA gene
also allows us to take advantage of the curated taxonomy from the
[SILVA](https://www.arb-silva.de/) project, whereas the taxonomic annotation of
genome sequences deposited in public databases is typically user-reported and
not always consistent.

# How can I use the results?

You can browse the results as a text file or HTML formatted report (with
embedded graphics) to look at the results from one sample.

To compare multiple samples you can produce a heatmap or barplot, which are
especially useful with large numbers of samples. [Several
utilities](utilities.html) are available for working with phyloFlash output.

Assembled and/or reconstructed full-length rRNA seqeuences can be used for
downstream analyses, for phylogenetic analysis.

Reads which map to SSU rRNA are extracted into Fastq formatted read files by
phyloFlash. These can also be used for further downstream analyses, e.g.
assembly to full-length sequences by other software, or mapping to other
reference databases.

# Should I run SPAdes and EMIRGE, or skip them?

phyloFlash incorporates two different tools for getting full-length rRNA
sequences from unassembled metagenomic/transcriptomic reads - one is the SPAdes
assembler (on by default), the other is EMIRGE (off by default) which uses an
expectation maximization method to reconstruct rRNA sequences.

If you don't need full-length sequences, you may want to turn off SPAdes
(`-skip_spades`) so that phyloFlash runs faster. In that case, only the initial
mapping to SILVA and the taxonomic summary will be produced. However
full-length sequences are useful to have because NTUs alone may give a
misleadingly high impression of taxonomic diversity.

EMIRGE in particular can be tricky to install because of its own software
dependencies, and you when you first try the pipeline you may want to skip
EMIRGE if you don't already have it installed on your system.

# How is the taxonomic summary produced?

The taxonomic affiliation of the SSU rRNA reads is taken from the reference
sequence to which it maps in the SILVA database. From v3.2b1 onwards, the
taxonomy reported is the last-common-ancestor consensus of the top hits, i.e.
if a read has more than one best-scoring hit, it will report the lowest
taxonomic level which they have in common. In previous versions, only the
single best hit was retained, but this would result in overly-preecise
taxonomic assignments especially for divergent sequences. To replicate the old
behavior, use the `--tophit` option.

# How can I use phyloFlash in a Snakemake pipeline?

[Snakemake](https://snakemake.readthedocs.io/en/stable/) is a workflow engine
for bioinformatics, which allows users to document workflows and automatically
detect dependencies that have to be updated. The concept is similar to [GNU
Make](https://www.gnu.org/software/make/).

Steps in a Snakemake pipeline are called "rules", and each rule definition
specifies an input file pattern, an output file pattern, and the command that
produces the output from the input. The input/output file patterns can be
defined e.g. by a list or regular expressions, so this allows plenty of
flexibility.

There are a few possible stumbling blocks when incorporating phyloFlash into a
Snakemake pipeline:
 * Snakemake takes files as input/output arguments, whereas the `-lib` option
   of phyloFlash is just the file prefix.
 * phyloFlash produces output in the current working folder, whereas a user
   might want to put all output from a given rule into a single folder. However
   phyloFlash will not accept a path prefix (with `/` character) in the `-lib`
   option.

Here is an example Snakemake rule for phyloFlash, assuming that you have a
[configuration
file](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html)
listing all the sample names under `samples` and the path to phyloFlash
database under `phyloFlash_db`.

```python
rule phyloFlash:
    input:
        "reads/{sample}.fastq.gz"
    output:
        expand("phyloFlash/{sample}.phyloFlash.tar.gz",sample=config['samples'])
    threads: 16
    log: "logs/{sample}_phyloFlash.log"
    conda: "envs/phyloflash.yml" # Specify a Conda environment containing phyloFlash
    params:
        dbhome=config['phyloFlash_db']
    shell:
        "phyloFlash.pl -lib {wildcards.sample} -dbhome {params.dbhome} -read1 {input} -interleaved -CPUs {threads} -almosteverything; mv {wildcards.sample}.phyloFlash.* phyloFlash/"
```

# How do I get help if something doesn't work?

First check if you have all the dependencies installed properly with
`phyloFlash.pl -check_env`.

See the built-in help messages and manuals for any script in the phyloFlash
folder by using either the `--help` or `--man` options.

To get in touch, please use the [phyloFlash Google
group](https://groups.google.com/forum/#!forum/phyloflash), or the [issues
tracker](https://github.com/HRGV/phyloFlash/issues) on GitHub to submit bug
reports, and supply us with information that will help us diagnose the problem.
We'll try to help if we can but please understand that phyloFlash is for
academic, non-commercial use and we may not always be able to respond promptly.
