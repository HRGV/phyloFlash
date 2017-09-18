#!/usr/bin/env perl

# Given run accession number, download Fastq files from ENA and run phyloflash

use 5.010; # For the 'say' command
use strict;
use warnings;

use LWP::Simple;    # to read data from URLs
use FindBin;
use lib $FindBin::RealBin;
use Getopt::Long;
use POSIX;

# Input params
my $acc;
my $CPUs = 8; # Num processors - for nhmmer
my $download = 1; # By default, download
my $run_pF = 0; # By default, do not run phyloFlash
my $pF_dbhome;

# Get input from command line
GetOptions ("acc=s" => \$acc,
            "CPUs=i" => \$CPUs,
            "download!" => \$download,
            "phyloFlash!" => \$run_pF,
            "dbhome=s" => \$pF_dbhome,
            ) or die ("Error with input options");

# Variables
my @fastq_urls; # List of Fastq URLs
my @fastq_basenames;
my $phyloFlash = "phyloFlash.pl"; # assume phyloFlash script in path
my $wget = "wget --no-verbose";   # Wget binary and options

## MAIN ########################################################################

@fastq_urls = get_fastq_urls($acc);

# Check if any URLs were retrieved
if (defined $fastq_urls[0]) {
    # Open log file to record details on this file
    # Header line for log file
    foreach my $fastq (@fastq_urls) {
        system (join " ", ($wget, $fastq)) if ($download == 1);
        # Strip URL dirs from filename, save to an array
        my $filename = $1 if $fastq =~ m/([^\/]+)$/;
        push @fastq_basenames, $filename;
    }
}

# Run phyloFlash (v3.0beta1) to extract SSU reads
my @pF_args = ("-lib pF_$acc",
               "-CPUs $CPUs",
               "-emirge",               # Run both SPAdes and EMIRGE
               "-zip",                  # Output to archive
               "-log",
               );

# If specific dbhome specified, pass to phyloFlash, otherwise use default detected
push @pF_args, "-dbhome $pF_dbhome" if defined $pF_dbhome;

# Check how many Fastq files (paired or unpaired)
if (scalar (@fastq_basenames) == 1) {
    push @pF_args, "-read1 ".$fastq_basenames[0];
} elsif (scalar (@fastq_basenames) == 2) {
    push @pF_args, "-read1 ".$fastq_basenames[0];
    push @pF_args, "-read2 ".$fastq_basenames[1];
}

system (join " ", ($phyloFlash, @pF_args)) if $run_pF == 1;

## SUBS ########################################################################


sub get_fastq_urls {
    # Get the URL(s) of ENA Fastq file(s) for a given ENA entry

    # Input: Accession number of sample or run
    my ($acc) = @_;
    my @urls_arr;    # Output array containing URLs
    # Get report table using ENA REST query
    my $url = "http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=$acc&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes";
    my $tab = get $url;
    foreach my $line (split /\n/, $tab) {
        next if $line =~ m/^run_accession/; # Skip header
        my @splitline = split /\t/, $line;
        # Multiple URLs are separated by a semicolon
        my @spliturl = split /;/, $splitline[1];
        #my $run_accession = $splitline[0]; # run_accession, different from sample acc
        foreach my $fastq_url (@spliturl) {
            push @urls_arr, $fastq_url;
        }
    }
    # Return list of URLs
    return (@urls_arr);
}
