#!/usr/bin/env perl

# Given run accession number, download Fastq files from ENA and run phyloflash

=head1 NAME ENA_phyloFlash.pl - Download read files from ENA read archive and run phyloFlash

=head1 SYNOPSIS

./ENA_phyloFlash.pl -acc [ENA ACCESSION] -phyloFlash -http_proxy="http://..."

./ENA_phyloFlash.pl -help

=head1 DESCRIPTION

Download Fastq files associated with an ENA Read Archive run accession number,
and run phyloFlash on those files if option -phyloFlash is chosen.

Requires phyloFlash v3.0+ and wget

=cut

use 5.010; # For the 'say' command
use strict;
use warnings;

#use LWP::Simple;    # to read data from URLs
use LWP::UserAgent;
use FindBin;
use lib $FindBin::RealBin;
use Getopt::Long;
use POSIX;
use Pod::Usage;
use PhyloFlash qw(msg err @msg_log $VERSION);

# Input params
my $acc;
my $CPUs = 8; # Num processors - for nhmmer
my $download = 1; # By default, download
my $download_dir = ".";
my $cleanup;
my $run_pF; # By default, do not run phyloFlash
my $pF_dbhome;
my $http_proxy;

pod2usage(2) if (!@ARGV);

# Get input from command line
GetOptions ("acc=s" => \$acc,
            "CPUs=i" => \$CPUs,
            "download|d!" => \$download,
            "download_dir=s" => \$download_dir,
            "cleanup" => \$cleanup,
            "phyloFlash|p!" => \$run_pF,
            "dbhome=s" => \$pF_dbhome,
            "http_proxy=s" => \$http_proxy,
            'help|h' => sub { pod2usage(1) },
            'man' => sub { pod2usage(-exit=>0, -verbose=>2) },
            'version|v'  => sub { welcome(); exit(); },
            ) or pod2usage(2);

=head1 ARGUMENTS

=over 8

=item --acc STRING

Accession number of read run. Please ensure that this is the I<Run> accession
and not e.g. sample or study accession number. They are typically prefixed by
ERR..., SRR..., or DRR.... 

=item --CPUs I<INTEGER>

Number of processors to use (passed to phyloFlash) (Default: 8)

=item -d|--download

Download files? Turn off with -nodownload (Default: Yes)

=item --download_dir I<PATH>

Path to download files.

Default: Current folder

=item --cleanup

Delete downloaded read files

Default: No

=item -p|--phyloFlash 

Run phyloFlash? (Default: No, download only)

=item --dbhome I<PATH>

Path to phyloFlash database, otherwise use default, which is to look in folder
where phyloFlash script is located

=item --http_proxy I<URL>

URL for http proxy - be sure to specify the protocol with "http://..."
(Default: none)

=item --help

Show help message

=item --man

Show manual page in pager

=item --version

Print version number

=back

=cut

# Variables
my @fastq_urls; # List of Fastq URLs
my @fastq_basenames;
my @fastq_fullpaths;
my $phyloFlash = "$FindBin::RealBin/phyloFlash.pl"; # assume phyloFlash script in same folder
my $wget = "wget --no-verbose";   # Wget binary and options

## MAIN ########################################################################

welcome();

@fastq_urls = get_fastq_urls($acc);

if (@fastq_urls) {
    msg ("Found the following Fastq files associated with accession $acc");
    msg (join ("\n", @fastq_urls));
}

# Check download path to make sure is a directory and is writable
if (defined $download_dir) {
    if (! -d $download_dir) {
        err ("Path to place downloaded files $download_dir is not a valid directory");
    } elsif (! -w $download_dir) {
        err ("You do not have write permissions to path $download_dir");
    } else {
        msg ("Downloading read files to folder $download_dir");
        $wget .= " -P $download_dir"; # Add path option to wget command
    }
}

# Check if any URLs were retrieved
if (defined $fastq_urls[0]) {
    # Open log file to record details on this file
    # Header line for log file
    foreach my $fastq (@fastq_urls) {
        msg ("Running wget with command: $wget $fastq");
        system (join " ", ($wget, $fastq)) if ($download == 1);
        # Strip URL dirs from filename, save to an array
        my $filename = $1 if $fastq =~ m/([^\/]+)$/;
        push @fastq_basenames, $filename;
    }
}

# Run phyloFlash (v3.0beta1+) to extract SSU reads
my @pF_args = ("-lib pF_$acc",
               "-CPUs $CPUs",
               "-emirge",               # Run both SPAdes and EMIRGE
               "-zip",                  # Output to archive
               "-log",
               );

# If specific dbhome specified, pass to phyloFlash, otherwise use default detected
push @pF_args, "-dbhome $pF_dbhome" if defined $pF_dbhome;

# Add paths to Fastq file names

foreach my $file (@fastq_basenames) {
    push @fastq_fullpaths, "$download_dir/$file";
}

# Check how many Fastq files (paired or unpaired)
if (scalar (@fastq_fullpaths) == 1) {
    push @pF_args, "-read1 ".$fastq_fullpaths[0];
} elsif (scalar (@fastq_basenames) == 2) {
    push @pF_args, "-read1 ".$fastq_fullpaths[0];
    push @pF_args, "-read2 ".$fastq_fullpaths[1];
}

if (defined $run_pF) {
    my $pF_cmd = join " ", ($phyloFlash, @pF_args);
    msg ("Running phyloFlash with command: $pF_cmd");
    system (join " ", ($phyloFlash, @pF_args));
}

if (defined $cleanup) {
    msg ("Deleting read files");
    unlink @fastq_fullpaths if defined $cleanup; # Delete read files 
}

msg ("Thank you for using ENA_phyloFlash.pl"); # Be nice

## SUBS ########################################################################

sub welcome {
    my $message = "This is ENA_phyloFlash.pl from phyloFlash v$VERSION";
    print STDERR $message;
    print STDERR "\n";
}

sub get_fastq_urls {
    # Get the URL(s) of ENA Fastq file(s) for a given ENA entry

    # Input: Accession number of sample or run
    my ($acc) = @_;
    my @urls_arr;    # Output array containing URLs
    # Get report table using ENA REST query
    my $url = "http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=$acc&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes";
    msg ("Sending REST query to ENA at URL: $url");
    
    my $ua = LWP::UserAgent->new;
    if (defined $http_proxy) {
        # Use proxies if defined
        # Note that https proxy is also with http protocol;
        # see https://stackoverflow.com/questions/27852041/perl-lwpuseragent-https-proxy-to-specific-https-web-site-unknown-protocol-er
        $ua->proxy('http',$http_proxy);
        $ua->proxy('https',$http_proxy);
    }
    
    my $response = $ua->get($url);
    err ("Cannot get $url: ". $response->status_line) unless $response->is_success;
    my $tab = $response->content();
    
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


=head1 COPYRIGHT AND LICENSE

Copyright (C) 2017 - Brandon Seah <kbseah@mpi-bremen.de>

LICENCE

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

=cut
