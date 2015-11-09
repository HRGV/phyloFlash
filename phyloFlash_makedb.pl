#!/usr/bin/perl -w
=pod

=head1 NAME

phyloFlash_makedb.pl - prepares the phyloFlash dbdir

=head1 SYNOPSIS

### download databases via FTP

phyloFlash_makedb.pl --remote

### use local copies

phyloFlash makedb.pl --silva_file F<path/to/silva_db> --univec_file F<path/to/univec_db> 

=head1 DESCRIPTION

=over 3

=item -

Assumes the script is executed in a dir with write access and

=item -

Dependencies: wget, vsearch, bbmap and Emirge

=item -

This script generates approx. 5Gbytes of extracted databases

=back

=head1 ARGUMENTS

=over 3

=item --remote

Download databases via FTP

=item --silva_file F<path/to/silva_db>

Path to local copy of SILVA database file. Ignored if --remote flag is used.

=item --univec_file F<path/to/univec_db>

Path to local copy of Univec database file. Ignored if --remote flag is used.

=item --keep

Do not delete intermediary files

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2014- by Harald Gruber-Vodicka <hgruber@mpi-bremen.de>
                       Elmar Pruesse <elmar.pruesse@ucdenver.edu>

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

use strict;
use FindBin;
use lib $FindBin::Bin;
use PhyloFlash;
use Pod::Usage;
use Getopt::Long;
use Net::FTP;
use Digest::MD5;
use IO::Uncompress::AnyUncompress qw(anyuncompress $AnyUncompressError);
use Cwd;

# URLS
my $silva_url  = "ftp.arb-silva.de/current/Exports/*_SSURef_Nr99_tax_silva_trunc.fasta.gz";
my $univec_url = "ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec";

# constants
my $dbsource;
my $silva_file = "";
my $univec_file = "";
my $cwd = getcwd;
my $timer = new Timer;
my $cpus = get_cpus      # num cpus to use
my $keep = 0;

# commandline arguments
my $use_remote = 0;     # Download SILVA and Univec databases via FTP

# parse cmdline

GetOptions("remote|r" => \$use_remote,
           "silva_file=s" => \$silva_file,
           "univec_file=s" => \$univec_file,
           "CPUs=i" => \$cpus,  # override default CPU usage if necessary
           "keep|k" => \$keep,
           "help|h" => sub{pod2usage(0)})
or pod2usage(2);

if ($use_remote==0 && ($silva_file eq "" | $univec_file eq "")) {
    pod2usage("Please specify paths to local copies of SILVA and Univec databases.");
}

# binaries needed by phyloFlash
require_tools((
    barrnapHGV => "$FindBin::RealBin/barrnap-HGV/bin/barrnap_HGV",
    grep => "grep",
    bbmask => "bbmask.sh",
    bbduk => "bbduk.sh",
    bbmap => "bbmap.sh",
    bowtiebuild => "bowtie-build",
    ));

check_environment();

### MAIN ###

# check whether user has local copy of Univec and Silva files, or if
# they should be downloaded
if ($use_remote==1 && $univec_file eq "") {
    msg("downloading latest univec from ncbi");
    $univec_file = file_download($univec_url);
}
elsif ($use_remote==0 && $univec_file ne "") {
    msg("using local copy of univec: $univec_file");
}
else {
    err("No univec file found and downloading disabled.");
}

if ($use_remote==1 && $silva_file eq "") {
    msg("downloading latest SSU RefNR from www.arb-silva.de");
    $silva_file  = file_download($silva_url);
}
elsif ($use_remote==0 && $silva_file ne "") {
    msg("using local copy of Silva SSU RefNR: $silva_file");
}
else {
    err("No SILVA database found and downloading disabled.");
}

# extract SILVA version
my ($silva_release) = ($silva_file =~ m/SILVA_([^_]+)_/);

if (!$&) {
    err("Unable to extract version from SILVA database filename:",
        "Expected 'SILVA_<version>_...' in '$silva_file'.");
}

# create database directory
my $dbdir = "./".$silva_release."/";
if (! -e $dbdir) {
    mkdir $dbdir
        or err("Failed to create dir $silva_release", $!);
}

msg("unpacking SILVA database");
anyuncompress $silva_file => "$dbdir/SILVA_SSU.fasta"
    or err("unpacking failed:  $AnyUncompressError");

my @lsu_in_ssh = find_LSU("$dbdir/SILVA_SSU.fasta");


fasta_copy_except("$dbdir/SILVA_SSU.fasta",
                  "$dbdir/SILVA_SSU.noLSU.fasta",
                  @lsu_in_ssh);
unlink "$dbdir/SILVA_SSU.fasta" unless ($keep==1);

mask_repeats("$dbdir/SILVA_SSU.noLSU.fasta",
             "$dbdir/SILVA_SSU.noLSU.masked.fasta");
unlink "$dbdir/SILVA_SSU.noLSU.fasta" unless ($keep==1);

univec_trim($univec_file,
            "$dbdir/SILVA_SSU.noLSU.masked.fasta",
            "$dbdir/SILVA_SSU.noLSU.masked.trimmed.fasta");
unlink "$dbdir/SILVA_SSU.noLSU.masked.fasta" unless ($keep==1);

#the cleaned, masked and trimmed databases are clustered
# 1) at 99id with full labels for bbmap
# 2) at 96id for emirge

cluster("./$silva_release/SILVA_SSU.noLSU.masked.trimmed.fasta",
           "./$silva_release/SILVA_SSU.noLSU.masked.trimmed.NR99.fasta",
           "0.99",
           $cpus);
# no unlink -> trimmed.fasta needed for vsearch


fasta_copy_iupac_randomize("./$silva_release/SILVA_SSU.noLSU.masked.trimmed.NR99.fasta",
             "./$silva_release/SILVA_SSU.noLSU.masked.trimmed.NR99.fixed.fasta");


bbmap_db("./$silva_release/SILVA_SSU.noLSU.masked.trimmed.NR99.fixed.fasta", "./$silva_release/");


cluster("./$silva_release/SILVA_SSU.noLSU.masked.trimmed.NR99.fasta",
           "./$silva_release/SILVA_SSU.noLSU.masked.trimmed.NR96.fasta",
           "0.96",
           $cpus);
unlink "$dbdir/SILVA_SSU.noLSU.masked.trimmed.NR99.fasta" unless ($keep==1);

fasta_copy_iupac_randomize("./$silva_release/SILVA_SSU.noLSU.masked.trimmed.NR96.fasta",
              "./$silva_release/SILVA_SSU.noLSU.masked.trimmed.NR96.fixed.fasta");
unlink "$dbdir/SILVA_SSU.noLSU.masked.trimmed.NR96.fasta" unless ($keep==1);

bowtie_index("./$silva_release/SILVA_SSU.noLSU.masked.trimmed.NR96.fixed.fasta");

finish();

### subroutines ###

# searching for LSU genes in the SSU RefNR using a modified b
# barrnap version that only utilizes LSU hmm profiles
# barrnap was changed to use a different db folder that
# only contains the LSU profiles
sub find_LSU {
    msg("searching for LSU contamination in SSU RefNR");
    my @lsus;

    my $fh;
    foreach ('bac', 'arch', 'euk') {
        my $log = "tmp.barrnap_hits.$_.barrnap.out";
        my $res = "tmp.barrnap_hits.$_.gff";
        #1 or run_prog("barrnapHGV",
        run_prog("barrnapHGV",
                 "  --kingdom $_ "
                 . "--threads $cpus "
                 . "--evalue 1e-10 "
                 . "--gene lsu "
                 . "--reject 0.01 "
                 . "./$silva_release/SILVA_SSU.fasta ",
                 $res, $log);
        open_or_die(\$fh, "-|",
                    "grep -E '23S_rRNA\|28S_rRNA' $res");
        while (my $row = <$fh>) {
            $row =~ s/\t.*//;
            chomp($row);
            push @lsus, $row;
        }
    }
    return @lsus;
}

#bbmask masks low entropy regions and repeats in the fasta file
sub mask_repeats {
    my ($src, $dst) = @_;
    msg("masking low entropy regions in SSU RefNR");
    run_prog("bbmask",
             "  overwrite=t "
             . "-Xmx4g "
             . "threads=$cpus "
             . "in=$src "
             . "out=$dst "
             . "minkr=4 maxkr=8 mr=t minlen=20 minke=4 maxke=8 "
             . "fastawrap=0 ");
}

#db is screened against UniVec and hit bases are trimmed using bbduk
sub univec_trim {
    my ($univec, $src, $dst) = @_;
    msg("removing UniVec contamination in SSU RefNR");
    run_prog("bbduk",
             "ref=$univec "
             . "  overwrite=t "
             . "-Xmx4g "
             . "threads=$cpus "
             . "fastawrap=0 "
             . "ktrim=r ow=t minlength=800 mink=11 hdist=1 "
             . "in=$src "
             . "out=$dst "
             . "stats=$dst.UniVec_contamination_stats.txt");
}

sub iuppac_replace {
    my ($src, $dst) = @_;
    msg("replacing IUPAC coded ambiguous bases with randomly chosen bases");
    run_prog("fixchars",
             "<$src",
             $dst);
}

#the NR99 file is fixed and a reference file is generated for bbmap
sub bbmap_db {
    my ($ref, $path) = @_;
    msg("creating bbmap reference");
    run_prog("bbmap",
             "  -Xmx10g "   # Original 4 Gb limit was not enough
             . "threads=$cpus "
             . "ref=$ref "
             . "path=$path ");
}

#the NR96 file is fixed and a bowtie index is built
sub bowtie_index {
    my ($fasta) = @_;
    my $btidx = $fasta =~ s/\.fasta$/\.bt/r;
    msg("creating bowtie index (for emirge)");
    run_prog("bowtiebuild", "$fasta $btidx -q");
}

#--------------------------------------------- final timing stats and goodbye
sub finish {
  msg("Walltime used: ".$timer->minutes);;
  msg("The databases for Silva release $silva_release are ready for phyloFlash

    Please provide your location of
    the databases with -dbhome: /path/to/your/databases/
    or change script line XXX accordingly");#Fixme
}
