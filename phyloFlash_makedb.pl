#!/usr/bin/env perl
use strict;
use warnings;
=pod

=head1 NAME

phyloFlash_makedb.pl - prepares the phyloFlash dbdir

=head1 SYNOPSIS

### download databases via FTP

phyloFlash_makedb.pl --remote

### use local copies

phyloFlash_makedb.pl --silva_file F<path/to/silva_db> --univec_file F<path/to/univec_db>

## Get help

phyloFlash_makedb.pl --help

phyloFlash_makedb.pl --man

=head1 DESCRIPTION

=over 3

=item -

Assumes the script is executed in a dir with write access and

=item -

Dependencies: wget, vsearch, bbmap and Emirge

=item -

This script generates approx. 5 Gbytes of extracted databases (9 Gb if Vsearch
UDB database is also produced, with vsearch v2.5.0+)

=back

=head1 ARGUMENTS

=head2 INPUT FILES

=over 8

=item --remote

Download databases via FTP

=item --silva_file F<path/to/silva_db>

Path to local copy of SILVA database file. Ignored if --remote flag is used.

This should be the Fasta-formatted SILVA SSURef file, clustered at 99% identity,
with SILVA taxonomy strings in file header, and sequences truncated to SSU gene
boundaries. The file name should be in the form
I<SILVA_[Release]_SSURef_Nr99_tax_silva_trunc.fasta.gz>

=item --univec_file F<path/to/univec_db>

Path to local copy of Univec database file. Ignored if --remote flag is used.

=back

=head2 OPTIONAL TOOLS

=over 8

=item --emirge

Index database with Bowtie v1 for Emirge. Requires I<bowtie-build> to be in path.

Default: Yes (turn off with --noemirge)

=item --sortmerna

Index database for Sortmerna, if you wish to use it as an alternative to BBmap
for extracting rRNA reads from the read file. Requires I<indexdb_rna> to be in
path.

Default: No (--nosortmerna).

=back

=head2 CONFIGURATION AND HELP

=over 8

=item --keep

Do not delete intermediary files

=item --overwrite

Overwrite files if already present. Files are not overwritten by default,
allowing you to restart the DB build process if it was interrupted (but you will
have to do find and delete corrupted files manually).

Default: No ("--nooverwrite")

=item --CPUs I<N>

Number of processors to use

Default: All available

=item --mem I<N>

Memory limit (in Gb) for indexing tools. At least 10 is recommended.

Default: 10

=item --log I<FILE>

Write phyloFlash_makedb.pl log to a file.

Default: None

=item --check_env

Check that required dependencies are available in your path. If specifying
optional tools I<--sortmerna> and I<--emirge>, put the I<--check_env> argument
at the end of the command.

=item --help

Help message

=item --man

Full manual page

=item --version

Report version

=back

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2014-2018 by Harald Gruber-Vodicka <hgruber@mpi-bremen.de>
                           Elmar Pruesse <elmar.pruesse@ucdenver.edu>
                           Brandon Seah <kbseah@mpi-bremen.de>

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

use FindBin;
use lib $FindBin::RealBin;
use PhyloFlash;
use Pod::Usage;
use Getopt::Long;
use Net::FTP;
use Digest::MD5;
use IO::Uncompress::AnyUncompress qw(anyuncompress $AnyUncompressError);
use Cwd;
use Storable;
use File::Spec;

# URLS
my $silva_url  = "ftp.arb-silva.de/current/Exports/*_SSURef_N?99_tax_silva_trunc.fasta.gz";
my $univec_url = "ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec";

# constants
my $dbsource;
my $silva_file = "";
my $univec_file = "";
my $cwd = getcwd;
my $timer = new Timer;
my $cpus = get_cpus     # num cpus to use
my $keep = 0;           # Keep temp files
my $sortmerna = 0;      # Index sortmerna tool
my $emirge = 1;         # Index emirge
my $memlimitGb = 10;    # In Gb
my $makedb_log;

# commandline arguments
my $use_remote = 0;     # Download SILVA and Univec databases via FTP
my $overwrite = 0;

# parse cmdline
pod2usage(-verbose=>0) if !@ARGV;
GetOptions("remote|r" => \$use_remote,
           "silva_file=s" => \$silva_file,
           "univec_file=s" => \$univec_file,
           "emirge!" => \$emirge,
           "sortmerna!" => \$sortmerna,
           "CPUs=i" => \$cpus,  # override default CPU usage if necessary
           "mem=i" => \$memlimitGb,
           "check_env" => sub { process_dependencies();
                                check_environment();
                                exit; },
           "keep|k" => \$keep,
           "overwrite!" => \$overwrite,
           "log=s" => \$makedb_log,
           "help|h" => sub{pod2usage(-verbose=>1);},
           "man|m" => sub {pod2usage(-verbose=>2);},
           "version" => sub { welcome();
                              exit(); },
           )
or pod2usage(-verbose=>0,-exit=>2);

if ($use_remote==0 && ($silva_file eq "" | $univec_file eq "")) {
    pod2usage("Please specify paths to local copies of SILVA and Univec databases.");
}

process_dependencies();
check_environment();

## MAIN ########################################################################

welcome();

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
if (! -e "$dbdir/SILVA_SSU.fasta" || $overwrite == 1 ) {
    anyuncompress $silva_file => "$dbdir/SILVA_SSU.fasta"
        or err("unpacking failed:  $AnyUncompressError");
} else {
    msg ("File $dbdir/SILVA_SSU.fasta exists, not overwriting");
}

my @lsu_in_ssh;
if (! -e "$dbdir/SILVA_SSU.noLSU.fasta" || $overwrite == 1) {
    @lsu_in_ssh = find_LSU("$dbdir/SILVA_SSU.noLSU.fasta");
    msg ("Removing sequences with potential LSU contamination");
    fasta_copy_except("$dbdir/SILVA_SSU.fasta",
                      "$dbdir/SILVA_SSU.noLSU.fasta",
                      @lsu_in_ssh,
                      $overwrite);
    unlink "$dbdir/SILVA_SSU.fasta" unless ($keep==1);
    unlink glob "$dbdir/tmp.barrnap_hits.*" unless ($keep==1);
} else {
    msg ("LSU-filtered file found, not overwriting");
}

mask_repeats("$dbdir/SILVA_SSU.noLSU.fasta",
             "$dbdir/SILVA_SSU.noLSU.masked.fasta",
             $overwrite);
unlink "$dbdir/SILVA_SSU.noLSU.fasta" unless ($keep==1);

univec_trim($univec_file,
            "$dbdir/SILVA_SSU.noLSU.masked.fasta",
            "$dbdir/SILVA_SSU.noLSU.masked.trimmed.fasta",
            $overwrite);
unlink "$dbdir/SILVA_SSU.noLSU.masked.fasta" unless ($keep==1);

# Index database into UDB file, if Vsearch v2.5.0+
# Speeds up run time in search phase of phyloFlash as db can be directly read to mem
my $vsearch_ver_check = check_vsearch_version();
if (defined $vsearch_ver_check) {
    my $fasta_in = "$dbdir/SILVA_SSU.noLSU.masked.trimmed.fasta";
    my $udb_out = "$dbdir/SILVA_SSU.noLSU.masked.trimmed.udb";
    make_vsearch_udb($fasta_in,
                     $udb_out,
                     $overwrite);
}

#the cleaned, masked and trimmed databases are clustered
# 1) at 99id with full labels for bbmap
# 2) at 96id for emirge and sortmerna

cluster("./$silva_release/SILVA_SSU.noLSU.masked.trimmed.fasta",
        "./$silva_release/SILVA_SSU.noLSU.masked.trimmed.NR99.fasta",
        "0.99",
        $cpus,
        $overwrite);
# no unlink -> trimmed.fasta needed for vsearch

fasta_copy_iupac_randomize("./$silva_release/SILVA_SSU.noLSU.masked.trimmed.NR99.fasta",
                           "./$silva_release/SILVA_SSU.noLSU.masked.trimmed.NR99.fixed.fasta",
                           $overwrite);

bbmap_db("./$silva_release/SILVA_SSU.noLSU.masked.trimmed.NR99.fixed.fasta",
         "./$silva_release/",
         $overwrite);

cluster("./$silva_release/SILVA_SSU.noLSU.masked.trimmed.NR99.fasta",
        "./$silva_release/SILVA_SSU.noLSU.masked.trimmed.NR96.fasta",
        "0.96",
        $cpus,
        $overwrite);
unlink "$dbdir/SILVA_SSU.noLSU.masked.trimmed.NR99.fasta" unless ($keep==1);

fasta_copy_iupac_randomize("./$silva_release/SILVA_SSU.noLSU.masked.trimmed.NR96.fasta",
                           "./$silva_release/SILVA_SSU.noLSU.masked.trimmed.NR96.fixed.fasta",
                           $overwrite);
unlink "$dbdir/SILVA_SSU.noLSU.masked.trimmed.NR96.fasta" unless ($keep==1);

bowtie_index("./$silva_release/SILVA_SSU.noLSU.masked.trimmed.NR96.fixed.fasta") if $emirge == 1;

if ($sortmerna == 1) {
    sortmerna_index("./$silva_release/SILVA_SSU.noLSU.masked.trimmed.NR96.fixed.fasta",
                    $overwrite);
    hash_SILVA_acc_taxstrings_from_fasta ("./$silva_release/SILVA_SSU.noLSU.masked.trimmed.NR96.fixed.fasta",
                                          $overwrite);
}

finish();

write_logfile($makedb_log) if defined $makedb_log;


## SUBS ########################################################################

sub welcome {
    print STDERR "\nThis is phyloFlash_makedb.pl from phyloFlash.pl v$VERSION\n\n";
}

sub process_dependencies {
    # binaries needed by phyloFlash
    require_tools((
        barrnapHGV => "$FindBin::RealBin/barrnap-HGV/bin/barrnap_HGV",
        grep => "grep",
        bbmask => "bbmask.sh",
        bbduk => "bbduk.sh",
        bbmap => "bbmap.sh",
        vsearch => 'vsearch',
        ));
    # Binaries for optional tools
    if ($emirge == 1) {
        require_tools ((
            bowtiebuild => "bowtie-build",
        ));
    }
    if ($sortmerna == 1) {
        # If Sortmerna is required, add to dependencies list
        require_tools ((
            indexdb_rna => "indexdb_rna",
        ));
    }
}

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

sub make_vsearch_udb {
    my ($infile, $outfile, $overwrite) = @_;
    msg ("Indexing $infile to make UDB file $outfile with Vsearch");
    my $log = "tmp.vsearch_make_udb.log";
    if (! -e $outfile || $overwrite == 1) {
        my @vsearch_params = ("--threads $cpus",
                              '--notrunclabels', # Keep full header including taxstring
                              "--makeudb_usearch $infile",
                              "--output $outfile",
                              );
        run_prog("vsearch",
                 join (" ", @vsearch_params),
                 undef,
                 $log
                 );
    } else {
        msg ("WARNING: File $outfile already exists. Not overwriting");
    }
}

#bbmask masks low entropy regions and repeats in the fasta file
sub mask_repeats {
    my ($src, $dst, $overwrite) = @_;
    msg("masking low entropy regions in SSU RefNR");
    my $log = "tmp.bbmask_mask_repeats.log";
    if (! -e $dst || $overwrite == 1) {
        run_prog("bbmask",
             "  overwrite=t "
             . "-Xmx".$memlimitGb."g "
             . "threads=$cpus "
             . "in=$src "
             . "out=$dst "
             . "minkr=4 maxkr=8 mr=t minlen=20 minke=4 maxke=8 "
             . "fastawrap=0 ",
             undef,
             $log);
    } else {
        msg ("File $dst exists, not overwriting");
    }
}

#db is screened against UniVec and hit bases are trimmed using bbduk
sub univec_trim {
    my ($univec, $src, $dst,$overwrite) = @_;
    msg("removing UniVec contamination in SSU RefNR");
    my $log = "tmp.bbduk_remove_univec.log";
    if (! -e $dst || $overwrite == 1 ) {
        run_prog("bbduk",
                 "ref=$univec "
                 . "  overwrite=t "
                 . "-Xmx".$memlimitGb."g "
                 . "threads=$cpus "
                 . "fastawrap=0 "
                 . "ktrim=r ow=t minlength=800 mink=11 hdist=1 "
                 . "in=$src "
                 . "out=$dst "
                 . "stats=$dst.UniVec_contamination_stats.txt",
                 undef,
                 $log);
    } else {
        msg ("File $dst exists, not overwriting");
    }
}

#the NR99 file is fixed and a reference file is generated for bbmap
sub bbmap_db {
    my ($ref, $path, $overwrite) = @_;
    msg("creating bbmap reference");
    my $log = "tmp.bbmap_index.log";
    if (!-d "$path/ref" || $overwrite == 1) {
        run_prog("bbmap",
                 "  -Xmx".$memlimitGb."g "   # Original 4 Gb limit was not enough
                 . "threads=$cpus "
                 . "ref=$ref "
                 . "path=$path ",
                 undef,
                 $log);
    } else {
        msg ("WARNING: Folder $path exists, not overwriting");
    }
}

#the NR96 file is fixed and a bowtie index is built
sub bowtie_index {
    my ($fasta) = @_;
    my $btidx = $fasta =~ s/\.fasta$/\.bt/r;
    my $log = "tmp.bowtiebuild.log";
    msg("creating bowtie index (for emirge)");
    run_prog("bowtiebuild",
             "$fasta $btidx -q",
             undef,
             $log);
}

# Generate index for Sortmerna from the filtered fixed SSU clustered at 96% id
sub sortmerna_index {
    my ($fasta,$overwrite) = @_;
    my $memlimitMb = $memlimitGb * 1000;
    my $prefix = $fasta =~ s/\.fasta$//r;
    my $log = "tmp.indexdb_rna.log";
    msg ("creating sortmerna index");
    if (! -e "$prefix.bursttrie_0.dat" || $overwrite == 1) {
        my @indexdb_args = ("--ref $fasta,$prefix",
                            "-m $memlimitMb",
                            "-v");
        run_prog("indexdb_rna",
                 join (" ", @indexdb_args),
                 undef,
                 $log
                 )
    } else {
        msg ("WARNING: SortMeRNA indices for file prefix $prefix already exist, not overwriting");
    }
}

sub hash_SILVA_acc_taxstrings_from_fasta {
    # Hash of accession numbers and taxonomy strings from SILVA fasta headers
    # Store as a perl hash image with Storable
    # For later when wrangling sortmerna SAM output to bbmap-like format
    my ($fasta) = @_;
    my $prefix = $fasta =~ s/\.fasta$//r;
    msg ("Hashing SILVA accession numbers to taxonomy strings");
    if (! -e "$prefix.acc2taxstring.hashimage" || $overwrite == 1) {
        my %hash;
        open(my $fh, "<", $fasta) or die ("$!");
        while (my $line = <$fh>) {
            if ($line =~ m/^>(.+)/) {
                my $header = $1;
                chomp $header;
                my ($id, @taxsplit) = split / /, $header; # Split on first space
                $hash{$id} = join " ", @taxsplit; # Rejoin tax string if it has spaces
            }
        }
        close($fh);
        store \%hash, "$prefix.acc2taxstring.hashimage";
    } else {
        msg ("WARNING: File $prefix.acc2taxstring.hashimage already exists, not overwriting");
    }
}

#--------------------------------------------- final timing stats and goodbye
sub finish {
  my $dbdir_abs = File::Spec->rel2abs($dbdir); # Convert to absolute path
  msg("Walltime used: ".$timer->minutes);;
  msg("The databases for Silva release $silva_release are ready for phyloFlash

    When running phyloFlash, please provide the location of
    the databases with the following option:
      -dbhome $dbdir_abs

    or add the following line to your .bashrc or .bash_profile:
      export PHYLOFLASH_DBHOME=$dbdir_abs");
}
