#!/usr/bin/env perl
use strict;
use warnings;
=head1 NAME

phyloFlash - A script to rapidly estimate the phylogenetic composition of
an Illumina (meta)genomic dataset and reconstruct SSU rRNA genes

=head1 SYNOPSIS

B<phyloFlash.pl> [OPTIONS] -lib I<name> -read1 F<<file>> -read2 F<<file>>

B<phyloFlash.pl> -help

B<phyloFlash.pl> -man

B<phyloFlash.pl> -check_env

=head1 DESCRIPTION

This tool rapidly approximates the phylogenetic composition of a (meta)genomic
library by mapping reads to a reference SSU rRNA database, and reconstructing
full-length SSU rRNA sequences. The pipeline is intended for Illumina paired-end
HiSeq and MiSeq reads.

For more information, see the manual at L<http://hrgv.github.io/phyloFlash/>

=head1 ARGUMENTS

=head2 INPUT

=over 8

=item -lib I<NAME>

Library I<NAME> to use for output file. The name must be one word comprising
only letters, numbers and "_" or "-" (no whitespace or other punctuation).

=item -read1 F<FILE>

Forward reads in FASTA or FASTQ formats. May be compressed with Gzip (.gz
extension). If interleaved reads are provided, please use I<--interleaved> flag
in addition for paired-end processing.

=item -read2 F<FILE>

File containing reverse reads. If this option is omitted, B<phyloFlash>
will run in B<experimental> single-end mode.

=item -check_env

Invokes checking of working environment and dependencies without data input.
Use to test setup.

=item -help

Print brief help.

=item -man

Show manual.

=item -outfiles

Show detailed list of output and temporary files and exit.

=item -version

Report version

=back

=head2 LOCAL SETTINGS

=over 8

=item -CPUs I<N>

Set the number of threads to use. Defaults to all available CPU cores.

=item -crlf

Use CRLF as line terminator in CSV output (to become RFC4180 compliant).

=item -decimalcomma

Use decimal comma instead of decimal point to fix locale problems.

Default: Off

=item -dbhome F<DIR>

Directory containing phyloFlash reference databases.
Use F<phyloFlash_makedb.pl> to create an appropriate directory.

If not specified, phyloFlash will check for an environment variable
I<$PHYLOFLASH_DBHOME>, then look in the current dir, the home dir, and
the dir where the phyloFlash.pl script is located, for a suitable database dir.

=back

=head2 ANALYSIS TOOLS

=over 8

=item -emirge

Turn on EMIRGE reconstruction of SSU sequences

Default: Off ("-noemirge")

=item -skip_spades

Turn off SPAdes assembly of SSU sequences

Default: No (use SPAdes by default)

=item -sortmerna

Use Sortmerna instead of BBmap to extract SSU rRNA reads. Insert size and
%id to reference statistics will not be available.

Default: No

=back

=head2 INPUT AND ANALYSIS OPTIONS

=over 8

=item -interleaved

Interleaved readfile with R1 and R2 in a single file at read1

=item -readlength I<N>

Sets the expected readlength. If no value is supplied, auto-detect read length.
Must be within 50..500.

Default: -1 (auto-detect)

=item -readlimit I<N>

Limits processing to the first I<N> reads in each input file. Use this
for transcriptomes with a lot of rRNA reads (use values <1000000).

Default: unlimited

=item -amplimit I<N>

Sets the limit of SSU read pairs to switch from emirge.py to
emirge_amplicon.py. This feature is not reliable as emirge_amplicon.py
has been problematic to run (use values >100000).

Default: 500000

=item -id I<N>

Minimum allowed identity for BBmap read mapping process in %. Must be within
50..98. Set this to a lower value for very divergent taxa

Default: 70

=item -evalue_sortmerna I<N>

E-value cutoff to use for SortMeRNA (only if I<-sortmerna> option chosen). In
scientific exponent notation (see below)

Default: 1e-09

=item -clusterid I<N>

Identity threshold for clustering with vsearch in %.
Must be within 50..100.

Default: 97

=item -taxlevel I<N>

Level in the taxonomy string to summarize read counts per taxon.
Numeric and 1-based (i.e. "1" corresponds to "Domain").

Default: 4 ("Order")

=item -tophit

Report taxonomic summary based on the best hit per read only. Otherwise, the
consensus of all ambiguous mappings of a given read will be used to assign its
taxonomy. (This was the default behavior before v3.2)

Default: No (use consensus)

=item -maxinsert I<N>

Maximum insert size allowed for paired end read mapping. Must be within
0..1200.

Default: 1200

=item -sc

Turn on single cell MDA data mode for SPAdes assembly of SSU sequences

=item -trusted F<FILE>

User-supplied Fasta file of trusted contigs containing SSU rRNA sequences. The
SSU sequences will be extracted with Barrnap, and the input read files will be
screened against these extracted "trusted" SSU sequences

=item -poscov

Use Nhmmer to find positional coverage of reads across Barrnap's HMM model of
the 16S and 18S rRNA genes from a subsample of reads, as an estimate of
coverage evenness.

Default: Off ("-noposcov")

=item -everything

Turn on all the optional analyses and output options. Options without defaults
and any local settings must still be specified. Equivalent to "-emirge -poscov
-treemap -zip -log"

=item -almosteverything

Like I<-everything> except without running EMIRGE.

=back

=head2 OUTPUT OPTIONS

=over 8

=item -html

Generate output in HTML format.

Default: On. (Turn off with "-nohtml")

=item -treemap

Draw interactive treemap of taxonomic classification in html-formatted report.
This uses Google Visualization API, which requires an internet connection,
requires that you agree to their terms of service (see
L<https://developers.google.com/chart/terms>), and is not open source, although
it is free to use.

Default: Off ("-notreemap")

=item -zip

Compress output into a tar.gz archive file. Overridden by I<-almosteverything>
or I<-everything>.

Default: Off ("-nozip")

=item -keeptmp

Keep temporary/intermediate files

Default: No ("-nokeeptmp")

=item -log

Write status messages printed to STDERR also to a log file

Default: Off ("-nolog")

=back

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2014-2018 by Harald Gruber-Vodicka <hgruber@mpi-bremen.de>
                           Elmar A. Pruesse <elmar.pruesse@ucdenver.edu>
                           and Brandon Seah <kbseah@mpi-bremen.de>

LICENSE

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
use File::Basename;
use File::Temp;
use File::Path 'remove_tree';
use Storable;
use IPC::Cmd qw(can_run);
use Cwd;
#use Data::Dumper;

# change these to match your installation
my @dbhome_dirs = (".", $ENV{"HOME"}, $FindBin::RealBin);

# constants
my $version       = "phyloFlash v$VERSION";       # Current phyloFlash version
my $progname      = $FindBin::Script;               # Current script name
my $cwd           = getcwd;                         # Current working folder
my $progcmd       = join " ", ($progname, @ARGV) ; # How the script was called

# configuration variables / defaults
my $DBHOME      = undef;
my $readsf_full = undef;        # path to forward read file
my $readsr_full = undef;        # path to reverse read file
my $readsf      = undef;        # forward read filename, stripped of directory names
my $readsr      = undef;        # reverse read filename, stripped of directory names
my $SEmode      = 0;            # single ended mode
my $libraryNAME = undef;        # output basename
my $interleaved = 0;            # Flag - interleaved read data in read1 input (default = 0, no)
my $id          = 70;           # minimum %id for mapping with BBmap
my $evalue_sortmerna = 1e-09;   # E-value cutoff for sortmerna, ignored if -sortmerna not chosen
my $readlength  = -1 ;          # length of input reads, -1 for automatic read length detection
my $readlimit   = -1;           # max # of reads to use
my $amplimit    = 500000;       # number of SSU pairs at which to switch to emirge_amplicon
my $maxinsert   = 1200;         # max insert size for paired end read mapping
my $cpus        = get_cpus      # num cpus to use
my $clusterid   = 97;           # threshold for vsearch clustering
my $trusted_contigs;            # Filename for trusted contigs, if available

my $html_flag   = 1;            # generate HTML output? (default = 1, on)
my $treemap_flag = 0;           # generate interactive treemap (default = 0, off)
my $crlf        = 0;            # csv line terminator
my $decimalcomma= 0;            # Decimal separator (default = .)
my $use_sortmerna = 0;          # Flag - Use Sortmerna instead of BBmap? (Default = 0, no)
my $skip_emirge = 1;            # Flag - skip Emirge step? (default = 1, yes)
my $skip_spades = 0;            # Flag - skip SPAdes assembly? (default = 0, no)
my $poscov_flag = 0;            # Flag - use Nhmmer to estimate positional coverage? (Default = 0, no)
my $sc          = 0;            # Flag - single cell data? (default = 0, no)
my $zip         = 0;            # Flag - Compress output into archive? (default = 0, no)
my $save_log    = 0;            # Save STDERR messages to log (Default = 0, no)
my $keeptmp     = 0;            # Do not delete temporary files (Default = 0, do delete temporary files)
my $tophit_flag;                # Use "old" way of taxonomic summary, by taking top bbmap hit (no ambiguities) (Default: undef, no)
my @tools_list;                 # Array to store list of tools required
                                # (0 will be turned into "\n" in parsecmdline)
# default database names for EMIRGE and Vsearch
my $emirge_db   = "SILVA_SSU.noLSU.masked.trimmed.NR96.fixed";
my $vsearch_db  = "SILVA_SSU.noLSU.masked.trimmed.udb";
my $sortmerna_db = $emirge_db;

my $ins_used = "SE mode!"; # Report insert size used by EMIRGE

my %outfiles;           # Hash to keep track of output files

# variables for report generation
my $sam_fixed_href;     # Ref to hash storing initial mapping data (stored as array of refs to hashes keyed by SAM column name), keyed by RNAME, read orientation, QNAME
my $sam_fixed_aref;     # Ref to array storing initial mapping data (pointing to same refs as above) in original order of SAM input file
my %ssu_sam_mapstats;   # Statistics of mapping vs SSU database, parsed from from SAM file
# Fields for %ssu_sam_mapstats
#   total_reads         Total read (pairs) input
#   total_read_segments Total read segments input
#   ssu_fwd_seg_map     Total fwd read segments mapped
#   ssu_rev_seg_map     Total rev read segments mapped
#   ssu_tot_map         Total read segments mapped
#   ssu_tot_pair_map    Total read pairs mapped (one or both segments)
#   ssu_tot_pair_ratio  Ratio of read (pairs) mapping to SSU (one or both segments
#   ssu_tot_pair_ratio_pc   Above as percentage
#   assem_tot_map       Total segments mapping to full length assembled
#   ssu_unassem         Total segments not mapping to full length
#   assem_ratio         Ratio assembled/unassmebled
#   assem_ratio_pc      Above as percentage
#   spades_fwd_map      Fwd read segments mapping to SPAdes
#   spades_rev_map
#   spades_tot_map
#   emirge_fwd_map
#   emirge_rev_map
#   emirge_tot_map
#   trusted_fwd_map
#   trusted_rev_map
#   trusted_tot_map

# Hashes to keep count of NTUs
my %taxa_ambig_byread_hash;     # Hash of all taxonomy assignments per read
my $taxa_full_href;             # Summary of read counts for full-length taxonomy strings (for treemap) - hash ref
my $taxa_summary_href;          # Summary of read counts at specific taxonomic level (hash ref)
my $taxa_unassem_summary_href;  # Summary of read counts of UNASSEMBLED reads at specific taxonomic level (hash ref)

my $taxon_report_lvl = 4;       # Taxonomic level to report counts

my @ssuassem_results_sorted;    # Sorted list of SSU sequences for reporting
my @ssurecon_results_sorted;
my %ssufull_hash;

# mapping statistics parsed from BBmap output
my $readnr = 0;
my $readnr_pairs;
my $SSU_total_pairs;
my $SSU_ratio;
my $SSU_ratio_pc;

# Insert size stats from BBmap output
my $ins_me = 0;
my $ins_std = 0;

# Alpha-diversity statistics calculated from taxon counts
my $chao1 = 0;
my @xtons = (0,0,0);  # singleton, doubleton and tripleton count

my $runtime;

# display welcome message incl. version
sub welcome {
    print STDERR "\nThis is $version\n\n";
}

# checks whether a given directory is a valid dbhome
# (contains the necessary files for bbmap, emirge and vsearch)
sub check_dbhome {
    my $dbhome = shift;
    my @required_list = ('ref/genome/1/summary.txt',
                         $emirge_db.".fasta",
                         $vsearch_db);
    push @required_list, ("$sortmerna_db.bursttrie_0.dat","$sortmerna_db.acc2taxstring.hashimage") if ($use_sortmerna == 1);
    foreach (@required_list) {
        return "${dbhome}/$_" unless -r "${dbhome}/$_"
    }
    return "";
}

# searches @dbhome_dirs for the dbdir with the highest version
# returns "" if none found
sub find_dbhome {
    my @dirs;
    my @dbdirs;
    foreach (@dbhome_dirs) {
        push(@dirs, get_subdirs($_));
    }
    foreach (@dirs) {
        if (check_dbhome($_) eq "") {
            push(@dbdirs, $_);
        }
    }
    if (scalar @dbdirs > 0) {
        @dbdirs = version_sort(@dbdirs);
        return $dbdirs[0];
    }
    return "";
}

# Specify list of required tools. Run this subroutine AFTER parse_cmdline()
sub process_required_tools {
  # binaries needed by phyloFlash
    require_tools(
        bbmap => "bbmap.sh",
        reformat => "reformat.sh",
        barrnap => "$FindBin::RealBin/barrnap-HGV/bin/barrnap_HGV",
        vsearch => "vsearch",
        mafft => "mafft",
        fastaFromBed => "fastaFromBed",
        sed => "sed",
        grep => "grep",
        awk => "awk",
        cat => "cat",
        plotscript_SVG => "$FindBin::RealBin/phyloFlash_plotscript_svg.pl",
    );
    if ($skip_spades == 0) {
        require_tools(spades => "spades.py");
    }
    if ($skip_emirge == 0) {
        require_tools(emirge => "emirge.py",
                      emirge_amp => "emirge_amplicon.py",
                      emirge_rename_fasta => "emirge_rename_fasta.py");
    }
    if ($use_sortmerna == 1) {
        require_tools (sortmerna => "sortmerna",
                       indexdb_rna => "indexdb_rna");
    }
    # Check operating system to decide which nhmmer to use
    my $opsys = $^O;
    msg ("Current operating system $opsys");
    if ($opsys eq "darwin") {
        require_tools (nhmmer => "$FindBin::RealBin/barrnap-HGV/binaries/darwin/nhmmer");
    } elsif ($opsys eq "linux") {
        require_tools (nhmmer => "$FindBin::RealBin/barrnap-HGV/binaries/linux/nhmmer");
    } else {
        msg ("Could not determine operating system, assuming linux");
        require_tools (nhmmer => "$FindBin::RealBin/barrnap-HGV/binaries/linux/nhmmer");
    }
}

sub verify_dbhome {
    # verify database present
    if (defined($DBHOME)) {
        if (my $file = check_dbhome($DBHOME)) {
            pod2usage("\nBroken dbhome directory: missing file \"$file\"");
        }
    } elsif (defined $ENV{"PHYLOFLASH_DBHOME"}) {
        # dbhome defined as environment variable
        if (my $file = check_dbhome($ENV{"PHYLOFLASH_DBHOME"})) {
            pod2usage("\nBroken dbhome directory $ENV{'PHYLOFLASH_DBHOME'}: missing file \"$file\"");
        } else {
            $DBHOME = $ENV{"PHYLOFLASH_DBHOME"};
        }
    } else {
        $DBHOME = find_dbhome();
        pod2usage("Failed to find suitable DBHOME. (Searched \""
                  .join("\", \"",@dbhome_dirs)."\".)\nPlease provide a path using -dbhome. "
                  ."You can build a reference database using phyloFlash_makedb.pl\n")
            if ($DBHOME eq "");
    }
    msg("Using dbhome '$DBHOME'");
}

sub auto_detect_readlength {
    # Get read lengths from input files instead of asking user for input
    my ($fwd,  # Fwd read file
        $rev,  # Rev read file, undef if not available
        $limit # Stop after this many reads
        ) = @_;
    require_tools("readlength" => "readlength.sh");
    my @readlength_cmd = ("readlength",
                       "reads=$limit",
                       "bin=1",
                       "in=$fwd");
    if (defined $rev) {
        push @readlength_cmd, "in2=$rev";
    }
    # Use readlength.sh from bbmap suite to count read lengths from Fasta or Fastq
    run_prog("readlength",
             join(" ", @readlength_cmd),
             $outfiles{"readlength_out"}{"filename"},
             $outfiles{"readlength_err"}{"filename"});
    $outfiles{"readlength_out"}{"made"} = 1;
    # Hash counts of reads per read length
    my @lens;
    my $fh;
    open_or_die(\$fh, "<", $outfiles{"readlength_out"}{"filename"});
    while (my $line = <$fh>) {
        next if $line =~ m/^#/; # Skip headers
        my @split = split /\t/, $line;
        push @lens, $split[0];
    }
    # Sort by highest counts
    if (scalar @lens > 1) {
        msg("Reads are not all the same length, using longest read length found...");
    }
    my @sort_lens = sort { $b <=> $a } @lens;
    # Return the longest read length
    msg("Auto-detected read length: $sort_lens[0]");
    return ($sort_lens[0]);
}

# parse arguments passed on commandline and do some
# sanity checks
sub parse_cmdline {
    pod2usage(-verbose=>0) if !@ARGV;
    my $emirge = 0;
    my $everything = 0 ;
    my $almosteverything = 0;
    GetOptions('read1=s' => \$readsf_full,
               'read2=s' => \$readsr_full,
               'lib=s' => \$libraryNAME,
               'dbhome=s' => \$DBHOME,
               'interleaved' => \$interleaved,
               'readlength=i' => \$readlength,
               'readlimit=i' => \$readlimit,
               'amplimit=i' => \$amplimit,
               'maxinsert=i' => \$maxinsert,
               'id=i' => \$id,
               'clusterid=i' => \$clusterid,
               'taxlevel=i' => \$taxon_report_lvl,
               'CPUs=i' => \$cpus,
               'html!' => \$html_flag,
               'treemap!' => \$treemap_flag,
               'crlf' => \$crlf,
               'decimalcomma' => \$decimalcomma,
               'sortmerna!' => \$use_sortmerna,
               'evalue_sortmerna=s' => \$evalue_sortmerna,
               'emirge!' => \$emirge,
               'skip_emirge' => \$skip_emirge,
               'skip_spades' => \$skip_spades,
               'trusted=s' => \$trusted_contigs,
               'poscov!' => \$poscov_flag,
               'sc' => \$sc,
               'zip!' => \$zip,
               'log!' => \$save_log,
               'keeptmp!' => \$keeptmp,
               'everything' => \$everything,
               'almosteverything' => \$almosteverything,
               'tophit!' => \$tophit_flag,
               'check_env' => \&check_env_wrapper,
               'outfiles' => \&output_description,
               'help|h' => sub { pod2usage(-verbose=>0); },
               'man' => sub { pod2usage(-exitval=>0, -verbose=>2); },
               'version' => sub { welcome(); exit; },
               ) or pod2usage(2);
    $skip_emirge = 0 if $emirge == 1; # ain't gonna not be less careful with no double negatives

    # Say hello
    welcome();

    # Verify DBHOME is a valid phyloFlash DB folder
    verify_dbhome();

    # verify valid lib name
    pod2usage("Please specify output file basename with -lib\n")
        if !defined($libraryNAME);
    pod2usage("\nERROR: Argument to -lib may not be empty\n")
        if length($libraryNAME) == 0;
    pod2usage("\nERROR: Argument to -lib may contain only alphanumeric characters,"
              ."'_' and '-'.\n")
        if $libraryNAME =~ m/[^a-zA-Z0-9_-]/ ;

    msg("working on library $libraryNAME");

    # verify read files
    pod2usage("\nPlease specify input forward read file with -read1")
        if !defined($readsf_full);
    pod2usage("\nUnable to open forward read file '$readsf_full'. ".
              "Make sure the file exists and is readable by you.")
        if ! -e $readsf_full;

    # strip directory paths from forward read filename
    if ($readsf_full =~ m/.*\/(.+)/) {
      $readsf = $1;
    } else {
        $readsf = $readsf_full;
    }

    if (defined($readsr_full)) {
        pod2usage("\nUnable to open reverse read file '$readsr_full'.".
                  "Make sure the file exists and is readable by you.")
            if ! -e $readsr_full;
        pod2usage("\nForward and reverse input read file need to be different")
            if ($readsr_full eq $readsf_full);
        if ($readsr_full =~ m/.*\/(.+)/) {    # strip directory paths from reverse read filename
          $readsr = $1;
        } else { $readsr = $readsr_full; }

    } elsif ( $interleaved == 1 ){
    msg("Using interleaved read data");

    } else {
        $SEmode = 1; # no reverse reads, we operate in single ended mode
        $readsr = "";
        $readsr_full = "";
    }

    msg("Forward reads $readsf_full");
    if ($SEmode == 0)  {
        if (defined($readsr_full)) {
            msg("Reverse reads $readsr_full");
        } else{
            msg("Reverse reads from interleaved read file $readsf_full");
        }
    } else {
        msg("Running in single ended mode");
    }

    # populate hash to keep track of output files
    my $outfiles_href = initialize_outfiles_hash($libraryNAME,$readsf);
    %outfiles = %$outfiles_href;
    # use hash to keep track of output file description, filename,
    # whether it should be deleted at cleanup, and if file was actually created

    if ($readlength < 0) {
        msg("Automatic read length detection");
        $readlength = auto_detect_readlength($readsf_full, $readsr_full, 10000);
    }

    # check lengths
    pod2usage("\nReadlength must be within 50...500")
        if ($readlength < 50 or $readlength > 500);
    pod2usage("\nMaxinsert must be within 0..1200")
        if ($maxinsert < 0 or $maxinsert > 1200);
    pod2usage("\nReadmapping identity (-id) must be within 50..98")
        if ($id < 50 or $id > 98);
    pod2usage("\nClustering identity (-clusterid) must be within 50..100")
        if ($clusterid < 50 or $clusterid > 100);

    $crlf = $crlf ? "\r\n":"\n";

    # check CPUS
    if ($cpus eq "all" or $cpus < 0) {
        $cpus = get_cpus();
    }

    # Warn user about Sortmerna requirements
    if ($use_sortmerna == 1) {
        msg ("Sortmerna will be used instead of BBmap for the initial SSU rRNA read extraction");
        msg ("Sortmerna requires uncompressed input files - Ensure that you have enough disk space");
        msg ("The following input parameters specific to BBMap will be ignored: --id, --maxinsert");
    }

    # check surplus arguments
    err("Command line contains extra words:", @ARGV)
        if ($ARGV[0]);

    # If trusted contigs supplied, check that file is at least a valid ASCII/UTF-8 text file
    if (defined $trusted_contigs) {
        msg ("Trusted contigs file $trusted_contigs supplied");
        if (! -T $trusted_contigs) {
            msg ("WARNING: Trusted contigs file $trusted_contigs does not appear to be plain text file. Ignoring...");
            $trusted_contigs = undef;
        }
    }

    # Activate all optional outputs if "everything" is asked for
    if ($everything + $almosteverything > 0) {
        ($poscov_flag,
         $html_flag,
         $treemap_flag,
         $zip,
         $save_log
         ) = (1,1,1,1,1);
        $skip_spades = 0; # Override any skip spades
        if ($everything == 1) {
            msg ("Running \"everything\" - overrides other command line options");
            $skip_emirge = 0;
        } else {
            # "almost everything" does not turn on EMIRGE
            msg ("Running \"everything\" except EMIRGE- overrides other command line options");
        }
    }
}

sub check_env_wrapper {
    # verify tools present
    welcome();
    process_required_tools();
    check_environment(); # will die on failure
    verify_dbhome();
    exit();
}

sub output_description {
    my $example_outfiles = initialize_outfiles_hash("LIBNAME","READS");

    print STDERR "Description of output files \n\n";
    print STDERR "FILENAME\tDESCRIPTION\n";
    foreach my $outfile (sort {$a cmp $b} keys %$example_outfiles) {
        my @out = ($example_outfiles->{$outfile}{"filename"},
                   $example_outfiles->{$outfile}{"description"},);
        print STDERR join ("\t", @out)."\n";
    }
    exit;
}

sub print_report {
    ## Write plaintext report file

    msg("writing final files...");
    my $fh;
    open_or_die(\$fh, '>', $outfiles{"report"}{"filename"});
    $outfiles{"report"}{"made"} = 1;

    print {$fh} qq ~
$version - high throughput phylogenetic screening using SSU rRNA gene(s) abundance(s)

Command:\t$progcmd
Library name:\t$libraryNAME
---
Forward read file\t$readsf_full~;
if ($SEmode == 0) {
  if ($interleaved == 1) { # reads are provided interleaved
        print {$fh} qq~
Reverse reads provided interleaved in file\t$readsf_full~;}
  else { # reads are provided in two separate files
    print {$fh} qq ~
Reverse read file\t$readsr_full~;
  }
}

print {$fh} qq~
---
Current working directory\t$cwd
---
Minimum mapping identity:\t$id%
~;

    if (defined($readsr_full) || $interleaved == 1) { # If in PE mode
        print {$fh} qq~
---
Input PE-reads:\t$readnr_pairs
Mapped SSU reads:\t$SSU_total_pairs
Mapping ratio:\t$SSU_ratio_pc%
Detected median insert size:\t$ins_me
Used insert size:\t$ins_used
Insert size standard deviation:\t$ins_std
~;
    } else { # Else if in SE mode
        print {$fh} qq~
---
Input SE-reads:\t$readnr_pairs
Mapped SSU reads:\t$SSU_total_pairs
Mapping ratio:\t$SSU_ratio_pc%
~;
    }

    if (defined $ssu_sam_mapstats{"assem_ratio"}) {
        print {$fh} "Ratio of assembled SSU reads:\t".$ssu_sam_mapstats{"assem_ratio"}."\n";
    }
    my $chao1_3dp = $chao1 eq "n.d." ? $chao1 : sprintf ("%.3f", $chao1);

    print {$fh} qq~
---
Runtime:\t$runtime
CPUs used:\t$cpus
---
Read mapping based higher taxa (NTUs) detection

NTUs observed once:\t$xtons[0]
NTUs observed twice:\t$xtons[1]
NTUs observed three or more times:\t$xtons[2]
NTU Chao1 richness estimate:\t$chao1_3dp

List of NTUs in order of abundance (min. 3 reads mapped):
NTU\treads
~;

    # sort keys numerically descending in hash of
    # mapping-based detected higher taxa
    my @keys = sort {${$taxa_summary_href}{$b} <=> ${$taxa_summary_href}{$a}} keys %$taxa_summary_href;
    foreach my $key (@keys) {
        if (${$taxa_summary_href}{$key} >= 3) {
            my @out = ($key, ${$taxa_summary_href}{$key});
            print {$fh} join("\t", @out)."\n";
        }
    }

    if (defined $outfiles{"spades_fasta"}{"made"}) {
        ## Print the table of SSU assembly-based taxa to report file
        print {$fh} "---\n";
        print {$fh} "SSU assembly based taxa:\n";
        print {$fh} "OTU\tread_cov\tcoverage\tdbHit\ttaxonomy\t%id\talnlen\tevalue\n";
        foreach my $seqid ( sort { $ssufull_hash{$b} <=> $ssufull_hash{$a} } keys %ssufull_hash) {
            next unless (defined $ssufull_hash{$seqid}{"source"} && $ssufull_hash{$seqid}{"source"} eq "SPAdes");
            my @out;
            push @out, $seqid;
            my @fields = qw(counts cov dbHit taxon pcid alnlen evalue);
            foreach my $field (@fields) {
                push @out, $ssufull_hash{$seqid}{$field};
            }
            print {$fh} join("\t", @out)."\n";
        }
    }

    if (defined $outfiles{"emirge_fasta"}{"made"}) {
        ## Print the table of SSU reconstruction-based taxa to report file
        print {$fh} "---\n";
        print {$fh} "SSU reconstruction based taxa:\n";
        print {$fh} "OTU\tread_cov\tratio\tdbHit\ttaxonomy\t%id\talnlen\tevalue\n";
        foreach my $seqid (sort { $ssufull_hash{$b} <=> $ssufull_hash{$a} } keys %ssufull_hash) {
            next unless (defined $ssufull_hash{$seqid}{"source"} && $ssufull_hash{$seqid}{"source"} eq "EMIRGE");
            my @out;
            push @out, $seqid;
            my @fields = qw(counts cov dbHit taxon pcid alnlen evalue);
            foreach my $field (@fields) {
                push @out, $ssufull_hash{$seqid}{$field};
            }
            print {$fh} join("\t", @out)."\n";
        }
    }

    if (defined $outfiles{"trusted_fasta"}{"made"}) {
        ## Print the table of SSU extracted from trusted contigs to report file
        print {$fh} "---\n";
        print {$fh} "SSU from trusted contigs:\n";
        print {$fh} "OTU\tread_cov\tratio\tdbHit\ttaxonomy\t%id\talnlen\tevalue\n";
        foreach my $seqid (sort { $ssufull_hash{$b} <=> $ssufull_hash{$a} } keys %ssufull_hash) {
            next unless (defined $ssufull_hash{$seqid}{"source"} && $ssufull_hash{$seqid}{"source"} eq "trusted");
            my @out;
            push @out, $seqid;
            my @fields = qw(counts cov dbHit taxon pcid alnlen evalue);
            foreach my $field (@fields) {
                push @out, $ssufull_hash{$seqid}{$field};
            }
            print {$fh} join("\t", @out)."\n";
        }
    }

    if (defined $outfiles{"emirge_fasta"}{"made"} || defined $outfiles{"spades_fasta"}{"made"} || defined $outfiles{"trusted_fasta"}{"made"}) {
        ## Print the table of taxonomic affiliations for unassembled SSU reads
        print {$fh} "---\n";
        print {$fh} "Taxonomic affiliation of unassembled reads (min. 3 reads mapped):\n";
        my @taxsort = sort {${$taxa_unassem_summary_href}{$b} <=> ${$taxa_unassem_summary_href}{$a}} keys %$taxa_unassem_summary_href;
        foreach my $uatax (@taxsort) {
            if (${$taxa_unassem_summary_href}{$uatax} >= 3) {
                my @out = ($uatax, ${$taxa_unassem_summary_href}{$uatax});
                print {$fh} join ("\t", @out)."\n";
            }
        }
    }

    close($fh);
}

sub write_csv {
    msg("exporting results to csv");

    # CSV file of phyloFlash run metadata
    my @report = csv_escape((
        "version",$version,
        "library name",$libraryNAME,
        "forward read file",$readsf_full,
        "reverse read file",$readsr_full,
        "cwd",$cwd,
        "minimum mapping identity",$id,
        "single ended mode",$SEmode,
        "input reads",$ssu_sam_mapstats{'total_reads'},
        "input read segments",$ssu_sam_mapstats{'total_read_segments'},
        "mapped SSU reads",$SSU_total_pairs,
        "mapped SSU read segments",$ssu_sam_mapstats{'ssu_tot_map'},
        "mapping ratio",$SSU_ratio_pc,
        "detected median insert size",$ins_me,
        "used insert size",$ins_used,
        "insert size stddev",$ins_std,
        "runtime",$runtime,
        "CPUs",$cpus,
        "NTUs observed once",$xtons[0],
        "NTUs observed twice",$xtons[1],
        "NTUs observed three or more times",$xtons[2],
        "NTU Chao1 richness estimate", $chao1 eq "n.d." ? $chao1 : sprintf ("%.3f", $chao1), # Round to 3 d.p.
        "program command",$progcmd,
        "database path",$DBHOME,
    ));
    my $fh;
    open_or_die(\$fh, ">", $outfiles{"report_csv"}{"filename"});
    $outfiles{"report_csv"}{"made"} = 1;
    while ($#report > 0) {
        print {$fh} shift(@report).",".shift(@report).$crlf;
    }
    close($fh);

    # CSV file of taxonomic units (NTUs) from mapping vs SILVA database
    # truncated to requested taxonomic level
    open_or_die(\$fh, ">", $outfiles{"ntu_csv"}{"filename"});
    $outfiles{"ntu_csv"}{"made"} = 1;
    # Sort results descending
    my @keys = sort {${$taxa_summary_href}{$b} <=> ${$taxa_summary_href}{$a}} keys %$taxa_summary_href;
    foreach my $key (@keys) {
        my @out = ($key, ${$taxa_summary_href}{$key});
        print {$fh} join(",",csv_escape(@out)).$crlf;
    }
    close($fh);

    if (!defined $tophit_flag) {
        #CSV file of taxonomic units (NTUs) from mapping vs SILVA database,
        # not truncated to requested taxonomic level
        open_or_die(\$fh, ">", $outfiles{"ntu_full_csv"}{"filename"});
        $outfiles{"ntu_full_csv"}{"made"} = 1;
        # Sort results descending
        my @keys_full = sort {${$taxa_full_href}{$b} <=> ${$taxa_full_href}{$a}} keys %$taxa_full_href;
        foreach my $key (@keys_full) {
            my @out = ($key, ${$taxa_full_href}{$key});
            print {$fh} join(",",csv_escape(@out)).$crlf;
        }
        close($fh);
    }

    # If full-length seqeunces were assembled or reconstructed
    if (defined $outfiles{"spades_fasta"}{"made"} || defined $outfiles{"emirge_fasta"}{"made"} || defined $outfiles{"trusted_fasta"}{"made"}) {

        # CSV file of assembled/reconstructed sequnces
        open_or_die(\$fh, ">", $outfiles{"full_len_class"}{"filename"});
        $outfiles{"full_len_class"}{"made"}++;
        print $fh "OTU,read_cov,coverage,dbHit,taxonomy,%id,alnlen,evalue\n"; # Header

        # Sort descending by read counts
        foreach my $seqid (sort {$ssufull_hash{$b}{"counts"} <=> $ssufull_hash{$a}{"counts"}} keys %ssufull_hash) {
            my @out;
            push @out, $seqid;
            foreach my $field (qw(counts cov dbHit taxon pcid alnlen evalue)) {
                if (defined $ssufull_hash{$seqid}{$field}) {
                    push @out, $ssufull_hash{$seqid}{$field};
                } else {
                    push @out, 'NA';
                }
            }
            print {$fh} join (",", csv_escape(@out)).$crlf;
        }
        close($fh);

        # CSV file of taxonomic affiliations for unassembled/unreconstructed reads
        my $fh2;
        open_or_die (\$fh2, ">", $outfiles{"unassem_csv"}{"filename"});
        my @taxsort = sort {${$taxa_unassem_summary_href}{$b} <=> ${$taxa_unassem_summary_href}{$a}} keys %$taxa_unassem_summary_href;
        foreach my $uatax (@taxsort) {
            my @out = ($uatax, ${$taxa_unassem_summary_href}{$uatax});
            print {$fh2} join (",", csv_escape(@out)).$crlf;
        }
        $outfiles{"unassem_csv"}{"made"}++;
        close ($fh2);
    }
}

sub bbmap_fast_filter_sam_run {
    # input: $readsf, $readsr
    # output: $libraryNAME.$readsf.SSU.sam
    #         $libraryNAME.$readsf.SSU.{1,2}.fq
    #         $libraryNAME.bbmap.out
    #         $libraryNAME.inserthistogram
    # tmp:    $libraryNAME.basecompositionhistogram
    msg("filtering reads with SSU db using minimum identity of $id%");
    if ($readlimit != -1) {
        msg("Only using the first $readlimit reads");
    }
    # Minimum mapping ID of 50%
    my $minID= $id / 100;
    if ($minID < 0.50) {
        $minID = 0.50
    }
    # Set up input arguments for BBmap
    my @bbmap_args = ("fast=t",
                      "minidentity=$minID",
                      "-Xmx20g reads=$readlimit",
                      "threads=$cpus",
                      "po=f",
                      "outputunmapped=f",
                      "path=$DBHOME",
                      "out=".$outfiles{"bbmap_sam"}{"filename"},
                      "outm=".$outfiles{"reads_mapped_f"}{"filename"},
                      'noheader=t', # Do not print header lines of SAM file
                      'ambiguous=all', # Report all ambiguous reads
                      "build=1",
                      "in=$readsf_full",
                      "bhist=".$outfiles{"basecompositionhist"}{"filename"},
                      "ihist=".$outfiles{"inserthistogram"}{"filename"},
                      "idhist=".$outfiles{"idhistogram"}{"filename"},
                      "scafstats=".$outfiles{"hitstats"}{"filename"},
                      "overwrite=t", # Overwrite existing files
                      );
    # Additional input arguments for paired-end input
    if ($SEmode == 0) {
        if ($interleaved == 1) {
            push @bbmap_args, ("outm2=".$outfiles{"reads_mapped_r"}{"filename"},
                               "pairlen=$maxinsert interleaved=t"
                               );
        } else {
            push @bbmap_args, ("outm2=".$outfiles{"reads_mapped_r"}{"filename"},
                               "pairlen=$maxinsert",
                               "in2=$readsr_full",
                               );
        }
        $outfiles{"reads_mapped_r"}{"made"}++;
    }
    # Run BBmap
    run_prog("bbmap",
             join (" ", @bbmap_args),
             undef,
             $outfiles{"bbmap_log"}{"filename"}
             );
    # Record which files were created
    foreach my $madekey (qw(bbmap_sam reads_mapped_f basecompositionhist inserthistogram idhistogram hitstats bbmap_log)) {
        $outfiles{$madekey}{"made"}++;
    }
    msg("done...");
}

sub sortmerna_filter_sam {
    msg ("Using Sortmerna to extract SSU rRNA reads instead of BBmap");

    # Reformat Fastq file
    msg ("Reformatting input files to uncompressed Fastq format for Sortmerna");
    my @reformat_args = ("in=$readsf_full",
                         "reads=$readlimit",
                         "threads=$cpus",
                         "out=".$outfiles{'reads_uncompressed'}{'filename'},
                         );
    if ($SEmode == 0) {
        # Additionally for paired-end input
        if ($interleaved == 1) {
            push @reformat_args, ("interleaved=t");
        } else {
            push @reformat_args, ("in2=$readsr_full");
        }
    }
    run_prog("reformat",
             join (" ", @reformat_args),
             );

    # Record that file was made
    $outfiles{'reads_uncompressed'}{'made'} ++;

    # Sortmerna arguments
    my @sortmerna_args = ("--ref $DBHOME/$sortmerna_db.fasta,$DBHOME/$sortmerna_db",
                          "--reads ".$outfiles{'reads_uncompressed'}{'filename'},
                          "--aligned $libraryNAME.sortmerna", # Output file prefix
                          "--fastx",        # Report extracted reads in Fastq format
                          "--paired_in",    # Report both segments in pair to Fastq format
                          "--sam",          # SAM alignment output
                          #"--blast 1",     # Blast-tabular alignment output
                          "--log",
                          #"--best 10",      # Take 10 best hits
                          "--min_lis 10",
                          "-e $evalue_sortmerna",    # E-value cutoff
                          "-a $cpus",
                          "-v",             # Verbose (turn off?)
                          );
    if (defined $tophit_flag) {
        msg ("Reporting only single best hit per read");
        push @sortmerna_args, "--best 1";
    } else {
        msg ("Reporting 10 best hits per read");
        push @sortmerna_args, "--best 10";
    }

    msg ("Running sortmerna:");
    # Run Sortmerna
    run_prog("sortmerna",
             join (" ", @sortmerna_args),
             );

    # Record which files created
    foreach my $madekey (qw (sortmerna_sam sortmerna_fastq sortmerna_log)) {
        $outfiles{$madekey}{'made'} ++;
    }

    msg ("Fixing bitflags and RNAME in SAM file output of Sortmerna");
    # Fix Sortmerna output file - bitflags and names of ref sequences
    my ($fixed_sam_href,
        $fixed_sam_aref,
        $fixed_sam_lines_aref
        ) = fix_hash_sortmerna_sam($outfiles{'sortmerna_fastq'}{'filename'},
                                   $outfiles{'sortmerna_sam'}{'filename'},
                                   "$DBHOME/$sortmerna_db.acc2taxstring.hashimage",
                                   $SEmode
                                   );
    # Write file
    if (defined $fixed_sam_lines_aref) {
        my $sam_fh;
        open_or_die(\$sam_fh,
                    ">",
                    $outfiles{"sam_map"}{"filename"});
        foreach my $line (@$fixed_sam_lines_aref) {
            print $sam_fh "$line\n";
        }
        close ($sam_fh);
        $outfiles{"sam_map"}{"made"} ++ ;
    }

    # Reformat fwd and rev mapped Fastq reads from Sortmerna
    if ($SEmode == 0) {
        msg ("Reformatting interleaved mapped Fastq from sortmerna to split files");
        # Split paired reads into fwd and rev files
        my @split_paired_mapped_args = ("in=".$outfiles{'sortmerna_fastq'}{'filename'},
                                        'interleaved=t', # Force interleaved input
                                        "out=".$outfiles{"reads_mapped_f"}{"filename"},
                                        "out2=".$outfiles{"reads_mapped_r"}{"filename"},
                                        "threads=$cpus",
                                        );
        run_prog('reformat',
                 join (" ", @split_paired_mapped_args),
                 );
        $outfiles{"reads_mapped_f"}{'made'} ++;
        $outfiles{"reads_mapped_r"}{'made'} ++;
    } else {
        system (join " ", ('mv',
                           $outfiles{'sortmerna_fastq'}{'filename'},
                           $outfiles{'reads_mapped_f'}{'filename'}));
        $outfiles{'sortmerna_fastq'}{'made'} = undef;
        $outfiles{'reads_mapped_f'}{'made'} ++;
    }

    msg ("Done");
    return ($fixed_sam_href, $fixed_sam_aref);
}

sub parse_mapstats_from_sam_arr {
    my ($aref,      # Array of hrefs to SAM mapping
        $readnr,    # Total reads (segments) input to mapper
        $statshref, # Ref to hash of mapping statistics
        $SEmode     # Flag for single-end mode
        ) = @_;
    #$readnr,$readnr_pairs,$SSU_total_pairs,$SSU_ratio,$SSU_ratio_pc

    # Variables used internally
    my $ssu_pairs     = 0;  # Pairs with both segments mapping
    my $ssu_bad_pairs = 0;  # Pairs where segments map to different references
    my $ssu_f_reads   = 0;  # Forward segments mapping
    my $ssu_r_reads   = 0;  # Reverse segments mapping
    my $SSU_total_pairs = 0;# Total read pairs with at least one segment mapping
    my $mapped_half = 0;    # Pairs where one segment maps but the other doesn't
    my $unmapped_segment = 0; # Don't use this for summary
    my $unmapped;

    foreach my $href (@$aref) {
        # Go through each SAM record and check bitflags to see if mapping or not
        next if $href->{'FLAG'} & 0x100; # Skip secondary alignments

        if ($SEmode == 0) {
            if ($href->{'FLAG'} & 0x40) { # Fwd read
                $ssu_f_reads ++ unless $href->{'FLAG'} & 0x4;   # Unless fwd read unmapped
            } elsif ($href->{'FLAG'} & 0x80) { # Rev read
                $ssu_r_reads ++ unless $href->{'FLAG'} & 0x4;   # Unless rev read unmapped
            }
        } else {
            # Single end reads
            if ($href->{'FLAG'} == 0 || $href->{'FLAG'} & 0x10) {
                $ssu_f_reads ++;
            }
        }

        if ($href->{'FLAG'} & 0x2) { # Each segment properly aligned
            $ssu_pairs ++; # This will have to be divided by two later
        } elsif ($href->{'FLAG'} & 0x8) { # Next segment unmapped
            $mapped_half ++;
        } elsif ($href->{'FLAG'} & 0x4) { # Unmapped segment
            $unmapped_segment ++;
        } else { # Both segments map but to different references (bbmap 'bad pairs')
            $ssu_bad_pairs ++; # Will have to be divided by two later
        }
    }

    # If paired reads, divide segment-based counts by two into pair-based counts
    $ssu_pairs /= 2 if $SEmode == 0;
    $ssu_bad_pairs /= 2 if $SEmode == 0;
    $SSU_total_pairs = $ssu_pairs + $ssu_bad_pairs + $mapped_half;

    # Report fwd and reverse reads mapping
    msg ("Forward read segments mapping: $ssu_f_reads");
    msg ("Reverse read segments mapping: $ssu_r_reads");

    # calculating mapping ratio
    $readnr_pairs = $readnr;
    # If paired reads, divide segment-based counts by two into pair-based counts
    $readnr_pairs /= 2 if $SEmode == 0;

    $SSU_ratio = $SSU_total_pairs / $readnr_pairs;
    $SSU_ratio_pc = sprintf ("%.3f", $SSU_ratio * 100);

    if ($SEmode == 0) {
        # Report summary stats
        msg("Reporting mapping statistics for paired end input");
        msg("Total read pairs with at least one segment mapping: $SSU_total_pairs");
        msg("=> both segments mapping to same reference: $ssu_pairs");
        msg("=> both segments mapping to different references: $ssu_bad_pairs");
        msg("Read segments where next segment unmapped: $mapped_half");
    } elsif ($SEmode == 1) {
        $unmapped = (1-$SSU_ratio)*$readnr;
    }

    # Ratios of mapped vs unmapped to report
    my @mapratio_csv;
    if ($SEmode == 1) { # TO DO: Add numerical values to text labels

        push @mapratio_csv, "Unmapped,".$unmapped;
        push @mapratio_csv, "Mapped,".$SSU_ratio*$readnr;
    } elsif ($SEmode == 0) {
        # For PE reads, do not include Unmapped in piechart because it will be
        # impossible to read the other slices anyway. Instead report mapping
        # ratio in title.

        push @mapratio_csv, "Mapped pair,".$ssu_pairs*2;
        push @mapratio_csv, "Mapped bad pair,".$ssu_bad_pairs*2;
        push @mapratio_csv, "Mapped single,".$mapped_half;
    }

    # Update stats hash
    $statshref->{'total_reads'} = $readnr_pairs;
    $statshref->{'total_read_segments'} = $readnr;
    $statshref->{'ssu_fwd_seg_map'} = $ssu_f_reads;
    $statshref->{'ssu_rev_seg_map'} = $ssu_r_reads;
    $statshref->{'ssu_tot_map'} = $ssu_f_reads + $ssu_r_reads;
    $statshref->{'ssu_tot_pair_map'} = $SSU_total_pairs;
    $statshref->{'ssu_tot_pair_ratio'} = $SSU_ratio;
    $statshref->{'ssu_tot_pair_ratio_pc'} = $SSU_ratio_pc;

    # CSV file to draw piechart
    my $fh_csv;
    open_or_die (\$fh_csv, ">", $outfiles{"mapratio_csv"}{"filename"});
    $outfiles{"mapratio_csv"}{"made"}++;
    print $fh_csv join ("\n", @mapratio_csv);
    close ($fh_csv);

    msg("mapping rate: $SSU_ratio_pc%");

    my $skip_assembly_flag;
    if ($SSU_total_pairs * 2 * $readlength < 1800) {
        msg("WARNING: mapping coverage lower than 1x,\nreconstruction with SPADES and Emirge disabled.");
        $skip_assembly_flag = 1;
    }

    my @output_array = ($readnr,$readnr_pairs,$SSU_total_pairs,$SSU_ratio,$SSU_ratio_pc); # TO DO: Use hash ssu_sam_mapstats to manage this instead
    return (\@output_array,$skip_assembly_flag);
}


sub get_total_reads_from_bbmap_log {
    my ($infile) = @_;
    my $fh;
    my $readnr;
    open_or_die(\$fh, "<", $infile);
    while (<$fh>) {
        if (/^Reads\ Used\:\ +\t+([0-9]*).*$/) {
            $readnr = $1;
            msg("Total read segments processed: $readnr");
        }
        if (/^insert\ median:\ +\t+\ +([0-9]*).*$/) {
            $ins_me = $1;
            msg("insert size median: $ins_me");
        }
        if (/^insert\ std\ dev:\ +\t+\ +([0-9]*).*$/) {
            $ins_std = $1;
            msg("insert size std deviation: $ins_std");
        }
    }
    return ($readnr,$ins_me,$ins_std);
}

sub get_total_reads_from_sortmerna_log {
    my ($infile) = @_;
    my $fh;
    my $readnr;
    open_or_die(\$fh, "<", $infile);
    while (<$fh>) {
        if (/^\s+Total reads = (\d+)$/) {
            $readnr = $1;
            msg("Total read segments processed: $readnr");
        }
    }
    return ($readnr);
}

sub fix_hash_bbmap_sam {
    # Read BBMap SAM file to two structures:
    #  - Array arranged by original order of entries in file
    #  - Hash keyed by the QNAME, read segment, and RNAME
    # Each contain refs to hashes of the SAM entries
    # This allows us to keep the same data referenced in different ways,
    # and when we modify the data (e.g. fix bitflags) it is accessible from both
    # e.g. there are some fixes which are dependent on the order in the file
    #
    # Also fix the following known issues in the input SAM file:
    #  - rev read name is not properly represented for secondary alignments
    #  - bitflag not correct for rev read secondary alignments

    my ($infile,        # Raw SAM file from BBmap
        $stats_href,    # Ref to hash storing stats from mapping
        ) = @_;

    # Hash refs to each SAM alignment by QNAME, segment, and RNAME
    my %sam_hash;
    # Store refs to each SAM alignment in an array, in order encountered in SAM file
    my @sam_arr;

    # Open sam file and hash each read entry
    msg ("Reading SAM file $infile into memory");
    my $fh;
    open_or_die(\$fh, "<", $infile);
    my $current_read; # Keep track of current read pair
    my $current_orientation;
    while (my $line = <$fh>) {
        next if $line =~ m/^@/; # Skip headers
        chomp $line;
        my $line_href = split_sam_line_to_hash($line); # Hash fields accordingly
        my @discard;
        ($line_href->{'QNAME'},@discard) = split /\s/, $line_href->{'QNAME'}; # Split on first whitespace
        # Update current read if this is the first read
        if (! defined $current_read ) {
            $current_read = $line_href->{'QNAME'};
        }
        if (! defined $current_orientation) {
            $current_orientation = 'fwd';
        }

        # Check if SE or PE mode, update current read field accordingly
        if ($SEmode == 0) {
            # PE input
            if ($line_href->{'FLAG'} & 0x40) { # First segment in template
                unless ($line_href->{'FLAG'} & 0x100) { # Ignore secondary alignments
                    $current_read = $line_href->{'QNAME'};  # Update name; split on whitespace
                    $current_orientation = 'fwd';           # Update orientation.
                    $stats_href->{'ssu_fwd_map'}++;         # Update stats counter
                }
            } elsif ($line_href->{'FLAG'} & 0x80) { # Last segment in template
                unless ($line_href->{'FLAG'} & 0x100) {
                    $current_orientation = 'rev';           # Update orientation
                    $stats_href->{'ssu_rev_map'}++;         # Update stats counter
                }
            }
        } else {
            # SE input
            unless ($line_href->{'FLAG'} & 0x100) {
                $current_read = $line_href->{'QNAME'};      # Ignore secondary alignments
                $stats_href->{'ssu_fwd_map'}++;             # Update stats counter
            }
        }

        # Fix wrong assignment of secondary alignments of reverse read
        if ($line_href->{'FLAG'} & 0x100) {
            if ($current_orientation eq 'rev' && $line_href->{'FLAG'} & 0x40) {
                # Fix bug in bbmap handling of BItflags (current as of v37.99)
                $line_href->{'QNAME'} = $current_read;
                $line_href->{'FLAG'} -= 64;     # Turn off flag 0x40
                $line_href->{'FLAG'} += 128;    # Turn on flag 0x80
            }
        }

        # Hash refs to each SAM alignment by QNAME, segment, and RNAME
        push @{$sam_hash{$current_read}{$current_orientation}{$line_href->{'RNAME'}}}, $line_href;
        # Store array of SAM alignments in order encountered in file
        push @sam_arr, $line_href;

    }
    close ($fh);

    # Write copy of fixed SAM file
    if (@sam_arr) {
        msg ("Writing fixed SAM file to ".$outfiles{'sam_map'}{'filename'});
        my $sam_fh;
        open_or_die(\$sam_fh,
                    ">",
                    $outfiles{"sam_map"}{"filename"});
        foreach my $href (@sam_arr) {
            my @outarr;
            my @fields = qw(QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL);
            foreach my $field (@fields) {
                push @outarr, $href->{$field};
            }
            push @outarr, $href->{'OPTIONAL'} if defined $href->{'OPTIONAL'};
            print $sam_fh join "\t", @outarr;
            print $sam_fh "\n";
        }
        $outfiles{"sam_map"}{"made"}++;
        close $sam_fh;
        msg ("Done");
    }

    return (\%sam_hash,
            \@sam_arr);
}

sub parse_stats_taxonomy_from_sam_array {
    # From array of fixed SAM records ...
    # Parse out taxonomy summary to the following global vars
    #  $taxa_ambig_byread_hash
    #  $taxa_full_href
    #  @taxa_full
    #  $taxa_summary_href
    #  $chao1
    #  @xtons

    my ($aref) = @_;

    my @taxa_full;

    msg ("Summarizing taxonomy from mapping hits to SILVA database");
    msg ("Using best hit only") if defined $tophit_flag;
    foreach my $href (@$aref) {
        # Hash taxa hits by read
        if ($href->{'RNAME'} =~ m/\w+\.\d+\.\d+\s(.+)/) { # Match to SILVA header: Accession Taxstring
            my $taxonlongstring = $1;
            if (defined $tophit_flag) { # Global var $tophit_flag
                # "old" way - just take the top hit of each read as the taxon
                # Therefore ignore secondary alignments
                # Count full-length taxon for treemap
                $taxa_full_href->{$taxonlongstring}++ unless $href->{'FLAG'} & 0x100;
                # Save full taxonomy string
                push @taxa_full, $taxonlongstring unless $href->{'FLAG'} & 0x100;
            } else {
                # Use consensus of ambiguous mapping hits to get the consensus taxon per read
                my @taxonlongarr = split /;/, $taxonlongstring;
                taxstring2hash(\%{$taxa_ambig_byread_hash{$href->{'QNAME'}}}, \@taxonlongarr);
            }
        } else {
            msg ("Warning: malformed database entry for ".$href->{'RNAME'});
        }
    }

    if (defined $tophit_flag) {
        # Using only top alignment per read to get taxonomy
        $taxa_summary_href = summarize_taxonomy (\@taxa_full,$taxon_report_lvl); #
    } else {
        # Use consensus of all ambiguous alignments per read to get taxonomy

        # List of all the readcounter IDs
        my @allreadcounters = (keys %taxa_ambig_byread_hash);
        ## Count taxonomy in tree from consensus taxstrings, truncated to taxonomic level requested
        $taxa_summary_href = consensus_taxon_counter(\%taxa_ambig_byread_hash, \@allreadcounters, $taxon_report_lvl);
        ## Count taxonomy in tree from consensus taxstrings, to full taxonomic level (7)
        $taxa_full_href= consensus_taxon_counter(\%taxa_ambig_byread_hash, \@allreadcounters, 7);
    }

    # Calculate Chao1 index from taxon counts
    ($chao1, @xtons) = diversity_stats_from_taxonomy($taxa_summary_href);

    msg("done...");
}

sub diversity_stats_from_taxonomy {
    # From hash of counts keyed by taxon (global var $taxa_summary_href)
    # Count 1-tons, 2-tons, and 3+-tons and calculate Chao1 statistic
    my ($href) = @_;
    my ($chao);
    my ($oneton, $twoton, $moreton) = (0,0,0);
    foreach my $taxon (keys %$href) {
        if ($href->{$taxon} == 1) {        # 1-tons
            $oneton++;
        } elsif ($href->{$taxon} == 2) {   # 2-tons
            $twoton++;
        } elsif ($href->{$taxon} >= 3) {   # 3+-tons
            $moreton++;
        }
    }
    if ($twoton > 0) {
        $chao =
          $moreton +                               # Is there an error here? Should be sum of all spp. observed
          ($oneton * $oneton) / 2 / $twoton;
    } else {
        $chao = 'n.d.';
    }
    return ($chao, $oneton, $twoton, $moreton);
}

sub truncate_taxonstring {
    my ($in, $level) = @_;
    my $out;
    my $lvl = $level - 1;
    my @arr = split (";", $in);
    if ( $#arr > $lvl ) {
        $out = join (";", @arr[0..$lvl]);
    } else {
        $out = $in;
    }
    return ($out);
}

sub summarize_taxonomy {
    # Given list of taxon strings and desired taxonomic level:
    # Count number of occurrences of given taxon substring and output counts
    my ($in_aref,       # Ref to input array of taxon strings
        $lvl            # Level to output summary; 1-based
        ) = @_;
    my @input = @$in_aref; # Dereference array input
    my %taxhash;        # Hash to store counts per taxon at each taxonomic level

    foreach my $taxstring (@input) {
        my $taxshort = truncate_taxonstring ($taxstring, $lvl);
        $taxhash{$taxshort}++;
    }

    return (\%taxhash); # Return reference to output array
}

sub spades_run {
    # running SPADES on bbmap output
    msg("creating phylotypes with SPAdes");

    my $kmer;
    if ($readlength >= 134) {
        $kmer = "99,111,127";
    } else {
        my $spades_rl = $readlength - $readlength % 2; # drop to even number
        $kmer = ($spades_rl - 27).",".($spades_rl - 17).",".($spades_rl - 7);
    }
    msg ("kmers for SPAdes are ".$kmer);

    my $args = "";
    if ($sc == 1) {
      $args = '--sc ';
    }
    if ($SEmode == 1) {
        $args = $args."-s ".$outfiles{"reads_mapped_f"}{"filename"};
    } else {
        $args = $args."-1 ".$outfiles{"reads_mapped_f"}{"filename"}
                ." -2 ".$outfiles{"reads_mapped_r"}{"filename"};
    }

    # Limit number of SPAdes processors to 24 - if run on a server with all 64
    # processors SPAdes will exceed max memory and crash
    my $cpus_spades = $cpus > 24 ? 24 : $cpus;

    my $return = run_prog_nodie("spades",
                                "-o $libraryNAME.spades -t $cpus_spades -m 20 -k $kmer "
                                . $args,
                                $outfiles{"spades_log"}{"filename"},"&1"
                                );
    $outfiles{"spades_log"}{"made"}++;

    msg("done...");
    return ($return);
}

sub spades_parse {
    # getting spades output and reformatting it...
    msg("getting SSU phylotypes and their coverages...");

    #spades scaffolds file is empty, settting skip_spades variable to avoid
    #further processing
    if (! -s "$libraryNAME.spades/scaffolds.fasta") {
        msg("no phylotypes assembled with SPAdes");
        msg("possible causes include low or uneven coverage of SSU sequences");
        $skip_spades = 1;
        remove_tree("$libraryNAME.spades");
        return 1;
    }

    else {
        # run barrnap once for each domain
        # if single cell data - accept partial rRNAs down to 0.1
        # and use lower e-value setting
        my $b_args = " --evalue 1e-100 --reject 0.6 " ;
        if ($sc == 1) {
            $b_args = " --evalue 1e-20 --reject 0.1 ";
        }

        foreach ('bac', 'arch', 'euk') {
            run_prog("barrnap",
                     $b_args.
                     "--kingdom $_ --gene ssu --threads $cpus " .
                     "$libraryNAME.spades/scaffolds.fasta",
                     $outfiles{"gff_".$_}{"filename"},
                     $outfiles{"barrnap_log"}{"filename"});
            $outfiles{"gff_".$_}{"made"}++;
            $outfiles{"barrnap_log"}{"made"}++;
        }

        # now merge multi-hits on the same scaffold-and-strand by picking
        # the lowest start and highest stop position.

        my %ssus;
        # pre-filter with grep for speed
        my $fh;
        open_or_die(\$fh, "-|",
                    "grep -hE '16S_rRNA\|18S_rRNA' ".
                    $outfiles{"gff_bac"}{"filename"}." ".
                    $outfiles{"gff_arch"}{"filename"}." ".
                    $outfiles{"gff_euk"}{"filename"});
        while (my $row = <$fh>) {
            my @cols    = split("\t", $row);
            # gff format:
            # 0 seqname, 1 source, 2 feature, 3 start, 4 end,
            # 5 score, 6 strand, 7 frame, 8 attribute

            my $seqname = $cols[0];
            my $start   = $cols[3];
            my $stop    = $cols[4];
            my $strand  = $cols[6];

            # put our desired output fasta name into "feature" col 3
            # the may be "bug" using / mixing bed and gff formats
            # but it saves us messing with the fasta afterwards
            $seqname =~ m/NODE_([0-9]*)_.*cov_([0-9\\.]*)/;
            $cols[2] = "$libraryNAME.PFspades_$1_$2";

            # do the actual merging, left most start and right most stop wins
            if (exists $ssus{$seqname.$strand}) {
                my $old_start = $ssus{$seqname.$strand}[3];
                my $old_stop  = $ssus{$seqname.$strand}[4];
                $cols[3] = ($start < $old_start) ? $start : $old_start;
                $cols[4] = ($stop  > $old_stop)  ? $stop  : $old_stop;
            }
            $ssus{$seqname.$strand} = [@cols];
        }
        close($fh);

        if (scalar keys %ssus == 0) {
            msg("no contig in spades assembly found to contain rRNA");
            return 1;
        } else {
            open_or_die(\$fh, ">", $outfiles{"gff_all"}{"filename"});
            for my $key (sort keys %ssus) {
                print $fh join("\t",@{$ssus{$key}});
            }
            close($fh);
            $outfiles{"gff_all"}{"made"}++;

            # fastaFromBed will build a .fai index from the source .fasta
            # However, it does not notice if the .fasta changed. So we
            # delete the .fai if it is older than the .fasta.
            if ( -e "$libraryNAME.spades/scaffolds.fasta.fai" &&
                 file_is_newer("$libraryNAME.spades/scaffolds.fasta",
                               "$libraryNAME.spades/scaffolds.fasta.fai")) {
                unlink("$libraryNAME.spades/scaffolds.fasta.fai");
            }

            # extract rrna fragments from spades scaffolds accoding to gff
            my @fastaFromBed_args = ("-fi $libraryNAME.spades/scaffolds.fasta",
                                     "-bed",$outfiles{"gff_all"}{"filename"},
                                     "-fo",$outfiles{"spades_fastaFromBed"}{"filename"},
                                     "-s","-name",
                                     );
            run_prog("fastaFromBed",
                     join (" ",@fastaFromBed_args),
                     $outfiles{"fastaFromBed_out"}{"filename"},
                     "&1");
            # Fix bedtools headers from bedtools v2.26 and v2.27+
            fix_bedtools_header ($outfiles{"spades_fastaFromBed"}{"filename"},
                                 $outfiles{"spades_fasta"}{"filename"}
                                 );
            $outfiles{"spades_fastaFromBed"}{"made"}++;
            $outfiles{"spades_fasta"}{"made"}++;
            $outfiles{"fastaFromBed_out"}{"made"}++;
            msg("done...");
            return 0;
        }
    }
}

sub fix_bedtools_header {
    # Bedtools getfasta v2.26+ produces Fasta headers with chromosome, position,
    # and strand when options -s -name are used, in addition to the sequence
    # name field from GFF/BED file. However, this breaks the sequence name
    # used for downstream processing.
    # Bedtools v.2.27.0+ changes this behavior once again, giving chromosome,
    # position and strand with new -name+ option, but using -s -name still
    # gives the strand information.
    # This sub uses a regex to fix the header names
    my ($infile, $outfile) = @_;
    my ($fh_in, $fh_out);
    open_or_die (\$fh_in,
                 '<',
                 $infile);
    open_or_die(\$fh_out,
                '>',
                $outfile);
    while (my $line = <$fh_in>) {
        chomp $line;
        if ($line =~ m/^>/) {
            # If Fasta header line
            if ($line =~ m/^(>.*)::.*/) {
                # bedtools v2.26 behavior
                print $fh_out $1;
                print $fh_out "\n";
            } elsif ($line =~ m/^(>.*)\([+-]\)/) {
                # bedtools v2.27 behavior
                print $fh_out $1;
                print $fh_out "\n";
            } else {
                # bedtools v2.25 and earlier
                print $fh_out $line;
                print $fh_out "\n";
            }
        } else {
            # If sequence line, print unchanged
            print $fh_out $line;
            print $fh_out "\n";
        }
    }
    close($fh_in);
    close($fh_out);
}

sub trusted_contigs_parse {
    # Extract SSU rRNA sequences from "trusted contigs" Fasta file
    msg("Extracting SSU rRNA from trusted contigs $trusted_contigs...");

    # run barrnap once for each domain
    # if single cell data - accept partial rRNAs down to 0.1
    # and use lower e-value setting
    my $b_args = " --evalue 1e-100 --reject 0.6 " ;
    if ($sc == 1) {
        $b_args = " --evalue 1e-20 --reject 0.1 ";
    }

    foreach ('bac', 'arch', 'euk') {
        run_prog("barrnap",
                 $b_args.
                 "--kingdom $_ --gene ssu --threads $cpus " .
                 "$trusted_contigs",
                 $outfiles{"trusted_gff_".$_}{"filename"},
                 $outfiles{"trusted_barrnap_log"}{"filename"});
        $outfiles{"trusted_gff_".$_}{"made"}++;
        $outfiles{"trusted_barrnap_log"}{"made"}++;
    }

    # now merge multi-hits on the same scaffold-and-strand by picking
    # the lowest start and highest stop position.

    my %ssus;
    # pre-filter with grep for speed
    my $fh;
    open_or_die(\$fh, "-|",
                "grep -hE '16S_rRNA\|18S_rRNA' ".
                $outfiles{"trusted_gff_bac"}{"filename"}." ".
                $outfiles{"trusted_gff_arch"}{"filename"}." ".
                $outfiles{"trusted_gff_euk"}{"filename"});
    my $counter = 0;
    while (my $row = <$fh>) {
        $counter++;
        my @cols    = split("\t", $row);
        # gff format:
        # 0 seqname, 1 source, 2 feature, 3 start, 4 end,
        # 5 score, 6 strand, 7 frame, 8 attribute

        my $seqname = $cols[0];
        my $start   = $cols[3];
        my $stop    = $cols[4];
        my $strand  = $cols[6];

        # put our desired output fasta name into "feature" col 3
        # the may be "bug" using / mixing bed and gff formats
        # but it saves us messing with the fasta afterwards
        $cols[2] = "$libraryNAME.PFtrusted_$counter"."_0";

        # do the actual merging, left most start and right most stop wins
        if (exists $ssus{$seqname.$strand}) {
            my $old_start = $ssus{$seqname.$strand}[3];
            my $old_stop  = $ssus{$seqname.$strand}[4];
            $cols[3] = ($start < $old_start) ? $start : $old_start;
            $cols[4] = ($stop  > $old_stop)  ? $stop  : $old_stop;
        }
        $ssus{$seqname.$strand} = [@cols];
    }
    close($fh);

    if (scalar keys %ssus == 0) {
        msg("no SSU rRNA sequences found in trusted contigs by Barrnap");
        return 1;
    } else {
        open_or_die(\$fh, ">", $outfiles{"trusted_gff_all"}{"filename"});
        for my $key (sort keys %ssus) {
            print $fh join("\t",@{$ssus{$key}});
        }
        close($fh);
        $outfiles{"trusted_gff_all"}{"made"}++;

        # fastaFromBed will build a .fai index from the source .fasta
        # However, it does not notice if the .fasta changed. So we
        # delete the .fai if it is older than the .fasta.
        if ( -e "$trusted_contigs.fai" &&
             file_is_newer("$trusted_contigs",
                           "$trusted_contigs.fai")) {
            unlink("$trusted_contigs.fai");
        }

        # extract rrna fragments from trusted contigs accoding to gff
        my @fastaFromBed_args = ("-fi $trusted_contigs",
                                 "-bed",$outfiles{"trusted_gff_all"}{"filename"},
                                 "-fo",$outfiles{"trusted_fastaFromBed"}{"filename"},
                                 "-s","-name",
                                 );
        run_prog("fastaFromBed",
                 join (" ",@fastaFromBed_args),
                 $outfiles{"trusted_fastaFromBed_out"}{"filename"},
                 "&1");
        # Fix fasta headers from bedtools v2.26 and v2.27+
        fix_bedtools_header ($outfiles{"trusted_fastaFromBed"}{"filename"},
                             $outfiles{"trusted_fasta"}{"filename"}
                             );
        $outfiles{"trusted_fastaFromBed"}{"made"}++;
        $outfiles{"trusted_fasta"}{"made"}++;
        $outfiles{"trusted_fastaFromBed_out"}{"made"}++;
        msg("done...");
        return 0;
    }
}

sub bbmap_remap {
    # Map extracted reads back to the SPAdes or EMIRGE output or trusted contigs to see what
    # proportion of reads can be explained by the assembled sequences
    my ($which,     # SPAdes or EMIRGE output, or trusted contigs SSU ?
        ) = @_;
    my ($ref,$outsam,$outlog);
    if ($which eq "SPAdes") {
        $ref = $outfiles{"spades_fasta"}{"filename"};
        $outsam = $outfiles{"sam_remap_spades"}{"filename"};
        $outlog = $outfiles{"bbmap_remap_log_spades"}{"filename"};
        # Record that output files were made
        $outfiles{"sam_remap_spades"}{"made"}++;
        $outfiles{"bbmap_remap_log_spades"}{"made"}++;
    } elsif ($which eq "EMIRGE") {
        $ref = $outfiles{"emirge_fasta"}{"filename"};
        $outsam = $outfiles{"sam_remap_emirge"}{"filename"};
        $outlog = $outfiles{"bbmap_remap_log_emirge"}{"filename"};
        # Record that output files were made
        $outfiles{"sam_remap_emirge"}{"made"}++;
        $outfiles{"bbmap_remap_log_emirge"}{"made"}++;
    } elsif ($which eq 'trusted') {
        $ref = $outfiles{'trusted_fasta'}{'filename'};
        $outsam = $outfiles{'sam_remap_trusted'}{'filename'};
        $outlog = $outfiles{'bbmap_remap_log_trusted'}{'filename'};
        # Record which output files made
        $outfiles{'sam_remap_trusted'}{'made'}++;
        $outfiles{'bbmap_remap_log_trusted'}{'made'}++;
    }
    # Report for log
    msg("mapping extracted SSU reads back on $which SSU sequences");
    # Define input arguments to BBmap
    my @remap_args = ("fast=t",
                      "minidentity=0.98", # Note high identity
                      "-Xmx20g",
                      "threads=$cpus",
                      "po=f",
                      "outputunmapped=t", # This is important
                      "ref=$ref",
                      "nodisk",
                      "in=".$outfiles{"reads_mapped_f"}{"filename"},
                      "out=$outsam",
                      'noheader=t',
                      'overwrite=t', # Overwrite existing files
                      );
    # If running in PE mode, include reverse reads
    if ($SEmode == 0) {
        if ($interleaved == 1) {
            push @remap_args, ("interleaved=t",
                               "pairlen=$maxinsert",
                               );
        } else {
            push @remap_args, ("in2=".$outfiles{"reads_mapped_r"}{"filename"},
                               "pairlen=$maxinsert",
                               );
        }
    }
    # If mapping vs trusted contigs, separately output reads that don't map to trusted SSU
    if ($which eq 'trusted') {
        push @remap_args, "outu=$outfiles{'reads_mapped_notrusted_f'}{'filename'}";
        $outfiles{'reads_mapped_notrusted_f'}{'made'}++; # Record file as made
        if ($SEmode == 0) { # If paired reads, also create reverse
            push @remap_args, "outu2=$outfiles{'reads_mapped_notrusted_r'}{'filename'}";
            $outfiles{'reads_mapped_notrusted_r'}{'made'}++; # Record file as made
        }
    }

    # Run BBmap
    run_prog("bbmap",
             join (" ", @remap_args),
             undef,
             $outlog,
             );

    msg("done...");
}

sub screen_remappings {
    # Read SAM file from re-mapping and identify whether given read has mapped
    # to a full length sequence from either SPAdes or EMIRGE

    msg ("Reading remappings to summarize taxonomy of unmapped reads"); # TK

    # first SAM in %$sam_fixed_href
    my %sam_spades; # too good to be true
    my %sam_emirge;
    my %sam_trusted;
    # If SPAdes results were mapped, read into memory
    if (defined $outfiles{"spades_fasta"}{"made"}) {
        flag_unmapped_sam("SPAdes");
    }
    # If EMIRGE results were mapped, read into memory
    if (defined $outfiles{"emirge_fasta"}{"made"}) {
        flag_unmapped_sam("EMIRGE");
    }
    # If trusted contigs were remapped, read into memory
    if (defined $outfiles{'trusted_fasta'}{'made'}) {
        flag_unmapped_sam('trusted');
    }

    my @unassem_taxa;
    my @unassem_readcounters; # Readcounter IDs for unassembled reads
    my $ssu_fwd_map = 0;
    my $ssu_rev_map = 0;
    my $total_assembled = 0;    # Count total read segments that map to a full-length SSU

    # Go through hash of initial mapping
    # Count total reads, total assembled/reconstructed, and total unassembled/reconstructed
    foreach my $read (keys %$sam_fixed_href) {
        # Iterate through all segments of read pairs
        my @pairs;
        if ($SEmode == 1) {
            push @pairs, 'fwd';
        } else {
            push @pairs, qw(fwd rev);
        }

        foreach my $pair (@pairs) {
            # Check if read segment exists and add to total reads
            if (defined $sam_fixed_href->{$read}{$pair}) {
                # Check whether read segment been flagged as mapping to SPAdes or Emirge or trusted contig
                if (defined $sam_fixed_href->{$read}{$pair}{"mapped2spades"} | defined $sam_fixed_href->{$read}{$pair}{"mapped2emirge"} | defined $sam_fixed_href->{$read}{$pair}{'mapped2trusted'}) {
                    # Add to total of read segments mapping to a full-length seq
                    $total_assembled ++;
                } else {
                    # If not mapping to full-length seq, check if it has mapped to a SILVA ref sequence
                    foreach my $ref (keys %{$sam_fixed_href->{$read}{$pair}}) {
                        foreach my $href (@{$sam_fixed_href->{$read}{$pair}{$ref}}) {
                            if (ref($href) eq 'HASH' && defined $href->{'FLAG'}) { # If this is a hashed SAM record
                                unless ($href->{'FLAG'} & 0x4 || $href->{'FLAG'} & 0x100) { # Unless segment unmapped or is a secondary alignment
                                    push @unassem_readcounters, $read;
                                    # If using top hit
                                    if ($href->{'RNAME'} =~ m/\w+\.\d+\.\d+\s(.+)/) {
                                        my $taxonlongstring = $1;
                                        push @unassem_taxa, $taxonlongstring;
                                    }
                                }
                            }
                        }


                    }
                }
            }
        }
    }

    # Summarize counts per NTU of unassembled reads
    if (defined $tophit_flag) {
        $taxa_unassem_summary_href = summarize_taxonomy(\@unassem_taxa, $taxon_report_lvl);
    } else {
        $taxa_unassem_summary_href = consensus_taxon_counter(\%taxa_ambig_byread_hash, \@unassem_readcounters, $taxon_report_lvl);
    }

    # Calculate ratio of reads segments assembled and store in hash
    my $total_unassembled = $ssu_sam_mapstats{'ssu_tot_map'} - $total_assembled;
    #$ssu_sam_mapstats{"ssu_tot_map"} = $ssu_tot_map;            # Total reads mapped
    $ssu_sam_mapstats{"assem_tot_map"} = $total_assembled;      # Total mapping to full-length
    $ssu_sam_mapstats{"ssu_unassem"} = $total_unassembled;      # Total not mapping to a full-length
    $ssu_sam_mapstats{"assem_ratio"} = $total_assembled/$ssu_sam_mapstats{'ssu_tot_map'};   # Ratio assembled/reconstructed
    $ssu_sam_mapstats{"assem_ratio_pc"} = sprintf ("%.3f", $ssu_sam_mapstats{"assem_ratio"} * 100); # As a percentage to 3 dp

    # Calculate totals for each tool used
    if (defined $ssu_sam_mapstats{"spades_fwd_map"}) {
        my $spades_tot_map = $ssu_sam_mapstats{"spades_fwd_map"};
        if (defined $ssu_sam_mapstats{"spades_rev_map"}) {
            $spades_tot_map += $ssu_sam_mapstats{"spades_rev_map"};
        }
        $ssu_sam_mapstats{"spades_tot_map"} = $spades_tot_map;
    }
    if (defined $ssu_sam_mapstats{"emirge_fwd_map"}) {
        my $emirge_tot_map = $ssu_sam_mapstats{"emirge_fwd_map"};
        if (defined $ssu_sam_mapstats{"emirge_rev_map"}) {
            $emirge_tot_map +=  $ssu_sam_mapstats{"emirge_rev_map"};
        }
        $ssu_sam_mapstats{"emirge_tot_map"} = $emirge_tot_map;
    }
    if (defined $ssu_sam_mapstats{'trusted_fwd_map'}) {
        my $trusted_tot_map = $ssu_sam_mapstats{'trusted_fwd_map'};
        if (defined $ssu_sam_mapstats{'trusted_rev_map'}) {
            $trusted_tot_map += $ssu_sam_mapstats{'trusted_rev_map'};
        }
        $ssu_sam_mapstats{'trusted_tot_map'} = $trusted_tot_map;
    }

    # Write CSV of reads assembled for piechart
    my @assemratio_csv;
    push @assemratio_csv, "Unassembled,".$total_unassembled;
    push @assemratio_csv, "Assembled,".$total_assembled;
    my $fh_csv;
    open_or_die (\$fh_csv, ">", $outfiles{"assemratio_csv"}{"filename"});
    $outfiles{"assemratio_csv"}{"made"}++;
    print $fh_csv join ("\n", @assemratio_csv);
    close($fh_csv);

    msg ("Done");
}

sub flag_unmapped_sam {
    # Given a SAM file from remapping, and hash of original mapping, it will:
    # 1. Flag reads in hash of original mapping to whether mapped to spades or emirge
    # 2. Count number of reads mapping to each reference to %ssufull_hash
    my ($which) = @_;
    my ($sam_in,$mappedname,$mapped_fwd,$mapped_rev);
    if ($which eq "SPAdes") {
        $sam_in = $outfiles{"sam_remap_spades"}{"filename"};
        $mappedname = "mapped2spades";
        $mapped_fwd = "spades_fwd_map";
        $mapped_rev = "spades_rev_map";
    } elsif ($which eq "EMIRGE") {
        $sam_in = $outfiles{"sam_remap_emirge"}{"filename"};
        $mappedname = "mapped2emirge";
        $mapped_fwd = "emirge_fwd_map";
        $mapped_rev = "emirge_rev_map";
    } elsif ($which eq 'trusted') {
        $sam_in = $outfiles{'sam_remap_trusted'}{'filename'};
        $mappedname = 'mapped2trusted';
        $mapped_fwd = 'trusted_fwd_map';
        $mapped_rev = 'trusted_rev_map';
    }
    msg ("Reading results of remapping to $which results");
    my $fh_in;
    open_or_die(\$fh_in, "<", $sam_in);
    while (my $line = <$fh_in>) {
        next if ($line =~ m/^@/); # Skip headers
        my ($readfull, $bitflag, $ref, @discard) = split /\t/, $line;
        my ($read, @discard2) = split /\s/, $readfull; # Split on first whitespace
        # Check whether fwd or rev
        my $pair;
        if ($bitflag & 0x1) { # If PE read
            if ($SEmode != 0) { # Sanity check
                msg ("Error: bitflag in SAM file conflicts with SE mode flag");
            }
            if ($bitflag & 0x40) {
                $pair = "fwd";
            } elsif ($bitflag & 0x80) {
                $pair = "rev";
            }
        } else {
            $pair = "fwd";
        }
        # Check whether read has mapped to reference
        unless ($bitflag & 0x4) { # negate condition because bitflag 0x4 means "segment unmapped"
            if (defined $sam_fixed_href->{$read}{$pair}) {
                # Mark reads that map to SPAdes or EMIRGE or trusted
                $sam_fixed_href->{$read}{$pair}{$mappedname}++;
                # Add to read count for that reference sequence
                my ($refshort) = $ref =~ /($libraryNAME\.PF\w+)_[\d\.]+/;
                $ssufull_hash{$refshort}{"counts"}++;
                # Count how many fwd and rev reads map to SPAdes
                if ($pair eq "fwd") {
                    $ssu_sam_mapstats{$mapped_fwd}++;
                } elsif ($pair eq "rev") {
                    $ssu_sam_mapstats{$mapped_rev}++;
                }
            } else {
                msg ("Warning: The $pair segment of read $read is not in initial mapping result");
            }
        }
    }
    msg ("total fwd reads remapping:".$ssu_sam_mapstats{$mapped_fwd});
    msg ("total rev reads remapping:".$ssu_sam_mapstats{$mapped_rev}) if $SEmode == 0;
    close($fh_in);
}

sub emirge_run {
    # running Emirge on bbmap output
    msg("creating phylotypes with Emirge");

    my $cmd = "emirge";
    my @emirge_args = ("-1",$outfiles{"reads_mapped_f"}{"filename"});

    if ($SEmode == 1) {
        msg("only one read file provided - running in single end mode");
    } elsif ($readlength < 152) {
        #long reads will hardly map with bowtie in paired end mode, using
        #forward+reverse reads combined as single read data when > 151bp
        msg("reads < 152 bp - running in paired end mode");

        # checking for too phyphysmall insert size
        msg("setting library insert size according to insert size distribution");
        # Calculate minimum insert size for emirge.
        my $min_ins = int(2.2 * $readlength + 0.5);
        if ($ins_me < $min_ins) {
            msg("Warning: estimated insert size very small: $ins_me, "
                . "using $min_ins instead");
            $ins_used = $min_ins;
        } else {
            $ins_used = $ins_me;
        }
        if ($ins_std == 0) {
            # insert size std dev has to be nonzero or Emirge will refuse to run
            # SortMeRNA doesn't report insert size, so we make a guess
            $ins_std = $ins_used / 2;
            msg ("Warning: No insert size reported by mapper. Using initial guess $ins_std");
        }

        msg("the insert size used is $ins_used +- $ins_std");
        # FIXME: EMIRGE dies with too many SSU reads, the cutoff needs to be adjusted...
        if ($SSU_total_pairs <= $amplimit) {
            msg("Less than $amplimit SSU read pairs - using Emirge");
        } else {
            $cmd = "emirge_amp";
            msg("Warning: More than $amplimit SSU reads - using Emirge Amplicon");
        }

        @emirge_args = ("-1",$outfiles{"reads_mapped_f"}{"filename"},
                        "-2",$outfiles{"reads_mapped_r"}{"filename"},
                        "-i $ins_used",
                        "-s $ins_std",
                        );
    } else {
        msg("reads > 151 bp - running in single end mode");
        run_prog("cat",
                 $outfiles{"reads_mapped_f"}{"filename"}." ".
                 $outfiles{"reads_mapped_r"}{"filename"},
                 $outfiles{"reads_mapped_cat"}{"filename"});
        $outfiles{"reads_mapped_cat"}{"made"}++;
        # rename the reads with a running number to make emirge happy
        # using awk for speed, these files can be huge

        run_prog("awk",
                 "\'{print (NR%4 == 1) ? \"\@\" ++i  : (NR%4 == 3) ? \"+\" :\$0}\' "
                 .$outfiles{"reads_mapped_cat"}{"filename"},
                 $outfiles{"reads_mapped_cat_rename"}{"filename"});
        $outfiles{"reads_mapped_cat_rename"}{"made"}++;

        @emirge_args = ("-1",$outfiles{"reads_mapped_cat_rename"}{"filename"});
    }

    unshift @emirge_args, "$libraryNAME.emirge"; # Add library name to arguments
    push @emirge_args, ("-f ${DBHOME}/${emirge_db}.fasta",
                        "-b ${DBHOME}/${emirge_db}.bt",
                        "-l $readlength",
                        "-a $cpus",
                        "--phred33",
                        );
    my $return = run_prog_nodie($cmd,
                                join (" ", @emirge_args),
                                $outfiles{"emirge_log"}{"filename"},
                                "&1",
                                );
    $outfiles{"emirge_log"}{"made"}++;
    msg("done...");
    return ($return);
}

sub emirge_parse {
    # getting emirge output and reformatting it...
    msg("getting Emirge phylotypes and their abundances...");
    run_prog("emirge_rename_fasta",
             "./$libraryNAME.emirge/iter.40/",
             $outfiles{"emirge_raw_fasta"}{"filename"});
    $outfiles{"emirge_raw_fasta"}{"made"}++;

    my $fh_in;
    my $fh_out;
    open_or_die(\$fh_in, "<", $outfiles{"emirge_raw_fasta"}{"filename"});
    open_or_die(\$fh_out, ">", $outfiles{"emirge_fasta"}{"filename"});
    $outfiles{"emirge_fasta"}{"made"}++;
    while (<$fh_in>) {
        chomp;
        if ($_ =~ />(\d*)\|.*NormPrior=([\d.]*)/ ) {
            print {$fh_out} ">$libraryNAME.PFemirge_".$1."_".$2."\n";
        } else {
            print {$fh_out} $_."\n";
        }
    }
    close($fh_in);
    close($fh_out);

    msg("done...");
}

sub vsearch_best_match {
    # running vsearch to map the SPADES output to the taxonomy
    msg("searching for DB matches with vsearch");

    # join emirge and spades hits and trusted contig SSU into one file
    # (vsearch takes a while to load, one run saves us time)
    my @input_fasta;
    foreach my $ff ('emirge_fasta','spades_fasta','trusted_fasta') {
        if (defined $outfiles{$ff}{"made"}) {
           push @input_fasta, $outfiles{$ff}{"filename"};
        }
    }
    run_prog("cat",
             join(" ", @input_fasta),
             $outfiles{"all_final_fasta"}{"filename"});
    $outfiles{"all_final_fasta"}{"made"}++;

    if (-s $outfiles{"all_final_fasta"}{"filename"}) {
        my @vsearch_args = ("-usearch_global", $outfiles{"all_final_fasta"}{"filename"},
                            "-id 0.7",
                            "-userout", $outfiles{"vsearch_csv"}{"filename"},
                            "-userfields query+target+id+alnlen+evalue+id3+qs+pairs+gaps+mism+ids",
                            "-threads $cpus",
                            "--strand plus --notrunclabels",
                            "-notmatched", $outfiles{"notmatched_fasta"}{"filename"},
                            "-dbmatched", $outfiles{"dbhits_all_fasta"}{"filename"},
                            "-db ${DBHOME}/${vsearch_db}",
                            );
        # Run Vsearch
        run_prog("vsearch",
                 join(" ", @vsearch_args),
                 $outfiles{"vsearch_out"}{"filename"},
                 "&1");
        $outfiles{"vsearch_csv"}{"made"}++;
        $outfiles{"dbhits_all_fasta"}{"made"}++;
        $outfiles{"notmatched_fasta"}{"made"}++;
        $outfiles{"vsearch_out"}{"made"}++;

       # query, target: labels
       # id: "100* matching colums / (alignment length - terminal gaps)"
       # alnlen: "number of alignment columns"
       # id3: "MBL definition of %id - counting each extended gap as single difference"
       # qs: query segment length
       # pairs: # letter pair cols
       # gaps: # gap cols
       # mism: # mismatches
       # ids: # matches
    }
}

sub vsearch_parse {
    my @vsearch_matches;            # list of vsearch best-match results

    # Parse the output from Vsearch and store in the hash %SSU_assembly
    my $fh;
    open_or_die(\$fh, "<", $outfiles{"vsearch_csv"}{"filename"});
    while (<$fh>) {
        chomp;
        # lib.PFspades_1_1.23332\tAH12345.1.1490 Bacteria;...\t...
        s/PF(\w+)_([^_]+)_([^\t]+)\t(\w+\.\d+\.\d+)\s/PF$1_$2\t$3\t$4\t/;
        push @vsearch_matches, [split("\t",$_)];                                # To be replaced by ssufull_hash

        # Split into individual fields
        my @line = split /\t/, $_;
        my $seqid = $line[0];
        # Check whether hit is from SPAdes, EMIRGE, or trusted contig
        my $source;
        if ($line[0] =~ m/PFspades/) {
            $source = "SPAdes";
        } elsif ($line[0] =~ m/PFemirge/) {
            $source = "EMIRGE";
        } elsif ($line[0] =~ m/PFtrusted/) {
            $source = 'trusted';
        }

        #fields: qw(source cov dbHit taxon pcid alnlen evalue); counts added later
        $ssufull_hash{$seqid}{"source"} = $source;
        $ssufull_hash{$seqid}{"cov"} = $line[1];
        $ssufull_hash{$seqid}{"dbHit"} = $line[2];
        $ssufull_hash{$seqid}{"taxon"} = $line[3];
        $ssufull_hash{$seqid}{"pcid"} = $line[4];
        $ssufull_hash{$seqid}{"alnlen"} = $line[5];
        $ssufull_hash{$seqid}{"evalue"} = $line[6];

    }
    close($fh);

    msg("done...");
}

sub vsearch_cluster {
    my ($clusterid) = @_;
    msg("clustering DB matches at $clusterid%");
    run_prog("vsearch",
             "  --cluster_fast ".$outfiles{"dbhits_all_fasta"}{"filename"}
             . " -id ".($clusterid/100)
             . " -centroids ".$outfiles{"dhbits_nr97_fasta"}{"filename"}
             . " -notrunclabels"
             . " --threads $cpus",
             $outfiles{"vsearch_clusterdb_out"}{"filename"},
             "&1");
    $outfiles{"dbhits_all_fasta"}{"made"}++;
    $outfiles{"dhbits_nr97_fasta"}{"made"}++;
    $outfiles{"vsearch_clusterdb_out"}{"made"}++;
    msg("done...");
}

sub mafft_run {
    msg("creating alignment and tree...");

    my @input_fasta;
    # Add closest DB hits to assem/recon sequences
    push @input_fasta, $outfiles{"dhbits_nr97_fasta"}{"filename"};
    # Add assem/recon sequences if they were produced
    foreach my $ff ('emirge_fasta','spades_fasta','trusted_fasta') {
        if (defined $outfiles{$ff}{"made"}) {
           push @input_fasta, $outfiles{$ff}{"filename"};
        }
    }
    # Cat all to one file for alignment
    run_prog("cat",
             join(" ", @input_fasta),
             $outfiles{"ssu_coll_fasta"}{"filename"});
    $outfiles{"ssu_coll_fasta"}{"made"}++;

    run_prog("mafft",
             "--treeout ".$outfiles{"ssu_coll_fasta"}{"filename"},
             $outfiles{"ssu_coll_aln_fasta"}{"filename"},
             $outfiles{"ssu_coll_aln_mafftout"}{"filename"});
    $outfiles{"ssu_coll_aln_fasta"}{"made"}++;
    $outfiles{"ssu_coll_tree"}{"made"}++; # Tree filename is input filename with suffix
    $outfiles{"ssu_coll_aln_mafftout"}{"made"}++;

    # fix missing ; at and of MaFFT newick tree
    my $fh;
    open_or_die(\$fh, ">>", $outfiles{"ssu_coll_tree"}{"filename"});
    print {$fh} ";";
    close($fh);

    msg("done...");
}

sub nhmmer_model_pos {
    # Use nhmmer to find how evenly the SSU reads are distributed across the
    # length of the gene. Cannot use mapping position from SAM file because
    # SSU reads in DB are of different lengths and some are fragments.
    # Use the HMM models and nhmmer binary packaged (conveniently) with barrnap-HGV
    msg ("subsampling SSU reads and running nhmmer to check coverage evenness across gene");

    # Inputs
    # File containing HMM models for archaea, bacteria, and eukaryotes SSU
    my $hmm = "$FindBin::RealBin/barrnap-HGV/db/ssu/ssu_ABE.hmm";
    my $sam = $outfiles{"sam_map"}{"filename"}; # SAM file of mapping
    my $readsf = $outfiles{'reads_mapped_f'}{'filename'};
    my $readsr = $outfiles{'reads_mapped_r'}{'filename'};
    my $subsample = $outfiles{"readsf_subsample"}{"filename"};
    my $samplelimit = 10000; # Take sample of max this number of reads

    # Subsample reads with reformat.sh
    my @reformat_args = ("in=$readsf",
                         # "in2=$readsr",
                         #"path=$DBHOME",
                         "out=$subsample",
                         "srt=$samplelimit",
                         "ow=t", # Overwrite existing files
                         );

    # Check if SE reads or interleaved reads
    if ($SEmode == 0) {
        if ($interleaved == 1) {
            push @reformat_args, ("interleaved=t");
        } else {
            push @reformat_args, ("in2=$readsr");
        }
    }

    run_prog ("reformat",
              join (" ", @reformat_args),
              undef,
              $outfiles{"reformat_out"}{"filename"}
              );
    $outfiles{"readsf_subsample"}{"made"}++;
    $outfiles{"reformat_out"}{"made"}++;

    # Run nhmmer
    my @nhmmer_args = ("--cpu $cpus",
                       "--tblout", $outfiles{"nhmmer_tblout"}{"filename"},
                       $hmm, $subsample);
    run_prog ("nhmmer",
              join (" ", @nhmmer_args),
              "/dev/null", # Discard human-readable output
              undef
              );
    $outfiles{"nhmmer_tblout"}{"made"}++;

    # Parse nhmmer output
    my %hash;
    my $fh_tbl;
    open_or_die(\$fh_tbl,"<",$outfiles{"nhmmer_tblout"}{"filename"});
    while (<$fh_tbl>) {
        next if m/^#/;                  # Ignore comment lines
        my @spl = split /\s+/, $_;      # Split on whitespace
        my ($readname,                  # Name of read from Fasta file
            $pos,                       # Position of alignment on model
            $score                      # HMMer score
            ) = ($spl[0], $spl[4], $spl[13]);
        # Figure out model type
        my $type;
        if ($spl[2] eq "18S_rRNA") {
            $type = "euk";
        } elsif ($spl[2] eq "16S rRNA") {
            if ($spl[3] eq "RF01959") {
                $type = "arc";
            } elsif ($spl[3] eq "RF00177") {
                $type = "bac";
            }
        }
        if (defined $hash{$readname}) {
            if ($hash{$readname}{"score"} < $score) {
                # Update recorded values if new hit for same read has better score
                $hash{$readname}{"score"} = $score;
                $hash{$readname}{"pos"} = $pos;
                $hash{$readname}{"type"} = $type;
            }
        } else {
            # Record values if read not yet hashed
            $hash{$readname}{"score"} = $score;
            $hash{$readname}{"pos"} = $pos;
            $hash{$readname}{"type"} = $type;
        }
    }
    close($fh_tbl);

    # Count number of hits to make histogram
    my %prok_pos;
    my %euk_pos;
    foreach my $readname (keys %hash) {
        if (defined $hash{$readname}{"type"} && $hash{$readname}{"type"} eq "euk") {
            $euk_pos{$hash{$readname}{"pos"}}++;
        } else {
            $prok_pos{$hash{$readname}{"pos"}}++;
        }
    }

    # Total up number of hits per model. If fewer than 50 hits, do not draw
    # histogram
    my ($prok_total, $euk_total) = (0,0);
    foreach my $pos (keys %prok_pos) { $prok_total += $prok_pos{$pos}; }
    foreach my $pos (keys %euk_pos) { $euk_total += $euk_pos{$pos}; }

    # Write result to prokaryote histogram
    if ($prok_total >= 50) {
        my $fh_prok;
        open_or_die (\$fh_prok, ">", $outfiles{"nhmmer_prok_histogram"}{"filename"});
        foreach my $key (sort {$a <=> $b} keys %prok_pos) {
                print $fh_prok $key."\t". $prok_pos{$key}."\n";
        }
        close($fh_prok);
        $outfiles{"nhmmer_prok_histogram"}{"made"}++;
    } else {
        msg ("Fewer than 50 maps to prokaryotic SSU HMM model, skip plotting histogram...");
    }

    # Write results for eukaryote histogram
    if ($euk_total >= 50) {
        my $fh_euk;
        open_or_die (\$fh_euk, ">", $outfiles{"nhmmer_euk_histogram"}{"filename"});
        foreach my $key (sort {$a <=> $b} keys %euk_pos) {
            print $fh_euk $key."\t". $euk_pos{$key}."\n";
        }
        close($fh_euk);
        $outfiles{"nhmmer_euk_histogram"}{"made"}++;
    } else {
        msg ("Fewer than 50 maps to eukaryotic SSU HMM model, skip plotting histogram...");
    }

    msg ("done");
}

sub clean_up {
    # cleanup of intermediate files and folders
    if ($keeptmp == 1) {
        msg ("Retaining temp files and folders...");
    } elsif ($keeptmp == 0) {
        msg("Cleaning temp files...");

        # Delete SPAdes and/or EMIRGE output folders
        if ($skip_spades == 0) {
            remove_tree ("$libraryNAME.spades");
        }
        if ($skip_emirge == 0) {
            remove_tree("$libraryNAME.emirge");
        }

        # Remove files that are earmarked for destruction
        # (Change earmarking in PhyloFlash.pm)
        foreach my $key (keys %outfiles) {
            if (defined $outfiles{$key}{"made"} && $outfiles{$key}{"discard"} == 1) {
                unlink $outfiles{$key}{"filename"} or msg ("WARNING: Could not delete file ".$outfiles{$key}{"filename"});
            }
        }
    }
}

sub do_zip {
    # Compress output into tar.gz archive if requested
    my @filelist;
    my @todelete;
    my $tarfile = $outfiles{"phyloFlash_archive"}{"filename"};
    msg ("Compressing results into tar.gz archive $tarfile");
    foreach my $key (keys %outfiles) {
        if (defined $outfiles{$key}{"made"} && $outfiles{$key}{"discard"} ==0) {
            push @filelist, $outfiles{$key}{"filename"};
            push @todelete, $outfiles{$key}{'filename'} unless ($key eq 'report_html' || $key eq 'phyloFlash_log');
            # Retain report and log file after zip
        }
    }
    my $to_tar = join " ", @filelist;
    system ("tar -czf $tarfile $to_tar");
    unlink @todelete or msg ("WARNING: Could not delete files ".join(" ", @todelete));
    msg("Done...");
}

sub run_plotscript_SVG {
    msg ("generating graphics for report in SVG format");
    # Plot mapping ID histogram
    if (defined $outfiles{"idhistogram"}{"made"}) {
        # If using BBmap, this file is generated
        msg ("Plotting mapping ID histogram");
        my @idhist_args = ( "--hist",
                            $outfiles{"idhistogram"}{"filename"},
                            " --title=\"Mapping identity (%)\" ",
                          );
        if ($decimalcomma == 1) {
            push @idhist_args, "--decimalcomma";
        }
        run_prog("plotscript_SVG",
             join (" ", @idhist_args),
             $outfiles{"plotscript_out"}{"filename"},
             "&1");
        $outfiles{"idhistogram_svg"}{"made"}++;
    }

    # Piechart of proportion mapped
    msg("Plotting piechart of mapping ratios");
    my @map_args = ("-pie",
                    $outfiles{"mapratio_csv"}{"filename"},
                    #"-title=\"$SSU_ratio_pc % reads mapped\""
                    );
    if ($SEmode == 0) {
        push @map_args, "-title=\"$SSU_ratio_pc % pairs mapped\"";
    } else {
        push @map_args, "-title=\"$SSU_ratio_pc % reads mapped\"";
    }

    run_prog("plotscript_SVG",
             join(" ", @map_args),
             $outfiles{"plotscript_out"}{"filename"},
             "&1");
    $outfiles{"mapratio_svg"}{"made"}++;

    # Piechart of proportion assembled, if SPAdes or EMIRGE or trusted contigs output
    if (defined $outfiles{"assemratio_csv"}{"made"}) {
        msg ("Plotting piechart of assembled reads ratio");
        my $assem_ratio_pc = $ssu_sam_mapstats{"assem_ratio_pc"};
        my @pie_args = ("--pie",
                        $outfiles{"assemratio_csv"}{"filename"},
                        " --title=\"Reads assembled\"",
                        );
        run_prog("plotscript_SVG",
                 join (" ", @pie_args),
                 $outfiles{"plotscript_out"}{"filename"},
                 "&1");
        $outfiles{"assemratio_svg"}{"made"}++;
    }

    # Plot insert size histogram unless running in SE mode or using sortmerna
    if ($SEmode != 1 && defined $outfiles{"inserthistogram"}{"made"}) { # If not running in SE mode ...
        msg ("Plotting histogram of insert sizes");
        my @inshist_args = ("--hist ",
                            $outfiles{"inserthistogram"}{"filename"},
                            " --title=\"Insert size (bp)\" ",
                            );
        if ($decimalcomma == 1) {
            push @inshist_args, "--decimalcomma";
        }

        run_prog("plotscript_SVG",
                 join (" ", @inshist_args),
                 $outfiles{"plotscript_out"}{"filename"},
                 "&1");
        $outfiles{"inserthistogram_svg"}{"made"}++;
    }

    # Plot positional coverage plots for 18S and 16S rRNA genes unless skipped
    if ($poscov_flag == 1) {
        msg ("Plotting positional coverage across SSU rRNA hmm models");
        # Eukaryotic gene model
        my @args_euk = ("-hist",
                        $outfiles{"nhmmer_euk_histogram"}{"filename"},
                        "-height 150",
                        "-width 480",
                        "-title=\"Coverage on 18S model\"",
                        "-color=\"blue\"",
                        );
        if (-s $outfiles{"nhmmer_euk_histogram"}{"filename"}) {
            # Check that Nhmmer has output a nonzero-sized file (if no reads
            # found, it will still produce output)
            run_prog("plotscript_SVG",
                     join (" ", @args_euk),
                     $outfiles{"plotscript_out"}{"filename"},
                     "&1");
            $outfiles{"nhmmer_euk_histogram_svg"}{"made"}++;
        }

        # Prokaryotic gene model
        my @args_prok = ("-hist",
                         $outfiles{"nhmmer_prok_histogram"}{"filename"},
                         "-height 150",
                         "-width 480",
                         "-title=\"Coverage on 16S model\"",
                         "-color=\"blue\"",
                         );
        if (-s $outfiles{"nhmmer_prok_histogram"}{"filename"}) {
            run_prog("plotscript_SVG",
                     join (" ", @args_prok),
                     $outfiles{"plotscript_out"}{"filename"},
                     "&1");
            $outfiles{"nhmmer_prok_histogram_svg"}{"made"}++;
        }
    }

    # Plot tree if spades/emirge unless no alignment was output
    if (defined $outfiles{"ssu_coll_tree"}{"made"}) {
        msg ("Plotting tree of assembled + reference SSU sequences");
        my @treeplot_args = ("-tree",
                             $outfiles{"ssu_coll_tree"}{"filename"},
                             );
         # Add "bubbles" corresponding to read coverage
        if (defined $outfiles{"full_len_class"}{"made"}) {
            push @treeplot_args, ("-assemcov",
                                  $outfiles{"full_len_class"}{"filename"},
                                  "-treefasta",
                                  $outfiles{"ssu_coll_fasta"}{"filename"}
                                  );
            # Check that unassembled count is more than zero. Can be less than
            # zero in certain cases
            if ($ssu_sam_mapstats{"ssu_unassem"} >= 0) {
                push @treeplot_args, ("-unassemcount",
                                      $ssu_sam_mapstats{"ssu_unassem"});
            } else {
                push @treeplot_args, ("-unassemcount",
                                      '0');
            }

        }
        run_prog("plotscript_SVG",
                 join(" ",@treeplot_args),
                 $outfiles{"plotscript_out"}{"filename"},
                 "&1");
        $outfiles{"ssu_coll_tree_svg"}{"made"}++;
    }

    # Generate barplot of taxonomy at level XX
    msg "Plotting barplot of taxonomic summary";
    my @barplot_args = ("-bar",
                        $outfiles{"ntu_csv"}{"filename"},
                        "-title=\"Taxonomic summary from reads mapped\"",
                        );
    run_prog("plotscript_SVG",
             join(" ", @barplot_args),
             $outfiles{"plotscript_out"}{"filename"},
             "&1");
    $outfiles{"ntu_csv_svg"}{"made"}++;

    # Mark plotscript out tmp file as made
    $outfiles{"plotscript_out"}{"made"}++,

    msg("done");
}

sub generate_treemap_data_rows {
    # Generate data rows for drawing treemap chart
    my %parents_HoH; # Hash of count vals by parent and child taxa
    my @output; # Array to store output
    # Parse taxstrings into hash of child-parent relationships

    foreach my $taxstring (keys %$taxa_full_href) {
        my @taxsplit = split ";", $taxstring; # Get taxonomy string, split by semicolons
        my $taxlen = scalar @taxsplit; # Get length of tax string
        while (scalar @taxsplit > 1) {
            my $countval = 0; # Initialize count value as dummy "zero" by default
            if (scalar @taxsplit == $taxlen) { # If leaf node, set count value to real value
                $countval = $taxa_full_href->{$taxstring};
            }
            my $child_taxstring = join ";", @taxsplit;
            my $dummy = pop @taxsplit; # Get parent node
            my $parent_taxstring = join ";", @taxsplit;
            # Remove non-word and non-semicolon chars to avoid problems with Javascript
            $child_taxstring =~ s/[^\w;_]/_/g;
            $parent_taxstring =~ s/[^\w;_]/_/g;
            # Update the parent-child hash if this taxon not already represented
            if (!exists $parents_HoH{$parent_taxstring}{$child_taxstring}) {
                $parents_HoH{$parent_taxstring}{$child_taxstring} = $countval;
            }
        }
    }

    # Write output in dataRow format for treemap chart
    # Root and top-level taxa
    push @output, "[\'Cellular organisms\',\t,\t0],\n";
    push @output, "[\'Bacteria\',\t\'Cellular organisms\',\t0],\n";
    push @output, "[\'Archaea\',\t\'Cellular organisms\',\t0],\n";
    push @output, "[\'Eukaryota\',\t\'Cellular organisms\',\t0],\n";
    # Go through sorted hash and write parent-child data rows
    foreach my $parent (sort {$a cmp $b} keys %parents_HoH) {
        foreach my $child (sort {$a cmp $b} keys %{$parents_HoH{$parent}}) {
            # Concatenate output as string
            my $outstring = "[\'".$child."\',\t\'".$parent."\',\t".$parents_HoH{$parent}{$child}."],\n";
            push @output, $outstring;
        }
    }
    # Return array of output lines
    return @output;
}

sub write_report_html {
    # Get input parameters
    my ($version,$progcmd,$cwd,$DBHOME,
        $libraryNAME, $id, $SEmode,
        $readsf_full,$readsr_full,
        $readnr,$readnr_pairs,
        $ins_me,$ins_std,$ins_used,
        $SSU_ratio, $SSU_ratio_pc, $SSU_total_pairs,
        $skip_spades, $skip_emirge, $treemap_flag,
        $taxon_report_lvl,
        $mapstats_href,
        $outfiles_href,
        $xtons_aref,$chao1,
        $taxa_summary_href,
        $taxa_unassem_summary_href,
        ) = @_;
    my %outfiles = %$outfiles_href;
    my @xtons = @$xtons_aref;
    # Location of template file
    my $template = "$FindBin::RealBin/phyloFlash_report_template.html",
    msg ("Generating HTML-formatted report and graphics ... ");
    my %flags;  # Hash to store substitution flags in template
    # Hash parameters that are compulsory
    %flags = (
        "VERSION" => $version,
        "LIBNAME" => $libraryNAME,
        "PROGCMD" => $progcmd,
        "CWD" => $cwd,
        "DBHOME" => $DBHOME,
        "ID" => $id,
        "READSF_FULL" => $readsf_full,
        "READNR" => $readnr,
        "SSU_RATIO" => $SSU_ratio, # Not in report?
        "SSU_RATIO_PC" => $SSU_ratio_pc,
        "SSU_TOTAL_PAIRS" => $SSU_total_pairs,
        "TAXON_REPORT_LVL" => $taxon_report_lvl,
        "XTONS0" => $xtons[0],
        "XTONS1" => $xtons[1],
        "XTONS2" => $xtons[2],
        "CHAO1" => $chao1 eq "n.d." ? $chao1 : sprintf ("%.3f", $chao1), # Round to 3 decimal places
    );

    # Define suppress flags (which turn off writing of report)
    my %suppress_flags;
    my %suppress_end_flags;

    if ($skip_spades == 1) {
        $suppress_flags{"SUPPRESS_IF_SKIP_SPADES"} = 1;
        $suppress_end_flags{"SUPPRESS_IF_SKIP_SPADES_END"} = 1;
    }
    if ($skip_emirge == 1) {
        $suppress_flags{"SUPPRESS_IF_SKIP_EMIRGE"} = 1;
        $suppress_end_flags{"SUPPRESS_IF_SKIP_EMIRGE_END"} = 1;
    }
    if (!defined $outfiles{'trusted_fasta'}{'made'}) {
        $suppress_flags{"SUPPRESS_IF_SKIP_TRUSTED"} = 1;
        $suppress_end_flags{"SUPPRESS_IF_SKIP_TRUSTED"} = 1;
    }
    if (!defined $outfiles{'ssu_coll_tree_svg'}{'made'})  {
        $suppress_flags{"SUPPRESS_IF_NO_TREE"} = 1;
        $suppress_end_flags{"SUPPRESS_IF_NO_TREE_END"} = 1;
    }
    if ($SEmode == 1) {
        $suppress_flags{"SUPPRESS_IF_SE_READS"} = 1;
        $suppress_end_flags{"SUPPRESS_IF_SE_READS_END"} = 1;
    } elsif ($SEmode == 0) {
        $suppress_flags{"SUPPRESS_IF_PE_READS"} = 1;
        $suppress_end_flags{"SUPPRESS_IF_PE_READS_END"} = 1;
    }
    unless ($treemap_flag == 1) {
        $suppress_flags{"SUPPRESS_IF_NO_TREEMAP"} = 1 ;
        $suppress_end_flags{"SUPPRESS_IF_NO_TREEMAP_END"} = 1 ;
    }
    if ($poscov_flag == 0) {
        $suppress_flags{"SUPPRESS_IF_NO_POSCOV"} = 1;
        $suppress_end_flags{"SUPPRESS_IF_NO_POSCOV_END"} = 1;
    }

    # Slurp in the SVG plots to embed in HTML file
    $flags{"IDHISTOGRAM"} = slurpfile($outfiles{"idhistogram_svg"}{"filename"}) if $outfiles{"idhistogram_svg"}{"made"};
    $flags{"MAPRATIOPIE"} = slurpfile($outfiles{"mapratio_svg"}{"filename"}) if $outfiles{"mapratio_svg"}{"made"};
    $flags{"TAXONSUMMARYBAR"} = slurpfile($outfiles{"ntu_csv_svg"}{"filename"}) if $outfiles{"ntu_csv_svg"}{"made"};

    # Table of output files produced
    my @table_outfiles;
    my @fileskeys = sort {$outfiles{$a}{"filename"} cmp $outfiles{$b}{"filename"}} keys %outfiles;
    foreach my $key (@fileskeys) {
        if ($outfiles{$key}{"made"} && $outfiles{$key}{"intable"} == 1) {
            push @table_outfiles, "  <tr>\n";
            push @table_outfiles, "    <th>".$outfiles{$key}{"description"}."</th>\n";
            push @table_outfiles, "    <td>".$outfiles{$key}{"filename"}."</td>\n";
            push @table_outfiles, "  </tr>\n";
        }
    }
    $flags{"OUTPUT_FILES_TABLE"} = join "", @table_outfiles;

    # Table of taxonomic affiliations
    my @table_db_map;
    my @taxkeys = sort {${$taxa_summary_href}{$b} <=> ${$taxa_summary_href}{$a}} keys %$taxa_summary_href;
    foreach my $key (@taxkeys) {
        if ($taxa_summary_href->{$key} > 3) { # Only display those with at least three reads mapping
            push @table_db_map, "  <tr>\n";   # Start HTML table row
            push @table_db_map, "    <td>".$key."</td>\n"; # Taxon name
            push @table_db_map, "    <td>".$taxa_summary_href->{$key}."</td>\n"; # Count
            push @table_db_map, "  </tr>\n";
        }
    }
    $flags{"READ_MAPPING_NTU_TABLE"} = join "", @table_db_map;

    # Data rows for Treemap
    if ($treemap_flag == 1) {
        my @treemap_output_array = generate_treemap_data_rows(); # TO DO - scope this fn
        $flags{"TREEMAPDATAROWS"} = join "", @treemap_output_array;
    }

    # Params defined only for PE reads
    if ($SEmode == 0) {
        $flags{"INS_ME"} = $ins_me;
        $flags{"INS_STD"} = $ins_std;
        $flags{"READSR_FULL"} = $readsr_full;
        $flags{"READNR_PAIRS"} = $readnr_pairs;
        $flags{"INSERTHISTOGRAM"} = slurpfile($outfiles{"inserthistogram_svg"}{"filename"}) if $outfiles{"inserthistogram_svg"}{"made"};
    }

    # Slurp in positional coverage histograms, if defined
    if ($poscov_flag == 1) {
        $flags{"POSCOVHIST_PROK"} = slurpfile($outfiles{"nhmmer_prok_histogram_svg"}{"filename"}) if $outfiles{"nhmmer_prok_histogram_svg"}{"made"};
        $flags{"POSCOVHIST_EUK"} = slurpfile($outfiles{"nhmmer_euk_histogram_svg"}{"filename"}) if $outfiles{"nhmmer_euk_histogram_svg"}{"made"};
    }

    # Params defined only for assembled (SPAdes) reads
    if ($skip_spades == 0) {
        # Table of assembled SSU sequences
        my @table_assem_seq;

        foreach my $seqid (sort { $ssufull_hash{$b}{"counts"} <=> $ssufull_hash{$a}{"counts"} } keys %ssufull_hash) {
            # Check that sequence was assembled by spades
            next unless $ssufull_hash{$seqid}{"source"} eq "SPAdes";
            # Parse database entry number of hit to get Genbank accession
            my ($gbk, @discard) = split /\./, $ssufull_hash{$seqid}{"dbHit"};
            $ssufull_hash{$seqid}{"gbklink"} = "<a href=\"http://www.ncbi.nlm.nih.gov/nuccore/".
                                               $gbk.
                                               "\">".
                                               $ssufull_hash{$seqid}{"dbHit"}.
                                               "</a>";
            # Build an output line
            #fields: qw(source cov dbHit taxon pcid alnlen evalue); counts added later
            my @fields = qw(counts cov gbklink taxon pcid alnlen evalue);
            push @table_assem_seq, "  <tr>\n";                  # Start table row
            push @table_assem_seq, "    <td>".$seqid."</td>\n"; # Push sequence ID to line
            foreach my $field (@fields) {                       # Push individual fields to line
                push @table_assem_seq, "    <td>".$ssufull_hash{$seqid}{$field}."</td>\n";
            }
            push @table_assem_seq, "  </tr>\n";                 # End table row
        }
        $flags{"ASSEMBLED_SSU_TABLE"} = join "", @table_assem_seq;
    }

    # Params defined only for reconstructed (EMIRGE) reads
    if ($skip_emirge == 0) {
        $flags {"INS_USED"} = $ins_used;
        my @table_recon_seq;

        foreach my $seqid (sort {$ssufull_hash{$b}{"counts"} <=> $ssufull_hash{$a}{"counts"} } keys %ssufull_hash) {
            # Check that sequence was assembled by spades
            next unless $ssufull_hash{$seqid}{"source"} eq "EMIRGE";
            # Parse database entry number of hit to get Genbank accession
            my ($gbk, @discard) = split /\./, $ssufull_hash{$seqid}{"dbHit"};
            $ssufull_hash{$seqid}{"gbklink"} = "<a href=\"http://www.ncbi.nlm.nih.gov/nuccore/".
                                               $gbk.
                                               "\">".
                                               $ssufull_hash{$seqid}{"dbHit"}.
                                               "</a>";
            # Build an output line
            #fields: qw(source cov dbHit taxon pcid alnlen evalue); counts added later
            my @fields = qw(counts cov gbklink taxon pcid alnlen evalue);
            push @table_recon_seq, "  <tr>\n";                  # Start table row
            push @table_recon_seq, "    <td>".$seqid."</td>\n"; # Push sequence ID to line
            foreach my $field (@fields) {                       # Push individual fields to line
                push @table_recon_seq, "    <td>".$ssufull_hash{$seqid}{$field}."</td>\n";
            }
            push @table_recon_seq, "  </tr>\n";                 # End table row
        }
        $flags{"EMIRGE_TABLE"} = join "", @table_recon_seq;
    }

    # Params defined only for trusted reads
    if (defined $outfiles{'trusted_fasta'}{'made'}) {
        # Table of SSU sequences extracted from trusted contigs
        my @table_trusted_seq;

        foreach my $seqid (sort { $ssufull_hash{$b}{"counts"} <=> $ssufull_hash{$a}{"counts"} } keys %ssufull_hash) {
            # Check that sequence was assembled by spades
            next unless $ssufull_hash{$seqid}{"source"} eq "trusted";
            # Parse database entry number of hit to get Genbank accession
            my ($gbk, @discard) = split /\./, $ssufull_hash{$seqid}{"dbHit"};
            $ssufull_hash{$seqid}{"gbklink"} = "<a href=\"http://www.ncbi.nlm.nih.gov/nuccore/".
                                               $gbk.
                                               "\">".
                                               $ssufull_hash{$seqid}{"dbHit"}.
                                               "</a>";
            # Build an output line
            #fields: qw(source cov dbHit taxon pcid alnlen evalue); counts added later
            my @fields = qw(counts gbklink taxon pcid alnlen evalue);
            push @table_trusted_seq, "  <tr>\n";                  # Start table row
            push @table_trusted_seq, "    <td>".$seqid."</td>\n"; # Push sequence ID to line
            foreach my $field (@fields) {                       # Push individual fields to line
                push @table_trusted_seq, "    <td>".$ssufull_hash{$seqid}{$field}."</td>\n";
            }
            push @table_trusted_seq, "  </tr>\n";                 # End table row
        }
        $flags{"TRUSTED_SSU_TABLE"} = join "", @table_trusted_seq;
    }

    # Params defined if either SPAdes or EMIRGE are run
    if (defined $outfiles{'spades_fasta'}{'made'} || defined $outfiles{'emirge_fasta'}{'made'} || defined $outfiles{'trusted_fasta'}{'made'}) {
        # Assembly ratios
        $flags{"ASSEM_RATIO"} = $mapstats_href->{"assem_ratio_pc"};
        $flags{"ASSEMBLYRATIOPIE"} = slurpfile($outfiles{"assemratio_svg"}{"filename"}) if (-f $outfiles{"assemratio_svg"}{"filename"});

        # Sequence tree
        $flags{"SEQUENCES_TREE"} = slurpfile($outfiles{"ssu_coll_tree_svg"}{"filename"}) if (-f $outfiles{"ssu_coll_tree_svg"}{"filename"});

        # Table of unassembled SSU reads affiliations
        my @taxsort = sort {${$taxa_unassem_summary_href}{$b} <=> ${$taxa_unassem_summary_href}{$a}} keys %$taxa_unassem_summary_href;
        my @table_unassem_lines;
        foreach my $uatax (@taxsort) {
            if (${$taxa_unassem_summary_href}{$uatax} > 3) { # Only display those with at least three reads mapping
                push @table_unassem_lines, "  <tr>\n";                  # Start HTML table row
                push @table_unassem_lines, "    <td>".$uatax."</td>\n"; # Name of taxon
                push @table_unassem_lines, "    <td>".${$taxa_unassem_summary_href}{$uatax}."</td>\n";
                push @table_unassem_lines, "  </tr>\n";
            }
        } # Push into flags hash
        $flags{"UNASSEMBLED_READS_TABLE"} = join "", @table_unassem_lines;
    }

    # Open template and process output
    my ($fh_in, $fh_out);
    open_or_die(\$fh_in, "<", $template);
    open_or_die(\$fh_out, ">", $outfiles{"report_html"}{"filename"});
    $outfiles{"report_html"}{"made"}++;
    my $write_flag = 1;
    while (my $line = <$fh_in>) {
        chomp $line;
        if ($line =~ m/<!--(.+)-->/) {
            my $flag = $1;
            if (defined $flags{$flag}) {
                my $replace = $flags{$flag};
                $line =~ s/<!--$flag-->/$replace/;
            }
            # Turn off writing if option was not specified
            $write_flag = 0 if (defined $suppress_flags{$flag});
            # Turn writing back on if optional section is finished
            $write_flag = 1 if (defined $suppress_end_flags{$flag});
        }
        print $fh_out $line."\n" if ($write_flag == 1);
    }
    close($fh_out);
    close($fh_in);
}



######################### MAIN ###########################

parse_cmdline();
process_required_tools();
check_environment();

my $timer = new Timer;

# Initial mapping vs. SILVA database with either BBmap or Sortmerna
if ($use_sortmerna == 1 ) {
    # Run Sortmerna if explicitly called for
    ($sam_fixed_href, $sam_fixed_aref) = sortmerna_filter_sam();
    # Get total reads processed from Sortmerna log file
    $readnr = get_total_reads_from_sortmerna_log($outfiles{'sortmerna_log'}{'filename'});
    # Set insert size summary stats to 0 because not reported by sortmerna
    ($ins_me, $ins_std) = (0,0);

} else {
    # Run BBmap against the SILVA database
    bbmap_fast_filter_sam_run();
    # Fix bitflags and read names in SAM file if necessary
    ($sam_fixed_href, $sam_fixed_aref) = fix_hash_bbmap_sam($outfiles{'bbmap_sam'}{'filename'},
                                                            \%ssu_sam_mapstats);
    # Get total reads and insert size stats from BBmap log file
    ($readnr,$ins_me,$ins_std) = get_total_reads_from_bbmap_log($outfiles{'bbmap_log'}{'filename'});
}

if ($readnr == 0) {
    # Exit if no reads mapped
    err("ERROR: No reads were detected! Possible reasons: Coverage in library too low; option -readlimit too low; input is not a (meta)genome/transcriptome dataset.");
}

# Get taxonomic summary from hashed SAM records
parse_stats_taxonomy_from_sam_array($sam_fixed_aref);

# Get mapping statistics from hashed SAM records
my ($mapping_stats_aref, $skipflag) = parse_mapstats_from_sam_arr($sam_fixed_aref,$readnr,\%ssu_sam_mapstats,$SEmode);
# Dereference stats
($readnr,$readnr_pairs,$SSU_total_pairs,$SSU_ratio,$SSU_ratio_pc) = @$mapping_stats_aref;
if (defined $skipflag && $skipflag == 1) {
    # If coverage too low, skip assembly
    $skip_spades = 1;
    $skip_emirge = 1;
}

# Find positional coverage along prok and euk SSU rRNA models using nhmmer
# from a subsample of reads
nhmmer_model_pos() if ($poscov_flag == 1 );

# Extract rRNA from trusted contigs if assembly file is supplied
if (defined $trusted_contigs) {
    trusted_contigs_parse();
    bbmap_remap('trusted');
}

# Run SPAdes if not explicitly skipped
if ($skip_spades == 0) {
    my $spades_return = spades_run();
    if ($spades_return == 0) {
        my $spades_parse_return = spades_parse();
        if ($spades_parse_return == 0) {
            bbmap_remap("SPAdes");
        }
    } else {
        msg("SPAdes exited with error, this may happen if no sequences were assembled Possible causes include coverage per sequence too low or uneven");
    }
}
# Run Emirge if not explicitly skipped
if ($skip_emirge == 0) {
    my $emirge_return = emirge_run();
    # Parse output from emirge if it exits without errors
    if ($emirge_return == 0) {
        emirge_parse();
        bbmap_remap("EMIRGE");
    } else {
        msg("EMIRGE exited with error, this may happen if no sequences were reconstructed Possible causes include coverage per sequence too low or uneven");
    }
}
# If full length sequenced produced, match vs SILVA database with Vsearch and align with best hits
if (defined $outfiles{"spades_fasta"}{"made"} || defined $outfiles{"emirge_fasta"}{"made"} || defined $outfiles{'trusted_fasta'}{'made'}) {
    vsearch_best_match();
    vsearch_parse();
    vsearch_cluster($clusterid);
    mafft_run();
    screen_remappings();
}

$runtime = $timer->minutes; # Log run time - This should go behind print_report()

# Capture output parameters for reports
my @report_inputs = (
    $version, $progcmd, $cwd, $DBHOME, $libraryNAME, $id, $SEmode,              # input parameters
    $readsf_full, $readsr_full, $readnr, $readnr_pairs,                         # Read input stats
    $ins_me,$ins_std,$ins_used,                                                 # Insert size stats
    $SSU_ratio, $SSU_ratio_pc, $SSU_total_pairs,                                # Mapping stats
    $skip_spades, $skip_emirge, $treemap_flag,                                  # Option flags used
    $taxon_report_lvl,                                                          # Taxonomic summary
    \%ssu_sam_mapstats,\%outfiles,\@xtons,$chao1,
    $taxa_summary_href, $taxa_unassem_summary_href,
    #\@ssuassem_results_sorted, \@ssurecon_results_sorted,
    );

# Print report file and CSV output
print_report();
write_csv();

# SVG formatted graphics and HTML report
# Must be run after write_csv()
if ($html_flag) {
    run_plotscript_SVG();
    write_report_html(@report_inputs);
}

# Clean up temporary files
clean_up();

msg("Walltime used: $runtime with $cpus CPU cores");
msg("Thank you for using phyloFlash
You can find your results in $libraryNAME.*,
Your main result file is $libraryNAME.phyloFlash");

if ($save_log == 1) {
    write_logfile($outfiles{"phyloFlash_log"}{"filename"});
    $outfiles{'phyloFlash_log'}{'made'} ++;
}

do_zip() if ($zip == 1) ;
