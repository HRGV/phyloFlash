#!/usr/bin/env perl
use strict;
use warnings;
=head1 NAME

phyloFlash - A script to rapidly estimate the phylogenetic composition of
an illumina (meta)genomic dataset and reconstruct SSU rRNA genes.

=head1 SYNOPSIS

B<phyloFlash.pl> [OPTIONS] -lib I<name> -read1 F<<file>> -read2 F<<file>>

B<phyloFlash.pl> -help

B<phyloFlash.pl> -man

B<phyloFlash.pl> -check_env

=head1 DESCRIPTION

This tool rapidly approximates the phylogenetic composition of a
(meta)genomic read set based on SSU mapping and reconstruction.
Right now Illumina paired end or single HiSeq and MiSeq reads
are supported.

=head1 ARGUMENTS

=over 15

=item -lib I<name>

Library I<name> to use for output file. The name must be one word comprising
only letters, numbers and "_" or "-" (no whitespace or other punctuation).

=item -read1 F<file>

File containing forward reads. Both FASTA and FASTQ formats are understood.
The file may be compressed (.gz). If interleaved read data is provided, please
use --interleaved flag for paired-end processing

=item -read2 F<file>

File containing reverse reads. If this option is omitted, B<phyloFlash>
will run in B<experimental> single-end mode.

=item -check_env

Invokes checking of working environment and dependencies without data input.
Use to test setup.

=item -help

Print brief help.

=item -man

Show manual.

=back

=head1 OPTIONS

=over 15


=item -dbhome F<dir>

Directory containing phyloFlash reference databases.
Use F<phyloFlash_makedb.pl> to create an appropriate directory.

=item -interleaved

Interleaved readfile with R1 and R2 in a single file at read1

=item -readlength I<N>

Sets the expected readlength. Always use this if your read length
differs from 100 (the default). Must be within 50..500.

=item -CPUs I<N>

Set the number of threads to use. Defaults to all available CPU cores.

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

Minimum allowed identity for read mapping process in %. Must be within
63..98. Set this to a lower value for very divergent taxa
Default: 70

=item -clusterid I<N>

Identity threshold for clustering with vsearch in %.
Must be within 50..100. Default: 97

=item -maxinsert I<N>

Maximum insert size allowed for paired end read mapping. Must be within
0..1200. Default: 1200

=item -html

Also generate output in HTML format.

=item -treemap

Draw interactive treemap of taxonomic classification in html-formatted
report. This uses Google Visualization API, which requires an internet
connection, requires that you agree to their terms of service (see
https://developers.google.com/chart/terms), and is not open source,
although it is free to use.

=item -crlf

Use CRLF as line terminator in CVS output (to become RFC4180 compliant).

=item -decimalcomma

Use decimal comma instead of decimal point to fix locale problems
(default: off)

=item -skip_emirge

Turn off EMIRGE reconstruction of SSU sequences

=item -skip_spades

Turn off SPAdes assembly of SSU sequences

=item -sc

Turn on single cell MDA data mode for SPAdes assembly of SSU sequences

=back

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2014- by Harald Gruber-Vodicka <hgruber@mpi-bremen.de>
                    and Elmar A. Pruesse <elmar.pruesse@ucdenver.edu>
                    with help from Brandon Seah <kbseah@mpi-bremen.de>

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

use strict;
use warnings;

use FindBin;
use lib $FindBin::RealBin;

use PhyloFlash;

use Pod::Usage;
use Getopt::Long;
use File::Basename;
use IPC::Cmd qw(can_run);
use Cwd;
use FindBin;

# change these to match your installation
my @dbhome_dirs = (".", $ENV{"HOME"}, $FindBin::RealBin);

# constants
my $version       = 'phyloFlash v3.0 beta 1';
my $progname      = $FindBin::Script;
my $cwd = getcwd;

# configuration variables / defaults
my $DBHOME      = undef;
my $readsf_full = undef;        # path to forward read file
my $readsr_full = undef;        # path to reverse read file
my $readsf      = undef;        # forward read filename, stripped of directory names
my $readsr      = undef;        # reverse read filename, stripped of directory names
my $SEmode      = 0;            # single ended mode
my $libraryNAME = undef;        # output basename
my $interleaved = 0;            # Flag - interleaved read data in read1 input (default = 0, no)
my $id          = 70;           # minimum %id for mapping
my $readlength  = 100;          # length of input reads
my $readlimit   = -1;           # max # of reads to use
my $amplimit    = 500000;       # number of SSU pairs at which to switch to emirge_amplicon
my $maxinsert   = 1200;         # max insert size for paired end read mapping
my $cpus        = get_cpus      # num cpus to use
my $clusterid   = 97;           # threshold for vsearch clustering

my $html_flag   = 0;            # generate HTML output? (default = 0, off)
my $treemap_flag = 0;           # generate interactive treemap (default = 0, off)
my $crlf        = 0;            # csv line terminator
my $decimalcomma= 0;            # Decimal separator (default = .)
my $skip_emirge = 0;            # Flag - skip Emirge step? (default = 0, no)
my $skip_spades = 0;            # Flag - skip SPAdes assembly? (default = 0, no)
my $sc          = 0;            # Flag - single cell data? (default = 0, no)
my $check_env   = 0;            # Check environment (runs check_environment subroutine only)
my @tools_list;                 # Array to store list of tools required
                                # (0 will be turned into "\n" in parsecmdline)
# default database names for EMIRGE and Vsearch
my $emirge_db   = "SILVA_SSU.noLSU.masked.trimmed.NR96.fixed";
my $vsearch_db  = "SILVA_SSU.noLSU.masked.trimmed";

my $ins_used = "SE mode!"; # Report insert size used by EMIRGE

# variables for report generation
my %ssu_sam;            # Readin sam file from first mapping vs SSU database
my %ssu_sam_mapstats;   # Statistics of mapping vs SSU database, parsed from from SAM file

my %taxa_from_hitmaps;          # Hash of read counts from first mapping, keyed by taxon string         # To be replaced
my @taxa_from_hitmaps_sorted;   # Sorted list of read counts                                            # To be replaced
my %taxa_from_hitmaps_unassem;          # Hash of read counts for unassembled reads, keyed by taxon string  # To be replaced
my @taxa_from_hitmaps_unassem_sorted;   # Sorted list of read counts for unassembled reads                  # To be replaced

my $taxa_summary_href;          # Summary of read counts at specific taxonomic level (hash ref)
my $taxa_unassem_summary_href;  # Summary of read counts of UNASSEMBLED reads at specific taxonomic level (hash ref)

my $taxon_report_lvl = 4;       # Taxonomic level to report counts

my @ssuassem_results_sorted;    # Sorted list of SSU sequences for reporting
my %ssuassem_cov;               # Coverage of assembled SSU sequences
my @ssurecon_results_sorted;

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
    foreach ('ref/genome/1/summary.txt', $emirge_db.".fasta",
	     $vsearch_db.".fasta") {
	return "${dbhome}/$_"
	    unless -r "${dbhome}/$_"
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
        barrnap => "$FindBin::RealBin/barrnap-HGV/bin/barrnap_HGV",
        vsearch => "vsearch",
        mafft => "mafft",
        fastaFromBed => "fastaFromBed",
        sed => "sed",
        grep => "grep",
        awk => "awk",
        cat => "cat",
        plotscript => "$FindBin::RealBin/phyloFlash_plotscript.R",
        plotscript_SVG => "$FindBin::RealBin/phyloFlash_plotscript_svg.pl"
    );
    if ($skip_spades == 0) {
        require_tools(spades => "spades.py");
    }
    if ($skip_emirge == 0) {
        require_tools(emirge => "emirge.py",
                      emirge_amp => "emirge_amplicon.py",
                      emirge_rename_fasta => "emirge_rename_fasta.py");
    }
}

# parse arguments passed on commandline and do some
# sanity checks
sub parse_cmdline {
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
               'CPUs=i' => \$cpus,
               'html' => \$html_flag,
               'treemap' => \$treemap_flag,
               'crlf' => \$crlf,
               'decimalcomma' => \$decimalcomma,
               'skip_emirge' => \$skip_emirge,
               'skip_spades' => \$skip_spades,
               'sc' => \$sc,
               'check_env' => \$check_env,
               'help' => sub { pod2usage(1) },
               'man' => sub { pod2usage(-exitval=>0, -verbose=>2) },
           )
        or pod2usage(2);


    # verify tools present
    if ($check_env == 1) {
      process_required_tools();
      check_environment(); # will die on failure
    }


    # verify database present
    if (defined($DBHOME)) {
	if (my $file = check_dbhome($DBHOME)) {
	    pod2usage("\nBroken dbhome directory: missing file \"$file\"")
	}
    } else {
	$DBHOME = find_dbhome();
	pod2usage("Failed to find suitable DBHOME. (Searched \""
		  .join("\", \"",@dbhome_dirs)."\".)\nPlease provide a path using -dbhome. "
		  ."You can build a reference database using phyloflash_makedb.pl\n")
	    if ($DBHOME eq "");
    }
    msg("Using dbhome '$DBHOME'");


    # verify valid lib name
    pod2usage("Please specify output file basename with -lib")
        if !defined($libraryNAME);
    pod2usage("\nArgument to -lib may not be empty")
        if length($libraryNAME) == 0;
    pod2usage("\nArgument to -lib may contain only alphanumeric characters,"
              ."'_' and '-'.")
        if not $libraryNAME =~ m/[a-zA-Z0-9_-]/ ;

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
        $readsr = "<NONE>";
        $readsr_full = "<NONE>";
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


    # check lengths
    pod2usage("\nReadlength must be within 50...500")
        if ($readlength < 50 or $readlength > 500);
    pod2usage("\nMaxinsert must be within 0..1200")
        if ($maxinsert < 0 or $maxinsert > 1200);
    pod2usage("\nReadmapping identity (-id) must be within 63..98")
        if ($id < 63 or $id > 98);
    pod2usage("\nClustering identidy (-clusterid) must be within 50..100")
        if ($clusterid < 50 or $clusterid > 100);

    $crlf = $crlf ? "\r\n":"\n";

    # check CPUS
    if ($cpus eq "all" or $cpus < 0) {
        $cpus = get_cpus();
    }

    # check surplus arguments
    err("Command line contains extra words:", @ARGV)
        if ($ARGV[0]);
}

sub print_report {
    ## Write plaintext report file

    msg("writing final files...");
    my $fh;
    open_or_die(\$fh, '>', $libraryNAME.'.phyloFlash');

    print {$fh} qq ~
$version - high throughput phylogenetic screening using SSU rRNA gene(s) abundance(s)
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

    print {$fh} qq~
---
Runtime:\t$runtime
CPUs used:\t$cpus
---
Read mapping based higher taxa (NTUs) detection

NTUs observed once:\t$xtons[0]
NTUs observed twice:\t$xtons[1]
NTUs observed three or more times:\t$xtons[2]
NTU Chao1 richness estimate:\t$chao1

List of NTUs in order of abundance (min. 3 reads mapped):
NTU\treads
~;

    # sort keys numerically descending in hash of
    # mapping-based detected higher taxa
    foreach my $taxonshortstring ( @taxa_from_hitmaps_sorted ) {
        # Print the name of this higher taxon, and the
        # corresponding no. of unambig hits mapped
        print {$fh} join("\t",@$taxonshortstring)."\n";
    }

    if ($skip_spades == 0) {
        ## Print the table of SSU assembly-based taxa to report file
        print {$fh} "---\n";
        print {$fh} "SSU assembly based taxa:\n";
        print {$fh} "OTU\tread_cov\tcoverage\tdbHit\ttaxonomy\t%id\talnlen\tevalue\n";

        foreach my $arr (@ssuassem_results_sorted) {
            my @out = @$arr;
            my $spades_id = $out[0];
            splice @out, 1, 0, $ssuassem_cov{$spades_id};
            print {$fh} join("\t",  @out)."\n";
        }

        ## Print the table of taxonomic affiliations for unassembled SSU reads
        print {$fh} "---\n";
        print {$fh} "Taxonomic affiliation of unassembled reads (min. 3 reads mapped):\n";
        foreach my $taxonshortstring ( @taxa_from_hitmaps_unassem_sorted ) {
            # Print the name of this higher taxon, and the
            # corresponding no. of unambig hits mapped
            print {$fh} join("\t",@$taxonshortstring)."\n";
        }
    }

    if ($skip_emirge == 0) {
        ## Print the table of SSU reconstruction-based taxa to report file
        print {$fh} "---\n";
        print {$fh} "SSU reconstruction based taxa:\n";
        print {$fh} "OTU\tratio\tdbHit\ttaxonomy\t%id\talnlen\tevalue\n";
        foreach my $arr (@ssurecon_results_sorted) {
            print {$fh} join("\t", @$arr)."\n";
        }
    }

    close($fh);
}

sub write_csv {
    msg("exporting results to csv");
    my $fh;

    my @report = csv_escape((
        "version",$version,
        "library name",$libraryNAME,
        "forward read file",$readsf_full,
        "reverse read file",$readsr_full,
        "cwd",$cwd,
        "minimum mapping identity",$id,
        "single ended mode",$SEmode,
        "input reads",$readnr_pairs,
        "mapped SSU reads",$SSU_total_pairs,
        "mapping ratio",$SSU_ratio_pc,
        "detected median insert size",$ins_me,
        "used insert size",$ins_used,
        "insert size stddev",$ins_std,
        "runtime",$runtime,
        "CPUs",$cpus,
        "NTUs observed once",$xtons[0],
        "NTUs observed twice",$xtons[1],
        "NTUs observed three or more times",$xtons[2],
        "NTU Chao1 richness estimate",$chao1
    ));

    open_or_die(\$fh, ">", "$libraryNAME.phyloFlash.report.csv");
    while ($#report > 0) {
        print {$fh} shift(@report).",".shift(@report).$crlf;
    }
    close($fh);

    open_or_die(\$fh, ">", "$libraryNAME.phyloFlash.NTUabundance.csv");
    print $fh "NTU,reads\n";
    foreach my $arr ( @taxa_from_hitmaps_sorted ) {
        print {$fh} join(",",csv_escape(@$arr)).$crlf;
    }
    close($fh);

    open_or_die(\$fh, ">", "$libraryNAME.phyloFlash.taxonsummary.csv");
    # Sort results descending
    my @keys = sort {${$taxa_summary_href}{$b} <=> ${$taxa_summary_href}{$a}} keys %$taxa_summary_href;
    foreach my $key (@keys) {
        my @out = ($key, ${$taxa_summary_href}{$key});
        print {$fh} join(",",csv_escape(@out)).$crlf;
    }
    close($fh);

    if ($skip_spades + $skip_emirge < 2) {  # Check if SPAdes or Emirge were skipped
      open_or_die(\$fh, ">", "$libraryNAME.phyloFlash.extractedSSUclassifications.csv");
      print $fh "OTU,read_cov,coverage,dbHit,taxonomy,%id,alnlen,evalue\n";
      if ($skip_spades == 0) {
        foreach my $arr (@ssuassem_results_sorted) {
            my $otu = ${$arr}[0];
            my @out = @$arr;
            splice @out, 1, 0, $ssuassem_cov{$otu}; # Insert read coverage into output array
            print {$fh} join(",",csv_escape(@out)).$crlf;
        }
        my $fh2;
        open_or_die (\$fh2, ">", "$libraryNAME.phyloFlash.unassembled.NTUabundance.csv");
        print $fh2 "NTU,reads\n";
        foreach my $arr ( @taxa_from_hitmaps_unassem_sorted ) {
            print {$fh2} join (",",csv_escape(@$arr)).$crlf;
        }
      }
      if ($skip_emirge == 0) {
        foreach my $arr (@ssurecon_results_sorted) {
            print {$fh} join(",",csv_escape(@$arr)).$crlf;
        }
      }
      close($fh);
    }
}

sub bbmap_fast_filter_sam_run {
    # input: $readsf, $readsr
    # output: lib.$readsf.SSU.sam
    #         lib.$readsf.SSU.{1,2}.fq
    #         lib.bbmap.out
    #         lib.inserthistogram
    # tmp:    tmp.lib.basecompositionhistogram
    msg("filtering reads with SSU db using minimal identity of $id%");
    if ($readlimit != -1) {
        msg("Only using the first $readlimit reads");
    }

    my $minID= $id / 100;
    if ($minID < 0.63) {
        $minID = 0.63
    }

    my $args = "";
    if ($SEmode == 0) {
	if ($interleaved == 1) {
	    $args =
	    "  outm2=$libraryNAME.$readsf.SSU.2.fq "
	    . "pairlen=$maxinsert interleaved=t";
	} else {
	    $args =
	    "  outm2=$libraryNAME.$readsf.SSU.2.fq "
	    . "pairlen=$maxinsert in2=$readsr_full";
	}
    }
    run_prog("bbmap",
             "  fast=t "
             . "minidentity=$minID "
             . "-Xmx20g reads=$readlimit "
             . "threads=$cpus "
             . "po=f "
             . "outputunmapped=f "
             . "path=$DBHOME "
	     . "out=$libraryNAME.$readsf.SSU.sam " # Also skip SAM header?
             . "outm=$libraryNAME.$readsf.SSU.1.fq "
             . "build=1 "
             . "in=$readsf_full "
             . "bhist=tmp.$libraryNAME.basecompositionhistogram "
             . "ihist=$libraryNAME.inserthistogram "
             . "idhist=$libraryNAME.idhistogram "
             . "scafstats=$libraryNAME.hitstats "
             . $args,
             undef,
             "$libraryNAME.bbmap.out");

    msg("done...");
}

sub bbmap_fast_filter_parse {
    # parsing bbmap.out for used read numbers, mapped reads,
    # insert size median and standard deviation

    # input: lib.bbmap.out
    my ($infile,            # Input file
        $SEmode,            # Flag for single-end mode
        ) = @_;
    # Variables used internally
    my $ssu_pairs     = 0;
    my $ssu_bad_pairs = 0;
    my $ssu_f_reads   = 0;
    my $ssu_r_reads   = 0;
    my $forward_count = 0;
    # Variables for output
    my ($readnr,$readnr_pairs,$SSU_total_pairs,$SSU_ratio,$SSU_ratio_pc);

    my $fh;
    #open_or_die(\$fh, "<", $libraryNAME.'.bbmap.out');
    open_or_die(\$fh, "<", $infile);
    while (<$fh>) {
        if (/^Reads\ Used\:\ +\t+([0-9]*).*$/) {
            $readnr = $1;
            msg("single reads mapped: $readnr");
        }
        if (/^mated pairs:\s+\S+\s+([0-9]*).*$/) {
            $ssu_pairs = $1;
            msg("mapped SSU pairs: $ssu_pairs");
        }
        if (/^bad pairs:\s+\S+\s+([0-9]*).*$/) {
            $ssu_bad_pairs = $1;
            msg("mapped bad SSU pairs: $ssu_bad_pairs");
        }
        if (/^insert\ median:\ +\t+\ +([0-9]*).*$/) {
            $ins_me = $1;
            msg("insert size median: $ins_me");
        }
        if (/^insert\ std\ dev:\ +\t+\ +([0-9]*).*$/) {
            $ins_std = $1;
            msg("insert size std deviation: $ins_std");
        }
        if ($forward_count == 0) {
            if (/^mapped:\s+\S+\s+([0-9]*).*$/) {
                $ssu_f_reads = $1;
                msg("mapped forward SSU reads: $ssu_f_reads");
                $forward_count++;
                next;
            }
        }
        if ($forward_count == 1) {
            if (/^mapped:\s+\S+\s+([0-9]*).*$/) {
                $ssu_r_reads = $1;
                msg("mapped reverse SSU reads: $ssu_r_reads");
                last;
            }
        }
    }
    close($fh);

    # calculating mapping ratio
    $readnr_pairs = $readnr;
    if ($SEmode == 0) {
        $readnr_pairs /= 2;
    }

    # Total read pairs with at least one partner mapping
    $SSU_total_pairs =
        $ssu_f_reads + $ssu_r_reads
        - $ssu_pairs - $ssu_bad_pairs;
    if ($SEmode == 0) {
        msg("mapped pairs output: $SSU_total_pairs")
    };

    $SSU_ratio = $SSU_total_pairs / $readnr_pairs;
    $SSU_ratio_pc = $SSU_ratio * 100;

    msg("mapping rate: $SSU_ratio_pc%");

    my $skip_assembly_flag;
    if ($SSU_total_pairs * 2 * $readlength < 1800) {
      msg("WARNING: mapping coverage lower than 1x,\n
      reconstruction with SPADES and Emirge disabled.");
      $skip_assembly_flag = 1;
    }

    my @output_array = ($readnr,$readnr_pairs,$SSU_total_pairs,$SSU_ratio,$SSU_ratio_pc);
    return (\@output_array,$skip_assembly_flag);
}


sub readsam {
    # Read SAM file into memory

    # Input params
    my $infile = "$libraryNAME.$readsf.SSU.sam";    # Input SAM filename
    my $href = \%ssu_sam;                           # Reference to hash to store SAM data
    my $stats_href = \%ssu_sam_mapstats;            # Reference to hash to store mapping statistics

    # Internal vars
    my @taxa_full;                                  # Arr of taxon names from first mapping

    msg ("reading mapping into memory");
    my $fh;
    open_or_die(\$fh, "<", "$infile");
    while (my $line = <$fh>) {
        next if ($line =~ m/^@/); # Skip header lines
        my ($read, $bitflag, $ref, @discard) = split /\t/, $line;
        # If not mapped, skip entry, do NOT record into hash
        if ($bitflag & 0x4) {
            next;
        } else {
            # Check if reads are PE or SE
            my $pair;
            if ($bitflag & 0x1) { # If PE read
                if ($SEmode != 0) { # Sanity check
                    msg ("ERROR: bitflag in SAM file conflicts with SE mode flag ");
                }
                if ($bitflag & 0x40) {
                    $pair="F";
                    ${$stats_href}{"ssu_fwd_map"}++;
                } elsif ($bitflag & 0x80) {
                    $pair="R";
                    ${$stats_href}{"ssu_rev_map"}++;
                }
            } else {
                $pair = "U";
                ${$stats_href}{"ssu_fwd_map"}++;
            }
            # Record into hash
            ${$href}{$read}{$pair}{"ref"} = $ref;
            ${$href}{$read}{$pair}{"bitflag"} = $bitflag;

            # Shorten taxonomy string and save into NTU table
            if ($ref =~ m/\w+\.\d+\.\d+\s(.+)/) {
                my $taxonlongstring = $1;
                # Truncate to 6 levels
                my $taxonshortstring = truncate_taxonstring($taxonlongstring, 6); # To be superseded
                # Increment count for taxon
                $taxa_from_hitmaps{$taxonshortstring}++;                          # To be superseded
                # Save full taxonomy string
                push @taxa_full, $taxonlongstring;
            } else {
                msg ("Warning: malformed database entry $ref");
            }
        }
    }
    close($fh);

    # Sort counts of taxa from hits, minimum of three hits                      # To be superseded
    @taxa_from_hitmaps_sorted =
        sort { @$b[1] <=> @$a[1] }
        grep { if (@$_[1] < 3) { @xtons[@$_[1]-1]++;} @$_[1] > 2 }
        map  { [$_, $taxa_from_hitmaps{$_}] }
        keys %taxa_from_hitmaps;

    # Calculate Chao1 statistic
    $xtons[2] = $#taxa_from_hitmaps_sorted +1;
    if ($xtons[1] > 0) {
        $chao1 =
          $xtons[2] +
          ($xtons[0] * $xtons[0]) / 2 / $xtons[1];
    } else {
        $chao1 = 'n.d.';
    }

    # Summarize taxonomy
    $taxa_summary_href = summarize_taxonomy(\@taxa_full, $taxon_report_lvl); # Summarize

    msg("done...");
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
    my %output;         # Hash to store output

    foreach my $taxstring (@input) {
        my $taxshort = truncate_taxonstring ($taxstring, $lvl);
        $taxhash{$taxshort}++;
    }
    # Sort output in descending order
    foreach my $key (sort {$taxhash{$b} <=> $taxhash{$a}} keys %taxhash) {
        $output{$key} = $taxhash{$key};
    }
    return (\%output); # Return reference to output array
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

    my $args;
    if ($sc == 1) {
      $args = '--sc ';
    }
    if ($SEmode == 1) {
        $args = $args."-s $libraryNAME.$readsf.SSU.1.fq";
    } else {
        $args = $args."-1 $libraryNAME.$readsf.SSU.1.fq -2 $libraryNAME.$readsf.SSU.2.fq";
    }

    run_prog("spades",
             "-o $libraryNAME.spades -t $cpus -m 20 -k $kmer "
             . $args,
             "$libraryNAME.spades.out","&1"
         );

    msg("done...");
}

sub spades_parse {
    # getting spades output and reformatting it...
    msg("getting SSU phylotypes and their coverages...");

    #spades scaffolds file is empty, settting skip_spades variable to avoid
    #further processing
    if (! -s "$libraryNAME.spades/scaffolds.fasta") {
        msg("no phylotypes assembled with SPAdes");
        $skip_spades = 1;
        system ("rm ./$libraryNAME.spades -r");
    }

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
                 "$libraryNAME.scaffolds.$_.gff",
                 "$libraryNAME.barrnap.out");
    }

    # now merge multi-hits on the same scaffold-and-strand by picking
    # the lowest start and highest stop position.

    my %ssus;
    # pre-filter with grep for speed
    my $fh;
    open_or_die(\$fh, "-|",
                "grep -hE '16S_rRNA\|18S_rRNA' ".
                "$libraryNAME.scaffolds.bac.gff ".
                "$libraryNAME.scaffolds.arch.gff ".
                "$libraryNAME.scaffolds.euk.gff");
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
	return;
    }

    open_or_die(\$fh, ">", "tmp.$libraryNAME.scaffolds.final.gff");
    for my $key (sort keys %ssus) {
        print $fh join("\t",@{$ssus{$key}});
    }
    close($fh);

    # fastaFromBed will build a .fai index from the source .fasta
    # However, it does not notice if the .fasta changed. So we
    # delete the .fai if it is older than the .fasta.
    if ( -e "$libraryNAME.spades/scaffolds.fasta.fai" &&
         file_is_newer("$libraryNAME.spades/scaffolds.fasta",
                       "$libraryNAME.spades/scaffolds.fasta.fai")) {
        unlink("$libraryNAME.spades/scaffolds.fasta.fai");
    }

    # extract rrna fragments from spades scaffolds accoding to gff
    run_prog("fastaFromBed",
             "  -fi $libraryNAME.spades/scaffolds.fasta "
             . "-bed tmp.$libraryNAME.scaffolds.final.gff "
             . "-fo $libraryNAME.spades_rRNAs.final.fasta "
             . "-s -name",
             "tmp.$libraryNAME.fastaFromBed.out",
             "&1");

    msg("done...");
}

sub bbmap_spades_out {
    # Map extracted reads back to the SPAdes output to see what proportion
    # of reads can be explained by the assembled sequences
    msg("mapping extracted SSU reads back on assembled SSU sequences");
    my $args = ""; # Check whether running in PE or SE mode
    if ($SEmode == 0) {
	if ($interleaved == 1) {
	    $args =
	    "  interleaved=t "
	    . "pairlen=$maxinsert ";
	} else {
	    $args =
	    "  in2=$libraryNAME.$readsf.SSU.2.fq "
	    . "pairlen=$maxinsert ";
	}
    }
    run_prog("bbmap",
         "  fast=t "
         . "minidentity=0.98 "
         . "-Xmx20g "
         . "threads=$cpus "
         . "po=f "
         . "outputunmapped=t " # This is important
         . "ref=$libraryNAME.spades_rRNAs.final.fasta "
         . "nodisk "
         . "in=$libraryNAME.$readsf.SSU.1.fq "
         . "out=$libraryNAME.$readsf.SSU_assem.sam " # Also skip SAM header?
         . $args,
         undef,
         "$libraryNAME.remap.bbmap.out");
    msg("done...");
}

sub taxonomy_spades_unmapped {
    # Filter output of mapping to SPAdes assembled SSU sequences
    # Report taxonomy of "leftover" sequences
    my $in = "$libraryNAME.$readsf.SSU_assem.sam";  # mapping of extracted reads vs assembled SSU seq
    my $sam_href = \%ssu_sam;                       # data from first mapping vs. SILVA, read into memory
    my $stats_href = \%ssu_sam_mapstats;
    my $out_href = \%taxa_from_hitmaps_unassem;     # hash to store taxonomy results from unassembled reads
    my $out_sorted_aref = \@taxa_from_hitmaps_unassem_sorted; # array of sorted taxonomy results from unassembled reads
    my $cov_href = \%ssuassem_cov;                  # Hash to store read coverage of assembled SSU sequences
    my @taxa_full;
    msg ("extracting taxonomy of unassembled SSU reads");
    my $fh;
    open_or_die(\$fh, "<", $in);
    while (my $line = <$fh>) {
        next if ( $line =~ m/^@/ ); # Skip headers
        my ($read, $bitflag, $ref, @discard) = split /\t/, $line;
        # Check whether fwd or rev read
        my $pair;
        if ($bitflag & 0x1) { # If PE read
            if ($SEmode != 0) { # Sanity check
                msg ("ERROR: bitflag in SAM file conflicts with SE mode flag ");
            }
            if ($bitflag & 0x40) {
                $pair="F";
            } elsif ($bitflag & 0x80) {
                $pair="R";
            }
        } else {
            $pair = "U";
        }
        # Check whether read mapped to assembled SSU
        if ($bitflag & 0x4) { # If not mapped, check corresponding read in original mapping
            if (defined ${$sam_href}{$read}{$pair}{"ref"}) {
                my $ref = ${$sam_href}{$read}{$pair}{"ref"};
                # Shorten taxonomy string and save into NTU table
                if (${$sam_href}{$read}{$pair}{"ref"} =~ m/\w+\.\d+\.\d+\s(.+)/) {
                    my $taxonlongstring = $1;
                    push @taxa_full, $taxonlongstring;
                    # Truncate to 6 levels
                    my $taxonshortstring = truncate_taxonstring($taxonlongstring, 6);           # To be superseded
                    # Increment count for taxon
                    ${$out_href}{$taxonshortstring}++;                                          # To be superseded
                } else {
                    msg ("warning: malformed database entry $ref");
                }
            }
        } else {
            # Add to read count
            my ($refshort) = $ref =~ /($libraryNAME\.PF\w+)_[\d\.]+/;
            ${$cov_href}{$refshort}++;
            if ($pair eq "F" | $pair eq "U") {
                ${$stats_href}{"assem_fwd_map"}++;
            } elsif ($pair eq "R") {
                ${$stats_href}{"assem_rev_map"}++;
            }

        }
    }

    # Sort taxon strings for unassembled reads affiliation
    my @discard2;
    @{$out_sorted_aref} =
        sort { @$b[1] <=> @$a[1] }
        grep { if (@$_[1] < 3) { @discard2[@$_[1]-1]++;} @$_[1] > 2 }
        map  { [$_, ${$out_href}{$_}] }
        keys %{$out_href};

    close($fh);

    # Summarize taxonomy
    $taxa_unassem_summary_href = summarize_taxonomy(\@taxa_full, $taxon_report_lvl); # Summarize

    # Calculate ratio of reads assembled
    my $assem_tot_map = ${$stats_href}{"assem_fwd_map"};
    if (defined ${$stats_href}{"assem_rev_map"}) {
        $assem_tot_map += ${$stats_href}{"assem_rev_map"};
    }
    my $ssu_tot_map = ${$stats_href}{"ssu_fwd_map"};
    if (defined ${$stats_href}{"ssu_rev_map"}) {
        $ssu_tot_map += ${$stats_href}{"ssu_rev_map"};
    }
    ${$stats_href}{"assem_ratio"} = $assem_tot_map / $ssu_tot_map;

    msg ("done ... ");
}

sub emirge_run {
    # running Emirge on bbmap output
    msg("creating phylotypes with Emirge");

    my $cmd = "emirge";
    my $args = "-1 $libraryNAME.$readsf.SSU.1.fq ";

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

        msg("the insert size used is $ins_used +- $ins_std");
        # FIXME: EMIRGE dies with too many SSU reads, the cutoff needs to be adjusted...
        if ($SSU_total_pairs <= $amplimit) {
            msg("Less than $amplimit SSU read pairs - using Emirge");
        } else {
            $cmd = "emirge_amp";
            msg("Warning: More than $amplimit SSU reads - using Emirge Amplicon");
        }

        $args = "  -1 $libraryNAME.$readsf.SSU.1.fq "
                . "-2 $libraryNAME.$readsf.SSU.2.fq "
                . "-i $ins_used -s $ins_std ";
    } else {
        msg("reads > 151 bp - running in single end mode");
        run_prog("cat",
                 "$libraryNAME.$readsf.SSU.1.fq $libraryNAME.$readsf.SSU.2.fq",
                 "tmp.$libraryNAME.SSU_all.fq");
        # ename the reads with a running number to make emirge happy
        # using awk for speed, these files can be huge

        run_prog("awk",
                 "\'{print (NR%4 == 1) ? \"\@\" ++i  : (NR%4 == 3) ? \"+\" :\$0}\'"
                 . " tmp.$libraryNAME.SSU_all.fq",
                 "tmp.$libraryNAME.renamed.SSU_all.fq");

        $args = " -1 tmp.$libraryNAME.renamed.SSU_all.fq ";
    }

    run_prog($cmd,
             " $libraryNAME.emirge "
             . $args
             . " -f ${DBHOME}/${emirge_db}.fasta"
             . " -b ${DBHOME}/${emirge_db}.bt "
             . " -l $readlength -a $cpus --phred33 ",
             , "$libraryNAME.emirge.out", "&1");

    msg("done...");
}

sub emirge_parse {
    # getting emirge output and reformatting it...
    msg("getting Emirge phylotypes and their abundances...");
    run_prog("emirge_rename_fasta",
             "./$libraryNAME.emirge/iter.40/",
             "tmp.$libraryNAME.emirge.result.fasta");

    my $fh_in;
    my $fh_out;
    open_or_die(\$fh_in, "<","tmp.$libraryNAME.emirge.result.fasta");
    open_or_die(\$fh_out, ">", "$libraryNAME.emirge.final.fasta");
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

    # join emirge and spades hits into one file
    # (vsearch takes a while to load, one run saves us time)
    if ($skip_spades == 1 && $skip_emirge == 0) {
      run_prog("cat",
               "   $libraryNAME.emirge.final.fasta",
               "$libraryNAME.all.final.fasta");
    }
    elsif ($skip_emirge == 1 && $skip_spades == 0){
      run_prog("cat",
               "   $libraryNAME.spades_rRNAs.final.fasta",
               "$libraryNAME.all.final.fasta");
    }
    elsif ($skip_emirge == 0 && $skip_spades == 0){
      run_prog("cat",
             "   $libraryNAME.emirge.final.fasta"
             . " $libraryNAME.spades_rRNAs.final.fasta",
             "$libraryNAME.all.final.fasta");
    }

    if (-s "$libraryNAME.all.final.fasta") {
       run_prog("vsearch",
             "-usearch_global $libraryNAME.all.final.fasta "
             . "-db ${DBHOME}/${vsearch_db}.fasta "
             . "-id 0.7 "
             . "-userout $libraryNAME.all.vsearch.csv "
             . "-userfields query+target+id+alnlen+evalue+id3+qs+pairs+gaps+mism+ids "
             . "-threads $cpus --strand plus --notrunclabels "
             . "-notmatched $libraryNAME.all.final.phyloFlash.notmatched.fa "
             . "-dbmatched $libraryNAME.all.final.phyloFlash.dbhits.fa ",
             "tmp.$libraryNAME.all.vsearch.out",
             "&1");

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
    open_or_die(\$fh, "<", "$libraryNAME.all.vsearch.csv");
    while (<$fh>) {
        chomp;
        # lib.PFspades_1_1.23332\tAH12345.1.1490 Bacteria;...\t...
        s/PF(\w+)_([^_]+)_([^\t]+)\t(\w+\.\d+\.\d+)\s/PF$1_$2\t$3\t$4\t/;
        push @vsearch_matches, [split("\t",$_)];
    }
    close($fh);

    # Sort numerically descending the list of
    # assembled SSU sequences by coverage value

    @ssuassem_results_sorted =
        sort { @{$b}[1] <=> @{$a}[1] or @{$a}[3] cmp @{$b}[3] }
        grep { @{$_}[0] =~ /^$libraryNAME.PFspades/ }
        @vsearch_matches;

    # Sort numerically descending the list of reconstructed
    # SSU sequences by mapping ratio
    @ssurecon_results_sorted =
        sort { @{$b}[1] <=> @{$a}[1] or @{$a}[3] cmp @{$b}[3] }
        grep { @{$_}[0] =~ /^$libraryNAME.PFemirge/ }
        @vsearch_matches;

    msg("done...");
}

sub vsearch_cluster {
    my $clusterid = 97;
    msg("clustering DB matches at $clusterid%");
    run_prog("vsearch",
             "  --cluster_fast $libraryNAME.all.final.phyloFlash.dbhits.fa "
             . "-id ".($clusterid/100)." "
             . "-centroids $libraryNAME.all.dbhits.NR97.fa "
             . "-notrunclabels "
             . "--threads $cpus ",
             "tmp.$libraryNAME.clusterdbhits.out",
             "&1");

    msg("done...");
}

sub mafft_run {
    msg("creating alignment and tree...");

    if ($skip_spades == 0 && $skip_emirge == 1) {
        run_prog("cat",
                 "$libraryNAME.all.dbhits.NR97.fa "
                 . "$libraryNAME.spades_rRNAs.final.fasta ",
                 "$libraryNAME.SSU.collection.fasta");

    } elsif ($skip_spades == 1 && $skip_emirge == 0) {
        run_prog("cat",
                 "$libraryNAME.all.dbhits.NR97.fa "
                 . "$libraryNAME.emirge.final.fasta ",
                 "$libraryNAME.SSU.collection.fasta");

    } elsif ($skip_spades == 0 && $skip_emirge == 0) {
        run_prog("cat",
                 "$libraryNAME.all.dbhits.NR97.fa "
                 . "$libraryNAME.spades_rRNAs.final.fasta "
                 . "$libraryNAME.emirge.final.fasta ",
                 "$libraryNAME.SSU.collection.fasta");
    }

    run_prog("mafft",
             "--treeout $libraryNAME.SSU.collection.fasta ",
             "$libraryNAME.SSU.collection.alignment.fasta",
             "tmp.$libraryNAME.SSU.collection.alignment.mafftout");

    # fix missing ; at and of MaFFT newick tree
    my $fh;
    open_or_die(\$fh, ">>", "$libraryNAME.SSU.collection.fasta.tree");
    print {$fh} ";";
    close($fh);

    msg("done...");
}

sub clean_up {
    # cleanup of intermediate files and folders
    msg("cleaning temp files...");
    if ($skip_spades == 0) {
        system ("rm ./$libraryNAME.spades -r");
    }
    if ($skip_emirge == 0) {
        system ("rm ./$libraryNAME.emirge -r");
    }
    system ("rm tmp.$libraryNAME.* -r");
    msg("done...");
}

sub run_plotscript_SVG {
    msg ("generating histogram and tree graphics in SVG format");
    # Plot mapping ID histogram
    run_prog("plotscript_SVG",
         "--hist $libraryNAME.idhistogram "
         ."--title=\"Mapping identity (%)\" ",
         #. "$decimalcomma ",
         "tmp.$libraryNAME.plotscript.out",
         "&1");

    # Plot insert size histogram unless running in SE mode
    if ($SEmode != 1) { # If not running in SE mode ...
        run_prog("plotscript_SVG",
                 "--hist $libraryNAME.inserthistogram "
                 ."--title=\"Insert size (bp)\" ",
                 #. "$decimalcomma ",
                 "tmp.$libraryNAME.plotscript.out",
                 "&1");
    }

    # Plot tree if spades/emirge unless both skipped
    unless ($skip_spades + $skip_emirge == 2) {
        run_prog("plotscript_SVG",
                 "--tree $libraryNAME.SSU.collection.fasta.tree ",
                 #. "$decimalcomma ",
                 "tmp.$libraryNAME.plotscript.out",
                 "&1");
    }

    # Generate barplot of taxonomy at level XX
    run_prog("plotscript_SVG",
             "--bar $libraryNAME.phyloFlash.taxonsummary.csv "
             ."--title=\"Taxonomic summary from mapping to database\" ",
             "tmp.$libraryNAME.plotscript.out",
             "&1");
}

sub generate_treemap_data_rows {
  # Generate data rows for drawing treemap chart
  my %parents_HoH; # Hash of count vals by parent and child taxa
  my @output; # Array to store output
  # Parse taxstrings into hash of child-parent relationships

  foreach my $taxstring (keys %taxa_from_hitmaps) {
    my @taxsplit = split ";", $taxstring; # Get taxonomy string, split by semicolons
    my $taxlen = scalar @taxsplit; # Get length of tax string
    while (scalar @taxsplit > 1) {
      my $countval = 0; # Initialize count value as dummy "zero" by default
      if (scalar @taxsplit == $taxlen) { # If leaf node, set count value to real value
          $countval = $taxa_from_hitmaps{$taxstring};
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
  push @output,"[\'Eukaryota\',\t\'Cellular organisms\',\t0],\n";
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

sub write_report_html_new {
    # Get input parameters
    my ($version,
        $libraryNAME, $id, $SEmode,
        $readsf_full,$readsr_full,
        $readnr,$readnr_pairs,
        $ins_me,$ins_std,$ins_used,
        $SSU_ratio, $SSU_ratio_pc, $SSU_total_pairs,
        $skip_spades, $skip_emirge, $treemap_flag,
        $taxon_report_lvl,
        $mapstats_href,
        $taxa_summary_href,
        $taxa_unassem_summary_href,
        $ssuassem_results_sorted_aref,
        $ssurecon_results_sorted_aref,
        ) = @_;
    # Location of template file
    my $template = "$FindBin::RealBin/phyloFlash_report_template.html",
    msg ("Generating HTML-formatted report and graphics ... ");
    my %flags;  # Hash to store substitution flags in template
    # Hash parameters that are compulsory
    %flags = (
        "VERSION" => $version,
        "LIBNAME" => $libraryNAME,
        "ID" => $id,
        "READSF_FULL" => $readsf_full,
        "READNR" => $readnr,
        "IDHISTOGRAM" => "<embed width=240 height=240 src=\"".$libraryNAME.".idhistogram.svg\" />\n",
        #MAPRATIOPIE
        "TAXONSUMMARYBAR" => "<embed width=500 src=\"".$libraryNAME.".phyloFlash.taxonsummary.csv.svg\" />\n",
        "SSU_RATIO" => $SSU_ratio, # Not in report?
        "SSU_RATIO_PC" => $SSU_ratio_pc,
        "SSU_TOTAL_PAIRS" => $SSU_total_pairs,
        "TAX_REPORT_LVL" => $taxon_report_lvl,
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
        $flags{"INSERTHISTOGRAM"} = "<embed width=240 height=240 src=\"".$libraryNAME.".inserthistogram.svg\" />\n";
        $flags{"INS_ME"} = $ins_me;
        $flags{"INS_STD"} = $ins_std;
        $flags{"READSR_FULL"} = $readsr_full;
        $flags{"READNR_PAIRS"} = $readnr_pairs;
    }
    # Params defined only for assembled (SPAdes) reads
    if ($skip_spades == 0) {
        #$flags{"ASSEMBLYRATIOPIE"};
        $flags{"ASSEM_RATIO"} = $mapstats_href->{"assem_ratio"};
        $flags{"SEQUENCES_TREE"} = "<embed src=\"".$libraryNAME.".SSU.collection.fasta.tree.svg\" width=1200px />\n";

        # Table of assembled SSU sequences
        my @table_assem_seq;
        foreach (@$ssuassem_results_sorted_aref) {
            my @split_entry = @$_;
            # Parse the database entry number of the reference sequence to get Genbank accession
            my @get_genbank = split (/\./,$split_entry[2]);
            push @table_assem_seq, "  <tr>\n";
            push @table_assem_seq, "    <td>".$split_entry[0]."</td>\n";
            push @table_assem_seq, "    <td>".$ssuassem_cov{$split_entry[0]}."</td>\n";
            push @table_assem_seq, "    <td>".$split_entry[1]."</td>\n";
            push @table_assem_seq, "    <td><a href=\"http://www.ncbi.nlm.nih.gov/nuccore/".$get_genbank[0]."\">".$split_entry[2]."</a></td>\n";	# Link to Genbank entry using accession no.
            push @table_assem_seq, "    <td>".$split_entry[3]."</td>\n";
            push @table_assem_seq, "    <td>".$split_entry[4]."</td>\n";
            push @table_assem_seq, "    <td>".$split_entry[5]."</td>\n";
            push @table_assem_seq, "    <td>".$split_entry[6]."</td>\n";
            push @table_assem_seq, "  </tr>\n";
        }
        $flags{"ASSEMBLED_SSU_TABLE"} = join "", @table_assem_seq;

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
    # Params defined only for reconstructed (EMIRGE) reads
    if ($skip_emirge == 0) {
        $flags {"INS_USED"} = $ins_used;
        my @table_recon_seq;
        foreach (@ssurecon_results_sorted) {
            push @table_recon_seq, "  <tr>\n";
            my @split_entry = @$_;
            my $test = $split_entry[2];
            my @get_genbank = split (/\./,$test);
            push @table_recon_seq, "    <td>".$split_entry[0]."</td>\n";
            push @table_recon_seq, "    <td>".$split_entry[1]."</td>\n";
            push @table_recon_seq, "    <td><a href=\"http://www.ncbi.nlm.nih.gov/nuccore/".$get_genbank[0]."\">".$split_entry[2]."</a></td>\n";
            push @table_recon_seq, "    <td>".$split_entry[3]."</td>\n";
            push @table_recon_seq, "    <td>".$split_entry[4]."</td>\n";
            push @table_recon_seq, "    <td>".$split_entry[5]."</td>\n";
            push @table_recon_seq, "    <td>".$split_entry[6]."</td>\n";
            push @table_recon_seq, "  </tr>\n";
        }
        $flags{"EMIRGE_TABLE"} = join "", @table_recon_seq;
    }

    # Open template and process output
    my ($fh_in, $fh_out);
    open_or_die(\$fh_in, "<", $template);
    open_or_die(\$fh_out, ">", "$libraryNAME.phyloFlash.newreport.html");
    my $write_flag = 1;
    while (my $line = <$fh_in>) {
        chomp $line;
        if ($line =~ m/<!--(.+)-->/) {
            my $flag = $1;
            if (defined $flags{$flag}) {
                my $replace = $flags{$flag};
                $line =~ s/<!--$flag-->/$replace/;
            }
            # Turn off writing if option not specified
            $write_flag = 0 if (defined $suppress_flags{$flag});
            # Turn writing back on if optional section is finished
            $write_flag = 1 if (defined $suppress_end_flags{$flag});
        }
        print $fh_out $line."\n" if ($write_flag == 1);
    }
    close($fh_out);
    close($fh_in);

}

sub write_report_html {
    # Generate HTML-formatted report file -- lots of blocks of verbatim HTML
    # in this section
    msg("Generating HTML-formatted report and graphics...");

    my $fh;
    open_or_die(\$fh, ">", "$libraryNAME.phyloFlash.html");

    print {$fh} <<ENDHTML;
<!DOCTYPE html>
<html>
<head>
ENDHTML

print {$fh} "  <title>phyloFlash results summary for library ".$libraryNAME."</title>\n";

print {$fh} <<ENDHTML;
  <script language="javascript" type="text/javascript">
  // adapted from http://www.cssnewbie.com/example/showhide-content/
  function showHide(shID) {
    if (document.getElementById(shID)) {
      if (document.getElementById(shID).style.display == 'none') {
        document.getElementById(shID).style.display = 'block';
      }
      else {
        document.getElementById(shID).style.display = 'none';
      }
    }
  }
  </script>
ENDHTML

if ($treemap_flag == 1) { # Script for interactive treemap if flag is on
print {$fh} <<ENDHTML;
  <script type="text/javascript" src="https://www.gstatic.com/charts/loader.js"></script>
  <script type="text/javascript">
  google.charts.load('current', {'packages':['treemap']});
  google.charts.setOnLoadCallback(drawChart);
  function drawChart() {
    var data = new google.visualization.DataTable();
    data.addColumn('string','ID');
    data.addColumn('string','Parent');
    data.addColumn('number','NumReads');
    data.addRows([
ENDHTML
my @output_array = generate_treemap_data_rows();
while (scalar @output_array > 0) {
  print {$fh} shift @output_array;
}
print {$fh} <<ENDHTML;
    ]);

    tree = new google.visualization.TreeMap(document.getElementById('chart_div'));

    tree.draw(data, {
      minColor: '#009688',
      midColor: '#f7f7f7',
      maxColor: '#ee8100',
      highlightOnMouseOver: true,
      minHighlightColor: '#ee8100',
      midHighlightColor: '#f7f7f7',
      maxHighlightColor: '#009688',
      maxPostDepth: 2,
      headerHeight: 15,
      fontColor: 'black',
      showScale: false,
      generateTooltip: showFullTooltip
    });

    function showFullTooltip(row,size,value) {
      return '<div style="background:white; padding:10px; border-style:none">' +
      '<b>' + data.getValue(row, 0) + '</b><br>' +
      'Total reads: ' + size + '</div>';
    }

  }
</script>
ENDHTML
}

print {$fh} <<ENDHTML;
  <style type="text/css">
    body {
      font-family: "Helvetica", "Gill Sans", "Gill Sans MT", sans-serif;}

    th {
      background-color: #dee;
      padding: .5em .5em .5em .5em;
      font-weight: bold;}
    td {
      padding: 0em .5em 0em .5em;}

    table.slimTable th {
      padding: 0em .5em 0em .5em;
      font-weight: normal;
      text-align: right;}

    .withHoverText {
      text-decoration: none;
      border-bottom: 1px grey dotted;}

    .more {
      display: none;}
    a.showLink {
      text-decoration: none;}
    a.showLink:link {
      color: grey;}
    a.showLink:visited {
      color: grey;}
    a.showLink:hover {
      color: white;
      background-color: grey;}

  </style>
</head>
<body>

<h3><a href="https://github.com/HRGV/phyloFlash">
ENDHTML
print {$fh} $version; # Print phyloFlash name and version number
print {$fh} <<ENDHTML;
 </a> by <a href="mailto:hgruber\@mpi-bremen.de">Harald Gruber-Vodicka</a> - high throughput phylogenetic screening using SSU rRNA gene(s) abundance(s)</h3>
 <p>Click on report section headers to expand, mouse-over underlined text to see explanations.</p>
ENDHTML

print {$fh} "<h1>Library name: ".$libraryNAME."</h1>\n";

print {$fh} <<ENDHTML;
<h2>Parameters</h2>

<h3>Input files</h3>
<table class="slimTable">
  <tr>
    <th>Forward read file</th>
ENDHTML

print {$fh} "    <td>".$readsf_full."</td>\n";

print {$fh} <<ENDHTML;
  </tr>
  <tr>
    <th>Reverse read file</th>
ENDHTML

print {$fh} "    <td>".$readsr_full."</td>\n";

print {$fh} <<ENDHTML;
  </tr>
</table>

<h3>Mapping parameters</h3>
<table class="slimTable">
  <tr>
    <th><span class="withHoverText" title="Minimum sequence identity for mapping to SSU database to be accepted">Minimum mapping identity</span></th>
ENDHTML

print {$fh} "    <td>".$id."\%</td>\n";
if ($SEmode == 0) {
print {$fh} <<ENDHTML;
  </tr>
  <tr>
    <th>Input PE-reads</th>
ENDHTML

print {$fh} "    <td>".$readnr_pairs."</td>\n";

print {$fh} <<ENDHTML;
  </tr>
  <tr>
    <th>Mapped SSU read pairs</th>
ENDHTML

print {$fh} "    <td>".$SSU_total_pairs."</td>\n";
} else
{
print {$fh} <<ENDHTML;
  </tr>
  <tr>
    <th>Input SE-reads</th>
ENDHTML

print {$fh} "    <td>".$readnr."</td>\n";

print {$fh} <<ENDHTML;
  </tr>
  <tr>
    <th>Mapped SSU SE-reads</th>
ENDHTML

print {$fh} "    <td>".$SSU_total_pairs."</td>\n";
}
print {$fh} <<ENDHTML;
  </tr>
  <tr>
    <th><span class="withHoverText" title="Fraction of library mapping to SSU references">Mapping ratio</span></th>
ENDHTML

print {$fh} "    <td>".$SSU_ratio_pc."\%</td>\n";
print {$fh} <<ENDHTML;
  </tr>
ENDHTML

if (defined $ssu_sam_mapstats{"assem_ratio"}) {
print {$fh} <<ENDHTML;
  <tr>
    <th><span class="withHoverText" title="Fraction of extracted sequences assembled">Fraction assembled</span></th>
ENDHTML
    print {$fh} "    <td>".$ssu_sam_mapstats{"assem_ratio"}."</td>\n";
    print {$fh} <<ENDHTML;
  </tr>
ENDHTML
}

if ($SEmode==0) { # Only print insert size information if in paired-end mode
print {$fh} <<ENDHTML;
  <tr>
    <th>Detected median insert size</th>
ENDHTML

print {$fh} "    <td>".$ins_me."</td>\n";
print {$fh} <<ENDHTML;
  </tr>
ENDHTML

if ($skip_emirge == 0) {
print {$fh} <<ENDHTML;
  <tr>
    <th>Used insert size</th>
ENDHTML

print {$fh} "     <td>".$ins_used."</td>\n";
print {$fh} <<ENDHTML;
  </tr>
ENDHTML
}

print {$fh} <<ENDHTML;
  <tr>
    <th>Insert size standard deviation</th>
ENDHTML

print {$fh} "    <td>".$ins_std."</td>\n";
print {$fh} <<ENDHTML;
  </tr>
ENDHTML
}
print {$fh} <<ENDHTML;
</table>
ENDHTML

if ($SEmode==0) { # Only display insert size histogram if in paired-end mode
print {$fh} <<ENDHTML;
<h3><a href="#" id="histo-show" class="showLink" onclick="showHide('histo');return false;" title="Click to expand">Insert size histogram</a></h3>
<div id="histo" class="more">
<p>Insert sizes for read pairs. Distribution should generally be unimodal; more than one peak may indicate contamination from other libraries.</p>
ENDHTML

print {$fh} "<embed width=240 height=240 src=\"".$libraryNAME.".inserthistogram.svg\" />\n";
print {$fh} <<ENDHTML;
</div>
ENDHTML
}

print {$fh} <<ENDHTML;

<h3><a href="#" id="idhisto-show" class="showLink" onclick="showHide('idhisto');return false;" title="Click to expand">Mapping identity histogram</a></h3>
<div id="idhisto" class="more">
<p>Read-mapping %identity of reads vs. reference database. Lower %identity hits may indicate presence of divergent taxa not represented in the database.</p>
ENDHTML

print {$fh} "<embed width=240 height=240 src=\"".$libraryNAME.".idhistogram.svg\" />\n";
print {$fh} <<ENDHTML;
</div>

<h3>Output files</h3>
<table class="slimTable">
ENDHTML

print {$fh} "  <tr>\n";
print {$fh} "    <th>PhyloFlash report (plaintext)</th>\n";
print {$fh} "    <td>".$libraryNAME.".phyloFlash</td>\n";
print {$fh} "  </tr>\n";

print {$fh} "  <tr>\n";
print {$fh} "    <th>PhyloFlash report (HTML)</th>\n";
print {$fh} "    <td>".$libraryNAME.".phyloFlash.html</td>\n";
print {$fh} "  </tr>\n";

if ($SEmode==0) { # Only dislpay insert size information if in paired-end mode
print {$fh} "  <tr>\n";
print {$fh} "    <th>Insert size histogram (plaintext)</th>\n";
print {$fh} "    <td>".$libraryNAME.".inserthistogram</td>\n";
print {$fh} "  </tr>\n";

print {$fh} "  <tr>\n";
print {$fh} "    <th>Insert size histogram (graphics)</th>\n";
print {$fh} "    <td>".$libraryNAME.".inserthistogram.svg</td>\n";
print {$fh} "  </tr>\n";
}

print {$fh} "  <tr>\n";
print {$fh} "    <th>Mapping identity histogram (plaintext)</th>\n";
print {$fh} "    <td>".$libraryNAME.".idhistogram</td>\n";
print {$fh} "  </tr>\n";

print {$fh} "  <tr>\n";
print {$fh} "    <th>Mapping identity histogram (graphics)</th>\n";
print {$fh} "    <td>".$libraryNAME.".idhistogram.svg</td>\n";
print {$fh} "  </tr>\n";

if ($skip_spades + $skip_emirge < 2) {

print {$fh} "  <tr>\n";
print {$fh} "    <th>Tree of recovered SSU sequences (Newick)</th>\n";
print {$fh} "    <td>".$libraryNAME.".SSU.collection.fasta.tree</td>\n";
print {$fh} "  </tr>\n";

print {$fh} "  <tr>\n";
print {$fh} "    <th>Tree of recovered SSU sequences (graphics)</th>\n";
print {$fh} "    <td>".$libraryNAME.".SSU.collection.fasta.tree.svg</td>\n";
print {$fh} "  </tr>\n";

}

print {$fh} <<ENDHTML;
</table>

<h2>Results</h2>
ENDHTML

print {$fh} <<ENDHTML;
<h3><a href="#" id="taxa-show" class="showLink" onclick="showHide('taxa');return false;">Taxonomic affiliation of SSU rRNA reads in library</a></h3>
<div id="taxa" class="more">
ENDHTML
print {$fh} "<p>Approximate overview of taxonomic composition, by mapping reads to database sequences and summarizing their taxonomy strings to taxon level $taxon_report_lvl.</p>\n";
print {$fh} "<embed width=500 src=\"".$libraryNAME.".phyloFlash.taxonsummary.csv.svg\" />\n";
print {$fh} "</div>\n";

if ($skip_spades == 0) {
## write list of assembled SSU sequences
print {$fh} <<ENDHTML;
<h3><a href="#" id="spades-show" class="showLink" onClick="showHide('spades');return false;">SSU assembly based taxa</a></h3>
<div id="spades" class="more">
<table>
  <tr>
    <th><span class="withHoverText" title="Name assigned to assembled SSU sequence">OTU</span></th>
    <th><span class="withHoverText" title="Read coverage of assembled SSU sequence, from mapping">read_cov</span></th>
    <th><span class="withHoverText" title="Per-base kmer coverage of assembled SSU sequence, from assembler">coverage</span></th>
    <th><span class="withHoverText" title="Closest hit in SSU reference database">dbHit</span></th>
    <th>taxonomy</th>
    <th>\%id</th>
    <th><span class="withHoverText" title="Length of pairwise alignment between assembled SSU sequence and reference">alnlen</span></th>
    <th>evalue</th>
  </tr>
ENDHTML

foreach (@ssuassem_results_sorted) {
    my @split_entry = @$_;
    my @get_genbank = split (/\./,$split_entry[2]);
        # Parse the database entry number of the reference sequence to get Genbank accession
    print {$fh} "  <tr>\n";
    print {$fh} "    <td>".$split_entry[0]."</td>\n";
    print {$fh} "    <td>".$ssuassem_cov{$split_entry[0]}."</td>\n";
    print {$fh} "    <td>".$split_entry[1]."</td>\n";
    print {$fh} "    <td><a href=\"http://www.ncbi.nlm.nih.gov/nuccore/".$get_genbank[0]."\">".$split_entry[2]."</a></td>\n";	# Link to Genbank entry using accession no.
    #	print {$fh} "    <td>".$split_entry[2]."</td>\n";
    print {$fh} "    <td>".$split_entry[3]."</td>\n";
    print {$fh} "    <td>".$split_entry[4]."</td>\n";
    print {$fh} "    <td>".$split_entry[5]."</td>\n";
    print {$fh} "    <td>".$split_entry[6]."</td>\n";
    print {$fh} "  </tr>\n";
}

## write table of taxonomic affiliations for unassembled reads
print {$fh} <<ENDHTML;
</table>
</div>
<h3><a href="#" id="unassemtaxa-show" class="showLink" onclick="showHide('unassemtaxa');return false;">Taxonomic affiliation of unassembled SSU rRNA reads</a></h3>
<div id="unassemtaxa" class="more">
<p>Approximate overview of taxonomic composition for reads that did NOT assemble into full-length sequences. Only displaying taxa with > 3 reads mapped.</p>
<table>
  <tr>
    <th><span class="withHoverText" title="Higher taxon found by SSU mapping to SILVA reference database">Taxon</span></th>
    <th><span class="withHoverText" title="No. reads (single) mapped to this taxonomic group">Reads</span></th>
  </tr>
ENDHTML

## write list of higher taxa found
my @taxsort = sort {${$taxa_unassem_summary_href}{$b} <=> ${$taxa_unassem_summary_href}{$a}} keys %$taxa_unassem_summary_href;
foreach my $uatax (@taxsort) {
    if (${$taxa_unassem_summary_href}{$uatax} > 3) { # Only display those with at least three reads mapping
        # For each of these higher taxa
        print {$fh} "  <tr>\n";
        # print name of taxon
        print {$fh} "    <td>".$uatax."</td>\n";
        # print no . of reads mapping
        print {$fh} "    <td>".${$taxa_unassem_summary_href}{$uatax}."</td>\n";
        print {$fh} "  </tr>\n";
    }
}

print {$fh} <<ENDHTML;
</table>
</div>
ENDHTML

}

if ($skip_emirge == 0) {
## write list of reconstructed SSU sequences
print {$fh} <<ENDHTML;
</table>
</div>

<h3><a href="#" id="emirge-show" class="showLink" onclick="showHide('emirge');return false;">SSU reconstruction based taxa</a></h3>
<div id="emirge" class="more">
<table>
  <tr>
    <th>OTU</th>
    <th>ratio</th>
    <th>dbHit</th>
    <th>taxonomy</th>
    <th>\%id</th>
    <th>alnlen</th>
    <th>evalue</th>
  </tr>
ENDHTML

foreach (@ssurecon_results_sorted) {
    print {$fh} "  <tr>\n";
    my @split_entry = @$_;
    my $test = $split_entry[2];
    my @get_genbank = split (/\./,$test);
    print {$fh} "    <td>".$split_entry[0]."</td>\n";
    print {$fh} "    <td>".$split_entry[1]."</td>\n";
    print {$fh} "    <td><a href=\"http://www.ncbi.nlm.nih.gov/nuccore/".$get_genbank[0]."\">".$split_entry[2]."</a></td>\n";
#	print {$fh} "    <td>".$split_entry[2]."</td>\n";
    print {$fh} "    <td>".$split_entry[3]."</td>\n";
    print {$fh} "    <td>".$split_entry[4]."</td>\n";
    print {$fh} "    <td>".$split_entry[5]."</td>\n";
    print {$fh} "    <td>".$split_entry[6]."</td>\n";
    print {$fh} "  </tr>\n";
}
}

if ($skip_spades + $skip_emirge < 2) {
print {$fh} <<ENDHTML;
</table>
</div>
<h3><a href="#" id="tree-show" class="showLink" onclick="showHide('tree');return false;">Combined tree of sequences</a></h3>
<div id="tree" class="more">
ENDHTML

print {$fh} "<embed src=\"".$libraryNAME.".SSU.collection.fasta.tree.svg\" width=1200px />\n";
}

print {$fh} <<ENDHTML;
</div>
ENDHTML

if ($treemap_flag == 1) { # Display interactive treemap if flag is on
print {$fh}<<ENDHTML;
<h3>Interactive treemap of mapping-based taxonomic read classification</h3>
    <p><b>Navigation</b>: Left-click to go down, right-click to go up in taxonomic hierarchy, hover to see counts.</p>
    <p>Based on read-mapping hits to reference database, provides an approximate overview of taxonomic composition.</p>
    <div id="chart_div" style="width: 900px; height: 500px;"></div>
    <p>Drawn with Google Visualization API (<a href="https://developers.google.com/chart/terms">terms of service</a>)</p>
ENDHTML
}

print {$fh} <<ENDHTML;
<h4>Please cite...</h4>
<p>List of citations including dependencies</p>
</body>
</html>
ENDHTML

close($fh);
}

######################### MAIN ###########################
welcome();
parse_cmdline();
process_required_tools();
check_environment();

my $timer = new Timer;

bbmap_fast_filter_sam_run();
## Replaced: bbmap_hitstats_parse();

# Parse statistics from BBmap initial mapping
my ($bbmap_stats_aref, $skipflag) = bbmap_fast_filter_parse($libraryNAME.".bbmap.out", $SEmode);
# Dereference stats
($readnr,$readnr_pairs,$SSU_total_pairs,$SSU_ratio,$SSU_ratio_pc) = @$bbmap_stats_aref;
if (defined $skipflag && $skipflag == 1) {
    # If coverage too low, skip assembly
    $skip_spades = 1;
    $skip_emirge = 1;
}

# Parse sam file
readsam();

# Run SPAdes if not explicitly skipped
if ($skip_spades == 0) {
    spades_run();
    spades_parse();
}
# Run Emirge if not explicitly skipped
if ($skip_emirge == 0) {
    emirge_run();
    emirge_parse();
}
# If at least one of either SPAdes or Emirge is activated, parse results
if ($skip_spades + $skip_emirge < 2) {
    bbmap_spades_out();
    taxonomy_spades_unmapped();
    vsearch_best_match();
    vsearch_parse();
    vsearch_cluster();
    mafft_run();
}

$runtime = $timer->minutes;

print_report();
write_csv();

run_plotscript_SVG() if $html_flag;
#write_report_html()  if ($html_flag);
write_report_html_new(
    $version,
    $libraryNAME, $id, $SEmode,
    $readsf_full,$readsr_full,
    $readnr,$readnr_pairs,
    $ins_me,$ins_std,$ins_used,
    $SSU_ratio, $SSU_ratio_pc, $SSU_total_pairs,
    $skip_spades, $skip_emirge, $treemap_flag, # flags
    $taxon_report_lvl,
    \%ssu_sam_mapstats, # $mapstats_href
    $taxa_summary_href,
    $taxa_unassem_summary_href,
    \@ssuassem_results_sorted, # ssuassem_results_sorted_aref
    \@ssurecon_results_sorted, # ssurecon_results_sorted_aref
    ) if $html_flag;
#clean_up();

msg("Walltime used: $runtime with $cpus CPU cores");
msg("Thank you for using phyloFlash
You can find your results in $libraryNAME.*,
Your main result file is $libraryNAME.phyloFlash");
