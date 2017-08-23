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

=head1 INPUT ARGUMENTS

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

=head1 LOCAL SETTINGS

=over 15

=item -CPUs I<N>

Set the number of threads to use. Defaults to all available CPU cores.

=item -crlf

Use CRLF as line terminator in CVS output (to become RFC4180 compliant).

=item -decimalcomma

Use decimal comma instead of decimal point to fix locale problems.
Default: Off

=item -dbhome F<dir>

Directory containing phyloFlash reference databases.
Use F<phyloFlash_makedb.pl> to create an appropriate directory.

=back

=head1 INPUT AND ANALYSIS OPTIONS

=over 15

=item -interleaved

Interleaved readfile with R1 and R2 in a single file at read1

=item -readlength I<N>

Sets the expected readlength. Always use this if your read length
differs from 100 (the default). Must be within 50..500.

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

=item -taxlevel I<N>

Level in the taxonomy string to summarize read counts per taxon.
Numeric and 1-based (i.e. "1" corresponds to "Domain").
Default: 4 ("Order")

=item -maxinsert I<N>

Maximum insert size allowed for paired end read mapping. Must be within
0..1200. Default: 1200

=item -emirge

Turn on EMIRGE reconstruction of SSU sequences
Default: Off ("-noemirge")

=item -skip_spades

Turn off SPAdes assembly of SSU sequences

=item -sc

Turn on single cell MDA data mode for SPAdes assembly of SSU sequences

=item -poscov

Use Nhmmer to find positional coverage of reads across Barrnap's HMM model of
the 16S and 18S rRNA genes from a subsample of reads, as an estimate of
coverage evenness.
Default: Off ("-noposcov")

=back

=head1 OUTPUT OPTIONS

=over 15

=item -html

Generate output in HTML format.
Default: On.
(Turn off with "-nohtml")

=item -treemap

Draw interactive treemap of taxonomic classification in html-formatted
report. This uses Google Visualization API, which requires an internet
connection, requires that you agree to their terms of service (see
https://developers.google.com/chart/terms), and is not open source,
although it is free to use.
Default: Off ("-notreemap")

=item -zip

Compress output into a tar.gz archive file
Default: Off ("-nozip")

=item -keeptmp

Keep temporary/intermediate files
Default: No ("-nokeeptmp")

=item -log

Write status messages printed to STDERR also to a log file
Default: Off ("-nolog")

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


use FindBin;
use lib $FindBin::RealBin;

use PhyloFlash;

use Pod::Usage;
use Getopt::Long;
use File::Basename;
use IPC::Cmd qw(can_run);
use Cwd;

# change these to match your installation
my @dbhome_dirs = (".", $ENV{"HOME"}, $FindBin::RealBin);

# constants
my $version       = 'phyloFlash v3.0 beta 1';       # Current phyloFlash version
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
my $id          = 70;           # minimum %id for mapping
my $readlength  = 100;          # length of input reads
my $readlimit   = -1;           # max # of reads to use
my $amplimit    = 500000;       # number of SSU pairs at which to switch to emirge_amplicon
my $maxinsert   = 1200;         # max insert size for paired end read mapping
my $cpus        = get_cpus      # num cpus to use
my $clusterid   = 97;           # threshold for vsearch clustering

my $html_flag   = 1;            # generate HTML output? (default = 1, on)
my $treemap_flag = 0;           # generate interactive treemap (default = 0, off)
my $crlf        = 0;            # csv line terminator
my $decimalcomma= 0;            # Decimal separator (default = .)
my $skip_emirge = 1;            # Flag - skip Emirge step? (default = 1, yes)
my $skip_spades = 0;            # Flag - skip SPAdes assembly? (default = 0, no)
my $poscov_flag = 0;            # Flag - use Nhmmer to estimate positional coverage? (Default = 0, no)
my $sc          = 0;            # Flag - single cell data? (default = 0, no)
my $zip         = 0;            # Flag - Compress output into archive? (default = 0, no)
my $check_env   = 0;            # Check environment (runs check_environment subroutine only)
my $save_log    = 0;            # Save STDERR messages to log (Default = 0, no)
my $keeptmp     = 0;            # Do not delete temporary files (Default = 0, do delete temporary files)
my @tools_list;                 # Array to store list of tools required
                                # (0 will be turned into "\n" in parsecmdline)
# default database names for EMIRGE and Vsearch
my $emirge_db   = "SILVA_SSU.noLSU.masked.trimmed.NR96.fixed";
my $vsearch_db  = "SILVA_SSU.noLSU.masked.trimmed";

my $ins_used = "SE mode!"; # Report insert size used by EMIRGE

my %outfiles;           # Hash to keep track of output files

# variables for report generation
my %ssu_sam;            # Readin sam file from first mapping vs SSU database
my %ssu_sam_mapstats;   # Statistics of mapping vs SSU database, parsed from from SAM file

# Hashes to keep count of NTUs
my $taxa_full_href;             # Summary of read counts for full-length taxonomy strings (for treemap) - hash ref
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
        reformat => "reformat.sh",
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

# parse arguments passed on commandline and do some
# sanity checks
sub parse_cmdline {
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
               'emirge!' => \$emirge,
               'skip_spades' => \$skip_spades,
               'poscov!' => \$poscov_flag,
               'sc' => \$sc,
               'zip!' => \$zip,
               'log!' => \$save_log,
               'keeptmp!' => \$keeptmp,
               'everything' => \$everything,
               'almosteverything' => \$almosteverything,
               'check_env' => \$check_env,
               'help' => sub { pod2usage(1) },
               'man' => sub { pod2usage(-exitval=>0, -verbose=>2) },
           )
        or pod2usage(2);
    $skip_emirge = 0 if $emirge == 1; # ain't gonna not be less careful with no double negatives

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

    # populate hash to keep track of output files
    my $outfiles_href = initialize_infiles_hash($libraryNAME,$readsf);
    %outfiles = %$outfiles_href;
    # use hash to keep track of output file description, filename,
    # whether it should be deleted at cleanup, and if file was actually created

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
            $emirge = 1;
        } else {
            # "almost everything" does not turn on EMIRGE
            msg ("Running \"everything\" except EMIRGE- overrides other command line options");
        }
    }
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
    my @keys = sort {${$taxa_summary_href}{$b} <=> ${$taxa_summary_href}{$a}} keys %$taxa_summary_href;
    foreach my $key (@keys) {
        if (${$taxa_summary_href}{$key} >= 3) {
            my @out = ($key, ${$taxa_summary_href}{$key});
            print {$fh} join("\t", @out)."\n";
        }
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
            print {$fh} join("\t", @out)."\n";
        }

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
    
    # CSV file of phyloFlash run metadata
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
    my $fh;
    open_or_die(\$fh, ">", $outfiles{"report_csv"}{"filename"});
    $outfiles{"report_csv"}{"made"} = 1;
    while ($#report > 0) {
        print {$fh} shift(@report).",".shift(@report).$crlf;
    }
    close($fh);
    
    # CSV file of taxonomic units (NTUs) from mapping vs SILVA database
    open_or_die(\$fh, ">", $outfiles{"ntu_csv"}{"filename"});
    $outfiles{"ntu_csv"}{"made"} = 1;
    # Sort results descending
    my @keys = sort {${$taxa_summary_href}{$b} <=> ${$taxa_summary_href}{$a}} keys %$taxa_summary_href;
    foreach my $key (@keys) {
        my @out = ($key, ${$taxa_summary_href}{$key});
        print {$fh} join(",",csv_escape(@out)).$crlf;
    }
    close($fh);
    
    # If full-length seqeunces were assembled or reconstructed
    if ($skip_spades + $skip_emirge < 2) {
        
        # CSV file of assembled/reconstructed sequnces 
        open_or_die(\$fh, ">", $outfiles{"full_len_class"}{"filename"});
        $outfiles{"full_len_class"}{"made"}++;
        print $fh "OTU,read_cov,coverage,dbHit,taxonomy,%id,alnlen,evalue\n";
        if ($skip_spades == 0) {
            foreach my $arr (@ssuassem_results_sorted) {
                my $otu = ${$arr}[0];
                my @out = @$arr;
                splice @out, 1, 0, $ssuassem_cov{$otu}; # Insert read coverage into output array
                print {$fh} join(",",csv_escape(@out)).$crlf;
            }
        
            # CSV file of taxonomic affiliations for unassembled reads
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
        # Add sequences from EMIRGE, if available
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
    # output: $libraryNAME.$readsf.SSU.sam
    #         $libraryNAME.$readsf.SSU.{1,2}.fq
    #         $libraryNAME.bbmap.out
    #         $libraryNAME.inserthistogram
    # tmp:    $libraryNAME.basecompositionhistogram
    msg("filtering reads with SSU db using minimum identity of $id%");
    if ($readlimit != -1) {
        msg("Only using the first $readlimit reads");
    }
    # Minimum mapping ID of 63%
    my $minID= $id / 100;
    if ($minID < 0.63) {
        $minID = 0.63
    }
    # Set up input arguments for BBmap
    my @bbmap_args = ("fast=t",
                      "minidentity=$minID",
                      "-Xmx20g reads=$readlimit",
                      "threads=$cpus",
                      "po=f",
                      "outputunmapped=f",
                      "path=$DBHOME",
                      "out=".$outfiles{"sam_map"}{"filename"}, # Also skip SAM header?
                      "outm=".$outfiles{"reads_mapped_f"}{"filename"},
                      "build=1",
                      "in=$readsf_full",
                      "bhist=".$outfiles{"basecompositionhist"}{"filename"},
                      "ihist=".$outfiles{"inserthistogram"}{"filename"},
                      "idhist=".$outfiles{"idhistogram"}{"filename"},
                      "scafstats=".$outfiles{"hitstats"}{"filename"},
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
    foreach my $madekey (qw(sam_map reads_mapped_f basecompositionhist inserthistogram idhistogram hitstats bbmap_log)) {
        $outfiles{$madekey}{"made"}++;
    }
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
    # Read reported statistics from BBmap log file
    my $fh;
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
    $SSU_ratio_pc = sprintf ("%.3f", $SSU_ratio * 100);

    # Ratios of mapped vs unmapped to report
    my @mapratio_csv;
    if ($SEmode == 1) { # TO DO: Add numerical values to text labels
        my $unmapped = 1-$SSU_ratio;
        push @mapratio_csv, "Unmapped,".$unmapped;
        push @mapratio_csv, "Mapped,".$SSU_ratio;
    } elsif ($SEmode == 0) {
        # For PE reads, do not include Unmapped in piechart because it will be
        # impossible to read the other slices anyway. Instead report mapping
        # ratio in title.
        my $mapped_half = $ssu_f_reads + $ssu_r_reads - 2 * ($ssu_pairs + $ssu_bad_pairs);
        my $mapped_pairs = $ssu_pairs + $ssu_bad_pairs;
        my $unmapped_pairs = $readnr_pairs - $mapped_half - $mapped_pairs;
        push @mapratio_csv, "Mapped pair,".$ssu_pairs;
        push @mapratio_csv, "Mapped single,".$mapped_half;
    }

    # CSV file to draw piechart
    my $fh_csv;
    open_or_die (\$fh_csv, ">", $outfiles{"mapratio_csv"}{"filename"});
    $outfiles{"mapratio_csv"}{"made"}++;
    print $fh_csv join ("\n", @mapratio_csv);
    close ($fh_csv);

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
    my $infile = $outfiles{"sam_map"}{"filename"};  # Input is SAM from first mapping step
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
                $taxa_full_href->{$taxonlongstring}++; # Count full-length taxon for treemap
                # Save full taxonomy string
                push @taxa_full, $taxonlongstring;
            } else {
                msg ("Warning: malformed database entry $ref");
            }
        }
    }
    close($fh);

    # Summarize taxonomy
    $taxa_summary_href = summarize_taxonomy(\@taxa_full, $taxon_report_lvl); # Summarize

    # Count 1-tons, 2-tons, and 3+-tons and calculate Chao1 statistic
    foreach my $taxon (keys %$taxa_summary_href) {
        if ($taxa_summary_href->{$taxon} == 1) {        # 1-tons
            $xtons[0]++;
        } elsif ($taxa_summary_href->{$taxon} == 2) {   # 2-tons
            $xtons[1]++;
        } elsif ($taxa_summary_href->{$taxon} >= 3) {   # 3+-tons
            $xtons[2]++;
        }
    }
    if ($xtons[1] > 0) {
        $chao1 =
          $xtons[2] +                               # Is there an error here? Should be sum of all spp. observed
          ($xtons[0] * $xtons[0]) / 2 / $xtons[1];
    } else {
        $chao1 = 'n.d.';
    }
    
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

    run_prog("spades",
             "-o $libraryNAME.spades -t $cpus_spades -m 20 -k $kmer "
             . $args,
             $outfiles{"spades_log"}{"filename"},"&1"
         );
    $outfiles{"spades_log"}{"made"}++;

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
        return;
    }

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
                             "-fo",$outfiles{"spades_fasta"}{"filename"},
                             "-s","-name",
                             );
    run_prog("fastaFromBed",
             join (" ",@fastaFromBed_args),
             $outfiles{"fastaFromBed_out"}{"filename"},
             "&1");
    $outfiles{"spades_fasta"}{"made"}++;
    $outfiles{"fastaFromBed_out"}{"made"}++;
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
            "  in2=".$outfiles{"reads_mapped_r"}{"filename"}
            . " pairlen=$maxinsert ";
        }
    }
    run_prog("bbmap",
         "  fast=t"
         . " minidentity=0.98"
         . " -Xmx20g"
         . " threads=$cpus"
         . " po=f"
         . " outputunmapped=t" # This is important
         . " ref=".$outfiles{"spades_fasta"}{"filename"}
         . " nodisk"
         . " in=".$outfiles{"reads_mapped_f"}{"filename"}
         . " out=".$outfiles{"sam_remap"}{"filename"} # Also skip SAM header?
         . $args,
         undef,
         $outfiles{"bbmap_remap_log"}{"filename"});
    $outfiles{"sam_remap"}{"made"}++;
    $outfiles{"bbmap_remap_log"}{"made"}++;
    msg("done...");
}

sub taxonomy_spades_unmapped {
    # Filter output of mapping to SPAdes assembled SSU sequences
    # Report taxonomy of "leftover" sequences
    my $in = $outfiles{"sam_remap"}{"filename"};    # mapping of extracted reads vs assembled SSU seq
    my $sam_href = \%ssu_sam;                       # data from first mapping vs. SILVA, read into memory
    my $stats_href = \%ssu_sam_mapstats;            # mapping statistics
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
    close($fh);

    # Summarize counts per NTU of unassembled reads
    $taxa_unassem_summary_href = summarize_taxonomy(\@taxa_full, $taxon_report_lvl);

    # Calculate ratio of reads assembled and store in hash
    my $assem_tot_map = ${$stats_href}{"assem_fwd_map"};
    if (defined ${$stats_href}{"assem_rev_map"}) {
        $assem_tot_map += ${$stats_href}{"assem_rev_map"};
    }
    my $ssu_tot_map = ${$stats_href}{"ssu_fwd_map"};
    if (defined ${$stats_href}{"ssu_rev_map"}) {
        $ssu_tot_map += ${$stats_href}{"ssu_rev_map"};
    }
    my $ssu_unassem = $ssu_tot_map - $assem_tot_map;
    ${$stats_href}{"assem_tot_map"} = $assem_tot_map;
    ${$stats_href}{"ssu_tot_map"} = $ssu_tot_map;
    ${$stats_href}{"ssu_unassem"} = $ssu_unassem;
    ${$stats_href}{"assem_ratio"} = $assem_tot_map / $ssu_tot_map;
    ${$stats_href}{"assem_ratio_pc"} = sprintf ("%.3f", ${$stats_href}{"assem_ratio"} * 100) ;

    # Write CSV of reads assembled for piechart
    my @assemratio_csv;
    push @assemratio_csv, "Unassembled,".$ssu_unassem;
    push @assemratio_csv, "Assembled,".$assem_tot_map;
    my $fh_csv;
    open_or_die (\$fh_csv, ">", $outfiles{"assemratio_csv"}{"filename"});
    $outfiles{"assemratio_csv"}{"made"}++;
    print $fh_csv join ("\n", @assemratio_csv);
    close($fh_csv);

    msg ("done ... ");
}

sub emirge_run {
    # running Emirge on bbmap output
    msg("creating phylotypes with Emirge");

    my $cmd = "emirge";
    my @emirge_args = ("-1",$outfiles{"reads_mapped_f"}{"filename"});
    my $args = "-1 ".$outfiles{"reads_mapped_f"}{"filename"}." ";

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

        @emirge_args = ("-1",$outfiles{"reads_mapped_f"}{"filename"},
                        "-2",$outfiles{"reads_mapped_r"}{"filename"},
                        "-i $ins_used",
                        "-s $ins_std",
                        );
        $args = "  -1 ".$outfiles{"reads_mapped_f"}{"filename"}
                . " -2 ".$outfiles{"reads_mapped_r"}{"filename"}
                . " -i $ins_used -s $ins_std ";
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
                 "\'{print (NR%4 == 1) ? \"\@\" ++i  : (NR%4 == 3) ? \"+\" :\$0}\'"
                 .$outfiles{"reads_mapped_cat"}{"filename"},
                 $outfiles{"reads_mapped_cat_rename"}{"filename"});
        $outfiles{"reads_mapped_cat_rename"}{"made"}++;

        @emirge_args = ("-1",$outfiles{"reads_mapped_cat_rename"}{"filename"});
        $args = " -1 ".$outfiles{"reads_mapped_cat_rename"}{"filename"}." ";
    }

    unshift @emirge_args, "$libraryNAME.emirge"; # Add library name to arguments
    push @emirge_args, ("-f ${DBHOME}/${emirge_db}.fasta",
                        "-b ${DBHOME}/${emirge_db}.bt",
                        "-l $readlength",
                        "-a $cpus",
                        "--phred33",
                        );
    run_prog($cmd,
             join (" ", @emirge_args),
             $outfiles{"emirge_log"}{"filename"},
             "&1",
             );
    $outfiles{"emirge_log"}{"made"}++;
    msg("done...");
}

sub emirge_parse {
    # getting emirge output and reformatting it...
    msg("getting Emirge phylotypes and their abundances...");
    run_prog("emirge_rename_fasta",
             "./$libraryNAME.emirge/iter.40/",
             $outfiles{"emirge_raw_fasta"}{"filename"});
    $outfiles{"emirge_result_fasta"}{"made"}++;

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

    # join emirge and spades hits into one file
    # (vsearch takes a while to load, one run saves us time)
    if ($skip_spades == 1 && $skip_emirge == 0) {
        run_prog("cat",
                 "   ".$outfiles{"emirge_fasta"}{"filename"},
                 $outfiles{"all_final_fasta"}{"filename"});
        $outfiles{"all_final_fasta"}{"made"}++;
    }
    elsif ($skip_emirge == 1 && $skip_spades == 0){
        run_prog("cat",
                 "   ".$outfiles{"spades_fasta"}{"filename"},
                 $outfiles{"all_final_fasta"}{"filename"});
        $outfiles{"all_final_fasta"}{"made"}++;
    }
    elsif ($skip_emirge == 0 && $skip_spades == 0){
        run_prog("cat",
                 "   ".$outfiles{"emirge_fasta"}{"filename"}
                 ." ".$outfiles{"spades_fasta"}{"filename"},
                 $outfiles{"all_final_fasta"}{"filename"});
        $outfiles{"all_final_fasta"}{"made"}++;
    }

    if (-s $outfiles{"all_final_fasta"}{"filename"}) {
       run_prog("vsearch",
             "-usearch_global ".$outfiles{"all_final_fasta"}{"filename"}
             . " -db ${DBHOME}/${vsearch_db}.fasta"
             . " -id 0.7"
             . " -userout ".$outfiles{"vsearch_csv"}{"filename"}
             . " -userfields query+target+id+alnlen+evalue+id3+qs+pairs+gaps+mism+ids"
             . " -threads $cpus --strand plus --notrunclabels"
             . " -notmatched ".$outfiles{"notmatched_fasta"}{"filename"}
             . " -dbmatched ".$outfiles{"dbhits_all_fasta"}{"filename"},
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

    if ($skip_spades == 0 && $skip_emirge == 1) {
        run_prog("cat",
                 $outfiles{"dhbits_nr97_fasta"}{"filename"}
                 . " ".$outfiles{"spades_fasta"}{"filename"},
                 $outfiles{"ssu_coll_fasta"}{"filename"});
        $outfiles{"ssu_coll_fasta"}{"made"}++;

    } elsif ($skip_spades == 1 && $skip_emirge == 0) {
        run_prog("cat",
                 $outfiles{"dhbits_nr97_fasta"}{"filename"}
                 ." ".$outfiles{"emirge_fasta"}{"filename"},
                 $outfiles{"ssu_coll_fasta"}{"filename"});
         $outfiles{"ssu_coll_fasta"}{"made"}++;

    } elsif ($skip_spades == 0 && $skip_emirge == 0) {
        run_prog("cat",
                 $outfiles{"dhbits_nr97_fasta"}{"filename"}
                 ." ".$outfiles{"spades_fasta"}{"filename"}
                 ." ".$outfiles{"emirge_fasta"}{"filename"},
                 $outfiles{"ssu_coll_fasta"}{"filename"});
         $outfiles{"ssu_coll_fasta"}{"made"}++;
    }

    run_prog("mafft",
             "--treeout ".$outfiles{"ssu_coll_fasta"}{"filename"},
             $outfiles{"ssu_coll_aln_fasta"}{"filename"},
             $outfiles{"ssu_coll_aln_mafftout"}{"filename"});
    $outfiles{"ssu_coll_aln_fasta"}{"made"}++;
    $outfiles{"ssu_coll_tree"}{"made"}++;
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
    my $subsample = $outfiles{"readsf_subsample"}{"filename"};
    my $samplelimit = 10000; # Take sample of max this number of reads

    # Subsample reads with reformat.sh
    my @reformat_args = ("in=$sam","out=$subsample","srt=$samplelimit");
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
    # Write result to prokaryote histogram
    my $fh_prok;
    open_or_die (\$fh_prok, ">", $outfiles{"nhmmer_prok_histogram"}{"filename"});
    foreach my $key (sort {$a <=> $b} keys %prok_pos) {
        print $fh_prok $key."\t". $prok_pos{$key}."\n";
    }
    close($fh_prok);
    $outfiles{"nhmmer_prok_histogram"}{"made"}++;
    # Write results for eukaryote histogram
    my $fh_euk;
    open_or_die (\$fh_euk, ">", $outfiles{"nhmmer_euk_histogram"}{"filename"});
    foreach my $key (sort {$a <=> $b} keys %euk_pos) {
        print $fh_euk $key."\t". $euk_pos{$key}."\n";
    }
    close($fh_euk);
    $outfiles{"nhmmer_euk_histogram"}{"made"}++;

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
            system ("rm $libraryNAME.spades -r");
        }
        if ($skip_emirge == 0) {
            system ("rm $libraryNAME.emirge -r");
        }
        
        # Remove tmp files marked by "tmp" filename prefix                      # To be superseded
        system ("rm tmp.$libraryNAME.* -r");
    
        # Remove files that are earmarked for destruction
        # (Change earmarking in PhyloFlash.pm)
        foreach my $key (keys %outfiles) {
            if (defined $outfiles{$key}{"made"} && $outfiles{$key}{"discard"} == 1) {
                system ("rm -r ".$outfiles{$key}{"filename"});
            }
        }
    }
    
    # Compress output into tar.gz archive if requested
    if ($zip == 1) {
        my @filelist;
        my $tarfile = $outfiles{"phyloFlash_archive"}{"filename"};
        msg ("Compressing results into tar.gz archive $tarfile");
        foreach my $key (keys %outfiles) {
            if (defined $outfiles{$key}{"made"} && $outfiles{$key}{"discard"} ==0) {
                push @filelist, $outfiles{$key}{"filename"};
            }
        }
        my $to_tar = join " ", @filelist;
        system ("tar -czf $tarfile $to_tar");
        system ("rm -r $to_tar");
    }

    msg("Done...");
}

sub run_plotscript_SVG {
    msg ("generating graphics for report in SVG format");
    # Plot mapping ID histogram
    my @idhist_args = ( "--hist",
                        $outfiles{"idhistogram"}{"filename"},
                        " --title=\"Mapping identity (%)\" ",
                      );
    if (defined $decimalcomma) {
        push @idhist_args, "--decimalcomma";
    }
    run_prog("plotscript_SVG",
         join (" ", @idhist_args),
         $outfiles{"plotscript_out"}{"filename"},
         "&1");
    $outfiles{"idhistogram_svg"}{"made"}++;

    # Piechart of proportion mapped
    my @map_args = ("-pie",
                    $outfiles{"mapratio_csv"}{"filename"},
                    "-title=\"$SSU_ratio_pc % reads mapped\"");
    run_prog("plotscript_SVG",
             join(" ", @map_args),
             $outfiles{"plotscript_out"}{"filename"},
             "&1");
    $outfiles{"mapratio_svg"}{"made"}++;

    # Piechart of proportion assembled, if SPAdes run
    unless ($skip_spades == 1) {
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

    # Plot insert size histogram unless running in SE mode
    if ($SEmode != 1) { # If not running in SE mode ...
        my @inshist_args = ("--hist ",
                            $outfiles{"inserthistogram"}{"filename"},
                            " --title=\"Insert size (bp)\" ",
                            );
        if (defined $decimalcomma) {
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
        # Eukaryotic gene model
        my @args_euk = ("-hist",
                        $outfiles{"nhmmer_euk_histogram"}{"filename"},
                        "-height 150",
                        "-width 480",
                        "-title=\"Coverage on 18S model\"",
                        "-color=\"blue\"",
                        );
        run_prog("plotscript_SVG",
                 join (" ", @args_euk),
                 $outfiles{"plotscript_out"}{"filename"},
                 "&1");
        $outfiles{"nhmmer_euk_histogram_svg"}{"made"}++;
        # Prokaryotic gene model
        my @args_prok = ("-hist",
                         $outfiles{"nhmmer_prok_histogram"}{"filename"},
                         "-height 150",
                         "-width 480",
                         "-title=\"Coverage on 16S model\"",
                         "-color=\"blue\"",
                         );
        run_prog("plotscript_SVG",
                 join (" ", @args_prok),
                 $outfiles{"plotscript_out"}{"filename"},
                 "&1");
        $outfiles{"nhmmer_prok_histogram_svg"}{"made"}++;
    }

    # Plot tree if spades/emirge unless both skipped
    unless ($skip_spades + $skip_emirge == 2) {
        my @treeplot_args = ("-tree",
                             $outfiles{"ssu_coll_tree"}{"filename"},
                             );
         # Add "bubbles" corresponding to read coverage
        if (defined $outfiles{"full_len_class"}{"made"}) {
            push @treeplot_args, ("-assemcov",
                                  $outfiles{"full_len_class"}{"filename"},
                                  "-unassemcount",
                                  $ssu_sam_mapstats{"ssu_unassem"},
                                  );
        }
        run_prog("plotscript_SVG",
                 join(" ",@treeplot_args),
                 $outfiles{"plotscript_out"}{"filename"},
                 "&1");
        $outfiles{"ssu_coll_tree_svg"}{"made"}++;
    }

    # Generate barplot of taxonomy at level XX
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
        $ssuassem_results_sorted_aref,
        $ssurecon_results_sorted_aref,
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
        "CHAO1" => $chao1,
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
    if ($skip_spades + $skip_emirge == 2)  {
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
    $flags{"IDHISTOGRAM"} = slurpfile($outfiles{"idhistogram_svg"}{"filename"});
    $flags{"MAPRATIOPIE"} = slurpfile($outfiles{"mapratio_svg"}{"filename"});
    $flags{"TAXONSUMMARYBAR"} = slurpfile($outfiles{"ntu_csv_svg"}{"filename"});

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
        $flags{"INSERTHISTOGRAM"} = slurpfile($outfiles{"inserthistogram_svg"}{"filename"});
    }

    # Slurp in positional coverage histograms, if defined
    if ($poscov_flag == 1) {
        $flags{"POSCOVHIST_PROK"} = slurpfile($outfiles{"nhmmer_prok_histogram_svg"}{"filename"});
        $flags{"POSCOVHIST_EUK"} = slurpfile($outfiles{"nhmmer_euk_histogram_svg"}{"filename"});
    }

    # Params defined only for assembled (SPAdes) reads
    if ($skip_spades == 0) {
        $flags{"ASSEM_RATIO"} = $mapstats_href->{"assem_ratio_pc"};
        $flags{"ASSEMBLYRATIOPIE"} = slurpfile($outfiles{"assemratio_svg"}{"filename"});
        $flags{"SEQUENCES_TREE"} = slurpfile($outfiles{"ssu_coll_tree_svg"}{"filename"});

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
            push @table_assem_seq, "    <td><a href=\"http://www.ncbi.nlm.nih.gov/nuccore/".$get_genbank[0]."\">".$split_entry[2]."</a></td>\n";    # Link to Genbank entry using accession no.
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

sub write_logfile {
    # This should always be run last
    msg("Saving log to file");
    my $fh;
    open_or_die(\$fh, ">>", $outfiles{"phyloFlash_log"}{"filename"});
    print $fh join "\n", @PhyloFlash::msg_log;
    close($fh);
}

######################### MAIN ###########################
welcome();
parse_cmdline();
process_required_tools();
check_environment();

my $timer = new Timer;

bbmap_fast_filter_sam_run();

# Parse statistics from BBmap initial mapping
my ($bbmap_stats_aref, $skipflag) = bbmap_fast_filter_parse($outfiles{"bbmap_log"}{"filename"}, $SEmode);

# Dereference stats
($readnr,$readnr_pairs,$SSU_total_pairs,$SSU_ratio,$SSU_ratio_pc) = @$bbmap_stats_aref;
if (defined $skipflag && $skipflag == 1) {
    # If coverage too low, skip assembly
    $skip_spades = 1;
    $skip_emirge = 1;
}

# Parse sam file
readsam();

# Find positional coverage along prok and euk SSU rRNA models using nhmmer
# from a subsample of reads
nhmmer_model_pos() if ($poscov_flag == 1 );

# Run SPAdes if not explicitly skipped
if ($skip_spades == 0) {
    spades_run();
    spades_parse();
    bbmap_spades_out();
    taxonomy_spades_unmapped();
}
# Run Emirge if not explicitly skipped
if ($skip_emirge == 0) {
    emirge_run();
    emirge_parse();
}
# If at least one of either SPAdes or Emirge is activated, parse results
if ($skip_spades + $skip_emirge < 2) {
    vsearch_best_match();
    vsearch_parse();
    vsearch_cluster();
    mafft_run();
}

$runtime = $timer->minutes; # Log run time

# Capture output parameters for reports
my @report_inputs = (
    $version, $progcmd, $cwd, $DBHOME, $libraryNAME, $id, $SEmode,
    $readsf_full, $readsr_full, $readnr, $readnr_pairs,
    $ins_me,$ins_std,$ins_used,
    $SSU_ratio, $SSU_ratio_pc, $SSU_total_pairs,
    $skip_spades, $skip_emirge, $treemap_flag, # flags
    $taxon_report_lvl,
    \%ssu_sam_mapstats,\%outfiles,\@xtons,$chao1,
    $taxa_summary_href, $taxa_unassem_summary_href,
    \@ssuassem_results_sorted, \@ssurecon_results_sorted,
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

write_logfile() if ($save_log == 1);
