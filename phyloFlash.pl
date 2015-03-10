#!/usr/bin/perl -w
=head1 NAME

phyloFlash - A script to rapidly estimate the phylogenetic composition of
             an illumina (meta)genomic dataset.


=head1 SYNOPSIS

B<phyloFlash.pl> -dbhome F<dir> -lib B<name> -read1 F<<file>> -read2 F<<file>>

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
The file may be compressed (.gz).

=item -read2 F<file>

File containing reverse reads. If this option is omitted, B<phyloFlash>
will run in B<experimental> single-end mode.

=back

=head1 OPTIONS

=over 15

=item -dbhome F<dir>

Directory containing phyloFlash reference databases.
Use F<phyloFlash_makedb.pl> to create an appropriate directory.

=item -readlength I<N>

Sets the expected readlength. Always use this if your read length
differs from 100 (the default). Must be within 50..500.

=item -CPUs I<N>

Set the number of threads to use. Defaults to all available CPU cores.

=item -readlimit I<N>

Limits processing to the first I<N> reads in each input file. Use this
for transcriptomes with a lot of rRNA reads (use values <1000000).
Default: unlimited

=item -id I<N>

Minimum allowed identity for read mapping process in %. Must be within
63..95. Set this to a lower value for very divergent taxa
Default: 70

=item -clusterid I<N>

Identity threshold for clustering with vsearch in %.
Must be within 50..100. Default: 97

=item -maxinsert I<N>

Maximum insert size allowed for paired end read mapping. Must be within
0..1200. Default: 1200

=item -html

Also generate output in HTML format.

=item -crlf

Use CRLF as line terminator in CVS output (to become RFC4180 compliant).

=item -help

Print brief help.

=item -man

Show manual.

=back

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2014- by Harald Gruber-Vodicka <hgruber@mpi-bremen.de>
                    and Elmar A. Pruesse <elmar.pruesse@ucdenver.edu>
                    with help from Brandon Seah mail:


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

my $DBHOME = "$FindBin::RealBin/119";#HGV: edited to point to locally created db dir release version 119

# binaries needed by phyloFlash
require_tools((
    bbmap => "bbmap.sh",
    spades => "spades.py",
    barrnap => "$FindBin::RealBin/barrnap-HGV/bin/barrnap_HGV",
    emirge => "emirge.py",
    emirge_amp => "emirge_amplicon.py",
    emirge_rename_fasta => "emirge_rename_fasta.py",
    vsearch => "vsearch",
    mafft => "mafft",
    fastaFromBed => "fastaFromBed",
    sed => "sed",
    grep => "grep",
    awk => "awk",
    cat => "cat",
    plotscript => "$FindBin::RealBin/phyloFlash_plotscript.R"
    ));

# constants
my $version       = 'phyloFlash v2.0';
my $progname      = $FindBin::Script;
my $cwd = getcwd;

# configuration variables / defaults
my $readsf      = undef;        # forward read file
my $readsr      = undef;        # reverse read file
my $SEmode      = 0;            # single ended mode
my $libraryNAME = undef;        # output basename
my $id          = 70;           # minimum %id for mapping
my $readlength  = 100;          # length of input reads
my $readlimit   = -1;           # max # of reads to use
my $maxinsert   = 1200;         # max insert size for paired end read mapping
my $cpus        = get_cpus      # num cpus to use
my $clusterid   = 97;           # threshold for vsearch clustering

my $html_flag   = 0;            # generate HTML output? (default = 0, off)
my $crlf        = 0;            # csv line terminator
                                # (0 will be turned into "\n" in parsecmdline)

my $emirge_db   = "SILVA_SSU.noLSU.masked.trimmed.NR96.fixed";
my $vsearch_db  = "SILVA_SSU.noLSU.masked.trimmed";

my $ins_used = "SE mode!";


# variables for report generation
my @taxa_from_hitmaps_sorted;
my @ssuassem_results_sorted; #sorted list of SSU sequences for reporting
my @ssurecon_results_sorted;
my $readnr_pairs;
my $SSU_total_pairs;
my $SSU_ratio;
my $SSU_ratio_pc;
my $chao1 = 0;
my @xtons = (0,0);  # singleton and doubleton count
my $readnr = 0;
my $ins_me = 0;
my $ins_std =0;
my $runtime;


# display welcome message incl. version
sub welcome {
  print STDERR "\nThis is $version\n\n";
}


# parse arguments passed on commandline and do some
# sanity checks
sub parse_cmdline {
    GetOptions('read1=s' => \$readsf,
               'read2=s' => \$readsr,
               'lib=s' => \$libraryNAME,
               'dbhome=s' => \$DBHOME,
               'readlength=i' => \$readlength,
               'readlimit=i' => \$readlimit,
               'maxinsert=i' => \$maxinsert,
               'id=i' => \$id,
               'clusterid=i' => \$clusterid,
               'CPUs=i' => \$cpus,
               'html' => \$html_flag,
               'crlf' => \$crlf,
               'help' => sub { pod2usage(1) },
               'man' => sub { pod2usage(-exitval=>0, -verbose=>2) },
           )
        or pod2usage(2);

    # check correct $libraryNAME
    pod2usage("Please specify output file basename with -lib")
        if !defined($libraryNAME);
    pod2usage("Argument to -lib may not be empty")
        if length($libraryNAME) == 0;
    pod2usage("Argument to -lib may contain only alphanumeric characters,"
              ."'_' and '-'.")
        if not $libraryNAME =~ m/[a-zA-Z0-9_-]/ ;

    msg("working on library $libraryNAME");

    # check correct read file(s)
    pod2usage("Please specify input forward read file with -read1")
        if !defined($readsf);
    pod2usage("Unable to open forward read file '$readsf'. ".
              "Make sure the file exists and is readable by you.")
        if ! -e $readsf;

    if (defined($readsr)) {
        pod2usage("Unable to open reverse read file '$readsr'.".
                  "Make sure the file exists and is readable by you.")
            if ! -e $readsr;
        pod2usage("Forward and reverse input read file need to be different")
            if ($readsr eq $readsf);
    } else {
        $SEmode = 1; # no reverse reads, we operate in single ended mode
        $readsr = "<NONE>";
    }

    msg("Forward reads $readsf");
    if ($SEmode == 0)  {
        msg("Reverse reads $readsr");
    } else {
        msg("Running in single ended mode");
    }

    # check dbhome
    ## FIXME

    msg("Using dbhome '$DBHOME'");

    # check lengths
    pod2usage("Readlength must be within 50...500")
        if ($readlength < 50 or $readlength > 500);
    pod2usage("Maxinsert must be within 0..1200")
        if ($maxinsert < 0 or $maxinsert > 1200);
    pod2usage("Readmapping identity (-id) must be within 63..95")
        if ($id < 63 or $id > 95);
    pod2usage("Clustering identidy (-clusterid) must be within 50..100")
        if ($clusterid < 50 or $clusterid > 100);

    $crlf = $crlf ? "\r\n":"\n";

    # check CPUS
    if ($cpus eq "all" or $cpus < 0) {
        $cpus = get_cpus();
    }

    # check surplus arguments
    if ($ARGV[0]) {
        print "Commandline contains extra words:";
        foreach (@ARGV) {
            print "$_\n";
        }
        die("aborting");
    }
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
Forward read file\t$readsf
Reverse read file\t$readsr
Current working directory\t$cwd
---
Minimum mapping identity:\t$id%
~;

    if ($SEmode == 1) {
        print {$fh} qq~
---
Input PE-reads:\t$readnr_pairs
Mapped SSU read pairs:\t$SSU_total_pairs
~;
    } else {
        print {$fh} qq~
---
Input SE-reads:\t$readnr_pairs
Mapped SSU reads:\t$SSU_total_pairs
~;
    }

    print {$fh} qq~Mapping ratio:\t$SSU_ratio_pc%
Detected median insert size:\t$ins_me
Used insert size:\t$ins_used
Insert size standard deviation:\t$ins_std
---
Runtime:\t$runtime
CPUs used:\t$cpus
---
Read mapping based higher taxa (NTUs) detection

NTUs observed once:\t$xtons[0]
NTUs observed twice:\t$xtons[1]
NTUs observed three or more times:\t$#taxa_from_hitmaps_sorted
NTU Chao1 richness estimate:\t$chao1

List of NTUs in order of abundance:
NTU\treads
~;

    # sort keys numerically descending in hash of
    # mapping-based detected higher taxa
    foreach my $taxonshortstring ( @taxa_from_hitmaps_sorted ) {
        # Print the name of this higher taxon, and the
        # corresponding no. of unambig hits mapped
        print {$fh} join("\t",@$taxonshortstring)."\n";
    }

    ## Print the table of SSU assembly-based taxa to report file
    print {$fh} "---\n";
    print {$fh} "SSU assembly based taxa:\n";
    print {$fh} "OTU\tcoverage\tdbHit\ttaxonomy\t%id\talnlen\tevalue\n";

    foreach my $arr (@ssuassem_results_sorted) {
        print {$fh} join("\t",  @$arr)."\n";
    }

    ## Print the table of SSU reconstruction-based taxa to report file
    print {$fh} "---\n";
    print {$fh} "SSU reconstruction based taxa:\n";
    print {$fh} "OTU\tratio\tdbHit\ttaxonomy\t%id\talnlen\tevalue\n";
    foreach my $arr (@ssurecon_results_sorted) {
        print {$fh} join("\t", @$arr)."\n";
    }

    close($fh);
}

sub write_csv {
    msg("exporting results to csv");
    my $fh;

    my @report = csv_escape((
        "version",$version,
        "library, name",$libraryNAME,
        "forward\" read file",$readsf,
        "reverse \"read\" file",$readsr,
        "cwd",$cwd,
        "minimum mapping identity",$id,
        "single ended mode",$SEmode,
        "input reads",$readnr_pairs,
        "mapped SSU reads",$SSU_total_pairs,
        "mapping ration",$SSU_ratio_pc,
        "detected median insert size",$ins_me,
        "used insert size",$ins_used,
        "insert size stddev",$ins_std,
        "runtime",$runtime,
        "CPUs",$cpus,
        "NTUs observed once",$xtons[0],
        "NTUs observed twice",$xtons[1],
        "NTUs observed three or more times",$#taxa_from_hitmaps_sorted,
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

    open_or_die(\$fh, ">", "$libraryNAME.phyloFlash.extractedSSUclassifications.csv");
    print $fh "OTU,coverage,dbHit,taxonomy,%id,alnlen,evalue\n";
    foreach my $arr (@ssuassem_results_sorted) {
        print {$fh} join(",",csv_escape(@$arr)).$crlf;
    }
    foreach my $arr (@ssurecon_results_sorted) {
        print {$fh} join(",",csv_escape(@$arr)).$crlf;
    }
    close($fh);
}


sub bbmap_fast_filter_run {
    # input: $readsf, $readsr
    # output: lib.$readsf.SSU.{1,2}.fq
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
        $args =
        "  outm2=$libraryNAME.$readsf.SSU.2.fq "
        . "pairlen=$maxinsert in2=$readsr";
    }

    run_prog("bbmap",
             "  fast=t "
             . "minidentity=$minID "
             . "-Xmx20g reads=$readlimit "
             . "threads=$cpus "
             . "po=f "
             . "outputunmapped=f "
             . "path=$DBHOME "
             . "outm=$libraryNAME.$readsf.SSU.1.fq "
             . "build=1 "
             . "in=$readsf "
             . "bhist=tmp.$libraryNAME.basecompositionhistogram "
             . "ihist=$libraryNAME.inserthistogram "
             . "scafstats=$libraryNAME.hitstats "
             . $args,
             undef,
             "$libraryNAME.bbmap.out");

    msg("done...");
}


sub bbmap_fast_filter_parse() {
    # parsing bbmap.out for used read numbers, mapped reads,
    # insert size median and standard deviation

    # input: lib.bbmap.out

    my $ssu_pairs     = 0;
    my $ssu_bad_pairs = 0;
    my $ssu_f_reads   = 0;
    my $ssu_r_reads   = 0;
    my $forward_count = 0;

    my $fh;
    open_or_die(\$fh, "<", $libraryNAME.'.bbmap.out');
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

    $SSU_total_pairs =
        $ssu_f_reads + $ssu_r_reads
        - $ssu_pairs - $ssu_bad_pairs;
    if ($SEmode == 0) {
        msg("mapped pairs output: $SSU_total_pairs")
    };

    $SSU_ratio = $SSU_total_pairs / $readnr_pairs;
    $SSU_ratio_pc = $SSU_ratio * 100;

    msg("mapping rate: $SSU_ratio_pc%");
}

sub bbmap_hitstats_parse {
    # create taxa list from hitmaps with at least 3 unambiguously mapping reads
    # input: lib.hitstats
    msg("creating taxon list from read mappings");

    my $fh;
    open_or_die(\$fh, "<", "$libraryNAME.hitstats");

    # hitstats format/example
    #name %unambiguousReads unambiguousMB %ambiguousReads ambiguousMB unambiguousReads ambiguousReads
    #X03680.934.2693 Eukaryota;Opistho[...]abditis elegans\t4,72750\t0,57297\t52,23750\t6,33119\t5673\t62685

    my %taxa_from_hitmaps;
    while (<$fh>) {
        chomp;

        # skip comments
        next if ($_ =~ /^#/);

        # space to tab in name: "<acc>.<start>.<stop> <taxpath>"
        s/(\w+\.\d+\.\d+)\s/$1\t/;

        my @hitstats_line = split ("\t", $_);

        # skip 0 unambiguousReads
        next if ($hitstats_line[6] == 0);

        # truncate materialized path taxonomy string to 6 levels
        my $taxonshortstring;
        my @taxonstringarray = split (";", $hitstats_line[1]);
        if ( $#taxonstringarray > 5 ) {
            $taxonshortstring = join (";", @taxonstringarray[0..5]);
        } else {
            $taxonshortstring = $hitstats_line[1];
        }

        # add to hash or increment existing count
        if ( !exists $taxa_from_hitmaps{$taxonshortstring} ) {
            $taxa_from_hitmaps{$taxonshortstring} = $hitstats_line[6];
        } else {
            $taxa_from_hitmaps{$taxonshortstring} += $hitstats_line[6];
        }

    }
    close($fh);

    @taxa_from_hitmaps_sorted =
        sort { @$b[1] <=> @$a[1] }
        grep { if (@$_[1] < 3) { @xtons[@$_[1]-1]++;} @$_[1] > 2 }
        map  { [$_, $taxa_from_hitmaps{$_}] }
        keys %taxa_from_hitmaps;

    $chao1 =
        $#taxa_from_hitmaps_sorted +
        ($xtons[0] * $xtons[0]) / 2 / $xtons[1];

    msg("done...");
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
    if ($SEmode == 1) {
        $args = "-s $libraryNAME.$readsf.SSU.1.fq";
    } else {
        $args = "-1 $libraryNAME.$readsf.SSU.1.fq -2 $libraryNAME.$readsf.SSU.2.fq";
    }

    run_prog("spades",
             "-o $libraryNAME.spades -t $cpus -m 6 -k $kmer "
             . $args,
             "$libraryNAME.spades.out","&1"
         );

    msg("done...");
}

sub spades_parse {
    # getting spades output and reformatting it...
    msg("getting 16S phylotypes and their coverages...");

    # run barrnap once for each domain
    foreach ('bac', 'arch', 'euk') {
        run_prog("barrnap",
                 "--kingdom $_ --gene ssu --threads $cpus " .
		 "--evalue 1e-200 " .
                 "--reject 0.6 $libraryNAME.spades/scaffolds.fasta",
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
        $seqname =~ m/NODE_([0-9]*)_.*cov_([0-9\\.]*)_/;
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
        if ($SSU_total_pairs < 150000) {
            msg("Less than 300k SSU reads - using Emirge");
        } else {
            $cmd = "emirge_amp";
            msg("Warning: More than 25k SSU reads - using Emirge Amplicon");
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
                 . "tmp.$libraryNAME.SSU_all.fq",
                 "tmp.renamed.$libraryNAME.SSU_all.fq");

        $args = " -1 tmp.renamed.$libraryNAME.SSU_all.fq ";
    }

    run_prog($cmd,
             " $libraryNAME "
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
             "./$libraryNAME/iter.40/",
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

    run_prog("cat",
             "   $libraryNAME.emirge.final.fasta"
             . " $libraryNAME.spades_rRNAs.final.fasta",
             "$libraryNAME.all.final.fasta");

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
             . "-notrunclabels ",
             "tmp.$libraryNAME.clusterdbhits.out",
             "&1");

    msg("done...");
}

sub mafft_run {
    msg("creating alignment and tree...");

    run_prog("cat",
             "$libraryNAME.all.dbhits.NR97.fa "
             . "$libraryNAME.spades_rRNAs.final.fasta "
             . "$libraryNAME.emirge.final.fasta ",
             "$libraryNAME.SSU.collection.fasta");

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

sub run_plotscript {
    msg("generating histogram and tree graphics");

    #FIXME: crashes if .inserthistorgram contains no lines

    run_prog("plotscript",
             "--args $libraryNAME.SSU.collection.fasta.tree "
             . "$libraryNAME.inserthistogram ",
             "tmp.$libraryNAME.plotscript.out",
             "&1");
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

<h3>phyloFlash v1.5 by <a href="mailto:hgruber\@mpi-bremen.de">Harald Gruber-Vodicka</a> - high throughput phylogenetic screening using SSU rRNA gene(s) abundance(s)</h3>
ENDHTML

print {$fh} "<h1>Library name: ".$libraryNAME."</h1>\n";

print {$fh} <<ENDHTML;
<h2>Parameters</h2>

<h3>Input files</h3>
<table class="slimTable">
  <tr>
    <th>Forward read file</th>
ENDHTML

print {$fh} "    <td>".$readsf."</td>\n";

print {$fh} <<ENDHTML;
  </tr>
  <tr>
    <th>Reverse read file</th>
ENDHTML

print {$fh} "    <td>".$readsr."</td>\n";

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
  <tr>
    <th>Detected median insert size</th>
ENDHTML

print {$fh} "    <td>".$ins_used."</td>\n";
print {$fh} <<ENDHTML;
  </tr>
  <tr>
    <th>Used insert size</th>
ENDHTML

print {$fh} "     <td>".$ins_used."</td>\n";
print {$fh} <<ENDHTML;
  </tr>
  <tr>
    <th>Insert size standard deviation</th>
ENDHTML

print {$fh} "    <td>".$ins_std."</td>\n";
print {$fh} <<ENDHTML;
  </tr>
</table>

<h3><a href="#" id="histo-show" class="showLink" onclick="showHide('histo');return false;" title="Click to expand">Insert size histogram</a></h3>
<div id="histo" class="more">
ENDHTML

print {$fh} "<img src=\"".$libraryNAME.".inserthistogram.png\" />\n";
print {$fh} "<p><a href=\"".$libraryNAME.".inserthistogram.pdf\">PDF version</a></p>\n";
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

print {$fh} "  <tr>\n";
print {$fh} "    <th>Insert size histogram (plaintext)</th>\n";
print {$fh} "    <td>".$libraryNAME.".inserthistogram</td>\n";
print {$fh} "  </tr>\n";

print {$fh} "  <tr>\n";
print {$fh} "    <th>Insert size histogram (graphics)</th>\n";
print {$fh} "    <td>".$libraryNAME.".inserthistogram.png or .pdf</td>\n";
print {$fh} "  </tr>\n";

print {$fh} "  <tr>\n";
print {$fh} "    <th>Tree of recovered SSU sequences (Newick)</th>\n";
print {$fh} "    <td>".$libraryNAME.".SSU.collection.fasta.tree</td>\n";
print {$fh} "  </tr>\n";

print {$fh} "  <tr>\n";
print {$fh} "    <th>Tree of recovered SSU sequences (graphics)</th>\n";
print {$fh} "    <td>".$libraryNAME.".SSU.collection.fasta.tree.png or .pdf</td>\n";
print {$fh} "  </tr>\n";

print {$fh} <<ENDHTML;
</table>

<h2>Results</h2>

<h3><a href="#" id="taxa-show" class="showLink" onclick="showHide('taxa');return false;">Read mapping based detected higher taxa in order of appearance</a></h3>
<div id="taxa" class="more">
<table>
  <tr>
    <th><span class="withHoverText" title="Higher taxon found by SSU mapping to reference database">Taxon</span></th>
    <th><span class="withHoverText" title="No. reads mapped unambiguously to this taxonomic group">Unambig maps</span></th>
  </tr>
ENDHTML

## write list of higher taxa found
foreach my $taxonshortstring ( @taxa_from_hitmaps_sorted ) {
    # For each of these higher taxa
    print {$fh} "  <tr>\n";
    # print name of taxon
    print {$fh} "    <td>".@$taxonshortstring[0]."</td>\n";
    # print no . of reads mapping unambiguously
    print {$fh} "    <td>".@$taxonshortstring[1]."</td>\n";
    print {$fh} "  </tr>\n";
}

print {$fh} <<ENDHTML;
</table>
</div>
ENDHTML

## write list of assembled SSU sequences
print {$fh} <<ENDHTML;
<h3><a href="#" id="spades-show" class="showLink" onClick="showHide('spades');return false;">SSU assembly based taxa</a></h3>
<div id="spades" class="more">
<table>
  <tr>
    <th><span class="withHoverText" title="Name assigned to assembled SSU sequence">OTU</span></th>
    <th><span class="withHoverText" title="Per-base coverage of assembled SSU sequence, from assembler">coverage</span></th>
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
    print {$fh} "    <td>".$split_entry[1]."</td>\n";
    print {$fh} "    <td><a href=\"http://www.ncbi.nlm.nih.gov/nuccore/".$get_genbank[0]."\">".$split_entry[2]."</a></td>\n";	# Link to Genbank entry using accession no.
    #	print {$fh} "    <td>".$split_entry[2]."</td>\n";
    print {$fh} "    <td>".$split_entry[3]."</td>\n";
    print {$fh} "    <td>".$split_entry[4]."</td>\n";
    print {$fh} "    <td>".$split_entry[5]."</td>\n";
    print {$fh} "    <td>".$split_entry[6]."</td>\n";
    print {$fh} "  </tr>\n";
}

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

print {$fh} <<ENDHTML;
</table>
</div>
<h3><a href="#" id="tree-show" class="showLink" onclick="showHide('tree');return false;">Combined tree of sequences</a></h3>
<div id="tree" class="more">
ENDHTML

print {$fh} "<img src=\"".$libraryNAME.".SSU.collection.fasta.tree.png\" />\n";
print {$fh} "<p><a href=\"".$libraryNAME.".SSU.collection.fasta.tree.pdf\">PDF version</a></p>\n";

print {$fh} <<ENDHTML;
</div>
<h4>Please cite...</h4>
<p>List of citations including dependencies</p>
</body>
</html>
ENDHTML

close($fh);
}

######################### MAIN ###########################
welcome();
check_environment();
parse_cmdline();


my $timer = new Timer;

bbmap_fast_filter_run();
bbmap_fast_filter_parse();
bbmap_hitstats_parse();
spades_run();
spades_parse();
emirge_run();
emirge_parse();
vsearch_best_match();
vsearch_parse();
vsearch_cluster();
mafft_run();

$runtime = $timer->minutes;

print_report();
write_csv();
run_plotscript()    if ($html_flag);
write_report_html() if ($html_flag);


# cleanup of intermediate files
#msg("cleaning temp files...");
#system ("rm ./$libraryNAME.spades -r");
#system ("rm ./$libraryNAME -r");
#system ("rm tmp.* -r");

msg("Walltime used: $runtime with $cpus CPU cores");
msg("Thank you for using phyloFlash
You can find your results in $libraryNAME.*,
Your main result file is $libraryNAME.phyloFlash");

