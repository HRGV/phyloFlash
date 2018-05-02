#!/usr/bin/env perl
use strict;
use warnings;
=head1 NAME

phyloFlash_compare.pl - Compare phyloFlash NTU results for multiple libraries

=head1 SYNOPSIS

# NTU abundance table CSV files as input

phyloFlash_compare.pl --csv LIB1.phyloFlash.NTUabundance.csv,LIB2.phyloFlash.NTUabundance.csv -task barplot

phyloFlash_compare.pl --csv LIB1.phyloFlash.NTUabundance.csv,LIB2.phyloFlash.NTUabundance.csv -task heatmap

# phyloFlash tar.gz archives as input

phyloFlash_compare.pl --zip LIB1.phyloFlash.tar.gz,LIB2.phyloFlash.tar.gz -task barplot

phyloFlash_compare.pl --zip LIB1.phyloFlash.tar.gz,LIB2.phyloFlash.tar.gz -task heatmap

# Help/manual pages

phyloFlash_compare.pl --help

phyloFlash_compare.pl --man

=head1 DESCRIPTION

Compare the taxonomic composition multiple metagenomic/transcriptomic libraries
using the phyloFlash NTU abundance results.

This script is a convenient wrapper for the R scripts I<phyloFlash_heatmap.R> and
I<phyloFlash_barplot.R>. More options are available when running those scripts
separately (see help messages by running the commands without arguments).

=cut

# This is a file reading and data-munging wrapper for the R scripts that are
# actually producing the plots

use Pod::Usage;
use Getopt::Long;
use FindBin qw($Bin);
use lib $FindBin::RealBin;
use PhyloFlash;
use Cwd;
use File::Basename;
use File::Temp qw(tempfile tempdir);
use Data::Dumper;
use Archive::Tar;

# Input files
my ($samfiles_str, $csvfiles_str, $tarfiles_str); # Raw comma-separated input string
my ($csvfiles_aref, $tarfiles_aref);              # Refs to arrays of input file paths
my @samfiles_arr;
my @csvfiles_arr;
my @tarfiles_arr;

my $tempdir; # Temp folder to put extracted files, if using -zip option

my $task_opt;
my %task_hash;
my $useSAM;
my $taxlevel = 4;
my $barplot_display = 5;
my $out_prefix = 'test.phyloFlash_compare';
my $outfmt = "pdf";

my $heatmap_clustersamples = 'ward.D';
my $heatmap_clustertaxa = 'ward.D';
my $heatmap_longtaxnames;
my $heatmap_minntucount = 50;

my %ntuhash; # Hash of counts per taxon (primary key) and sample (secondary key)
my %ntubysamplehash; # Hash of counts per sample (primary key) and taxon (secondary key)
my %totalreads_per_sample; # Hash of total reads counted per sample

=head1 ARGUMENTS

=head2 INPUT FILES

Specify the output files from phyloFlash runs for the libraries that you wish
to compare.

The options I<-csv> and I<-zip> are mutually exclusive. If phyloFlash is run with the
I<-zip> option, the archives containing the results of each separate run can be
specified with the -zip option below. If the results are not compressed, you can
specify the NTU classification tables in CSV format.

=over 8

=item --csv I<FILES>

Comma-separated list of NTU abundance tables from phyloFlash runs. The files
should be named [LIBNAME].phyloFlash.NTUabundance.csv or
[LIBNAME].phyloFlash.NTUfull_abundance.csv

If running the I<heatmap> option, this assumes that the [LIBNAME].phyloFlash.report.csv
files for the corresponding NTU abundance tables are also available in the same
folder(s).

=item --zip I<FILES>

Comma-separated list of tar.gz archives from phyloFlash runs. These will be
parsed to search for the [LIBNAME].phyloFlash.NTUabundance.csv files
within the archive, to extract the NTU classifications. This assumes that the
archive filenames are named [LIBNAME].phyloFlash.tar.gz, and that the LIBNAME
matches the contents of the archive.

=item --recalculate_NTU_from_SAM

Ignore NTU abundance tables in CSV format, and recalculate the NTU abundances
from SAM files in the compressed tar.gz phyloFlash archives. Useful if e.g.
phyloFlash was originally called to summarize the taxonomy at a higher level than
you want to use for the comparison.

Only works if the tar.gz archives from phyloFlash runs are specified with the
I<--zip> option above.

Default: No.

=back

=head2 ANALYSIS OPTIONS

=over 8

=item --task I<STRING>

Type of analysis to be run. Options: "heatmap", "barplot", "matrix", or a
recognizable substring thereof. Supply more than one option as comma-separated
list.

For option "heatmap", the [LIBNAME].phyloFlash.report.csv files must be either
in the same location as the NTU abundance CSV files or in the tar.gz archives,
as they contain metadata which is parsed by the script.

Default: None

=item --out I<STRING>

Prefix for output files.

Default: "test.phyloFlash_compare"

=item --outfmt I<STRING>

Format for plots (tasks 'barplot' and 'heatmap' only). Options: "pdf", "png"

Default: "pdf"

=item --level I<INTEGER>

Taxonomic level to perform the comparison. Must be an integer between 1 and 7.

Default: 4 ('Order')

=back

=head2 ARGUMENTS FOR BARPLOT

The R script I<phyloFlash_barplot.R> can be run directly; run the script without
arguments to see the built-in help message. However, the input file to the
barplot script is produced by I<phyloFlash_compare.pl> (i.e. this script).

=over 8

=item --displaytaxa I<INTEGER>

Number of top taxa to display in barplot. Integer between 1 and 12 is preferable.

Default: 5

=back

=head2 ARGUMENTS FOR HEATMAP

More options are available by using the R script I<phyloFlash_heatmap.R> directly,
or by sourcing it in the R environment. Run the R script without arguments to
see the built-in help message.

=over 8

=item --cluster-samples I<STRING>

Clustering method to use for clustering/sorting samples in heatmap. Options:
"alpha", "ward.D", "single", "complete", "average", "mcquitty", "median", "centroid",
or "custom".

"custom" will use the Unifrac-like abundance weighted taxonomic distances (the
distance matrix can be separately output with I<--task matrix>)

Default: "ward.D"

=item --cluster-taxa I<STRING>

Clustering method to use for clustering/sorting taxa. Options: "alpha", "ward",
"single", "complete", "average", "mcquitty", "median", "centroid".

Default: "ward.D"

=item --long-taxnames

Do not shorten taxa names to two last groups

=item --min-ntu-count I<INTEGER>

Sum up NTUs with fewer counts into a pseudo-NTU "Other".

Default: 50

=back

=head1 OUTPUT

=over 8

=item "heatmap"

=item "barplot"

Produces a PDF plot of top N taxa (default N = 5) in the libraries being
compared, as bar plots showing proportional abundances in each library.

=item "matrix"

Outputs a tab-separated table containing pairwise distances of all possible pairs
of samples. The distances are abundance-weighted Unifrac-like, using the
taxonomy tree in place of a phylogenetic tree, and treating all branches in the
taxonomy tree as having length 1.

Output columns are "library 1", "library 2", "distance"

=back

=head2 HELP MESSAGES

=over 8

=item --help|-h

Help message

=item --man

Manual page in pager

=item --version|-v

Report phyloFlash version

=back

=cut

pod2usage (-verbose=>0, -exit=>1) if (!@ARGV);

GetOptions ("csv=s" => \$csvfiles_str,
            "sam=s" => \$samfiles_str,
            "zip=s" => \$tarfiles_str,
            "task=s" => \$task_opt,
            "level=i" => \$taxlevel,
            "recalculate_NTU_from_SAM" => \$useSAM,                             ## TODO
            "displaytaxa=i" => \$barplot_display,
            "cluster-samples=s" => \$heatmap_clustersamples,
            "cluster-taxa=s" => \$heatmap_clustertaxa,
            "long-taxnames" => \$heatmap_longtaxnames,
            "min-ntu-count=i" => \$heatmap_minntucount,
            "out=s" => \$out_prefix,
            "outfmt=s" => \$outfmt,                                             # TODO for barplot script
            "help|h" => sub { pod2usage(-verbose=>1); },
            "man" => sub { pod2usage(-verbose=>2); },
            "version|v" => sub { welcome(); exit; }, 
            ) or pod2usage(-verbose=>0, -exit=>1);

# Paths to R scripts
my $barplot_script = "$Bin/phyloFlash_barplot.R";
my $heatmap_script = "$Bin/phyloFlash_heatmap.R";
# msg ("Using barplot script: $barplot_script");
# msg ("Using heatmap script: $heatmap_script");

## MAIN ########################################################################

## Catch exceptions ############################################################

if (!defined $task_opt) {
    pod2usage ("ERROR: Please specify tasks [barplot, heatmap, matrix] to option --task");
    pod2usage(-verbose=>1,-exit=>1);
}
if ($taxlevel > 7) {
    msg ("Taxonomic level is > 7, this is unlikely to provide a meaningful result");
    exit;
}
if ($barplot_display > 20) {
    msg ("Number of taxa to display in barplot is > 20, this is unlikely to be legible, but continuing anyway...");
}
if ($outfmt ne 'pdf' && $outfmt ne 'png') {
    msg ("WARNING: Invalid output format $outfmt specified. Should be either \"pdf\" or \"png\". Using \"pdf\"...");
    $outfmt = 'pdf';
}

## Read in data ################################################################

welcome();

parse_task_options($task_opt,\%task_hash);

if (defined $csvfiles_str) {
    # Read from CSV files
    if (defined $samfiles_str || defined $tarfiles_str) {
        msg ("CSV files specified, ignoring SAM files and Tar archives");
    }
    $csvfiles_aref = filestr2arr ($csvfiles_str);
    foreach my $csv (@$csvfiles_aref) {
        if ($csv =~ m/^(.+)\.phyloFlash.NTU.*\.csv/) {
            my $samplename = $1;
            ntu_csv_file_to_hash($csv, $samplename, $taxlevel, \%ntuhash);
        } else {
            msg ("Filename of $csv does not match standard name of a phyloFlash NTU abundance file");
        }
    }
} elsif (defined $tarfiles_str) {
    # Read from Tar archives
    if (defined $samfiles_str) {
        msg ("Tar archives specified, ignoring SAM files");
    }
    $tarfiles_aref = filestr2arr($tarfiles_str);
    # my $tempdir = tempdir (CLEANUP=>1);
    $tempdir = tempdir (DIR => "."); # Testing version
    foreach my $tar (@$tarfiles_aref) {
        if ($tar =~ m/^(.+)\.phyloFlash\.tar\.gz/) {
            my $samplename = $1;
            # Create Archive::Tar object
            my $tarhandle = Archive::Tar->new;
            $tarhandle -> read($tar);
            my $ntufilename = "$samplename.phyloFlash.NTUabundance.csv";
            my $pFreportcsvname = "$samplename.phyloFlash.report.csv";
            my $ntufull_filename = "$samplename.phyloFlash.NTUfull_abundance.csv";  # TO DO: Implement check for full NTU table
            # Check that NTU abundance table is in archive
            if ($tarhandle->contains_file($ntufilename)) {
                # Extract NTU abundance table to a temporary file
                $tarhandle->extract_file($ntufilename, "$tempdir/$ntufilename");
                ntu_csv_file_to_hash("$tempdir/$ntufilename", $samplename, $taxlevel, \%ntuhash);
                msg ("Extracting NTU abundance table $ntufilename to temporary folder $tempdir");
            } elsif ($tarhandle->contains_file($ntufull_filename)) {
                # Extract NTU full abundance table to a temporary file
                $tarhandle->extract_file($ntufull_filename, "$tempdir/$ntufull_filename");
                ntu_csv_file_to_hash("$tempdir/$ntufull_filename", $samplename, $taxlevel, \%ntuhash);
                msg ("Extracting NTU abundance table $ntufull_filename to temporary folder $tempdir");
            } else {
                msg ("Expected NTU abundance file $ntufilename not found in tar archive $tar");
            }
            if ($tarhandle->contains_file($pFreportcsvname)) {
                # Extract pF CSV report file required by heatmap tool
                $tarhandle->extract_file($pFreportcsvname, "$tempdir/$pFreportcsvname") if (defined $task_hash{'heatmap'});
            } else {
                msg ("Expected phyloFlash report CSV file $pFreportcsvname not found in tar archive $tar");
            }
        } else {
            msg ("Filename of archive $tar does not match standard name of a phyloFlash output archive");
        }
    }
}

#print Dumper \%ntuhash;

## Re-hash NTU hash by sample rather than taxon first
foreach my $taxon (keys %ntuhash) {
    foreach my $sample (keys %{$ntuhash{$taxon}}) {
        $ntubysamplehash{$sample}{$taxon} = $ntuhash{$taxon}{$sample};
        $totalreads_per_sample{$sample} += $ntuhash{$taxon}{$sample}; # Add to running total per sample
    }
}

#print Dumper \%ntubysamplehash;
#print Dumper \%totalreads_per_sample;

## Perform plotting ############################################################

if (defined $task_hash{'barplot'} ) {
    ## Barplot #################################################################
    my $out_aref = abundance_hash_to_array(\%ntuhash);
    #my ($ntuall_fh, $ntuall_filename) = tempfile(DIR=>"."); # Temporarily put the temp file here for troubleshooting
    my ($ntuall_fh, $ntuall_filename) = tempfile();
    open ($ntuall_fh, ">", $ntuall_filename) or die ("Cannot open file $ntuall_filename for writing");
    foreach my $line (@$out_aref) {
        print $ntuall_fh $line."\n";
    }
    msg ("NTU abundance by sample table written to temporary file: $ntuall_filename");
    close ($ntuall_fh);
    my $outfile_name = "$out_prefix.barplot.$outfmt";
    my @barplot_args = ("-f $ntuall_filename",
                        "-t $barplot_display",
                        "-o $outfile_name");
    my $barplot_cmd = join " ", ('Rscript', $barplot_script, @barplot_args);
    msg ("Plotting barplot: $barplot_cmd");
    system ($barplot_cmd);
    msg ("Barplot written to file: $outfile_name");
}

if (defined $task_hash{'matrix'} || $heatmap_clustersamples eq 'custom') {
    ## Matrix of taxonomic weighted Unifrac-like distances #####################
    
    # This distance matrix is also produced if used by heatmap option 'custom'
    
    # For each sample, re-parse the taxon strings and counts and put them in a
    # tree structure with counts per sample stored on each taxon node

    my %taxon_tree_with_counts;
    my @outarr;

    foreach my $sample (keys %ntubysamplehash) {
        my $countername = "_COUNT_$sample";
        encode_persample_counts_on_tree(\%taxon_tree_with_counts,
                                        \%{$ntubysamplehash{$sample}},
                                        $countername);
    }

    #print Dumper \%taxon_tree_with_counts;
    #print Dumper \%totalreads_per_sample;

    # For each pair of samples, output weighted taxonomic unifrac distance by
    # walking through taxonomic tree and comparing the counts for that given pair
    foreach my $sample1 (keys %totalreads_per_sample) {
        foreach my $sample2 (keys %totalreads_per_sample) {
            my $dist = calc_weight_tax_unifrac_pair (\%taxon_tree_with_counts,
                                                     "_COUNT_$sample1",
                                                     "_COUNT_$sample2",
                                                     $totalreads_per_sample{$sample1},
                                                     $totalreads_per_sample{$sample2},
                                                     );
             push @outarr, join "\t", ($sample1, $sample2, $dist);
        }
    }
    # Write matrix file
    open (my $fhmatrix, ">", "$out_prefix.matrix.tsv") or die ("Cannot open file $out_prefix.matrix.tsv for writing");
    print $fhmatrix join "\n", @outarr;
    close ($fhmatrix);
    
    msg ("Matrix of Unifrac-like abundance-weighted taxonomic distances written to file $out_prefix.matrix.tsv");
}

if (defined $task_hash{'heatmap'}) {
    ## Heatmap of samples vs. taxa #############################################
    
    # Define input CSV files for heatmap R script
    my $heatmap_csv_input;
    if (defined $csvfiles_aref) {
        # If user gave CSV files as input
        $heatmap_csv_input = join " ", @$csvfiles_aref;
    } elsif (defined $tarfiles_aref && defined $tempdir) {
        # If user supplied tar.gz archives as input, the relevant CSV files have
        # been extracted to the temp dir
        $heatmap_csv_input = "$tempdir/*.csv";
    }
    
    # Define name of output plot file
    my $outfile_name = "$out_prefix.heatmap.$outfmt";

    my @heatmap_args = ("-o $outfile_name",
                        "--library-name-from-file",
                        "--cluster-samples=$heatmap_clustersamples",
                        "--cluster-taxa=$heatmap_clustertaxa",
                        "--min-ntu-count=$heatmap_minntucount");
    
    
    # If using custom distance matrix...
    if ($heatmap_clustersamples eq 'custom') {
        push @heatmap_args, "--custom-distance-matrix-sample=$out_prefix.matrix.tsv";
    }
    
    # Push input file names/path last
    push @heatmap_args, $heatmap_csv_input;
    
    my $heatmap_cmd = join " ", ($heatmap_script, @heatmap_args);
    
    msg ("Plotting heatmap: $heatmap_cmd");
    system ($heatmap_cmd);
    msg ("Heatmap written to file: $outfile_name");
}



## SUBS ########################################################################

sub parse_task_options {
    # Parse comma-separated string of task names (or substrings thereof) and
    # hash to options hash
    my ($task_string,
        $href) = @_;
    my @task_arr = split /,/, $task_string;
    foreach my $task (@task_arr) {
        $href->{'heatmap'} ++ if (index ('heatmap',$task) != -1);
        $href->{'barplot'} ++ if (index ('barplot',$task) != -1);
        $href->{'matrix'} ++ if (index ('matrix',$task) != -1);
    }
}

sub ntu_csv_file_to_hash {
    # Wrapper to directly read CSV file into hash of abundances vs taxon * samples
    my ($file,
        $samplename,
        $taxlevel,
        $href) = @_;
    my $counts_aref = ntu_csv_file_to_arr($file);
    my $counts_refactor_aref;
    if (defined $taxlevel) {
        # Refactor taxonstring to lower taxonomic level if requested
        $counts_refactor_aref = refactor_ntu_abundance_to_taxlevel($counts_aref,$taxlevel);
    } else {
        $counts_refactor_aref = $counts_aref;
    }
    ntu_csv_arr_to_hash($counts_refactor_aref,
                        $samplename,
                        $href);
}

sub filestr2arr {
    # Convert comma-separated string of filenames to array
    my ($str) = @_;
    my @arr = split /,/, $str;
    return \@arr;
}

sub refactor_ntu_abundance_to_taxlevel {
    # Refactor NTU abundance table to higher taxonomic level, given an NTU
    # abundance array, return array of refactored CSV
    my ($aref,  # Ref to array of NTU taxon strings (semicolon-separated) and respective abundances (comma-separated)
        $level  # Taxonomic level to summarize, should be higher than what is in array
        ) = @_;
    my %hash;
    my @arr;
    foreach my $entry (@$aref) {
        my ($taxstr, $count) = split /,/, $entry;
        my @taxarr = split /;/, $taxstr;
        if ($level > scalar @taxarr) {
            msg ("Specified taxonomic level for summary is higher than available in taxstring $taxstr; padding to lowest level...");
            my $diff = $level - scalar @taxarr;
            push @taxarr, ("($taxarr[$#taxarr])") x $diff;
        }
        my $shortstr = join ";", @taxarr[0..$level-1]; # -1 to convert from 1-based to 0-based numbering
        $hash{$shortstr} += $count;
    }
    foreach my $key (keys %hash) {
        push @arr, join ",", ($key, $hash{$key});
    }
    return \@arr;
}

sub ntu_csv_file_to_arr {
    # Read NTU abundance CSV file and convert to an array
    # Return ref to that array
    my ($csvfile) = @_;
    my @arr;
    open (my $fh, "<", $csvfile) or die ("Cannot open CSV file $csvfile");
    while (my $line = <$fh>) {
        chomp $line;
        push @arr, $line;
    }
    close ($fh);
    return \@arr;
}

sub ntu_csv_arr_to_hash {
    # Convert CSV array of NTU abundances to hash of abundances keyed by
    # taxon (primary key) and sample ID (secondary key)
    my ($csv_aref,   # CSV file converted to array, ref thereto
        $name,       # Sample name
        $href,       # Ref to hash containing counts for all taxa and samples
        ) = @_;
    foreach my $rec (@$csv_aref) {
        my ($taxon, $count) = split /,/, $rec;
        $href->{$taxon}{$name} = $count;
    }
}

sub abundance_hash_to_array {
    my ($href) = @_;
    my @arr;
    foreach my $taxon (sort keys %$href) {
        foreach my $sample (sort keys %{$href->{$taxon}}) {
            my $count = $href->{$taxon}{$sample};
            push @arr, join ",", ($taxon, $sample, $count);
        }
    }
    return \@arr;
}


sub encode_persample_counts_on_tree {
    my ($href,          # Output hash tree
        $href_in,       # Input hash for a given sample
        $countername,   # Name for the counter ref for this sample
        ) = @_;
    foreach my $taxstring (keys %$href_in) {
        my $display_taxstring = $taxstring;
        $display_taxstring =~ s/^;ROOT;//;
        my $count = $href_in->{$taxstring};
        my @taxarr = split /;/, $display_taxstring;
        taxstring2hash_count($href,\@taxarr,$count,$countername);
    }
}

sub calc_weight_tax_unifrac_pair_raw {
    # Calculate taxonomic weighted raw Unifrac for a given pair of counts on
    # the tree
    my ($href, # Ref to hash containing the tree
        $counts1, # Name of the key storing counts of sample 1
        $counts2, # Name of the key storing counts of sample 2
        $total1,  # Total count for sample 1
        $total2,  # Total count for sample 2
        $raw_unifrac_sref,
        ) = @_;

    # Calculate the contribution to the raw unifrac value for this branch
    my ($node1count, $node2count);
    # Account for missing taxa in one sample (assign count = 0)
    if (defined $href->{$counts1}) { $node1count = $href->{$counts1}; } else { $node1count = 0; }
    if (defined $href->{$counts2}) { $node2count = $href->{$counts2}; } else { $node2count = 0; }
    $$raw_unifrac_sref += abs($node1count/$total1 - $node2count/$total2);

    # Recursion
    foreach my $key (keys %$href) {
        if (ref ($href->{$key}) eq 'HASH') {
            calc_weight_tax_unifrac_pair_raw(\%{$href->{$key}},$counts1,$counts2,$total1,$total2,$raw_unifrac_sref);
        }
    }
}

sub calc_weight_tax_unifrac_pair {
    my ($href,
        $countname1,
        $countname2,
        $in1_total,
        $in2_total,
        ) = @_;
    my $raw_unifrac;
    calc_weight_tax_unifrac_pair_raw($href,
                                     $countname1,
                                     $countname2,
                                     $in1_total,
                                     $in2_total,
                                     \$raw_unifrac,
                                     );

    # The equivalent of the scaling factor D would be the total taxonomic rank depth
    # of the tax tree, e.g. 7 for species, 4 for order. The branch length between
    # each rank is implicitly = 1 in the calc_weight_tax_unifrac_pair() procedure
    return $raw_unifrac / $taxlevel;
}


sub taxstring2hash_count {
    # Convert taxonomy string to a nested hash structure recursively
    my ($href, # Reference to hash
        $aref, # Reference to taxonomy string as array
        $count, # Taxon abundance
        $countername, # Name for the taxon total count key
        ) = @_;
    my $taxon = shift @$aref;
    if (@$aref) { # Recursion
        taxstring2hash_count  (\%{$href->{$taxon}}, $aref, $count,$countername);
        $href->{$taxon}{$countername} += $count;
    } else { # End condition - count number of occurrences of this taxon
        $href->{$taxon}{$countername} += $count;
    }
}

sub welcome {
    my $message = "This is phyloFlash_compare.pl from phyloFlash v$VERSION";
    print STDERR $message;
    print STDERR "\n";
}
