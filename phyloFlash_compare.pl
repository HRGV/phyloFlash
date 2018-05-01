#!/usr/bin/env perl
use strict;
use warnings;
=head1 NAME

phyloFlash_compare.pl - Compare phyloFlash NTU results for multiple libraries

=head1 SYNOPSIS

=head1 DESCRIPTION

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
use File::Temp qw(tempfile);
use Data::Dumper;
use Archive::Tar;

my ($samfiles_str, $csvfiles_str, $tarfiles_str);
my @samfiles_arr;
my @csvfiles_arr;
my @tarfiles_arr;

my $task_opt = 'heatmap';
my $taxlevel = 4;
my $barplot_display = 10;
my $out_prefix = 'test.phyloFlash_compare';

my %ntuhash; # Hash of counts per taxon (primary key) and sample (secondary key)

=head1 INPUT FILES

The following options are mutually exclusive. If phyloFlash is run with the
I<-zip> option, the archives containing the results of each separate run can be
specified with the -zip option below. If the results are not compressed, you can
specify the NTU classification tables in CSV format.

=over 8

=item --csv I<FILES>

Comma-separated list of NTU abundance tables from phyloFlash runs. The files
should be named [LIBNAME].phyloFlash.NTUabundance.csv or
[LIBNAME].phyloFlash.NTUfull_abundance.csv

=item --zip I<FILES>

Comma-separated list of tar.gz archives from phyloFlash runs. These will be
parsed to search for the [LIBNAME].phyloFlash.NTUfull_abundance.csv files
within the archive, to extract the NTU classifications. This assumes that the
archive filenames are named [LIBNAME].phyloFlash.tar.gz, and that the LIBNAME
matches the contents of the archive.

=back

=head1 ARGUMENTS

=over 8

=item --task I<STRING>

Type of analysis to be run. Options: "heatmap", "barplot", "matrix", or a
recognizable substring thereof.

Default: None

=item --out I<STRING>

Prefix for output files.

Default: "test.phyloFlash_compare"

=item --level I<INTEGER>

Taxonomic level to perform the comparison. Must be an integer between 1 and 7.

Default: 4 ('Order')

=back

=head1 ARGUMENTS FOR BARPLOT

=over 8

=item --displaytaxa I<INTEGER>

Number of top taxa to display in barplot. Integer between 1 and 12.

Default: 10

=back

=cut

pod2usage (-verbose=>1, -exit=>1) if (!@ARGV);

GetOptions ("csv=s" => \$csvfiles_str,
            "sam=s" => \$samfiles_str,
            "zip=s" => \$tarfiles_str,
            "task=s" => \$task_opt,
            "level=i" => \$taxlevel,
            "displaytaxa=i" => \$barplot_display,
            "out=s" => \$out_prefix,
            "help|h" => sub { pod2usage(-verbose=>1); },
            "man" => sub { pod2usage(-verbose=>2); },
            ) or pod2usage(-verbose=>1, -exit=>1);

# Paths to R scripts
my $barplot_script = "$Bin/phyloFlash_barplot.R";
my $heatmap_script = "$Bin/phyloFlash_heatmap.R";
# msg ("Using barplot script: $barplot_script");
# msg ("Using heatmap script: $heatmap_script");

## MAIN ########################################################################

## Read in data ################################################################
if (defined $csvfiles_str) {
    # Read from CSV files 
    if (defined $samfiles_str || defined $tarfiles_str) {
        msg ("CSV files specified, ignoring SAM files and Tar archives");
    }
    my $csvfiles_aref = filestr2arr ($csvfiles_str);
    @csvfiles_arr = @$csvfiles_aref;
    foreach my $csv (@$csvfiles_aref) {
        if ($csv =~ m/^(.+)\.phyloFlash.NTU.*\.csv/) {
            my $samplename = $1;
            my $counts_aref = ntu_csv_file_to_arr($csv);
            my $counts_refactor_aref;
            if (defined $taxlevel) {
                # Refactor taxonstring to lower taxonomic level if requested
                $counts_refactor_aref = refactor_ntu_abundance_to_taxlevel($counts_aref,$taxlevel);
            } else {
                $counts_refactor_aref = $counts_aref;
            }
            ntu_csv_arr_to_hash($counts_refactor_aref,
                                $samplename,
                                \%ntuhash);
        } else {
            msg ("Filename of $csv does match standard name of a phyloFlash NTU abundance file");
        }
    }
}

#print Dumper \%ntuhash;

if (index ('barplot', $task_opt) != -1 ) {
    ## Barplot #################################################################
    my $out_aref = abundance_hash_to_array(\%ntuhash);
    my ($ntuall_fh, $ntuall_filename) = tempfile(DIR=>"."); # Temporarily put the temp file here for troubleshooting
    #my ($ntuall_fh, $ntuall_filename) = tempfile();
    open ($ntuall_fh, ">", $ntuall_filename) or die ("Cannot open file $ntuall_filename for writing");
    foreach my $line (@$out_aref) {
        print $ntuall_fh $line."\n";
    }
    msg ("NTU abundance by sample table written to temporary file: $ntuall_filename");
    close ($ntuall_fh);
    my @barplot_args = ("-f $ntuall_filename",
                        "-t $barplot_display",
                        "-o $out_prefix.barplot.pdf");
    my $barplot_cmd = join " ", ('Rscript', $barplot_script, @barplot_args);
    msg ("Plotting barplot: $barplot_cmd");
    system ($barplot_cmd);
    msg ("Barplot written to file: $out_prefix.barplot.pdf");
} else {
    msg ("No task specified, or task $task_opt not recognized");
}

## SUBS ########################################################################

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
