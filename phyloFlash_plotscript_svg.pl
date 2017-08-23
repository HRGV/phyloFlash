#!/usr/bin/env perl

=head1 NAME

phyloFlash_plotscript_svg.pl - Produce SVG-formatted plots for phyloFlash
pipeline

=head1 SYNOPSIS

B<phyloFlash_plotscript_svg.pl> -tree F<<FILE>> -hist F<<FILE>> -bar F<<FILE>> -pie F<<FILE>>

B<phyloFlash_plotscript_svg.pl> -help

B<phyloFlash_plotscript_svg.pl> -man

=head1 DESCRIPTION

Internal script used by B<phyloFlash.pl> to produce SVG-formatted plots for
HTML report file.

Output files are simply input filename with .svg suffix, silently overwritten.

=cut

use strict;
use warnings;
use Pod::Usage;
use POSIX qw (ceil floor);
use List::Util qw(min max);
use Getopt::Long;
use Math::Trig qw(pi cylindrical_to_cartesian);

# Plot insert size histogram and guide tree in SVG format without additional dependencies on R packages

# Input arguments
my ($treefile, $histofile, $barfile, $piefile, $title, $decimalcomma, $pipemode, $assemcov, $unassem_count);
my ($nbreaks, $barminprop) = (undef, 0.2); # Default values for params
my ($plotheight, $plotwidth, $plotcolor);
if (!@ARGV) { # Help msg if no arguments given
    pod2usage (2);
    exit();
}
GetOptions("tree|t=s" => \$treefile,        # Guide tree from MAFFT
           "assemcov=s" => \$assemcov,      # CSV file libNAME.phyloFlash.extractedSSUclassifications.csv
           "unassemcount=i" => \$unassem_count, # Number of unassembled reads
           "hist|h=s" => \$histofile,       # Insert size histogram from BBmap (PE reads only)
           "bar|r=s" => \$barfile,          # Table of counts to make barplot
           "pie|p=s" => \$piefile,          # Table of counts to make donut/piechart
           "pipe=s" => \$pipemode,          # Pipe mode - take input from STDIN and write to STDOUT - specify type of output
           "title=s" => \$title,            # Title for plot
           "height=i" => \$plotheight,      # Optional height for plot
           "width=i" => \$plotwidth,        # Optional width for plot in pixels
           "color=s" => \$plotcolor,        # Optional fill color for plot (currently implemented only for histogram)
           "decimalcomma" => \$decimalcomma,# BBmap is locale-aware and may produce histogram files with decimal comma!
                                            # Perl does not use locales unless requested so the other inputs should be safe
           "breakpoints|b=i" => \$nbreaks,  # Optional: manually specify number of breakpoints in histogram (e.g. 30)
           "help" => sub { pod2usage(1) },
           "man" => sub { pod2usage(-exitval=>0, -verbose=>2) },
           ) or pod2usage(2);

=head1 OPTIONS

=over 15

=item -tree F<<FILE>>

Phylogenetic tree plot from Newick-formatted tree. Does not support node or
branch labels, branch lengths required. Tree is oriented with root on left
and leaf labels on right. Height of plot scales with number of leaves.

=item -assemcov F<<FILE>>

CSV file containing coverage stats of assembled SSU sequences, from phyloFlash
output, for drawing bubbles representing read coverage per assembled sequence
on the tree.

When supplied, beneath the node for each leaf representing an assembled SSU
sequence will be a circle whose area represents the abundance of that sequence
(i.e. number of reads mapping to it in the re-mapping step) in the read library.

=item -unassemcount I<INT>

Number of reads that are not mapping to any of the assembled SSU sequences.
This is to draw a bubble representing coverage of unassembled reads, when
coverage of assembled sequences is also supplied via -assemcov parameter.

=item -hist F<<FILE>>

Histogram plot from TAB-separated histogram output file, e.g. those produced by
BBmap. Column 1: bin values, column 2: counts/frequencies. The counts must
already be binned into counts or frequencies! The counts can be re-binned into
new bins for plotting, but the user is responsible for making sure that the
number of new bins < number of original bins.

=item -bar F<<FILE>>

Interactive bar plot from CSV file. Column 1: label, column 2: counts. The
text labels are adjacent to the bars and their size and opacity are scaled to
the bar height to avoid overlapping text labels. However on mouseover (when
viewed in web browser) the box will be highlighted and the corresponding text
label will enlarge to legible size.

=item -pie F<<FILE>>

Pie chart from CSV file. Column 1: Label, column 2: counts/percentage/ratio.
Will automatically take total of the numbers provided. I.e. if using
percentages, ensure that they add up to 100%.

=item -pipe I<OUTPUT TYPE>

Use "pipe" mode. Input is read from STDIN and written to STDOUT. Use this
option to specify type of plot to produce: "tree", "hist", "bar", or "pie".
Naturally this can only read and write one file at a time.

=item -title="I<STRING>"

Optional title for plot. Enclose string in quotation marks if title has spaces.
Default: (empty)

=item -height I<INT>

=item -width I<INT>

Optional width and height of plot in pixels.
Defaults: Built-in defaults for each plot type.

=item -color="I<STRING>"

Optional fill color for plot. Currently implemented only for histogram.
Default: Built-in defaults for each plot type

=item -breakpoints I<INT>

Specify number of breakpoints in histogram.
Default: Sturges algorithm to calculate optimal breakpoints

=item -decimalcomma

Use comma as decimal separator (only used for histogram input). BBMap is locale-
aware and uses decimal comma in certain locales (e.g. Germany and France).
However, Perl does not unless explicitly requested. This option will replace
decimal commas with decimal periods in histogram input files only.

=item -help

This help message

=item -man

Manual page (identical to this help message)

=back

=cut


## MAIN ########################################################################

my $delim = "comma";
my $dl = defined $delim && $delim eq "tab" ? "\t" : ","; # Check if TSV

if (defined $pipemode) {
    # If in pipe mode, take input from STDIN and write to STDOUT
    # Get input into array
    my @input_arr;
    while (<>) {
        chomp;
        push @input_arr, $_;
    }
    # Process input depending on specified graphic type
    if ($pipemode eq "pie") {
        do_pie_or_bar_chart (\@input_arr, "array", $dl, "pie", $title);
    } elsif ($pipemode eq "bar") {
        do_pie_or_bar_chart (\@input_arr, "array", $dl, "bar", $title);
    }
} else {
    # Otherwise in "file" mode read and write to specified files
    if (defined $histofile) {
        do_histogram_plots($histofile, $title, $plotwidth, $plotheight, $plotcolor);
        #if (defined $plotheight || defined $plotwidth) {
        #    # If custom plot height and/or width specified
        #    do_histogram_plots($histofile, $title, $plotwidth, $plotheight);
        #} else {
        #    # Else use defaults (currently 240 x 240)
        #    do_histogram_plots($histofile, $title, 240, 240);
        #}
    }
    if (defined $treefile) {
        do_phylog_tree($treefile);
    }
    if (defined $barfile) {
        do_pie_or_bar_chart($barfile, "file", $dl, "bar", $title);
    }
    if (defined $piefile) {
        do_pie_or_bar_chart($piefile, "file", $dl, "pie", $title);
    }
}

## SUBROUTINES FOR PIECHART ###################################################

sub do_pie_or_bar_chart {
    my ($input, $mode, $dl, $type, $title) = @_;
    my $in_href;
    my $outname;
    my $outfh;
    # Array mode or infile mode
    if ($mode eq "file") {
        $in_href = csv2hash($input, $dl);
        $outname = $input.".svg";
        open ($outfh, ">", $outname) or die ("Cannot open file $outname for writing: $!");
    } elsif ($mode eq "array") { # pipe mode
        $in_href = csv_arr2hash ($input, $dl);
        $outfh = *STDOUT;
    }
    if ($type eq "pie") {
        hash2pie ($in_href, $outfh, $title);
    } elsif ($type eq "bar") {
        hash2barchart ($in_href, $outfh, $title);
    }
    close ($outfh) if $mode eq "file";
}

sub pc2xy {
    # Convert percentages along a circle to x-y coordinates
    my ($pc,    # Fractional (0 to 1) position along the circle
        $cx,    # Center of circle
        $cy,
        $rad    # Radius
        ) = @_;
    # Reminder - offset is counter-clockwise and starts from 3 o'clock
    $pc = 1 - $pc;      # Account for counter-clockwiseness
    # Use functions from Math::Trig
    my $theta = 2 * pi * $pc;
    my ($xoff, $yoff, $discard) = cylindrical_to_cartesian ($rad, $theta, 0);
    my ($x, $y) = ($xoff + $cx, $yoff + $cy); # Offset by circle center
    return ($x, $y);
}

sub hash2pie {
    my ($csv_href,  # Name of input hash
        $fh,        # Filehandle for output print
        $title,     # Title text, if available
        ) = @_;

    # Calculate cumulative percentages of input sorted by counts
    my $cumul_href = counthash_cumul_sum ($csv_href, 1, "counts");

    # Dereference output into arrays for plotting
    my @labels_arr = @{$cumul_href->{"labels"}};     # Array of text labels
    my @cumul_pc_arr = @{$cumul_href->{"cumul_pc"}}; # Array of cumulative percentages
    my @pc_arr = @{$cumul_href->{"pc"}};
    my @counts_arr = @{$cumul_href->{"counts"}};     # Array of counts
    my @colors_arr = @{$cumul_href->{"color"}} if defined $cumul_href->{"color"}; # Array of colors, if defined

    # SVG plot preferences
    my $viewBox_width = 240;     # width
    my $viewBox_height = 240;    # height
    my $margin = 60;
    my $font_size = 14;

    # Take shorter dimension, and calculate the pie diameter and circumference
    my $diam = $viewBox_width < $viewBox_height ? $viewBox_width : $viewBox_height;
    $diam = $diam - 2 * $margin;
    my $rad = $diam / 2;
    my $cx = $viewBox_width / 2;
    my $cy = $viewBox_height / 2;
    my $circum = $diam * 3.14159265; # probably precise enough
    #my $circum = $diam * 22 / 7; # probably not precise enough
    my $stroke_width = 0.75 * $rad; # stroke width is proportion of radius

    # Calculate stroke-dasharray and stroke-dashoffset
    # Tips: http://openstudio.redhat.com/scratch-made-svg-donut-pie-charts-in-html5/
    my @dasharray_arr;
    my @dashoffset_arr;
    for (my $i=0; $i <= $#labels_arr; $i++) {
        push @dashoffset_arr, $cumul_pc_arr[$i] * $circum;
        my $stroke = $pc_arr[$i] * $circum;
        my $space = $circum - $stroke;
        push @dasharray_arr, ("$stroke $space");
    }

    # SVG parameters
    my $viewBox = "0 0 $viewBox_width $viewBox_height";
    my $svg_open = "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"$viewBox\" height=\"100%\" >\n";
    my $donut_std_params = "cx=\"$cx\" ".
                            "cy=\"$cy\" ".
                            "r=\"$rad\" ".
                            "fill=\"transparent\" ".
                            "stroke-width=\"$stroke_width\" ";
    # Start SVG file
    print $fh $svg_open;
    # Print title
    if (defined $title) {
        print $fh "<text ".
                  "style=\"fill:black;font-size:$font_size;text-anchor:middle;font-weight:bold;\" ".
                  "x=\"".($viewBox_width / 2)."\" ".
                  "y=\"18\" ".
                  ">".
                  $title.
                  "</text>\n";
    }
    for (my $i=0; $i <= $#labels_arr; $i++) {
        my $color;
        if (defined $colors_arr[$i]) {
            $color = $colors_arr[$i];
        } else {
            my @rand_colors = (int(rand(256)),int(rand(256)),int(rand(256)));
            $color = join (",", @rand_colors);
            $color = "rgb(".$color.")";
        }
        print $fh "<circle class=\"donut-segment\" ".
                  $donut_std_params.
                  "stroke=\"$color\" ".
                  "stroke-dasharray=\"".$dasharray_arr[$i]."\" ".
                  "stroke-dashoffset=\"".$dashoffset_arr[$i]."\" ".
                  "></circle>\n";
    }
    # Print labels in separate loop because they must be on top of donut segments
    my $text_style_base = "fill:black;font-size:".$font_size."px;";
    for (my $i=0; $i <= $#labels_arr; $i++) {
        my ($x_text, $y_text) = pc2xy($cumul_pc_arr[$i] - $pc_arr[$i] / 2,
                                      $cx,
                                      $cy,
                                      $rad);
        print $fh "<text ".
                "x=\"$x_text\" ".
                "y=\"$y_text\" ".
                "style=\"".$text_style_base."text-anchor:middle;\" ".
                ">".
                $labels_arr[$i].
                "</text>\n";
    }
    print $fh "</svg>\n";
    # Mysteries:
    # why does stroke-dashoffset move counterclockwise????!?
    # why does SVG y coordinate start from the top????!??
}

### SUBROUTINES FOR BARCHART ##################################################

sub csv2hash {
    # CSV file to hash
    my ($infile, $delim) = @_;
    my %hash;
    open (my $fhin, "<", $infile) or die ("Cannot open file $infile for reading: $!");
    while (<$fhin>) {
        chomp;
        my @splitline = split "$delim";
        $hash{$splitline[0]} = $splitline[1];
    }
    close($fhin);
    return (\%hash);
}

sub csv_arr2hash {
    # CSV array to hash
    my ($aref, $delim) = @_;
    my %hash;
    foreach my $line (@$aref) {
        my @splitline = split "$delim", $line;
        $hash{$splitline[0]} = $splitline[1];
    }
    return (\%hash);
}

sub counthash_cumul_sum {
    my ($inhref,    # Ref to hash of raw counts (keyed by label)
        $maxprop,   # Maximum cumulative percentile to display taxonomic breakdown (<= 1)
        $sortby,    # Sort by what? either "counts" or "name"
        ) = @_;

    # Initialize variables
    my $minprop = 1 - $maxprop;
    my $total = 0;
    my %valhash;
    # Sort input by either counts or names
    my $running_total = 0;
    my @sortkey;
    if ($sortby eq "counts") {
        @sortkey = sort {$inhref->{$b} <=> $inhref->{$a}} keys %$inhref;
    }  elsif ($sortby eq "name") {
        @sortkey = sort {$a cmp $b} keys %$inhref;
    }
    # Add up total
    foreach my $key (@sortkey) {
        $total += $inhref->{$key};
    }
    # Add up cumulative totals
    my $index = 0;
    my $last_index;
    my @labels_arr;     # Array of text labels
    my @pc_arr;         # Array of percentages
    my @cumul_pc_arr;   # Array of cumulative percentages
    my @counts_arr;     # Array of raw counts
    my $other_count = 0;
    for (my $i=0; $i <= $#sortkey; $i++) {
        $valhash{$sortkey[$i]}{"raw"} = $inhref->{$sortkey[$i]};
        $running_total += $valhash{$sortkey[$i]}{"raw"}; # Update cumulative total
        $valhash{$sortkey[$i]}{"cumul"}  = $running_total;
        # Convert values to fractions of total
        $valhash{$sortkey[$i]}{"pc"} = $valhash{$sortkey[$i]}{"raw"} / $total;
        $valhash{$sortkey[$i]}{"cumul_pc"} = $valhash{$sortkey[$i]}{"cumul"} / $total;
        my $leftover = 1;
        $leftover = 1 - $valhash{$sortkey[$i-1]}{"cumul_pc"} if $i > 0; # Some gymnastics
        unless ($leftover < $minprop) {
            push @labels_arr, $sortkey[$i];
            push @pc_arr, $valhash{$sortkey[$i]}{"pc"};
            push @cumul_pc_arr, $valhash{$sortkey[$i]}{"cumul_pc"};
            push @counts_arr, $valhash{$sortkey[$i]}{"raw"};
        } else {
            # If there are taxa below minimum cumulative count, add counts to group "other"
            $other_count += $valhash{$sortkey[$i]}{"raw"};
        }
    }
    # Add last value for "Other" if it is defined
    if ($other_count > 0) {
        push @labels_arr, "Other taxa (below threshold)";
        push @cumul_pc_arr, 1;
        push @pc_arr, 1 - $#cumul_pc_arr;
        push @counts_arr, $other_count;
    }
    # Output is a hash of array references - this is necessary to preserve sort order
    my %outhash = ("labels" => \@labels_arr,
                   "cumul_pc" => \@cumul_pc_arr,
                   "pc" => \@pc_arr,
                   "counts" => \@counts_arr);

    return (\%outhash);
}

sub hash2barchart {
    # Read input
    my ($csv_href,       # Input CSV file
        $fh,            # Filehandle for output print
        $title,         # Optional title
        ) = @_;

    # Set preferences
    my $maxprop = 1;              # Maximum cumulative percentile to display taxonomic breakdown
    my $minprop = 1 - $maxprop;

    # SVG plot preferences
    my $orientation = "v";          # Horizontal (h) or vertical (v) alignment of figure long axis
    my $box_proportion = 0.10;      # Proportion of figure viewbox occupied by bar vs. text
    my $viewBox_longaxis = 300;     # Long dimension of the viewbox (parallel to main axis)
    my $viewBox_shortaxis = 350;    # Short dimension of the viewbox (perpendicular to main axis)
    my $margin = 5;

    # Calculate cumulative percentages of input sorted by counts
    my $cumul_href = counthash_cumul_sum ($csv_href, $maxprop, "counts");

    # Dereference output into arrays for plotting
    my @labels_arr = @{$cumul_href->{"labels"}};     # Array of text labels
    my @x1_arr = @{$cumul_href->{"cumul_pc"}};         # Array of values for right side of bars
    my @counts_arr = @{$cumul_href->{"counts"}};     # Array of counts

    # Generate array of values for left side of bars
    my @x0_arr = @x1_arr;
    pop @x0_arr;
    unshift @x0_arr, 0;
    my @widths_arr;
    my @y0_arr = (0) x scalar @x0_arr;
    my @y1_arr = (1) x scalar @x1_arr;

    # SVG plot parameters
    my @boxval = (0, 1, 0, 1);
    my $viewBox;
    my @box;

    if ($orientation eq "h") { # for horizontal bars
        $viewBox = "0 0 $viewBox_longaxis $viewBox_shortaxis"; # x y width height
        @box = ($margin,
                $viewBox_longaxis - $margin,
                $viewBox_shortaxis - $viewBox_shortaxis*$box_proportion + $margin,
                $viewBox_shortaxis - $margin
                ); # left right bottom top coordinates
    } else {
    # for vertical bars
        $viewBox = "0 0 $viewBox_shortaxis $viewBox_longaxis"; # x y width height
        @box = ($margin,
                $viewBox_shortaxis*$box_proportion + $margin,
                $margin,
                $viewBox_longaxis - $margin
                ); # left right bottom top coordinates
    }
    my @viewBox_arr = split " ", $viewBox;

    # Convert to coordinates
    my ($x0_rescale_aref, $y0_rescale_aref, $x1_rescale_aref, $y1_rescale_aref);
    if ($orientation eq "h") {
        ($x0_rescale_aref, $y0_rescale_aref) = val2coord ($viewBox, \@box, \@boxval, \@x0_arr, \@y0_arr);
        ($x1_rescale_aref, $y1_rescale_aref) = val2coord ($viewBox, \@box, \@boxval, \@x1_arr, \@y1_arr);
    } else {
        @x0_arr = map { 1 - $_ } @x0_arr; # Flip coordinates so that most abundant taxon appears on top
        @x1_arr = map { 1 - $_ } @x1_arr;
        ($x0_rescale_aref, $y0_rescale_aref) = val2coord ($viewBox, \@box, \@boxval, \@y0_arr, \@x1_arr);
        ($x1_rescale_aref, $y1_rescale_aref) = val2coord ($viewBox, \@box, \@boxval, \@y1_arr, \@x0_arr);
    }

    # Create rectangle values
    my %rect_vals;
    for (my $i=0; $i <= $#labels_arr; $i++) {
        $rect_vals{$i}{"label"} = $labels_arr[$i];
        $rect_vals{$i}{"x"} = $x0_rescale_aref->[$i];
        $rect_vals{$i}{"width"} = abs ($x1_rescale_aref->[$i] - $x0_rescale_aref->[$i]);
        $rect_vals{$i}{"y"} = $y1_rescale_aref->[$i];
        $rect_vals{$i}{"height"} = abs($y0_rescale_aref->[$i] - $y1_rescale_aref->[$i]);
        $rect_vals{$i}{"counts"} = $counts_arr[$i];
        $rect_vals{$i}{"id"} = "rect$i";    # ID of rect object in SVG file
    }

    my $svg_open = "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"$viewBox\" width=\"100%\" height=\"100%\">\n";
    #my $style_base = "fill-opacity:0.5;stroke:rgb(0,0,0);stroke-width:1;"; # Style for histogram bars
    my $style_base = "fill-opacity:0.5;".
                     #"stroke:rgb(0,0,0);".
                     "stroke-width:1;";

    print $fh $svg_open;
    # Print a title
    print $fh "<text ".
              "style=\"fill:black;text-anchor:middle;font-weight:bold;font-size:10px;\" ".
              "x=\"".($viewBox_arr[2]/2)."\" ".
              "y=\""."14"."\" ".
              ">".
              $title.
              "</text>\n";

    foreach my $rect (sort {$a <=> $b} keys %rect_vals) {
        my @bar_colors;
        if ($rect_vals{$rect}{"label"} eq "Other taxa (below threshold)") {
            @bar_colors = (155,155,155); # "Other taxa" should appear grey
        } else {
            # Named taxa have nice colors
            @bar_colors = (int(rand(256)),int(rand(256)),int(rand(256)));
        }
        my $bar_color="fill:rgb(".join(",",@bar_colors).");";
        my $style=$style_base.$bar_color;
        print $fh "<rect ".
                  "id=\"".$rect_vals{$rect}{"id"}."\" ".
                  "x=\"".$rect_vals{$rect}{"x"}."\" ".
                  "y=\"".$rect_vals{$rect}{"y"}."\" ".
                  "width=\"".$rect_vals{$rect}{"width"}."\" ".
                  "height=\"".$rect_vals{$rect}{"height"}."\" ".
                  "style=\"$style\" ".
                  "onmouseover=\"evt.target.setAttribute('stroke','red');\" ". # Surrounded by red border on mouseover
                  "onmouseout=\"evt.target.setAttribute('stroke','none');\" ".
                  "/>\n";
        # Print taxonomy string label
        my ($x_text, $y_text);
        my ($x_text_counts, $y_text_counts);
        my $font_size =  8;
        #my $text_style_base = "fill:black;font-size:".$font_size."px;";
        my $text_style_base = "fill:black;";
        my $text_style;

        # If bars are too narrowly spaced, make text more transparent and
        # smaller, but expand size and grow opaque on mouseover of bars for clarity
        # Tips from: http://www.petercollingridge.co.uk/data-visualisation/mouseover-effects-svgs
        # Optional params to make text smaller font size if bars too narrowsly spaced
        my $mouseover = "<set attributeName=\"font-size\" ".
                        "to=\"$font_size\" ".
                        "begin=\"".$rect_vals{$rect}{"id"}.".mouseover\" ".
                        "end=\"".$rect_vals{$rect}{"id"}.".mouseout\" ".
                        "/>";
        # ALL text labels set to fully opaque on mouseover of bar
        my $attr_text = "<set attributeName=\"fill-opacity\" ".
                        "to=\"1\" ".
                        "begin=\"".$rect_vals{$rect}{"id"}.".mouseover\" ".
                        "end=\"".$rect_vals{$rect}{"id"}.".mouseout\" ".
                        "/>";
        my $text_test_dimension;

        if ($orientation eq "h") {
            $text_style = $text_style_base."writing-mode:tb;";
            $x_text = $margin + $rect_vals{$rect}{"x"} + $rect_vals{$rect}{"width"}/2;
            $y_text = $margin + $rect_vals{$rect}{"y"} + $rect_vals{$rect}{"height"} + $font_size;
            $x_text_counts = $x_text;
            $y_text_counts = $margin + $rect_vals{$rect}{"y"} + $rect_vals{$rect}{"height"}/2;
            $text_test_dimension = $rect_vals{$rect}{"width"};
        } elsif ($orientation eq "v") {
            $text_style = $text_style_base;
            $x_text = $margin + $rect_vals{$rect}{"width"} + $font_size;
            $y_text = $margin + $rect_vals{$rect}{"y"} + $rect_vals{$rect}{"height"} / 2;
            $x_text_counts = $margin + $rect_vals{$rect}{"width"} / 2;
            $y_text_counts = $y_text;
            $text_test_dimension = $rect_vals{$rect}{"height"};
        }

        # Check if bars are too narrow
        if ($text_test_dimension < $font_size) {
            # If bars are too narrow, set text labels to be shorter/thinner
            $attr_text = $attr_text.$mouseover;
            my $new_opacity = 0.5 * ($text_test_dimension / $font_size);
            $text_style = $text_style."fill-opacity:$new_opacity;";
            $text_style = $text_style."font-size:".$text_test_dimension."px;";
        } else {
            # All text labels are of 0.6 opacity to begin with, for animation
            # to be more consistent to the eye
            $text_style = $text_style."fill-opacity:0.6;";
            $text_style = $text_style."font-size:".$font_size."px;";
        }

        print $fh "<text ".
                  "x=\"$x_text\" ".
                  "y=\"$y_text\" ".
                  "style=\"$text_style"."text-anchor:start;\" ".
                  ">".
                  $attr_text. # Attribute effects on mouseover
                  $rect_vals{$rect}{"label"}.
                  "</text>\n";
        # Print counts labels
        print $fh "<text ".
                  "x=\"$x_text_counts\" ".
                  "y=\"$y_text_counts\" ".
                  "style=\"$text_style;text-anchor:middle;\" ".
                  ">".
                  $attr_text. # Attribute effects on mouseover
                  $rect_vals{$rect}{"counts"}.
                  "</text>\n";
    }
    print $fh "</svg>\n";
}

### SUBROUTINES FOR HISTOGRAM #################################################

sub do_histogram_plots {
    # Get input parameters
    my ($infile,    # Input filename
        $title,     # Optional title
        $width,     # Optional SVG plot width - passed to viewBox
        $height,    # Optional SVG plot height -passed to viewBox
        $color,     # Optional SVG fill color
       ) = @_;

    # Check if options defined, else substitute defaults
    $width = defined $width ? $width : 240;
    $height = defined $height ? $height : 240;
    $color = defined $color ? $color : "rgb(155,155,155)";

    # SVG and Plot parameters for histograms
    my @viewBox_arr = (0, 0, $width, $height);         # Viewbox parameter for SVG header - x y width height
    my $viewBox = join " ", @viewBox_arr;
    my $svg_open = "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"$viewBox\" width=\"100%\" height=\"100%\">\n";
    # Bounding box coordinates for plot area
    my @box_coords = (40,               # Left margin -  leave space for labels
                      $width - 20,      # Right margin
                      20,               # Bottom margin
                      $height - 40,     # Top margin - leave space for title
                      ); #
    my $fill_style = "fill:$color;fill-opacity:0.5;stroke:none"; # Style for histogram bars

    # Plot histogram
    my $infile_out = "$infile.svg";       # Append .svg to get output file name
    my %histo_hash;                             # Define hash to hold data
    read_hist($infile, \%histo_hash);        # Read histogram file into memory
    open (my $histo_fh, ">", $infile_out)    # Open file for printing
        or die ("Cannot write to output file $infile_out: $!");
    print $histo_fh $svg_open;                  # Print SVG header
    if (defined $title) {                       # Print title if defined
        print $histo_fh "<text ".
                        "style=\"fill:black;font-size:14px;text-anchor:middle;font-weight:bold;\" ".
                        "x=\"".($viewBox_arr[2]/2)."\" ".
                        "y=\"18\" ".
                        ">".
                        $title.
                        "</text>\n";
    }
    draw_histogram ($viewBox, \@box_coords, \%histo_hash, $fill_style, $histo_fh, $nbreaks);
    print $histo_fh "</svg>\n";                 # Closing SVG tag
    close ($histo_fh);                          # Close file
}

sub draw_histogram {
    my ($viewBox,   # Viewbox string for SVG tag
        $box_aref,  # Ref to array (left right bottom top) of bounding box coords
        $href,      # Ref to histogram to draw
        $style,     # Style string for rect element
        $handle,    # Print output handle
        $nbreaks    # Number of breaks (if undefined, default to Sturges)
        ) = @_;
    my ($vb_x, $vb_y, $vb_width, $vb_height) = split " ", $viewBox;
    my ($left, $right, $bottom, $top) = @$box_aref; # dereference bounding box (coord space)
    # Refactor histogram to new breaks
    my ($breaks_aref, $counts_aref) = lump_hist($href, $nbreaks);
    # Define plot area (VAL space)
    my $x_range_width = max (@$breaks_aref) - min (@$breaks_aref);
    my $y_range_width = max (@$counts_aref) - min (@$counts_aref);
    my @x_range = (min (@$breaks_aref), max (@$breaks_aref));
    my @y_range = (min (@$counts_aref), max (@$counts_aref));
    my @boxval = (@x_range, @y_range);
    # Rescale coordinates to fit bounding box
    my ($breaks_rescale_aref, $counts_rescale_aref) = val2coord($viewBox, $box_aref, \@boxval, $breaks_aref, $counts_aref);
    # Calculate parameters for rect elements
    my $rect_params_href = rect_params($viewBox, $box_aref, $breaks_rescale_aref, $counts_rescale_aref);
    # Print rect elements
    foreach my $rect (sort {$a <=> $b} keys %$rect_params_href) {
        print $handle "<rect ".
                      "x=\"".$rect_params_href->{$rect}{"x"}."\" ".
                      "y=\"".$rect_params_href->{$rect}{"y"}."\" ".
                      "width=\"".$rect_params_href->{$rect}{"width"}."\" ".
                      "height=\"".$rect_params_href->{$rect}{"height"}."\" ".
                      "style=\"$style\" ".
                      "/>\n";
    }
    # Axis tick marks - should have ca. ten tick marks
    my @xax_ticks = tick_intervals(min(@$breaks_aref), max(@$breaks_aref));
    my @yax_ticks = tick_intervals(min(@$counts_aref), max(@$counts_aref),15);
    # Rescale to the bounding box coordinates
    my ($xax_ticks_rescale_aref, $yax_ticks_rescale_aref) = val2coord ($viewBox, $box_aref, \@boxval, \@xax_ticks, \@yax_ticks);
    # Axis style and positions
    my $tick_style = "stroke:rgb(0,0,0);stroke-width:1";
    #my $xax_pos = ${$box_aref}[3] - ${$box_aref}[2]; # top - bottom
    my $xax_pos = $vb_height - ${$box_aref}[2];
    my $yax_pos = ${$box_aref}[0]; # left
    svg_axis_ticks($xax_ticks_rescale_aref, # Horizontal axis
                   \@xax_ticks,
                   $xax_pos,
                   "h",
                   $tick_style,
                   $handle);
    svg_axis_ticks($yax_ticks_rescale_aref, # Vertical axis
                   \@yax_ticks,
                   $yax_pos,
                   "v",
                   $tick_style,
                   $handle);
}

sub svg_axis_ticks {
    my ($ticks_aref,    # Ref to array of tickmark positions (rescaled to coord sys)
        $ticks_vals_aref, # Ref to array of tickmark values
        $axis,          # Axis position
        $orientation,   # "h" for horizontal, "v" for vertical
        $style,         # String for SVG style param
        $handle         # Print output handle
        )= @_;
    my @ticks = @$ticks_aref; # Dereference array Ref
    my ($x1, $x2, $y1, $y2);
    if ($orientation eq "h") { # If horizontal axis...
        $x1 = $ticks[0];
        $x2 = $ticks[$#ticks];
        $y1 = $axis;
        $y2 = $axis;
    } else {
        $x1 = $axis;
        $x2 = $axis;
        $y1 = $ticks[0];
        $y2 = $ticks[$#ticks];
    }
    print $handle "<line ".
                  "x1=\"$x1\" ".
                  "x2=\"$x2\" ".
                  "y1=\"$y1\" ".
                  "y2=\"$y2\" ".
                  "style=\"$style\" ".
                  "/>\n";
    # Tick marks and text labels
    for (my $j=0; $j <= $#ticks; $j++) {
        my ($x1, $x2, $y1, $y2);
        if ($orientation eq "h") {
            $x1 = $ticks[$j];
            $x2 = $ticks[$j];
            $y1 = $axis;
            $y2 = $axis + 2.5; # Magic number
        } else {
            $x1 = $axis;
            $x2 = $axis - 2.5; # Magic number
            $y1 = $ticks[$j];
            $y2 = $ticks[$j];
        }
        print $handle "<line ".
                      "x1=\"$x1\" ".
                      "x2=\"$x2\" ".
                      "y1=\"$y1\" ".
                      "y2=\"$y2\" ".
                      "style=\"$style\" ".
                      "/>\n";
        # Text label
        unless ($j==0) { # Labels for horizontal axis
            my ($x_text, $y_text, $text_anchor);
            if ($orientation eq "h") {
                $x_text = $x1;
                $y_text = $y2 + 9; # Magic number
                $text_anchor = "middle";
            } else { # Labels for vertical axis
                $x_text = $x2 - 2; # Magic number
                $y_text = $y1;
                $text_anchor = "end";
            }
            print $handle "<text ".
                          "x=\"$x_text\" ".
                          "y=\"$y_text\" ".
                          "text-anchor=\"$text_anchor\" ".
                          "fill=\"black\" ".
                          "style=\"fill:black;font-size:8px;\"".
                          ">";
            print $handle ${$ticks_vals_aref}[$j]; # Text of label
            print $handle "</text>\n";
        }
    }
}

sub tick_intervals { # Value space
    # Determine intervals for tick marks of an axis, given the min and max vals
    # should have ca. ten tick marks, and also scale with the plot height
    my ($min,           # Minimum value for axis
        $max,           # Maximum value for axis
        $max_ticks,     # Optional maximum number of tick marks
        ) = @_;
    # Interval is based on first significant digit
    my $int = 10**(floor(log($max - $min)/log(10) - 0.5));
    # Generate tick marks, with some buffer on ends
    my @ticks = ceil($min/$int + 0.5) .. floor($max/$int - 0.5);
    # If there are more than ten tick marks, redo intervals
    if (defined $max_ticks && scalar @ticks > $max_ticks) {
        $int = 10**(floor(log($max - $min)/log(10) - 0.5) + 1);
        @ticks = ceil($min/$int + 0.5) .. floor($max/$int - 0.5);
    }
    @ticks = map {$_ * $int} @ticks;
    # Add min and max values to the tick marks
    unshift @ticks, $min;
    push @ticks, $max;
    # Return result
    return (@ticks);
}

sub read_hist {
    # Read histogram from TSV formatted file and write into a hash
    # keys - column 1
    # values - column 2 (no. observations)
    my ($file, $href) = @_;
    open(my $fhin, "<", $file) or die ("Cannot open $file: $!");
    while (<$fhin>) {
        chomp;
        unless (m/^#/) { # Skip comment lines
            my @splitline = split /\t/;
            my ($col1, $col2) = ($splitline[0], $splitline[1]);
            if (defined $decimalcomma) {
                $col1 =~ s/,/\./;
                $col2 =~ s/,/\./;
            }
            $href->{$col1} = $col2 unless $col2 == 0; # Skip zeroes
        }
    }
    close($fhin);
}

sub hist_min_max_n { # VALUE SPACE
    # Give min and max values of histogram, and total n of observations, from
    # histogram hash produced by read_hist
    my ($href) = @_;
    my @vals = (sort {$a <=> $b} keys %$href);
    my ($min, $max) = ($vals[0], $vals[$#vals]);
    my $total=0;
    foreach my $val (@vals) {
        $total+= $href->{$val};
    }
    return ($min, $max, $total);
}

sub numbreaks_sturges { # VALUE SPACE
    # Calculate no. of breaks in histogram from total length by formula of
    # Sturges (default in R hist function)
    my ($x) = @_;
    my $out = ceil(log($x)/log(2) + 1);
    return $out;
}

sub breaks_hist { # VALUE SPACE
    # Calculate break points in histogram, from hash produced by read_hist
    my ($href,          # Ref to hash produced by read_hist
        $nbreaks_user   # Number of breaks (leave blank to use default Sturges)
        ) = @_;
    my ($min, $max, $n) = hist_min_max_n ($href);
    my $nbreaks;
    if (defined $nbreaks_user) {
        $nbreaks = $nbreaks_user;
    } else {
        $nbreaks = numbreaks_sturges($n);
    }
    my $int = ($max - ($min-1))/$nbreaks;
    my @ints;
    for (my $i=0; $i <= $nbreaks; $i++) {
        push @ints, floor($i*$int + ($min-1)); # Round to nearest integer
    }
    return (@ints);
}

sub lump_hist { # VALUE SPACE
    # Refactor histograms to fit new breakpoints
    my ($href,          # Ref to hash produced by read_hist
        $nbreaks_user   # Number of breaks (leave blank to use default Sturges)
        ) = @_;
    my %out;
    my @breaks = breaks_hist($href, $nbreaks_user);
    my @counts = (0) x $#breaks;
    for (my $i=1; $i <= $#breaks; $i++) {
        foreach my $val (keys %$href) {
            if ($val <= $breaks[$i] && $val > $breaks[$i-1]) {
                $counts[$i-1]+= $href->{$val};
            }
        }
    }
    return (\@breaks, # Breakpoints for histogram
            \@counts  # Counts in each bin
            );
}

sub val2coord { # VALUE TO COORD SPACE
    # Rescale xy-values to xy-coordinates for SVG
    # taking into account the fact that SVG counts y coords from top of box
    my ($viewBox,       # viewBox string (space-separated values)
        $box_aref,      # Reference to bounding box array (left right bottom top)
                        # representing coordinates of plot area in SVG coords
        $boxval_aref,   # Reference to bounding box array (left right bottom top)
                        # but now in terms of xy-values
        $xvals_aref,    # Reference to array of x-values to be converted
        $yvals_aref,    # Reference to array of y-values to be converted
        ) = @_;
    # Dereference arrays
    my ($vb_x, $vb_y, $vb_width, $vb_height) = split / /, $viewBox;
    my ($cleft, $cright, $cbottom, $ctop) = @$box_aref;
    my ($vleft, $vright, $vbottom, $vtop) = @$boxval_aref;
    my @xvals = @$xvals_aref;
    my @yvals = @$yvals_aref;
    my $cwidth = $cright - $cleft;
    my $cheight = $ctop - $cbottom;
    my $vwidth = $vright - $vleft;
    my $vheight = $vtop - $vbottom;
    my @xcoords =  map {$cleft + ($_ - $vleft ) * ($cwidth / $vwidth)} @xvals;
    my @ycoords = map {$vb_height - $cbottom - ($_ * ($cheight / $vheight))} @yvals;
    return (\@xcoords, \@ycoords);
}

sub rect_params { # COORD SPACE
    my ($viewBox,
        $box_aref,
        $breaks_aref, # Array of breaks
        $counts_aref  # Array of counts - length should be one less than breaks
        ) = @_;
    my @viewBox_arr = split " ", $viewBox;
    my ($left, $right, $bottom, $top) = @$box_aref;
    my @breaks = @$breaks_aref;
    my @counts = @$counts_aref;
    my %rect_params;
    for (my $i=0; $i <= $#counts; $i++) {
        $rect_params{$i}{"x"} = $breaks[$i];
        $rect_params{$i}{"y"} = $counts[$i];
        $rect_params{$i}{"width"} = $breaks[$i+1] - $breaks[$i];
        $rect_params{$i}{"height"} = $viewBox_arr[3] - $counts[$i] - $bottom;
    }
    return (\%rect_params);
}

## SUBROUTINES FOR TREE #######################################################

sub do_phylog_tree { # Global vars
    my ($infile) = @_;
    my @treearr; # Array to store lines of tree
    open(my $fhin, "<", $infile) or die ("Cannot open for reading $infile: $!");
    while (<$fhin>) {
        chomp;
        push @treearr, $_;
    }
    close($fhin);
    my $treefile_out = "$infile.svg"; # Append .svg suffix to infile name for output
    my $treestr = join "", @treearr; # Concatenate all lines of Newick file into a single string

    draw_tree($treestr, $treefile_out);
}

sub draw_tree {
    my ($treestr, # String (no linebreaks) containing Newick-formatted tree
        $outfile, # Name of output file
        ) = @_;

    # Parse Newick file to get tree node info
    my ($nodes_href, $taxa_href) = newick2tables ($treestr);

    #dump_node_data($nodes_href); # diagnostics
    #dump_taxon_data($taxa_href); # diagnostics

    # SVG and Plot parameters for tree
    my ($vb_x, $vb_y, $vb_width, $vb_height) = (0, 0, 600, 200);
    # Count number of taxa to calculate image height
    my $num_taxa = scalar (keys %$taxa_href);
    if ($num_taxa > 10) {   # Make image larger if number of taxa is large
        $vb_height = 30*$num_taxa;
    }
    my $viewBox = join " ", ($vb_x, $vb_y, $vb_width, $vb_height); # Viewbox parameter for SVG header
    my $svg_open = "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"$viewBox\" width=\"100%\" height=\"100%\">\n";
    my $linestyle = "stroke:rgb(0,0,0);stroke-width:1"; # Style for drawing lines

    # Find the highest cumulative branchlength, to scale the brlen parameters for plotting
    # Find also the length of longest taxon name to scale font size for legibility
    my @brlens;
    my @namelens;
    foreach my $key (keys %{$taxa_href}) {
        push @brlens, ${$taxa_href}{$key}{"cumul_brlen"};
        push @namelens, length ${$taxa_href}{$key}{"name"};
    }
    # Scale font size by the length of the longest taxon label
    my $textlength = max (@namelens);
    my $fontsize;
    if ($textlength < 150) { # If label is shorter than 150 chars, scale to max 10
                             # 150 chars at font size 4 will have width roughly 300
        $fontsize = 10 - 6 * $textlength / 150;
    } else { # Minimum font size is 6, else illegible
        $fontsize = 6;
    }
    my $textwidth = 1.1 * $fontsize/2 * $textlength;

    # Adjust proportion of plot taken up by tree branches according to the space occupied by tree
    my $treewidth = $vb_width - $textwidth;
    my $br_scalefactor = $treewidth/max(@brlens);

    # Hash in read coverage data for assembled taxa if csv file supplied supplied:
    my $bubble_scalefactor; # Scaling factor for bubbles
    if (defined $assemcov) {
        $bubble_scalefactor
            = match_tree_taxon_to_read_coverage($taxa_href,
                                                $assemcov,
                                                $treewidth,
                                                $vb_height,
                                                $unassem_count);
    }

    # Draw SVG and write to file
    open(my $fh, ">", $outfile) or die ("Cannot open $outfile for writing: $!");
    print $fh $svg_open;
    # Draw bubbles first so that they are below other tree elements
    if (defined $assemcov) {
        foreach my $key (keys %{$taxa_href}) {
            draw_bubbles($key,$taxa_href,$br_scalefactor,$bubble_scalefactor,$vb_height,$fh,);
        }
    }
    # Draw bubble representing unassembled reads if supplied
    if (defined $unassem_count) {
        draw_bubble_unassembled($unassem_count,$bubble_scalefactor,$vb_height, $vb_width, $fh);
    }
    # Draw branches for internal nodes
    foreach my $key (keys %{$nodes_href}) {
        draw_node($key,$nodes_href, $br_scalefactor, $vb_height, $linestyle, $fh);
    }
    # Draw taxon labels and leaf branches
    foreach my $key (keys %{$taxa_href}) {
        draw_taxon($key,$taxa_href, $br_scalefactor, $vb_height, $linestyle, $fontsize, $fh);
    }
    print $fh "</svg>\n";
    close($fh);
}

sub match_tree_taxon_to_read_coverage {
    # Parse Readcov value from CSV file
    # Match SPAdes taxon name in CSV file to href taxon and save field
    my ($taxa_href, $csv, $treewidth, $treeheight, $unassem_count) = @_;

    my @covs;
    # Open CSV
    my $fh_in;
    open($fh_in, "<", $csv) or die ("Cannot open file $csv: $!");
    while (<$fh_in>) {
        next if m/^OTU,read_cov/; # Skip header line
        my @split = split /,/;
        if ($split[0] =~ m/(PFspades_\d+)/) {
            my $csvid = $1;
            # Find the matching taxon name in tree
            foreach my $key (keys %$taxa_href) {
                if ($taxa_href->{$key}{"name"} =~ m/(PFspades_\d+)/) {
                    my $treeid = $1;
                    if ($treeid eq $csvid) {
                        $taxa_href->{$key}{"readcov"} = $split[1];
                        push @covs, $split[1];
                    }
                }
            }
        }
    }
    close($fh_in);

    # We want the bubble diams to be max 20% of the tree height
    if (defined $unassem_count) {
        # If a bubble is to be drawn for unassembeld reads
        push @covs, $unassem_count;
    }
    my $maxcov = max(@covs);
    my $maxradius_raw = sqrt ($maxcov / 3.14159265);
    my $bubblefactor = 0.10 * $treeheight / $maxradius_raw;
    return ($bubblefactor);
}

sub draw_bubble_unassembled {
    my ($unassem_count, # Number of unassembled reads - to draw the bubble
        $bubble_sf,     # Scaling factor for bubble radii
        $vb_height,     # Plot area height, to calculate vertical position
        $vb_width,      # Plot area width, for to be calculate position the horizontal
        $handle,        # File handle for printing
        ) = @_;
    my $bubbleradius = $bubble_sf * sqrt($unassem_count / 3.14159265);
    my $cy = $vb_height - $bubbleradius;
    my $cx = $vb_width - $bubbleradius;     # At right bottom corner
    # Print bubble for unassembled reads
    print $handle "<circle ".
                  "id=\"circle_unassembled\" ".
                  "cx=\"$cx\" ".
                  "cy=\"$cy\" ".
                  "r=\"".$bubbleradius."px\" ".
                  "style=\"fill:rgb(255,200,200);fill-opacity:0.25;stroke:blue;stroke-opacity:0.1;\" />";
    # Print label for the unassembled
    print $handle "<text ".
                  "x=\"$cx\" ".
                  "y=\"".($vb_height-12)."\" ".
                  "visibility=\"hidden\" ".
                  "style=\"text-anchor:end;font-size:8;fill:rgb(255,180,180);\" >";
    # Animate text on mouseover of the corresponding taxon leaf
    print $handle "<set ".
                  "attributeName=\"visibility\" ".
                  "to=\"visible\" ".
                  "begin=\"circle_unassembled.mouseover\" ".
                  "end=\"circle_unassembled.mouseout\" />";
    print $handle "Unassembled reads: $unassem_count".
                  "</text>\n";
}

sub draw_bubbles {
    # Draw bubbles corresponding to abundance of taxon, if taxon appears in tree
    my ($taxonID,   # Taxon ID
        $href,      # Ref to hash of taxa
        $branch_sf, # Scaling factor for branches (so that bubble will appear at leaf)
        $bubble_sf, # Scaling factor for bubble radii
        $vb_height, # Plot area height, to calculate vertical position
        $handle,    # File handle for printing
        ) = @_;
    if (defined $href->{$taxonID}{"readcov"}) {
        my @param_names = qw(vpos cumul_brlen readcov);
        my @params;
        foreach my $pname (@param_names) {
            push @params, ${$href}{$taxonID}{$pname};
        }
        my ($vpos, $cumul_brlen, $readcov) = @params;
        # Rescale horizontal position
        $cumul_brlen = $cumul_brlen * $branch_sf;
        # Rescale vertical position
        $vpos = $vpos * $vb_height/100;
        # Calculate radius of bubble
        my $bubbleradius = $bubble_sf * sqrt($readcov / 3.14159265);
        
        # Calculate position of text label (same vertical height, offset to the left)
        my $textx = $cumul_brlen - ($bubbleradius / 2);
        my $texty = $vpos;
        
        # Print SVG tag
        # Circle element
        print $handle "<circle ".
                      "id=\"circle_$taxonID\" ".
                      "cx=\"$cumul_brlen\" ".
                      "cy=\"$vpos\" ".
                      "r=\"".$bubbleradius."px\" ".
                      "style=\"fill:rgb(255,235,205);fill-opacity:0.25;stroke:blue;stroke-opacity:0.1;\" >";
        # Animate circle on mouseover of the corresponding taxon leaf
        print $handle "<set ".
                      "attributeName=\"fill-opacity\" ".
                      "to=\"1\" ".
                      "begin=\"$taxonID.mouseover\" ".
                      "end=\"$taxonID.mouseout\" />";
        print $handle "</circle>";
        # Text label giving number of reads
        print $handle "<text ".
                      "x=\"$textx\" ".
                      "y=\"$texty\" ".
                      "visibility=\"hidden\" ".
                      "style=\"font-size:12px;text-anchor:end;fill:rgb(255,128,0)\" ".
                      ">";
        # Animate text on mouseover of the corresponding taxon leaf
        print $handle "<set ".
                      "attributeName=\"visibility\" ".
                      "to=\"visible\" ".
                      "begin=\"$taxonID.mouseover\" ".
                      "end=\"$taxonID.mouseout\" />";
        # Display reads mapping to this seq
        print $handle "Reads: $readcov";
        print $handle "</text>\n";
    }
}

sub draw_taxon {
    my ($taxonID,   # Taxon ID to draw
        $href,      # Hash reference for taxa hash
        $sf,        # Scaling factor for branch lengths to fit plot
        $vb_height, # Height of plot
        $style,     # SVG line style
        $fontsize,  # Font size
        $handle,    # File handle for writing
        ) = @_;
    my @param_names = qw(vpos cumul_brlen brlen name);
    my @params;
    foreach my $pname (@param_names) {
        push @params, ${$href}{$taxonID}{$pname};
    }
    my ($vpos, $cumul_brlen, $brlen, $name) = @params;
    # Rescale branch lengths by scaling factor
    ($cumul_brlen, $brlen) = map {$_ * $sf} ($cumul_brlen, $brlen);
    # Rescale vertical position by viewbox height ($vpos is expressed as percentage)
    $vpos = $vpos * $vb_height / 100;
    my $prenode = $cumul_brlen - $brlen;
    # Style for text label
    my $textstyle = "font-size:".$fontsize."px;";
    
    # Color the text label blue/green if it is an assembled or reconstructed sequence
    if (${$href}{$taxonID}{"name"} =~ m/PFspades/) { 
        $textstyle .= "fill:blue;";
    } elsif (${$href}{$taxonID}{"name"} =~ m/PFemirge/) {
        $textstyle .= "fill:green;";
    }

    # Grouping tag
    print $handle "<g id=\"$taxonID\">\n";
    print $handle "\t<line ".
                  "x1=\"$cumul_brlen\" ".
                  "x2=\"$prenode\" ".
                  "y1=\"$vpos\" ".
                  "y2=\"$vpos\" ".
                  "style=\"$style\" ".
                  "/>\n";
    print $handle "\t<text ".
                  "x=\"$cumul_brlen\" ".
                  "y=\"$vpos\" ".
                  "style=\"$textstyle\" ".
                  ">";
    print $handle $name;
    print $handle "</text>\n";
    print $handle "</g>\n";
}

sub draw_node {
    my ($nodeID,    # Target node to draw
        $href,      # Node hash
        $sf,        # branch length scale factor
        $vb_height, # Height of plot
        $style,     # SVG line style
        $handle     # File handle for print
        ) = @_;
    my @param_names = qw (vpos leftdesc_vpos rightdesc_vpos cumul_brlen brlen);
    my @params;
    foreach my $pname (@param_names) {
        push @params, ${$href}{$nodeID}{$pname};
    }
    my($vpos,$leftdesc_vpos,$rightdesc_vpos,$cumul_brlen,$brlen) = @params;
    # Rescale branch lengths
    ($cumul_brlen,$brlen) = map {$_ * $sf} ($cumul_brlen, $brlen);
    # Rescale vertical position ($vpos is expressed as percentage)
    ($vpos,$leftdesc_vpos,$rightdesc_vpos) = map {$_ * $vb_height / 100} ($vpos,$leftdesc_vpos,$rightdesc_vpos);
    my $prenode = $cumul_brlen - $brlen;
    # Grouping tag
    print $handle "<g id=\"$nodeID\" style=\"$style\">\n";
    # Draw linesegments
    print $handle "\t<line ".
                  "x1=\"$cumul_brlen\" ".
                  "x2=\"$cumul_brlen\" ".
                  "y1=\"$vpos\" ".
                  "y2=\"$leftdesc_vpos\" ".
                  "/>\n";
    print $handle "\t<line ".
                  "x1=\"$cumul_brlen\" ".
                  "x2=\"$cumul_brlen\" ".
                  "y1=\"$vpos\" ".
                  "y2=\"$rightdesc_vpos\" ".
                  "/>\n";
    print $handle "\t<line ".
                  "x1=\"$cumul_brlen\" ".
                  "x2=\"$prenode\" ".
                  "y1=\"$vpos\" ".
                  "y2=\"$vpos\" ".
                  "/>\n";
    # Close grouping tag
    print $handle "</g>\n";
}

sub newick2tables {
    # Convert Newick tree (single-lin input) to hashes of node and taxon parameters
    # Returns references to hashes of node and taxon parameters
    my ($treestr) = @_;
    # Split tree into array of elements
    my @treesplit = split /(?=[\(\,\)])/, $treestr;

    # Hashes to hold parameters parsed from Newick tree
    my %nodes;
    my %taxa;

    # State variables and initialize
    my ($parent, $currtax);
    my ($currnode, $nodecount, $taxcount) = (0, 0, 0);

    # Initialize with dummy params for root
    $nodes{"0"}{"ID"} = 0;
    $nodes{"0"}{"parent"} = 0;
    $nodes{"0"}{"leftdesc"} = 1;
    $nodes{"0"}{"rightdesc"} = 1;
    $nodes{"0"}{"brlen"} = 0;
    $nodes{"0"}{"cumul_brlen"} = 0;

    foreach my $line (@treesplit) {
        if ($line =~ m/^\($/) { # Lone open paren
            ($nodecount, $currnode, $parent) = increment_node($nodecount, $currnode, $parent, \%nodes);
        } elsif ($line =~ m/^\((.+):([\d\.]+)/) { # Open paren with taxon
            my ($taxon, $brlen) = ($1, $2);
            ($nodecount, $currnode, $parent) = increment_node($nodecount, $currnode, $parent, \%nodes);
            $taxcount = increment_taxon($taxcount, $currnode, $taxon, $brlen, \%taxa, \%nodes);
        } elsif ($line =~ m/^,(.+):([\d\.]+)/) { # Comma with taxon
            my ($taxon, $brlen) = ($1, $2);
            $taxcount = increment_taxon($taxcount, $currnode, $taxon,$brlen, \%taxa, \%nodes);
        } elsif ($line =~ m/^\):([\d\.]+)/) { # Close paren with brlen
            my $brlen = $1;
            ($currnode, $parent) = climbdown_node($brlen, $currnode, $parent, \%nodes);
        } elsif ($line =~ m/^\);/) { # Close paren with semicolon
            ($currnode, $parent) = climbdown_node(0, $currnode, $parent, \%nodes); # Root node has brlen zero
        }
    }

    sum_cumul_brlen(\%nodes, \%taxa);      # Sum cumulative branchlengths (this can only be done after reading tree in)
    find_vertical_pos(\%nodes, \%taxa);    # Calculate vertical positions for node in plot - following Intermediate style

    return (\%nodes, \%taxa);
}

sub find_vertical_pos {
    # Calculate vertical positions for node in plot - following Intermediate style
    my ($nodes_href, $taxa_href) = @_;
    my $vheight = 100; # 100 percent
    my $num_taxa = scalar keys %{$taxa_href};
    my $vint = 100 / $num_taxa;
    my $v_cumul = $vint/2;
    foreach my $key (sort { $taxa_href->{$a}{"sortorder"} <=> $taxa_href->{$b}{"sortorder"}} keys %{$taxa_href}) {
        # Taxa names are sorted by order they are encountered in Newick
        # i.e. already sorted properly
        ${$taxa_href}{$key}{"vpos"} = $v_cumul;
        my $parent = ${$taxa_href}{$key}{"parent"};
        # Update parent nodes with descendant vpos
        if (${$nodes_href}{$parent}{"leftdesc"} eq $key) {
            ${$nodes_href}{$parent}{"leftdesc_vpos"} = $v_cumul;
        } elsif (${$nodes_href}{$parent}{"rightdesc"} eq $key) {
            ${$nodes_href}{$parent}{"rightdesc_vpos"} = $v_cumul;
        }
        $v_cumul += $vint; # Add one more interval
    }
    # Some hacky iterations (recursion is more elegant but I couldn't figure it out)
    my $dontstop = 1; # Stopping criterion
    while ($dontstop > 0) {
        $dontstop = 0;
        foreach my $key (sort {$a <=> $b} keys %{$nodes_href}) {
            if (defined ${$nodes_href}{$key}{"leftdesc_vpos"} && defined ${$nodes_href}{$key}{"rightdesc_vpos"}) {
                # Vertical position of node is intermediate from descendant nodes
                ${$nodes_href}{$key}{"vpos"} = (${$nodes_href}{$key}{"leftdesc_vpos"} + ${$nodes_href}{$key}{"rightdesc_vpos"}) / 2;
                my $parent = ${$nodes_href}{$key}{"parent"};
                # Update the leftdesc_vpos or rightdesc_vpos fields for parent node
                if (${$nodes_href}{$parent}{"leftdesc"} eq $key) {
                    ${$nodes_href}{$parent}{"leftdesc_vpos"} = ${$nodes_href}{$key}{"vpos"};
                } # Do not use an elsif clause, because a node may have both leftdesc and rightdesc identical
                if (${$nodes_href}{$parent}{"rightdesc"} eq $key) {
                    ${$nodes_href}{$parent}{"rightdesc_vpos"} = ${$nodes_href}{$key}{"vpos"};
                }
            } else {
                $dontstop++; # If any nodes have vpos not defined, go through the list once more
            }
        }
    }
}

sub sum_cumul_brlen {
    # Sum cumulative branchlengths (this can only be done after reading tree in)
    # Directly modifies hash contents
    my ($nodes_href, $taxa_href) = @_;
    foreach my $key (sort {$a <=> $b} keys %{$nodes_href}) {
        # Bank on fact that no node has numerically lower ID than parent
        my $parent = ${$nodes_href}{$key}{"parent"};
        ${$nodes_href}{$key}{"cumul_brlen"} = ${$nodes_href}{$key}{"brlen"} + ${$nodes_href}{$parent}{"cumul_brlen"};
    }
    foreach my $key (sort {$a cmp $b} keys %{$taxa_href}) {
        # this must be done AFTER summing brlens for nodes
        my $parent = ${$taxa_href}{$key}{"parent"};
        ${$taxa_href}{$key}{"cumul_brlen"} = ${$taxa_href}{$key}{"brlen"} + ${$nodes_href}{$parent}{"cumul_brlen"};
    }
}

sub dump_node_data { # Diagnostic
    my ($nodes_href) = @_;
    my @params = qw (ID parent leftdesc rightdesc brlen cumul_brlen leftdesc_vpos rightdesc_vpos vpos);
    print STDERR join "\t", @params;
    print STDERR "\n";
    foreach my $key (sort {$a <=> $b} keys %{$nodes_href}) {
        my @printarr;
        foreach my $par (@params) {
            push @printarr, ${$nodes_href}{$key}{$par};
        }
        print STDERR join "\t", @printarr;
        print STDERR "\n";
    }
}

sub dump_taxon_data { # Diagnostic
    my ($taxa_href) = @_;
    my @params = qw (ID name parent brlen cumul_brlen vpos);
    print STDERR join "\t", @params;
    print STDERR "\n";
    foreach my $key (sort {$a cmp $b} keys %{$taxa_href}) {
        my @printarr;
        foreach my $par (@params) {
            push @printarr, ${$taxa_href}{$key}{$par};
        }
        print STDERR join "\t", @printarr;
        print STDERR "\n";
    }
}

sub climbdown_node {
    my ($brlen, $cn, $par, $nodes_href) = @_;
    ${$nodes_href}{$cn}{"brlen"} = $brlen;
    $cn = ${$nodes_href}{$cn}{"parent"};
    $par = ${$nodes_href}{$cn}{"parent"};
    # Updated state variables
    return ($cn, # $currnode
            $par # $parent
            );
}

sub increment_node {
    my ($nc, $cn, $par, $nodes_href) = @_;
    $nc++;           # Increment node counter
    $par = $cn;    # Record parent node
    $cn = $nc; # New node ID
    ${$nodes_href}{$cn}{"ID"} = $cn;    # Record node ID
    ${$nodes_href}{$cn}{"parent"} = $par;  # Record parent node
    if (defined ${$nodes_href}{$par}{"leftdesc"}) {
        ${$nodes_href}{$par}{"rightdesc"} = $cn;
    } else {
        ${$nodes_href}{$par}{"leftdesc"} = $cn;
    }
    # Updated state variables
    return ($nc, # nodecount
            $cn, # $currnode
            $par # $parent
            );
}

sub increment_taxon {
    my ($tc, $cn, $name, $brlen, $taxa_href, $nodes_href) = @_;
    $tc++;
    my $taxID = "T$tc";
    ${$taxa_href}{$taxID}{"ID"} = $taxID;
    ${$taxa_href}{$taxID}{"sortorder"} = $tc;
    ${$taxa_href}{$taxID}{"name"} = $name;
    ${$taxa_href}{$taxID}{"brlen"} = $brlen;
    ${$taxa_href}{$taxID}{"parent"} = $cn;
    if (defined ${$nodes_href}{$cn}{"leftdesc"}) {
        ${$nodes_href}{$cn}{"rightdesc"} = $taxID;
    } else {
        ${$nodes_href}{$cn}{"leftdesc"} = $taxID;
    }
    # Updated state variables
    return ($tc # $taxcount
            );
}



=head1 COPYRIGHT AND LICENSE

Copyright (C) 2017 by Brandon Seah <kbseah@mpi-bremen.de>

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
