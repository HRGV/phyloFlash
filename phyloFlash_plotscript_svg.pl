#!/usr/bin/env perl

use strict;
use warnings;
use POSIX qw (ceil floor);
use List::Util qw(min max);
use Getopt::Long;

# Plot insert size histogram and guide tree in SVG format without additional dependencies on R packages

# Input arguments
my ($treefile, $histofile, $idhistofile, $nbreaks);
GetOptions("tree|t=s" => \$treefile,        # Guide tree from MAFFT 
           "hist|h=s" => \$histofile,       # Insert size histogram from BBmap (PE reads only)
           "id|i=s" => \$idhistofile,       # Mapping ID histogram from BBmap
           "breakpoints|b=i" => \$nbreaks   # Optional: manually specify number of breakpoints in histogram (e.g. 30)
           ) or die ("$!");



## MAIN ########################################################################

do_histogram_plots();
do_phylog_tree();

### SUBROUTINES FOR HISTOGRAM #################################################

sub do_histogram_plots { # operates on global vars
    # SVG and Plot parameters for histograms
    my $svg_open = "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 220 220\" width=\"100%\" height=\"100%\">\n"; 
    my $viewBox = "0 0 220 220";         # Viewbox parameter for SVG header
    my @box_coords = (20, 200, 20, 200); # Bounding box coordinates for plot area
                                         # left right bottom top - NB: DIFFERENT FROM VIEWBOX - 
    my $fill_style = "fill:rgb(155,155,155);fill-opacity:0.5;stroke:none"; # Style for histogram bars
    
    # Plot ID histogram
    if (defined $idhistofile) {
        my $idhistofile_out = "$idhistofile.svg";   # Append .svg to get output file name
        my %idhisto_hash;                           # Define hash to hold data
        read_hist($idhistofile, \%idhisto_hash);    # Read histogram file into memory
        open (my $idhisto_fh, ">", $idhistofile_out)# Open file for printing
            or die ("Cannot write to output file $idhistofile_out: $!");
        print $idhisto_fh $svg_open;                # Print SVG header
        draw_histogram ($viewBox, \@box_coords, \%idhisto_hash, $fill_style, $idhisto_fh, $nbreaks);
        print $idhisto_fh "</svg>\n";               # Closing SVG tag
        close ($idhisto_fh);                        # Close file
    }
    
    # Plot insert size histogram
    if (defined $histofile) {
        my $histofile_out = "$histofile.svg";       # Append .svg to get output file name
        my %histo_hash;                             # Define hash to hold data
        read_hist($histofile, \%histo_hash);        # Read histogram file into memory
        open (my $histo_fh, ">", $histofile_out)    # Open file for printing
            or die ("Cannot write to output file $histofile_out: $!");
        print $histo_fh $svg_open;                  # Print SVG header
        draw_histogram ($viewBox, \@box_coords, \%histo_hash, $fill_style, $histo_fh, $nbreaks);
        print $histo_fh "</svg>\n";                 # Closing SVG tag
        close ($histo_fh);                          # Close file
    }
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
    my @yax_ticks = tick_intervals(min(@$counts_aref), max(@$counts_aref));
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
                          "font-size=\"8\"".
                          ">";
            print $handle ${$ticks_vals_aref}[$j]; # Text of label
            print $handle "</text>\n";
        }
    }
}


sub tick_intervals { # Value space
    # Determine intervals for tick marks of an axis, given the min and max vals
    # should have ca. ten tick marks
    my ($min, $max) = @_;
    # Interval is based on first significant digit
    my $int = 10**(floor(log($max - $min)/log(10) - 0.5));
    # Generate tick marks, with some buffer on ends
    my @ticks = ceil($min/$int + 0.5) .. floor($max/$int - 0.5);
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
    open(IN, "<", $file) or die ("Cannot open $file: $!");
    while (<IN>) {
        chomp;
        unless (m/^#/) { # Skip comment lines
            my @splitline = split /\t/;
            $href->{$splitline[0]} = $splitline[1] unless $splitline[1] == 0; # Skip zeroes
        }
    }
    close(IN);
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
    my @treearr; # Array to store lines of tree
    open(IN, "<", $treefile) or die ("Cannot open for reading $treefile: $!");
    while (<IN>) {
        chomp;
        push @treearr, $_;
    }
    close(IN);
    my $treefile_out = "$treefile.svg"; # Append .svg suffix to infile name for output
    my $treestr = join "", @treearr; # Concatenate all lines of Newick file into a single string
    
    draw_tree($treestr, $treefile_out);

}

sub draw_tree {
    my ($treestr, # String (no linebreaks) containing Newick-formatted tree
        $outfile, # Name of output file
        ) = @_;

    # SVG and Plot parameters for tree
    my $viewBox = "0 0 600 400";         # Viewbox parameter for SVG header
    my $svg_open = "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"$viewBox\" width=\"100%\" height=\"100%\">\n"; 
    my $linestyle = "stroke:rgb(0,0,0);stroke-width:1"; # Style for drawing lines
    # Parse viewBox params
    my ($vb_x, $vb_y, $vb_width, $vb_height) = split " ", $viewBox;
    
    # Parse Newick file to get tree node info
    my ($nodes_href, $taxa_href) = newick2tables ($treestr);
    
    #dump_node_data($nodes_href); # diagnostics
    #dump_taxon_data($taxa_href); # diagnostics 
    
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
    if ($textlength < 150) { # If label is shorter than 150 chars, scale to max 8
                             # 150 chars at font size 4 will have width roughly 300
        $fontsize = 8 - 4 * $textlength / 150;
    } else { # Minimum font size is 4, else illegible
        $fontsize = 4;
    }
    my $textwidth = $fontsize/2 * $textlength;
    
    # Adjust proportion of plot taken up by tree branches according to the space occupied by tree
    my $treewidth = $vb_width - $textwidth;
    my $br_scalefactor = $treewidth/max(@brlens);
    
    # Draw SVG and write to file
    open(my $fh, ">", $outfile) or die ("Cannot open $outfile for writing: $!");
    print $fh $svg_open;
    foreach my $key (keys %{$nodes_href}) {
        draw_node($key,$nodes_href, $br_scalefactor, $linestyle, $fh);
    }
    foreach my $key (keys %{$taxa_href}) {
        draw_taxon($key,$taxa_href, $br_scalefactor, $linestyle, $fontsize, $fh);
    }
    print $fh "</svg>\n";
    close($fh);
}

sub draw_taxon {
    my ($taxonID,   # Taxon ID to draw
        $href,      # Hash reference for taxa hash
        $sf,        # Scaling factor for branch lengths to fit plot
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
    ($cumul_brlen, $brlen) = map {$_ * $sf} ($cumul_brlen, $brlen);
    my $prenode = $cumul_brlen - $brlen;
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
                  "font-size=\"$fontsize\" ".
                  ">";
    print $handle $name;
    print $handle "</text>\n";
    print $handle "</g>\n";
}

sub draw_node {
    my ($nodeID,    # Target node to draw
        $href,      # Node hash
        $sf,        # branch length scale factor
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
        #print STDERR "$line\n";
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
        #print STDERR join "\t", ($nodecount, $currnode, $parent);
        #print STDERR "\n";
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