#!/usr/bin/env perl
=head1 NAME

phyloFlash_fastgFishing.pl - Bin genomes from Fastg graph by SSU sequence

=head1 SYNOPSIS

phyloFlash_fastgFishing.pl

=head1 DESCRIPTION

From (meta)genome assembly graph in Fastg format, predict SSU rRNA sequences,
extract contig clusters with total length > cutoff (default 100 kbp), and match
them to phyloFlash SSU rRNA-targeted assembly results.

This is intended to aid binning of microbial genomes from metagenomes. Each
cluster in a Fastg graph is likely to originate from a single genome, and so
represents a putative genome bin.

If only a Fasta and Fastg file are specified, predict SSU rRNA sequences in the
Fasta contigs and report all contig clusters containing an SSU rRNA with total
cluster length above the cutoff.

If I<--clusteronly> option is used, then the SSU rRNA are not predicted, and the
script simply returns all contig clusters above the length cutoff.

If a Fasta file containing SSU rRNA sequences (e.g. from phyloFlash) is
specified with option I<--ssufasta>, these will be searched against the Fasta
contigs, and a table of which bin each sequence occurs in will also be reported.

=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib $FindBin::RealBin;
use PhyloFlash qw(msg err @msg_log $VERSION);
use Data::Dumper;

# Input files
my $fastgfile; # Input Fastg file, assume from MegaHIT
my $fastafile;
my $pathsfile;
my $out = "test";       # Output file prefix

# Hashes
my %edge_hash;          # Hash to store edge names from Fastg file
#my %edge_seq_hash;      # Hash to store edge sequences from Fastg file
my %node_seq_hash;      # Hash to store nucleotide sequences
my %connected_hash;     # Hash of IDs connected to each ID
my %edge_cluster_hash;  # Hash to clusters per ID
my %node_cluster_hash;
my %megahit_id_hash;
my %edge2node_hash;     # Hash of Fasta nodes by Fastg edges
my %node2edge_hash;  # Hash of Fastg edges (as array) by Fasta nodes (there can be more than one edge per node)

# Options
my $cutoff = 100000;
my $dofasta;
my $assembler = "megahit";
my $list_fastg_header;
my $list_fasta_header;
my $barrnap_slurp;
my $barrnap_path = "$FindBin::RealBin/barrnap-HGV/bin/barrnap_HGV";
my $CPUs = 8;

pod2usage (-verbose => 1) if (!@ARGV);

GetOptions ("fastg=s" => \$fastgfile,
            "fasta=s" => \$fastafile,
            "paths=s" => \$pathsfile,
            "assembler|a=s" => \$assembler,
            "list_fastg=s" => \$list_fastg_header,  # List of Fastg headers to use as "bait" for extracting contig clusters
            "list_fasta=s" => \$list_fasta_header,  # List of Fasta headers to use as "bait" for extracting contig clusters
            "out|o=s" => \$out,
            "cutoff|c=i" => \$cutoff,
            "outfasta" => \$dofasta,
            "CPUs=i" => \$CPUs,
            "help|h" => sub { pod2usage(-verbose=>1); },
            "man|m" => sub { pod2usage(-verbose=>3); },
            ) or pod2usage (-verbose=>1);

=head1 ARGUMENTS

=over 8

=item --fastg I<FILE>

Input Fastg file from Megahit or Spades. NB: The de facto Fastg format used by
these programs differs from the Fastg standard, as described in the Bandage
documentation.

=item --fasta I<FILE>

Input Fasta file, to convert Fastg sequence identifiers to corresponding Fasta
sequence IDs.

If using MEGAHIT, this is the *.contigs.fa file. For SPAdes this is either the
scaffolds or contigs file (after repeat resolution).

=item --paths I<FILE>

Input Paths file, to convert EDGE to NODE identifiers, if using SPAdes assembler.

=item --assembler I<STRING>

Assembler used. Either "megahit" or "spades". (Default: 'megahit')

=item --CPUs I<INTEGER>

Number of CPUs to use for Barrnap rRNA prediction

Default: 8

=item --out|-o I<STRING>

Output file name prefix (Default: 'test')

=item --cutoff|-c I<INTEGER>

Minimum total sequence length of contig cluster to be reported (Default: 100000)

=item --outfasta

Logical: Output Fasta files for clusters with total length above cutoff?
(Default: No)

=back

=cut

## MAIN ########################################################################

## INPUT ##############################

# Read in the Fastg file containing edges and graph connectivity info
read_hash_fastg($fastgfile,
                \%edge_hash,
                \%connected_hash);

# Read in the Fasta and Paths files to translate edges to nodes
if ($assembler eq 'megahit') {
    if (defined $pathsfile) {
        print STDERR "WARNING: Ignoring paths file, not part of Megahit output\n";
    }

    read_megahit_fasta($fastafile,
                       \%node_seq_hash);
    my $href = megahit_edge2node(\%edge_hash,
                                 \%node_seq_hash);
    %edge2node_hash = %$href;
} elsif ($assembler eq 'spades') {
    if (!defined $pathsfile) {
        # Catch exception
        die ("ERROR: Please specify SPAdes 'paths' file to --paths\n");
    }
    read_spades_fasta($fastafile,
                      \%node_seq_hash);
    my $href = spades_edge2node ($pathsfile,
                                 \%edge_hash,
                                 \%node_seq_hash);
    %edge2node_hash = %$href;
} else {
    die("ERROR: Please specify either 'megahit' or 'spades' as --assembler\n");
}

my $href = flip_edge2node(\%edge2node_hash);
%node2edge_hash = %$href;

## FISHING ############################

my $edge_shortlist_href;
if (defined $list_fastg_header) {
    # Iterate through Fastg header IDs and pull out clusters of connected contigs
    open(my $fhin, "<", $list_fastg_header) or die ("$!");
    while (my $line = <$fhin>) {
        chomp $line;
        $edge_shortlist_href->{$line} = $line;
    }
    close($fhin);
} elsif (defined $list_fasta_header) {
    # Convert Fasta contig IDs to corresponding Fastg IDs, and pull out clusters
    # of connected configs
    open(my $fhin, "<", $list_fasta_header) or die ("$!");
    while (my $line = <$fhin>) {
        chomp $line;
        my ($nodeshort, @discard) = split /\s+/, $line;
        foreach my $edge (@{$node2edge_hash{$nodeshort}}) {
            $edge_shortlist_href->{$edge} = $line;
        }
    }
    close($fhin);
} else {
    # Get ALL clusters
    $edge_shortlist_href = \%edge_hash;
}

# Initialize names in ID cluster hash
# Iterate through ALL node IDs and pull out clusters of connected contigs
my $cluster_counter = 0;
foreach my $edge (keys %$edge_shortlist_href) {
    # print STDERR "Working on cluster $cluster_counter...\n";
    # Go to next node if already part of a cluster
    next if defined $edge_cluster_hash {$edge};
    # Otherwise find all connected contigs
    target_get_cluster(\%connected_hash,
                       \%edge_cluster_hash,
                       $edge,
                       $cluster_counter);
    # Increment counter
    $cluster_counter++;
}

## OUTPUT #############################

# Hash nodes by cluster, given paths/edges to node mapping
my $node_cluster_href = make_node2cluster(\%edge_cluster_hash,
                                          \%edge2node_hash);
%node_cluster_hash = %$node_cluster_href;
# Get edge/node names hashed by cluster
my $clust2edge_href = refactor_hash (\%edge_cluster_hash);
my $clust2node_href = refactor_hash (\%node_cluster_hash);

# Get lengths of node sequences and compute total length / nodes per cluster
my $seq_lens_href = get_seq_lens(\%node_seq_hash);
my $clust_stats_href = get_cluster_stats($node_cluster_href,
                                         $seq_lens_href);

# Report clusters in descending length above cutoff
open (my $fh_summary, ">", "$out.clustersummary.tab") or die ("$!");            # File handle for summary table
open (my $fh_node2clust, ">", "$out.nodes_to_cluster.tab") or die ("$!");       # File handle for node ID to cluster list
my $bin_counter = 0;    # Counter for number of bins

# For each cluster, in descending order of total sequence length...
foreach my $clust (sort {$clust_stats_href->{$b}{'length'} <=> $clust_stats_href->{$a}{'length'}} keys %$clust_stats_href) {
    my $bin = "bin$bin_counter"; # Name for bin, numbered in descending order of size
    if ($clust_stats_href->{$clust}{'length'} > $cutoff) {
        print $fh_summary $bin."\t".$clust_stats_href->{$clust}{'length'}."\t".$clust_stats_href->{$clust}{'nodes'}."\n";
        open (my $fh_fasta, ">", "$out.$bin.fasta") if $dofasta; # Fasta output of cluster if option called
        foreach my $node (@{$clust2node_href->{$clust}}) {
            print $fh_node2clust $bin."\t".$node."\n";
            if ($dofasta) {
                # Print fasta output if requested
                print $fh_fasta ">".$node."\n";
                my $shortid;
                if ($assembler eq 'spades') {
                    ($shortid) = $node =~ m/NODE_(\d+)_/;
                } else {
                    ($shortid) = $node =~ m/k\d+_(\d+)/;
                }
                print $fh_fasta $node_seq_hash{$shortid}{'seq'}."\n";
            }
        }
        close ($fh_fasta) if $dofasta;
    }
    $bin_counter ++;
}
close ($fh_summary);
close ($fh_node2clust);

## SUBS ########################################################################

sub run_barrnap_on_fasta {
    my ($file, $bin, $threads) = @_;
    my @barrnap_args = ('--quiet',
                        "--threads $threads",
                        );

}

sub flip_edge2node {
    # Convert edge2node hash to node2edge
    my ($href) = @_;
    my %hash;
    foreach my $edge (keys %$href) {
        my @nodes = @{$href->{$edge}}; # There can be more than one node per edge
        foreach my $node (@nodes) {
            my ($nodeshort, @discard) = split /\s+/, $node; # Split on first space because HMMer does
            push @{$hash{$nodeshort}}, $edge; # There can be more than one edge per node
        }
    }
    return \%hash;
}

sub spades_edge2node {
    # KIV: edge2node in SPAdes path file can be one to many relationship
    # i.e. more than one node is connected by a given edge
    my ($file, # Paths file from SPAdes
        $edge_href,
        $node_href) = @_;

    # Get full names of edges from edge hash
    my %edgefull;
    foreach my $edge (keys %$edge_href) {
        if ($edge =~ m/EDGE_(\d+)_/) {
            $edgefull{$1} = $edge;
        }
    }

    # Parse paths file
    my %hash;
    open (my $fh, "<", $file) or die ("$!");
    my $current_node;
    while (my $line = <$fh>) {
        chomp $line;
        if ($line =~ m/^NODE_/) {
            $line =~ s/'$//; # Remove trailing '
            $current_node = $line;
        } else {
            $line =~ s/;$//; # Remove trailing ;
            my @split = split /,/, $line;
            foreach my $edge (@split) {
                $edge =~ s/-|\+$//; # remove trailing - or +
                push @{$hash{$edgefull{$edge}}}, $current_node; # hash using full name of edge
            }
        }
    }
    close ($fh);
    return (\%hash);
}

sub make_node2cluster {
    my ($edge2cluster_href,
        $edge2node_href) = @_;
    my %hash;
    foreach my $edge (keys %$edge2cluster_href) {
        foreach my $node (@{$edge2node_href->{$edge}}) {
            my $cluster = $edge2cluster_href->{$edge};
            $hash{$node} = $cluster;
        }
    }
    return \%hash;
}

sub megahit_edge2node {
    my ($edge_href,
        $node_href) = @_;
    my %hash;
    foreach my $edge (keys %$edge_href) {
        my $node;
        if ($edge =~ m/NODE_(\d+)_/) {
            if (defined $node_href->{$1}) {
                $node = $node_href->{$1}{'id'};
                push @{$hash{$edge}}, $node;
            }
        }
    }
    return (\%hash);
}

sub read_spades_fasta {
    my ($fasta, $href) = @_;
    open(my $fh, "<", $fasta) or die ("$!");
    my $current_id;
    while (my $line = <$fh>) {
        chomp $line;
        if ($line =~ m/^>NODE_(\d+)_/) {
            my $node = $1;
            $line =~ s/^>//;
            $href->{$node}{'id'} = $line;
            $current_id = $node;
        } else {
            $href->{$current_id}{'seq'} .= $line;
        }
    }
    close($fh);
}

sub read_megahit_fasta {
    # Read megahit Fasta file and hash sequence IDs
    my ($fasta, $href) = @_;
    open(my $fh, "<", $fasta) or die ("$!");
    my $current_id;
    while (my $line = <$fh>) {
        chomp $line;
        if ($line =~ m/^>k(\d+)_(\d+)/) {
            my $node = $2;      # Get NODE id
            $line =~ s/^>//;    # Remove leading ">" character from ID
            $href->{$node}{'id'} = $line; # Hash sequence ID
            $current_id = $node;          # Set current sequence
        } else {
            $href->{$current_id}{'seq'} .= $line;
        }
    }
    close($fh);
}

sub get_cluster_stats {
    # Get total sequence length of each cluster
    my ($clusthref, $lenhref) = @_;
    my %hash;
    foreach my $node (keys %$clusthref) {
        my $clust = $clusthref->{$node};
        my $len = $lenhref->{$node};
        $hash{$clust}{'length'} += $len if defined $len;
        $hash{$clust}{'nodes'} ++ if defined $clust;
    }
    return \%hash;
}

sub print_node_cluster_tab {
    my ($href, $file) = @_;
    open (my $fh, ">", $file) or die ("Cannot open file $file: $!");
    foreach my $node (sort {$a cmp $b} keys %$href) {
        print $fh $node."\t".$href->{$node}."\n";
    }
    close ($fh);
}

sub refactor_hash {
    # Re-key hash by values, pushing multiple keys to array
    my ($href) = @_;
    my %hash;
    foreach my $key (sort {$a cmp $b} keys %$href) {
        push @{$hash{$href->{$key}}}, $key;
    }
    return \%hash;
}

sub target_get_cluster {
    # Given a target node, find the cluster of all nodes connected to it,
    # via an iterative approach
    my ($conhref, # Ref to hash of connected nodes
        $clusthref, # Ref to hash of cluster memberships
        $target_node, # Name of target node
        $number, # Current cluster number
        ) = @_;
    my %bait;
    $bait{$target_node} = $number;
    my ($init_count, $curr_count) = (0,1);
    while ($init_count != $curr_count) {
        ($init_count, $curr_count) = (0,0); # Reset counters
        # Count number of bait nodes
        $init_count = scalar keys %bait;
        # Mark contigs connected to the bait as bait themselves
        foreach my $baitnode (keys %bait) {
            foreach my $connnode (@{$conhref->{$baitnode}}) {
                $bait{$connnode} = $number;
            }
        }
        # Count number of bait nodes after fishing
        $curr_count = scalar keys %bait;
        # print "$init_count\t$curr_count\n";
    }
    # Hash result of fished nodes
    foreach my $node (keys %bait) {
        $clusthref->{$node} = $number;
    }
}

sub read_hash_fastg {
    # Read Fastg file
    #   store edge ID in hash
    #   store connectivity info into %connected_hash;
    my ($file, $edgehref, $conhref) = @_;
    open(my $fh, "<", $file) or die ("Cannot open Fastg file $file: $!");
    my $current_node;
    my $writeflag = 0;
    while (my $line = <$fh>) {
        chomp $line;
        if ($line =~ m/^>(.*);$/) {
            # Fastg ID lines start with > character and end with ;
            my $entry = $1;
            # : character separates current node id and list of connected nodes
            my ($id, $conn) = split /:/, $entry;
            my $id_only = $id =~ s/'$//r;

            if (defined $conn) {
                # List of connected nodes is separated by ,
                my @conns = split /,/, $conn;
                # Add to list of connected IDs for this node
                foreach my $currconn (@conns) {
                    $currconn =~ s/'$//; # Remove trailing ' character that indicates revcomp
                    push @{$conhref->{$id_only}}, $currconn;
                }
            }
            # Store Edge name
            $edgehref->{$id_only} = 1;
        }
    }
    close ($fh);
}

sub get_seq_lens {
    # Get lengths of sequences
    my ($href) = @_;
    my %hash;
    foreach my $key (keys %$href) {
        my $len = length $href->{$key}{'seq'};
        my $id = $href->{$key}{'id'};
        $hash{$id} = $len;
    }
    return \%hash;
}

sub welcome {
    print STDERR "This is phyloFlash_fastgFishing.pl from phyloFlash v$VERSION\n";
}


=head1 COPYRIGHT AND LICENSE

Code partially adapted from fastg_paths_fishing.pl script in the gbtools package
https://github.com/kbseah/genome-bin-tools

Copyright (C) 2018- Brandon Seah (kbseah@mpi-bremen.de)

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

=cut
