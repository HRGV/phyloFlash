#!/usr/bin/env perl

=head1 NAME

phyloFlash_fastgFishing.pl - Bin genomes from Fastg graph by SSU sequence

=head1 SYNOPSIS

### Megahit assembly

phyloFlash_fastgFishing.pl --fasta [Fasta] --fastg [Fastg] --out [PREFIX]

### SPAdes assembly

phyloFlash_fastgFishing.pl --fasta [Fasta] --fastg [Fastg] --paths [Paths] --out [PREFIX]

### Help

phyloFlash_fastgFishing.pl --help

phyloFlash_fastgFishing.pl --man

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
specified with option I<--compare-ssu>, these will be searched against the Fasta
contigs, and a table of which bin each sequence occurs in will also be reported.
If you have a phyloFlash (v3.1+) output tar.gz archive, you can specify it with
the option I<--compare-zip> and the assembled SSU file will be automatically
extracted for you.

=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib $FindBin::RealBin;
use PhyloFlash qw(msg err @msg_log $VERSION revcomp_DNA);
use File::Temp;
use File::Basename;
use Archive::Tar;
#use Data::Dumper;

# Input files
my $fastgfile; # Input Fastg file, assume from MegaHIT
my $fastafile;
my $pathsfile;
my $out = "test";       # Output file prefix
my $comparessu;
my $comparezip;

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
my $rejectthreshold = 0.2; # Minimum frac of full-length to reject as 16S (for barrnap)
my $dofasta;
my $assembler = "megahit";
my $list_fastg_header;
my $list_fasta_header;
my $barrnap_slurp;
my $barrnap_path = "$FindBin::RealBin/barrnap-HGV/bin/barrnap_HGV";
my $vsearch_path = "vsearch";
my $mafft_path = "mafft";
my $gff_aref;
my $contig_shortlist_href;
my $CPUs = 8;
my $clusteronly;

pod2usage (-verbose => 0) if (!@ARGV);

GetOptions ("fastg=s" => \$fastgfile,
            "fasta=s" => \$fastafile,
            "paths=s" => \$pathsfile,
            "compare-ssu=s" => \$comparessu,
            "compare-zip=s" => \$comparezip,
            "assembler|a=s" => \$assembler,
            "list_fastg=s" => \$list_fastg_header,  # List of Fastg headers to use as "bait" for extracting contig clusters
            "list_fasta=s" => \$list_fasta_header,  # List of Fasta headers to use as "bait" for extracting contig clusters
            "out|o=s" => \$out,
            "cutoff|c=i" => \$cutoff,
            "min-SSU-frac=f" => \$rejectthreshold,
            "outfasta" => \$dofasta,
            "clusteronly" => \$clusteronly,
            "CPUs=i" => \$CPUs,
            "barrnap-path=s" => \$barrnap_path,
            "help|h" => sub { pod2usage(-verbose=>1); },
            "man|m" => sub { pod2usage(-verbose=>3); },
            "version|v" => sub { welcome(); exit(); },
            ) or pod2usage (-verbose=>1);

=head1 ARGUMENTS

=head2 INPUT

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

Assembler used. Either "megahit" or "spades".

Default: 'megahit'

=item --compare-ssu I<FILE>

=item --compare-zip I<FILE>

If a targeted assembly of SSU rRNA has already been performed for this library,
compare the sequences to those extracted from the metagenome assembly by barrnap.

The sequences can be supplied either as a Fasta flie (option I<--compare-ssu>)
or as a phyloFlash tar.gz archive (option I<--compare-zip>).

Default: None.

=back

=head2 OUTPUT OPTIONS

=over 8

=item --out|-o I<STRING>

Output file name prefix

Default: 'test'

=item --clusteronly

Do not search for SSU rRNA sequences, instead only report all connected contig 
clusters above length threshold, regardless of whether they have SSU rRNA or not

Default: Off

=item --cutoff|-c I<INTEGER>

Minimum total sequence length of contig cluster to be reported (Default: 100000)

=item --min-SSU-frac I<NUMERIC>

Minimum fraction of SSU sequence to report (parameter passed to the '--reject'
option in Barrnap)

Default: 0.2

=item --outfasta

Logical: Output Fasta files for clusters with total length above cutoff?

Default: No

=back

=head2 USAGE OPTIONS

=over 8

=item --CPUs I<INTEGER>

Number of CPUs to use for Barrnap rRNA prediction

Default: 8

=item --barnap-path I<PATH>

Path to barrnap executable

Default: Use version bundled with phyloFlash

=item --help|-h

Brief help message

=item --man|-m

Manual page

=item --version|-v

Report phyloFlash version

=back

=cut

## MAIN ########################################################################

welcome();

## INPUT ##############################

# Run barrnap to predict 16S on Fasta file

unless ($clusteronly) {
    msg ("Running Barrnap to extract 16S rRNA sequences from Fasta file $fastafile");
    $gff_aref = run_barrnap_on_fasta($fastafile,$barrnap_path, $CPUs);
    if (defined $gff_aref) {
        $contig_shortlist_href = ssu_contigs_from_barrnap_gff($gff_aref);
        my $num_hits = scalar (keys %$contig_shortlist_href);
        msg ("Number of contigs containing SSU rRNA above $rejectthreshold of full-length found: $num_hits");
    }
}

# Read in the Fastg file containing edges and graph connectivity info
msg ("Reading contig graph from Fastg file $fastgfile");
read_hash_fastg($fastgfile,
                \%edge_hash,
                \%connected_hash);

# Read in the Fasta and Paths files to translate edges to nodes
msg ("Translating edge names to node names");
if ($assembler eq 'megahit') {
    if (defined $pathsfile) {
        msg "WARNING: Ignoring paths file, not part of Megahit output\n";
    }

    read_megahit_fasta($fastafile,
                       \%node_seq_hash);
    my $href = megahit_edge2node(\%edge_hash,
                                 \%node_seq_hash);
    %edge2node_hash = %$href;
} elsif ($assembler eq 'spades') {
    if (!defined $pathsfile) {
        # Catch exception
        err ("ERROR: Please specify SPAdes 'paths' file to --paths\n");
    }
    read_spades_fasta($fastafile,
                      \%node_seq_hash);
    my $href = spades_edge2node ($pathsfile,
                                 \%edge_hash,
                                 \%node_seq_hash);
    %edge2node_hash = %$href;
} else {
    err ("ERROR: Please specify either 'megahit' or 'spades' as --assembler\n");
}

my $href = flip_edge2node(\%edge2node_hash);
%node2edge_hash = %$href;

if (defined $contig_shortlist_href) {
    get_barrnap_hit_seqs($contig_shortlist_href, \%node_seq_hash);
}


## FISHING ############################

my $edge_shortlist_href;
if (defined $list_fastg_header) {
    msg ("Using shortlist of defined Fastg header names from file $list_fastg_header");
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
    msg ("Using shortlist of defined Fasta header names from file $list_fasta_header");
    open(my $fhin, "<", $list_fasta_header) or die ("$!");
    while (my $line = <$fhin>) {
        chomp $line;
        my ($nodeshort, @discard) = split /\s+/, $line;
        foreach my $edge (@{$node2edge_hash{$nodeshort}}) {
            $edge_shortlist_href->{$edge} = $line;
        }
    }
    close($fhin);
} elsif (defined $contig_shortlist_href && ref($contig_shortlist_href) eq 'HASH') {
    msg ("Using contigs identified with SSU rRNA sequences by barrnap");
    foreach my $contig (keys %$contig_shortlist_href) {
        foreach my $edge (@{$node2edge_hash{$contig}}) {
            $edge_shortlist_href->{$edge} = $contig;
        }
    }
} else {
    # Get ALL clusters
    msg ("Finding all clusters of connected contigs");
    $edge_shortlist_href = \%edge_hash;
}

msg ("Iteratively finding connected contigs from Fastg graph");

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

my ($clust_stats_href, $ssuid2clust_href) = get_cluster_stats($node_cluster_href,
                                                              $seq_lens_href,
                                                              $contig_shortlist_href);

msg ("Reporting contig cluster summary to $out.cluster_stats.tab");
msg ("Reporting table of extracted SSU IDs vs bins to $out.cluster_ssu_summary.tab") if defined $gff_aref;
report_cluster_summary($clust_stats_href,
                       $ssuid2clust_href,
                       "$out.cluster_stats.tab",
                       "$out.cluster_ssu_summary.tab");

msg ("Reporting list of contig names per cluster to $out.nodes_to_cluster.tab");
msg ("Writing Fasta files of each genome bin to $out.binXX.fasta") if $dofasta;
report_nodes_to_cluster($clust_stats_href,
                        $clust2node_href,
                        \%node_seq_hash,
                        "$out.nodes_to_cluster.tab");

if (defined $gff_aref) {
    msg ("Writing GFF file from Barrnap to $out.barrnap.gff");
    open(my $fh_gff, ">", "$out.barrnap.gff") or die ("$!");
    foreach my $entry (@$gff_aref) {
        print $fh_gff $entry;
        print $fh_gff "\n";
    }
    close ($fh_gff);
    
    if ($dofasta || defined $comparessu || defined $comparezip) {
        msg ("Writing sequences of extracted SSU rRNA to $out.barrnap_SSU.fasta");
        open(my $fh_ssu, ">", "$out.barrnap_SSU.fasta") or die ("$!");
        foreach my $contig (keys %$contig_shortlist_href) {
            foreach my $id (keys %{$contig_shortlist_href->{$contig}}) {
                print $fh_ssu ">$id\n";
                print $fh_ssu $contig_shortlist_href->{$contig}{$id}{'seq'};
                print $fh_ssu "\n";
            }
        }
        close ($fh_ssu);
    }
}

if (defined $comparezip && -f $comparezip) {
    msg ("Using assembled SSU rRNA from phyloFlash run archived in file $comparezip");
    my ($comparezip_base, $comparezip_dir) = fileparse $comparezip;
    if ($comparezip_base =~ m/(.+)\.phyloFlash\.tar\.gz/) {
        my $lib = $1;
        my $tarhandle = Archive::Tar->new();
        $tarhandle->read($comparezip);
        if ($tarhandle->contains_file("$lib.all.final.fasta")) {
            $tarhandle->extract("$lib.all.final.fasta");
            $comparessu = "$lib.all.final.fasta";
        } else {
            msg ("Fasta file of assembled SSU rRNA not found in phyloFlash archive. Perhaps archive file has been renamed, or no SSU was assembled?");
        }
    } else {
        msg ("Filename supplied to $comparezip does not match expected pattern for phyloFlash tar.gz archive");
    }
}

if (defined $comparessu && -f $comparessu) {
    msg ("Comparing extracted SSU rRNA from assembly to user-specified set of sequences in $comparessu");
    compare_ssu_vsearch($comparessu,
                        "$out.barrnap_SSU.fasta",
                        "$out.compare-ssu.alnout",
                        "$out.compare-ssu.samout",
                        );
    compare_ssu_mafft($comparessu,
                      "$out.barrnap_SSU.fasta");
}

## SUBS ########################################################################

sub compare_ssu_mafft {
    my ($fasta1,
        $fasta2,
       ) = @_;
    my $concatout = "$out.compare-ssu_concat.fasta";
    system ("cat $fasta1 $fasta2 > $concatout");
    my $treeout = "$concatout.tree";
    my @mafft_args = ('--quiet',
                      "--treeout $concatout",
                      ">",
                      "$concatout.mafft");
    my $mafft_cmd = join " ", ($mafft_path, @mafft_args);
    msg ("Running mafft with command: $mafft_cmd");
    system ($mafft_cmd);
}

sub compare_ssu_vsearch {
    my ($db,
        $query,
        $alnoutfile,
        $samoutfile,
        ) = @_;
    my @vsearch_args = ("--usearch_global $query",
                        "--db $db",
                        '--id 0.95',
                        '--maxhits 1',
                        "--threads $CPUs",
                        '--quiet',
#                        "--blast6out $blast6outfile",
                        "--alnout $alnoutfile",
                        "--samout $samoutfile",
                        );
    my $vsearch_cmd = join " ", ($vsearch_path, @vsearch_args);
    msg ("Running vsearch with command: $vsearch_cmd");
    my $returnval = system ($vsearch_cmd);
    return ($returnval);
}

sub get_barrnap_hit_seqs {
    my ($hit_href,
        $nodeseq_href) = @_;
    foreach my $contig (keys %$hit_href) {
        foreach my $id (keys %{$hit_href->{$contig}}) {
            my @splitgff = split /\t/, $hit_href->{$contig}{$id}{'gff'};
            my $start = $splitgff[3];
            my $stop = $splitgff[4];
            my $length = $stop - $start + 1;
            my $offset = $start - 1;
            my $shortid;
            if ($assembler eq 'spades') {
                ($shortid) = $contig =~ m/NODE_(\d+)_/;
            } else {
                ($shortid) = $contig =~ m/k\d+_(\d+)/;
            }
            my $seq = substr $nodeseq_href->{$shortid}{'seq'}, $offset, $length;
            $seq = revcomp_DNA($seq) if ($splitgff[6] eq '-');
            $hit_href->{$contig}{$id}{'seq'} = $seq if defined $seq;
            $hit_href->{$contig}{$id}{'length'} = $length if defined $seq;
            $hit_href->{$contig}{$id}{'offset'} = $offset if defined $seq;
        }
    }
}

sub report_nodes_to_cluster {
    my ($cluster_href,
        $clust2node_href,
        $node2seq_href,
        $outfile,
        ) = @_;
    # Report clusters in descending length above cutoff
    open (my $fh_node2clust, ">", $outfile) or die ("$!");       # File handle for node ID to cluster list
    my $bin_counter = 0;    # Counter for number of bins
    # For each cluster, in descending order of total sequence length...
    foreach my $clust (sort {$cluster_href->{$b}{'length'} <=> $cluster_href->{$a}{'length'}} keys %$cluster_href) {
        my $bin = "bin$bin_counter"; # Name for bin, numbered in descending order of size
        if ($cluster_href->{$clust}{'length'} > $cutoff) {
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
                    print $fh_fasta $node2seq_href->{$shortid}{'seq'}."\n";
                }
            }
            close ($fh_fasta) if $dofasta;
        }
        $bin_counter ++;
    }
    close ($fh_node2clust);
}

sub report_cluster_summary {
    my ($cluster_href,
        $ssuid2clust_href,
        $outfile,
        $outfile2,
       ) = @_;
    my %clust2bin;
    open (my $fh_summary, ">", $outfile) or die ("$!");            # File handle for summary table
    my @head_arr = qw(bin length contigs);
    push @head_arr, 'num_ssu' unless $clusteronly;
    print $fh_summary "#".join("\t", @head_arr)."\n"; # Header line
    my $bin_counter = 0;
    foreach my $clust (sort {$cluster_href->{$b}{'length'} <=> $cluster_href->{$a}{'length'}} keys %$cluster_href) {
        my $bin = "bin$bin_counter"; # Name for bin, numbered in descending order of size
        $clust2bin{$clust} = $bin;
        if ($cluster_href->{$clust}{'length'} > $cutoff) {
            my @out_arr = ($bin,
                           $cluster_href->{$clust}{'length'},
                           $cluster_href->{$clust}{'nodes'},
                           );
            push @out_arr, $cluster_href->{$clust}{'ssu'} unless $clusteronly;
            print $fh_summary join "\t", @out_arr;
            print $fh_summary "\n";
        }
        $bin_counter ++;
    }
    close ($fh_summary);
    if (defined $ssuid2clust_href && ref($ssuid2clust_href) eq 'HASH' && !defined $clusteronly) {
        open (my $fh_ssulist, ">", $outfile2) or die ("$!");
        print $fh_ssulist "#".join ("\t", qw(ssu_id bin note))."\n";
        foreach my $ssuid (sort keys %$ssuid2clust_href) {
            my $clust = $ssuid2clust_href->{$ssuid};
            my @outarr = ($ssuid, $clust2bin{$clust});
            if ($cluster_href->{$clust}{'length'} < $cutoff) {
                push @outarr, 'Below length cutoff';
            } else {
                push @outarr, '';
            }
            print $fh_ssulist join "\t", @outarr;
            print $fh_ssulist "\n";
        }
        close ($fh_ssulist);
    }
}

sub welcome {
    print STDERR "This is phyloFlash_fastgFishing.pl from phyloFlash v$VERSION\n";
}

sub ssu_contigs_from_barrnap_gff {
    # Get names of contigs containing 16S rRNA hits from barrnap GFF results
    # Return ref to array of contig names
    my ($aref,          # Ref to array of GFF file
        ) = @_;
    my %outhash;
    foreach my $line (@$aref) {
        next if $line =~ m/^#/; # Skip header
        my @splitline = split /\t/, $line;
        if (defined $splitline[8] && $splitline[8] =~ m/Name=16S_rRNA/) {
            #push @{$outhash{$splitline[0]}}, $line;
            my ($id) = ($splitline[8] =~ m/ID=(.+?);/);
            $outhash{$splitline[0]}{$id}{'gff'} = $line;
        }
    }
    return \%outhash;
}

sub run_barrnap_on_fasta {
    my ($file,          # Fasta file
        $bin,           # Path to Barrnap
        $threads        # Num threads for barrnap
        ) = @_;
    my @outarr;
    my @barrnap_args = ('--quiet',
                        "--threads $threads",
                        "--reject $rejectthreshold",
                        );
    my $barrnap_cmd = join " ", ($bin, @barrnap_args, $file);
    msg ("Running barrnap with the command: $barrnap_cmd");
    my $barrnap_gff_str = `$barrnap_cmd`;
    my @barrnap_gff_arr = split /\n/, $barrnap_gff_str;
    my %idcounter;
    foreach my $line (@barrnap_gff_arr) {
        my @splitline = split /\t/, $line;
        if (defined $splitline[8]) {
            $idcounter{$splitline[0]} ++;
            $splitline[8] = "ID=$splitline[0]_$idcounter{$splitline[0]};".$splitline[8];
        }
        push @outarr, join "\t", @splitline;
    }
    return \@outarr;
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
    my ($clusthref, $lenhref, $ssuhit_href) = @_;
    my %hash;       # Hash of cluster summary stats
    my %hash2;      # Hash of cluster hashed by ssu ID
    foreach my $node (keys %$clusthref) {
        my ($nodeshort,@discard) = split /\s+/, $node;
        my $clust = $clusthref->{$node};
        my $len = $lenhref->{$node};
        $hash{$clust}{'length'} += $len if defined $len;
        $hash{$clust}{'nodes'} ++ if defined $clust;
        if (defined $ssuhit_href->{$nodeshort}) {
            $hash{$clust}{'ssu'} += scalar(keys %{$ssuhit_href->{$nodeshort}});
            foreach my $id (keys %{$ssuhit_href->{$nodeshort}}) {
                $hash2{$id} = $clust;
            }
        }
    }
    return (\%hash, \%hash2);
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
