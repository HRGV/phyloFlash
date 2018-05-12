use strict;
# see bottom of file for copyright & license

package PhyloFlash;
use Exporter qw(import);
use Time::Piece;
use Text::Wrap;
use Config;

$Text::Wrap::huge = "overflow";

=head1 NAME

PhyloFlash Perl module - functions for phyloFlash

=head1 SYNOPSIS

   use PhyloFlash;

=head1 DESCRIPTION

This module contains helper functions shared by the phyloFlash scripts.

=head2 FUNCTIONS

=over

=cut

our $VERSION     = "3.2b1";
our @ISA         = qw(Exporter);
our @EXPORT      = qw(
  $VERSION
  get_cpus
  msg
  @msg_log
  err
  version_sort
  file_is_newer
  get_subdirs
  open_or_die
  slurpfile
  csv_escape
  require_tools
  check_environment
  check_vsearch_version
  run_prog
  run_prog_nodie
  file_download
  fasta_copy_except
  fasta_copy_iupac_randomize
  cluster
  hash2taxstring_counts
  hashtreeconsensus
  taxstring2hash
  consensus_taxon_counter
  revcomp_DNA
  fix_sortmerna_sam
  initialize_outfiles_hash
);

use IPC::Cmd qw(can_run);

=item get_cpus ()

Returns the number of CPUs as listed in F</proc/cpuinfo> (Linux)
or reported by F<sysctl> (Darwin/OSX)

=cut
sub get_cpus {
    if ($Config{"osname"} eq "linux") {
        open(my $fh, "<", "/proc/cpuinfo") or return 1;
        return scalar (map /^processor/, <$fh>)
    }
    elsif ($Config{"osname"} eq "darwin") {
        can_run("sysctl") or return 1;
        my $ncpu = `sysctl -n hw.ncpu`;
        chomp($ncpu);
        return $ncpu;
    }
    else {
        return 1;
    }
}


=item msg ($msg)

Logs a message to STDERR with time stamp prefix.

=cut
our @msg_log; # Global var to store log
sub msg {
    my $t = localtime;
    my $line = wrap ("[".$t->hms."] ", "           ",
                     join "\n", @_);
    push @msg_log, $line;
    print STDERR $line."\n";
}

=item err($msg)

Logs an error message to STDERR and dies

=cut
sub err {
    my @msg = (@_,"Aborting.");
    $msg[0] = "FATAL: ".$msg[0];
    msg(@msg);
    exit(3);
}

=item version_sort(@paths)

sorts paths by descending alphanumeric order ignoring directories and
the part of the filename preceding the first digit. E.g.:

/home/xx/silva_125
/var/lib/dbs/silva_ssu_123.1
/software/phyloFlash_2.0/ssu_123
./121

=cut
sub version_sort {
    # this will fail with Perl <5.13.2. Versions before that lack the
    # non-destructive subsitition flag /r used in the map line below.
    return  map { $_->[0] }
            sort { - ($a->[1] cmp $b->[1]) }
            map { [ $_, ($_ =~ s/.*\/[^\d]*//r) =~ s[(\d+)][pack "N", $1]ger ] }
            @_;
    # works like this:
    # last map creates a key->value map with the path the key and the value
    # 1. having leading directories and filename up to first digit removed
    # 2. all numbers converted to 4-byte strings (most significant first)
    # sort does descending sort by string comparison
    # first map picks just the key to return sorted array
}


=item file_is_newer ($file1, $file2)

Compares the time stamps of two files. Returns true if file1
is newer than file2.

=cut
sub file_is_newer {
    my $file1 = shift;
    my $file2 = shift;

    return (stat($file1))[9] > (stat($file2))[9];
}


=item get_subdirs ($parent)

Returns array of paths to subdirs of $parent (hidden excluded).

=cut
sub get_subdirs {
    my $parent = shift;
    opendir(my $dh, $parent) or return ();
    my @dirs = map { "$parent/" . $_} grep { !/^\./ && -d "$parent/$_" } readdir($dh);
    closedir($dh);
    return @dirs;
}


=item open_or_die (\$fh, $mode, $filename)

Opens a file handle and terminates with an error if the
open failed. $fh must be passed as reference!

=cut
sub open_or_die {
    my ($fh, $mode, $fname) = @_;
    my $msg;

    if ($mode eq ">") {
        $msg = "write to file";
    } elsif ($mode eq ">>") {
        $msg = "append to file";
    } elsif ($mode eq "<") {
        $msg = "read from file";
    } elsif ($mode eq "-|") {
        $msg = "run command";
    }

    open($$fh, $mode, $fname)
        or err("Failed to $msg '$fname': $!");
}

=item slurpfile ($filename)

Opens a file, slurps contents into string, and returns string.

=cut

sub slurpfile {
    my ($file) = @_;
    my $return;
    {
        my $fh_slurp;
        open_or_die(\$fh_slurp, "<", $file);
        local $/ = undef;
        $return = <$fh_slurp>;
        close($fh_slurp);
    }
    return ($return);
}

=item csv_escape ($var)

Escapes data for CSV output according to RFC4180.

=cut
sub csv_escape {
    return unless defined wantarray;

    my @parms = @_;

    for (@parms) {
        if (m/[",\r\n]/) {
            s/"/""/g;
            $_='"'.$_.'"';
        }
    }

    return wantarray ? @parms : $parms[0];
}

# we store a hash of tool_name=>tool_binary here
my %progs = (
    vsearch => "vsearch"
    );

=item require_tools (%hash)

Adds tools to the list of required tools. Argument must be
a hash. Use like so:

  require_tools(("tool" => "mytool.pl"))

=cut
sub require_tools {
    my %arg = @_;
    %progs = (%progs, %arg);
}

=item check_environment ()

Verifies availability of all tools requested using require_tools().

The required tools are located using can_run and the paths determined
printed using msg(). If any of the tools are missing, this function
will abort, asking the user to fix the prerequisites.

=cut
sub check_environment {
    my @missing;

    msg("Checking for required tools.");
    foreach my $prog (keys %progs) {
        my $progname = $progs{$prog};
        if ($progs{$prog} = can_run($progname)) {
            msg("Using $prog found at \"".$progs{$prog}."\".");
        } else {
            push @missing, "  $prog ($progname)";
        }
    }
    if (@missing) {
        msg("Unable to find all required tools. These are missing:");
        foreach my $prog (@missing) {
            msg($prog);
        }
        err("Please make sure these are installed and in your PATH.\n\n");
    } else {
        msg("All required tools found.");
    }
}

=item check_vsearch_version ()

Check that Vsearch version is at least 2.5.0. Return "1" if yes, otherwise undef

UDB file is only implemented from 2.5.0 onwards

=cut

sub check_vsearch_version {
    my $version_msg = `vsearch -v 2>&1`;
    my $return_value;
    if ($version_msg =~ m/vsearch v(\d+)\.(\d+)\.\d+/) {
        my ($maj,$min) = ($1, $2); # Major and minor version
        # Check if at least v2.5.0
        if ($maj >= 2 && $min >= 5) {
            $return_value = 1;
            msg ("Vsearch v2.5.0+ found, will index database to UDB file");
        }
    } else {
        msg ("Could not determine Vsearch version. Assuming to be < v2.5.0");
    }
    return $return_value;
}

=item run_prog ($progname, $args, $redir_stdout, $redir_sterr)

Runs a tool by the name given in require_tools.
Aborts the script if the tool returned with an error code.

=over 15

=item $progname

the name of the tool

=item $args

command line arguments to be passed

=item $redir_stdout

file or file descriptor to redirect stdout to

=item $redir_stdin

same for stderr

=back

=cut
sub run_prog {
    my ($prog, $args, $redir_stdout, $redir_stderr) = @_;

    if (not exists $progs{$prog}) {
        msg("trying to launch unknown tool \"$prog\". pls add to %progs");
        $progs{$prog} = can_run($prog) or err("Failed to find $prog");
    }
    my $cmd = $progs{$prog}." ".$args;
    $cmd .= " >".$redir_stdout if ($redir_stdout);
    $cmd .= " 2>".$redir_stderr if ($redir_stderr);

    msg("running subcommand:","$cmd");
    system($cmd) == 0
        or err("Tool execution failed!.",
               "Error was '$!' and return code '$?'");

    # FIXME: print tail of stderr if redirected
}

=item run_prog_nodie ($progname, $args, $redir_stdout, $redir_sterr)

Runs a tool by the name given in require_tools.
Does not abort the script if returned with error code

=over 15

=item $progname

the name of the tool

=item $args

command line arguments to be passed

=item $redir_stdout

file or file descriptor to redirect stdout to

=item $redir_stdin

same for stderr

=back

=cut
sub run_prog_nodie {
    my ($prog, $args, $redir_stdout, $redir_stderr) = @_;

    if (not exists $progs{$prog}) {
        msg("trying to launch unknown tool \"$prog\". pls add to %progs");
        $progs{$prog} = can_run($prog) or err("Failed to find $prog");
    }
    my $cmd = $progs{$prog}." ".$args;
    $cmd .= " >".$redir_stdout if ($redir_stdout);
    $cmd .= " 2>".$redir_stderr if ($redir_stderr);

    msg("running subcommand:","$cmd");
    my $return = system($cmd);
    if ($return != 0) {
        msg ("Tool execution failed!.",
             "Error was '$!' and return code '$?'");
    }
    return ($return);
    # FIXME: print tail of stderr if redirected
}

=item file_md5 ($filename)

Computes the (hex coded) MD5 sum for the contents of a local file ($filename).

=cut
sub file_md5 {
    my $file = shift;
    my $fh;
    open($fh, "<", $file)
        or err("Unable to open $file.");
    my $ctx = Digest::MD5->new;
    $ctx->addfile($fh);
    close($fh);
    return $ctx->hexdigest;
}

=item ftp_read_var ($ftp, $pattern)

Fetches the contents of a file from FTP. $ftp must be a connected
Net::FTP object and $pattern a file pattern relative to the
current path of $ftp. If the file exists, the contents are returned.
Otherwise returns an empty string

=cut
sub ftp_read_var {
    my $ftp = shift;
    my $file = shift;
    my $out = "";

    open(my $fh, '>', \$out);
    $ftp->get($file, $fh);
    close($fh);

    return $out;
}


=item file_download ($path)

Downloads a file from a FTP server.

$path should be of the form "ftp.somedomain.tld/some/dirs/pa.*ttern".

Returns the name of the downloaded file.

=over 3

=item -

Compares byte sizes before and after download.

=item -

Compares MD5 sum if <file>.md5 exists on server.

=item -

File download is skipped if both size and MD5 (if exists) are correct

=item -

Asks user to confirm license if the directory containing the target
file also contains a file named LICENSE.txt.

=back

=cut
sub file_download {
    my ($server,$path,$pat) = ($_[0] =~ /([^\/]*)(\/.*\/)([^\/]*)/);

    msg("  Connecting to $server");
    my $ftp = new Net::FTP($server,(
                               Passive => 1,
                               Debug => 0,
                               Timeout => 600
                           ))
        or err("Could not connect to $server");
    $ftp->login("anonymous", "-phyloFlash")
        or err("Could not login to $server:", $ftp->message);
    $ftp->binary();
    $ftp->pasv();

    msg("  Finding $path$pat");
    $ftp->cwd($path)
        or err("Could not enter path '\$path\': ", $ftp->message);
    my $files = $ftp->ls($pat)
        or err("Could not list files matching \'$pat'\ in \'$path\': ",
               $ftp->message);
    err("No files found?!")
        if (@$files == 0);
    msg("  Multiple files found?! Using first of ".join(@$files))
        if (@$files > 1);
    my $file = shift(@$files);
    my $file_size = $ftp->size($file)
        or err("Could not get file size:", $ftp->message);
    msg("  Found $file ($file_size bytes)");

    # try downloading md5
    my $file_md5 = ftp_read_var($ftp, $file.".md5");
    $file_md5 =~ s/ .*//;
    chomp($file_md5);

    # compare with potentially existing file
    if (-e $file) { # file exists
        msg("  Found existing $file");
        if (-s $file == $file_size) { # and has same size
            if ($file_md5 eq "") {
                msg("Sizes match. Skipping download");
                return $file;
            }
            msg("  Verifying MD5...");
            my $local_md5 = file_md5($file);
            if ($file_md5 eq $local_md5) { # and has same md5
                msg("  Verified. Skipping download");
                return $file;
            } else { # md5 mismatch
                msg("  MD5 sum mismatch: '$file_md5' != '$local_md5'");
            }
        } else {
            msg("  Size mismatch: ".$file_size." != ". -s $file);
        }
        msg("  -> re-downloading");
    }

    # verify license

    my $license = ftp_read_var($ftp, "LICENSE.txt");
    if (!$license eq "") {
        msg("The file you are about to download comes with a license:\n\n\n"
            .$license."\n\n");
        msg("Do you wish to continue downloading under the conditions");
        msg("specified above? [yes/no]: ");
        my $accept;
        while (<>) {
            chomp;
            if ($_ eq "yes") {
                $accept = 1;
                last;
            } elsif ($_ eq "no") {
                $accept = 0;
                last;
            } else {
                msg("[yes/no]: ");
            }
        }
        if ($accept == 0) {
            msg("Ok. Goodbye...");
            exit(0);
        }
    }

    # download file
    print STDERR "|" . "-" x 75 . "|\n";
    $ftp->hash(\*STDERR, $ftp->size($file)/76);
    $ftp->get($file)
        or err("Failed to download $file:", $ftp->message);
    print STDERR "\n";

    err("File size mismatch?!")
        if (-s $file != $file_size);

    if ($file_md5 eq "") { # had no md5
        return $file;
    }

    msg("  Verifying MD5...");
    my $local_md5 = file_md5($file);
    if ($local_md5 eq $file_md5) {
        msg("File ok");
        return $file;
    }

    err("  MD5 sum mismatch: '$file_md5' != '$local_md5'");
}

=item fasta_copy_except ($source, $dest, @accs)

Extracts all but the sequences listed in @accs from FASTA file
$source into FASTA file $dest.

=cut
sub fasta_copy_except {
    my ($source, $dest, @accs) = @_;
    my $sfh;
    my $dfh;
    open_or_die(\$sfh, "<", $source);
    open_or_die(\$dfh, ">", $dest);

    my %acc_hash = map { $_ => 1 } @accs;

    my $skip = 0;
    while(my $row = <$sfh>) {
        if (substr($row, 0, 1) eq '>') {
            my ($acc) = ($row =~ m/>([^ ]*) /);
            $skip = exists $acc_hash{$acc};
        }
        print $dfh $row if (!$skip);
    }
}

=item cluster ($source, $dest, $id)

Runs vsearch's cluster_fast algorithm to extract centroids
from $source into $dest at a cluster size of $id.

=cut
sub cluster {
    my ($src, $dst, $id, $cpus) = @_;
    msg("clustering database");
    run_prog("vsearch",
             "  --cluster_fast $src "
             . "--id $id "
             . "--centroids $dst "
             . "--notrunclabels "
             . "--threads $cpus ");
}

# hash of IUPAC characters coding for multiple possible bases
my %IUPAC_DECODE = (
    "X" => "ACGT",
    "R" => "AG",
    "Y" => "CT",
    "M" => "CA",
    "K" => "TG",
    "W" => "TA",
    "S" => "CG",
    "B" => "CTG",
    "D" => "ATG",
    "H" => "ATC",
    "V" => "ACG",
    "N" => "ACTG");


=item fasta_copy_iupac_randomize ($source, $dest)

Creates a normalized FASTA file $dest from FASTA file $source.

=over 3

=item -

removes alignment characters ("." and "-")

=item -

uppercases all bases

=item -

turns RNA into DNA (U->T)

=item -

replaces echo IUPAC coded ambiguous base with a base randomly chosen
from the set of options (i.e. replaces B with C, T or G).

=back

=cut
sub fasta_copy_iupac_randomize {
    my ($infile, $outfile) = @_;
    my ($ifh, $ofh);
    open_or_die(\$ifh, "<", $infile);
    open_or_die(\$ofh, ">", $outfile);

    # iterate over lines of FASTA file
    while(my $row = <$ifh>) {
        # pass through FASTA header
        if (substr($row, 0, 1) eq ">") {
            print $ofh $row;
            next;
        }

        # remove alignment, uppercase, turn U into T
        $row =~ tr/a-zA-Z.-/A-TTV-ZA-TTV-Z/d;

        # split into ok / not-ok letter segments
        for my $part ( split(/([^AGCT\n])/, $row) ) {
            if (not $part =~ m/[AGCT\n]/) {
                # segment in need of fix, iterate chars
                for my $i ( 0 .. (length($part)-1) ) {
                # get replacement character list
                    my $rpl = $IUPAC_DECODE{substr($part, $i, 1)};
                    substr($part, $i, 1) =
                        substr($rpl, rand(length($rpl)), 1)
                            if defined $rpl;
                }
            }
            print $ofh $part;
        }
    }
}

{
    package Timer;
    use strict;
    use Time::Piece;
    use Time::Seconds;

=item Timer->new ()

starts a new timer

=cut
    sub new {
        my ($class) = @_;
        my $self = bless {}, $class;
        $self->{time} = localtime;
        return $self;
    }

=item Timer->minutes ()

returns the runtime of the timer in minuts

=cut
    sub minutes {
        my ($self) = @_;
        my $endtime = localtime;
        my $diff = $endtime - $self->{time};
        return sprintf "%.2f minutes", $diff->minutes;
    }
}


=item hash2taxstring_counts

Collapse taxonomy tree into taxstrings and report counts. Taxa are encoded as
hash refs; counts as hash values.

=cut

sub hash2taxstring_counts {
    my ($href, $taxstring, $href2) = @_;
    foreach my $key (keys %$href) {
        if (ref ($href->{$key}) eq 'HASH') { # Recursion
            hash2taxstring_counts (\%{$href->{$key}}, "$taxstring;$key", $href2);
        } else { # End condition - have reached a count
            $href2->{"$taxstring;$key"} = $href->{$key};
        }
    }
}

=item hashtreeconsensus

Walk hash of taxonomy tree and report where it diverges from single branch. 
As a side effect, report taxstring of the congruent portion of the tree

=cut

sub hashtreeconsensus {
    my ($href, # Reference to hash of tax tree
        $aref  # Reference to array to store the congruent portion of tree
        ) = @_;
    my @keys = keys %$href;
    if (scalar @keys == 1) {
        push @$aref, $keys[0];
        if (ref($href->{$keys[0]}) eq 'HASH') { # Recursion
            hashtreeconsensus(\%{$href->{$keys[0]}}, $aref);
        } else {
            return '1'; # The tree is identical to the tips
        }
    } else {
        my @diverge = keys %$href;
        return \@diverge; # Return the level at which the tree diverges
    }
}

=item taxstring2hash

Convert taxonomy string to a nested hash structure recursively

=cut

sub taxstring2hash {
    my ($href, # Reference to hash
        $aref # Reference to taxonomy string as array
        ) = @_;
    my $taxon = shift @$aref;
    if (@$aref) { # Recursion
        taxstring2hash (\%{$href->{$taxon}}, $aref);
    } else { # End condition - count number of occurrences of this taxon
        $href->{$taxon} ++;
    }
}

=item consensus_taxon_counter

Take consensus portion of an array of taxon strings and hash into a taxonomy tree

=cut

sub consensus_taxon_counter {
    my ($href_in,   # Hash ref for taxonomy trees keyed by read
        $aref,      # Array ref of list of reads to summarize
        $taxlevel,  # Taxonomic level for summarizing
        ) = @_;
    
    my %hashout;
    my %taxhash;
    
    foreach my $read (@$aref) {
        next if (!defined $href_in->{$read}); # Skip if this read has no taxonomic assignment
        # Get consensus taxstring for a given read
        my @outarr;
        my $return = hashtreeconsensus(\%{$href_in->{$read}}, \@outarr);
        if (defined $return && @outarr) {
            if (scalar @outarr > $taxlevel) {
                # Trim to max taxonomic rank requested
                @outarr = @outarr[0..$taxlevel-1];
            } elsif (scalar @outarr < $taxlevel) {
                # If taxonomic rank of consensus doesn't reach to requested rank,
                # repeat the lowest taxonomic rank in brackets until requested 
                # rank, e.g. Bacteria;Proteobacteria;(Proteobacteria);(Proteobacteria); etc...
                my $diff = $taxlevel - scalar (@outarr);
                push @outarr, ("($outarr[$#outarr])") x $diff; # Parens before x operator make array
            }

            # Hash the consensus taxonomy into a taxonomic tree, with counts
            # as end values for all reads
            taxstring2hash(\%{$taxhash{'ROOT'}}, \@outarr);
        }
    }
    
    my %taxcounts;
    hash2taxstring_counts (\%taxhash, '', \%taxcounts);
    foreach my $key (keys %taxcounts) { # Replace ugly ;ROOT; from taxstring to be displayed
        my $display_taxstring = $key;
        $display_taxstring =~ s/^;ROOT;//; 
        $hashout{$display_taxstring} = $taxcounts{$key};
    }
    
    return \%hashout;
}

=item revcomp_DNA ($seq)

Return reverse complement of a DNA sequence

=cut 

sub revcomp_DNA {
    my ($seq) = @_;
    my $rev = reverse $seq;
    $rev =~ tr/ATCGatcgYRWSKMDVHByrwskmdvhb/TAGCtagcRYWSMKHBDVrywsmkhbdv/; # Include ambiguity codes
    return $rev;
}

=item fix_sortmerna_sam

Fix SAM file produced by Sortmerna. Returns ref to array of fixed SAM file lines

=cut

sub fix_sortmerna_sam {
    # Fix SAM file produced by Sortmerna v2.1b
    my ($fastq,
        $sam_original,
        $acc2tax,
        $semode) = @_;
    
    my $pe;
    if ($semode == 0) {
        $pe = 1;
    } else {
        $pe = 0;
    }
    
    # Read in Fastq file and hash read orientations and sequences
    my $fastq_href = read_interleaved_fastq ($fastq);
    #print Dumper $fastq_href;
    
    # Read in SAM file and fix bit flags 0x1, 0x40, 0x80, 0x100
    my ($sambyread_href, $sambyline_aref) = read_sortmerna_sam($sam_original, $fastq_href, $pe);
    
    # Fix flags related to read pairs: 0x2, 0x4, 0x8, 0x20
    fix_pairing_flags($sambyread_href, $fastq_href) if $pe == 1;
    #print Dumper $sambyread_href;
    
    # Add back the tax strings to RNAME field
    if (defined $acc2tax) {
        my $acc2tax_href = retrieve ($acc2tax);
        fix_rname_taxstr($sambyread_href, $acc2tax_href);
    }
    
    # Flatten hash back to array and print to output
    my $aref = samaref_to_lines($sambyline_aref, $fastq_href);
    #foreach my $line (@$aref) {
    #    print "$line\n";
    #}
    
    # Return ref to array containing the fixed SAM file
    return $aref;
}

sub fix_rname_taxstr {
    # Replace the accession number of RNAME with the original version
    # retrieved from hash of acc vs taxstring
    my ($sam_href,
        $acc_href
       ) = @_;
    foreach my $id (keys %$sam_href) {
        foreach my $segment (keys %{$sam_href->{$id}}) {
            foreach my $rname (keys %{$sam_href->{$id}{$segment}}) {
                my $new_rname = join " ", ($sam_href->{$id}{$segment}{$rname}{'RNAME'}, $acc_href->{$sam_href->{$id}{$segment}{$rname}{'RNAME'}});
                $sam_href->{$id}{$segment}{$rname}{'RNAME'} = $new_rname;
            }
        }
    }
}

sub samaref_to_lines {
    # Reconstruct SAM lines, and also insert dummy entries for unmapped read fwd segments
    # Input is an array of hash references and hash of Fastq sequences produced
    # by read_interleaved_fastq()
    my ($aref,
        $fastq_href) = @_;
    my @lines;
    foreach my $href (@$aref) {
        # Splice in first segment dummy entry if this is rev with fwd read unmapped
        if ($href->{'FLAG'} & 0x80 && $href->{'FLAG'} & 0x8) {
            my @splice = ($href->{'QNAME'}, # QNAME
                          0x1 + 0x4 + 0x40, # Paired read, segment unmapped, first segment
                          '*', # RNAME
                          '0', # POS
                          '255', # MAPQ
                          '*', # CIGAR
                          '*', # RNEXT
                          '0', # PNEXT
                          '0', # TLEN
                          $fastq_href->{$href->{'QNAME'}}{'fwd'}{'seq'}, # SEQ
                          $fastq_href->{$href->{'QNAME'}}{'fwd'}{'qual'} # QUAL
                          );
            push @lines, join "\t", @splice unless $href->{'FLAG'} & 0x100;
        }
        
        # Otherwise continue
        my @outarr;
        my @fields = qw(QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL);
        foreach my $field (@fields) {
            push @outarr, $href->{$field};
        }
        push @outarr, $href->{'OPTIONAL'} if defined $href->{'OPTIONAL'};
        push @lines, join "\t", @outarr;
    }
    return \@lines;
}

sub fix_pairing_flags {
    # Correct the following flags:
    #  - 0x2
    #  - 0x4
    #  - 0x8
    #  - 0x20
    my ($sam_href,
        $fastq_href
        ) = @_;
    foreach my $id (keys %$sam_href) {
        if (defined $sam_href->{$id}{'fwd'} && $sam_href->{$id}{'rev'}) {
            my %segment_revcomp;
            my %other_segment = ( 'rev' => 'fwd',
                                  'fwd' => 'rev' );
            # Flag 0x2 - both segments have alignments
            foreach my $segment (keys %{$sam_href->{$id}}) {
                foreach my $ref (keys %{$sam_href->{$id}{$segment}}) {
                    $sam_href->{$id}{$segment}{$ref}{'FLAG'} += 0x2 unless $sam_href->{$id}{$segment}{$ref}{'FLAG'} & 0x2;
                    # Record if flag 0x10 set for this segment
                    $segment_revcomp{$segment} ++ if $sam_href->{$id}{$segment}{$ref}{'FLAG'} & 0x10;
                }
            }
            # Sweep through once more
            # Flag 0x20 - next segment rev comp'ed
            foreach my $segment (keys %{$sam_href->{$id}}) {
                foreach my $ref (keys %{$sam_href->{$id}{$segment}}) {
                    if (defined $segment_revcomp{$other_segment{$segment}}) {
                        $sam_href->{$id}{$segment}{$ref}{'FLAG'} += 0x20 unless $sam_href->{$id}{$segment}{$ref}{'FLAG'} & 0x20;
                    }
                }
            }
        } else {
            # Unflag 0x2
            # Flag 0x8 - next segment unmapped
            foreach my $segment (keys %{$sam_href->{$id}}) {
                foreach my $ref (keys %{$sam_href->{$id}{$segment}}) {
                    $sam_href->{$id}{$segment}{$ref}{'FLAG'} -= 0x2 if $sam_href->{$id}{$segment}{$ref}{'FLAG'} & 0x2;
                    $sam_href->{$id}{$segment}{$ref}{'FLAG'} += 0x8 unless $sam_href->{$id}{$segment}{$ref}{'FLAG'} & 0x8;
                }
            }
        }
    }
    
    # No return value because it modifies hash in place
}

sub read_sortmerna_sam {
    # Read in SAM file,
    #  - hash entries by read and by line number
    #  - Correct the following bit flags:
    #     0x1
    #     0x40
    #     0x80
    #     0x100
    my ($file,
        $fastq_href,
        $pe,
        ) = @_;
    #my %samhash_by_line;
    my @samhref_arr;
    my %samhash_by_qname_segment_ref;
    open(my $fh, "<", $file) or die ("$!");
    my $counter = 0;
    while (my $line = <$fh>) {
        chomp $line;
        next if ($line =~ m/^@/); # Skip header lines
        
        # Split entry into fields
        my $href = split_sam_line_to_hash($line);
        
        # Hash SAM record by line number
        push @samhref_arr, $href;
        #$samhash_by_line {$counter} = $href;
        
        # Split read ID on first whitespace [not necessary as sortmerna already does this]
        #my ($id, @discard) = split / /, $href->{'QNAME'};
        my $id = $href->{'QNAME'};
        
        # Get sequence revcomp if bitflag 0x10 is set
        my $sequence_original;
        if ($href->{'FLAG'} & 0x16) {
            $sequence_original = revcomp_DNA($href->{'SEQ'});
        } else {
            $sequence_original = $href->{'SEQ'};
        }
        
        # Check whether fwd or rev segment if PE read
        my $segment;
        if ($pe == 1) {
            $segment = $fastq_href->{$id}{'byseq'}{$sequence_original};
            # Flag 0x1 (PE read) unless already set
            $href->{'FLAG'} += 0x1 unless ($href->{'FLAG'} & 0x1);
            # Flag 0x40 or 0x80 (fwd or rev) unless already set
            if ($segment eq 'fwd') {
                $href->{'FLAG'} += 0x40 unless $href->{'FLAG'} & 0x40;
            } elsif ($segment eq 'rev') {
                $href->{'FLAG'} += 0x80 unless $href->{'FLAG'} & 0x80;
            } else {
                # Diagnostics if segment not found
                print STDERR "Segment not found for this sequence under header $id\n";
                print STDERR "$sequence_original\n\n";
                #print Dumper $fastq_href->{$id};
                #print Dumper $href;
            }
        } else {
            $segment = 'fwd';
            # Unflag as 0x1 if erroneously set
            $href->{'FLAG'} -= 0x1 if ($href->{'FLAG'} & 0x1); 
        }
        
        # Check whether it is primary or secondary alignment for this read
        # Assume that supplementary alignments not reported
        if (defined $samhash_by_qname_segment_ref{$href->{'QNAME'}}{$segment}) {
            # Flag 0x100 as secondary alignment if entry for this read already encountered
            $href->{'FLAG'} += 0x100 unless $href->{'FLAG'} & 0x100;
        }
        $samhash_by_qname_segment_ref{$href->{'QNAME'}}{$segment}{$href->{'RNAME'}} = $href;
        $counter ++;
    }
    close($fh);
    
    # Diagnostics
    #print Dumper \%samhash_by_qname_segment_ref;
    #print Dumper \@samhref_arr;
    
    # Return both hashes
    return (\%samhash_by_qname_segment_ref, \@samhref_arr);
}

sub split_sam_line_to_hash {
    # According to SAM v1 specification 2018-04-27
    my $line = shift;
    my @split = split /\t/, $line;
    my %hash;
    # Mandatory fields (spec part 1.4)
    my @fields = qw(QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL);
    foreach my $field (@fields) {
        # Fields should occur in the correct order
        $hash{$field} = shift @split;
    }
    if (@split) {
        # If anything is left-  optional alignment section (spec part 1.5)
        $hash{'OPTIONAL'} = join "\t", @split;
    }
    return \%hash;
}

=item read_interleaved_fastq

=cut

sub read_interleaved_fastq {
    # Assume that interleaved Fastq is properly paired, i.e. produced with
    # "--paired_in" option to sortmerna
    my ($file, $limit) = @_;
    my $counter = 0;
    my %hash;
    my $currid;
    open(my $fh, "<", $file) or die ("$!");
    while (my $line = <$fh>){
        chomp $line;
        my $modulo = $counter % 8;
        last if defined $limit && ( $counter / 8 ) > $limit;
        if ($counter % 8 == 0 && $line =~ m/^@(.+)/) {
            # Header line of fwd read
            my $fullheader = $1;
            my ($id, @discard) = split / /, $fullheader; # Split on whitespace
            $currid = $id;
            if (defined $hash{$currid}) {
                # If another sequence of this name already defined, warn
                print STDERR "WARNING: Splitting Fastq header on whitespace yields non-unique seq ID: $currid\n";
            }
            $hash{$currid}{'fwd'}{'fullheader'} = $fullheader;
        } elsif ($counter % 8 == 1) {
            # Seq line of fwd read
            $hash{$currid}{'fwd'}{'seq'} = $line;
            $hash{$currid}{'byseq'}{$line} = 'fwd';
        } elsif ($counter % 8 == 3) {
            # Quality line of fwd read
            $hash{$currid}{'fwd'}{'qual'} = $line;
        } elsif ($counter % 8 == 4 && $line =~ m/^@(.+)/) {
            # Header line of rev read
            my $fullheader = $1;
            my ($id, @discard) = split / /, $fullheader; # Split on whitespace
            if ($id ne $currid) {
                print STDERR "WARNING: Fastq rev read header $id does not match fwd read $currid at line $counter\n";
                print STDERR "Check if interleaved file correctly formatted\n";
            }
            
            $hash{$currid}{'rev'}{'fullheader'} = $fullheader;
        } elsif ($counter % 8 == 5) {
            # Seq line of rev read
            $hash{$currid}{'rev'}{'seq'} = $line;
            $hash{$currid}{'byseq'}{$line} = 'rev';
        } elsif ($counter % 8 == 7) {
            # Qual line of rev read
            $hash{$currid}{'rev'}{'qual'} = $line;
        }
        $counter ++;
    }
    close ($fh);
    
    return \%hash;
}

=item initialize_outfiles_hash ($libraryNAME, $readsf)

Initialize a hash of the output filenames, descriptions, and flags when given
library name and name of read files. Returns hash reference.

=cut

sub initialize_outfiles_hash {
    my ($libraryNAME,$readsf) = @_;
    # field "made" will keep track of whether file was created
    my %hash = (
      "idhistogram_svg",
      {
        description => "SVG graphic of mapping identity histogram",
        discard     => 0,
        filename    => "$libraryNAME.idhistogram.svg",
        intable     => 0,
      },
      "gff_arch",
      {
        description => "GFF file for Barrnap run of Archaea model",
        discard     => 1,
        filename    => "$libraryNAME.scaffolds.arch.gff",
        intable     => 0,
      },
      "spades_fastaFromBed",
      {
        description => "FASTA file of sequences assembled by SPAdes retrieved by fastaFromBed",
        discard     => 1,
        filename    => "$libraryNAME.spades_rRNAs.fastaFromBed.fasta",
        intable     => 0,
      },
      "spades_fasta",
      {
        description => "FASTA file of sequences assembled by SPAdes",
        discard     => 0,
        filename    => "$libraryNAME.spades_rRNAs.final.fasta",
        intable     => 0,
      },
      "hitstats",
      {
        description => "Hitstats output from BBmap",
        discard     => 0,
        filename    => "$libraryNAME.hitstats",
        intable     => 0,
      },
      "inserthistogram",
      {
        description => "Insert size histogram from BBmap",
        discard     => 0,
        filename    => "$libraryNAME.inserthistogram",
        intable     => 1,
      },
      "ssu_coll_aln_fasta",
      {
        description => "FASTA file of alignment of all full-length sequences",
        discard     => 0,
        filename    => "$libraryNAME.SSU.collection.alignment.fasta",
        intable     => 1,
      },
      "report_csv",
      {
        description => "phyloFlash report in CSV format",
        discard     => 0,
        filename    => "$libraryNAME.phyloFlash.report.csv",
        intable     => 1,
      },
      "gff_euk",
      {
        description => "GFF file for Barrnap run of Eukaryote model",
        discard     => 1,
        filename    => "$libraryNAME.scaffolds.euk.gff",
        intable     => 0,
      },
      "reads_mapped_r",
      {
        description => "Reads (rev) mapping to SSU rRNA database",
        discard     => 0,
        filename    => "$libraryNAME.$readsf.SSU.2.fq",
        intable     => 1,
      },
      "ssu_coll_tree",
      {
        description => "Newick guide tree from MAFFT alignment of all full-length sequences and closest database hits",
        discard     => 0,
        filename    => "$libraryNAME.SSU.collection.fasta.tree",
        intable     => 1,
      },
      "spades_log",
      {
        description => "Log file from SPAdes assembler",
        discard     => 1,
        filename    => "$libraryNAME.spades.out",
        intable     => 0,
      },
      "full_len_class",
      {
        description => "Taxonomic classification of full-length sequences, in CSV format",
        discard     => 0,
        filename    => "$libraryNAME.phyloFlash.extractedSSUclassifications.csv",
        intable     => 1,
      },
      "report",
      {
        description => "phyloFlash report in plain text",
        discard     => 0,
        filename    => "$libraryNAME.phyloFlash",
        intable     => 1,
      },
      "inserthistogram_svg",
      {
        description => "SVG graphic of insert size histogram",
        discard     => 0,
        filename    => "$libraryNAME.inserthistogram.svg",
        intable     => 0,
      },
      "dhbits_nr97_fasta",
      {
        description => "FASTA file of sequences in database with hits to reconstructed sequences clustered at 0.97 identity",
        discard     => 0,
        filename    => "$libraryNAME.all.dbhits.NR97.fa",
        intable     => 0,
      },
      "vsearch_csv",
      {
        description => "CSV file of Vsearch output",
        discard     => 0,
        filename    => "$libraryNAME.all.vsearch.csv",
        intable     => 0,
      },
      "report_html",
      {
        description => "phyloFlash report in HTML format",
        discard     => 0,
        filename    => "$libraryNAME.phyloFlash.html",
        intable     => 1,
      },
      "unassem_csv",
      {
        description => "Taxonomic composition of unassembled SSU reads in CSV format",
        discard     => 0,
        filename    => "$libraryNAME.phyloFlash.unassembled.NTUabundance.csv",
        intable     => 1,
      },
      "mapratio_csv",
      {
        description => "Mapping ratio file in CSV format",
        discard     => 0,
        filename    => "$libraryNAME.mapratio.csv",
        intable     => 0,
      },
      "gff_bac",
      {
        description => "GFF file for Barrnap run of Bacteria model",
        discard     => 1,
        filename    => "$libraryNAME.scaffolds.bac.gff",
        intable     => 0,
      },
      "emirge_log",
      {
        description => "Log file from EMIRGE sequence reconstruction",
        discard     => 1,
        filename    => "$libraryNAME.emirge.out",
        intable     => 0,
      },
      "assemratio_csv",
      {
        description => "CSV file of ratio assembled to unassembled",
        discard     => 0,
        filename    => "$libraryNAME.assemratio.csv",
        intable     => 0,
      },
      "emirge_fasta",
      {
        description => "FASTA file of sequences reconstructed by EMIRGE",
        discard     => 0,
        filename    => "$libraryNAME.emirge.final.fasta",
        intable     => 0,
      },
      "ssu_coll_aln_mafftout",
      {
        description => "MAFFT output from aligning all full-length sequences",
        discard     => 1,
        filename    => "$libraryNAME.SSU.collection.alignment.mafftout",
        intable     => 0,
      },
      "ntu_csv",
      {
        description => "NTU abundances (truncated to requested taxonomic level) from initial mapping, in CSV format",
        discard     => 0,
        filename    => "$libraryNAME.phyloFlash.NTUabundance.csv",
        intable     => 1,
      },
      "ntu_full_csv",
      {
        description => "NTU abundances (untruncated) from initial mapping, in CSV format",
        discard     => 0,
        filename    => "$libraryNAME.phyloFlash.NTUfull_abundance.csv",
        intable     => 1,
      },
      "ntu_csv_svg",
      {
        description => "SVG graphic of taxonomic composition from initial read mapping",
        discard     => 0,
        filename    => "$libraryNAME.phyloFlash.NTUabundance.csv.svg",
        intable     => 1,
      },
      "idhistogram",
      {
        description => "Mapping identity histogram from BBmap",
        discard     => 0,
        filename    => "$libraryNAME.idhistogram",
        intable     => 1,
      },
      "sam_map",
      {
        description => "SAM file of initial read mapping to SSU rRNA database",
        discard     => 0,
        filename    => "$libraryNAME.$readsf.SSU.sam",
        intable     => 1,
      },
      "barrnap_log",
      {
        description => "Log file from Barrnap",
        discard     => 1,
        filename    => "$libraryNAME.barrnap.out",
        intable     => 0,
      },
      "mapratio_svg",
      {
        description => "SVG graphic of mapping ratio file in CSV format",
        discard     => 0,
        filename    => "$libraryNAME.mapratio.csv.svg",
        intable     => 0,
      },
      "ssu_coll_tree_svg",
      {
        description => "SVG graphic of guide tree from MAFFT alignment of all full-length sequences and closest database hits",
        discard     => 0,
        filename    => "$libraryNAME.SSU.collection.fasta.tree.svg",
        intable     => 1,
      },
      "dbhits_all_fasta",
      {
        description => "FASTA file of all sequences in database with hits to reconstructed sequences",
        discard     => 0,
        filename    => "$libraryNAME.all.final.phyloFlash.dbhits.fa",
        intable     => 0,
      },
      "notmatched_fasta",
      {
        description => "FASTA file of full-length sequences without any database hits",
        discard     => 0,
        filename    => "$libraryNAME.all.final.phyloFlash.notmatched.fa",
        intable     => 0,
      },
      "assemratio_svg",
      {
        description => "SVG graphic of assembly ratio",
        discard     => 0,
        filename    => "$libraryNAME.assemratio.csv.svg",
        intable     => 0,
      },
      "reads_mapped_f",
      {
        description => "Reads (fwd) mapping to SSU rRNA database",
        discard     => 0,
        filename    => "$libraryNAME.$readsf.SSU.1.fq",
        intable     => 1,
      },
      "bbmap_log",
      {
        description => "Log file from BBmap of initial mapping",
        discard     => 0,
        filename    => "$libraryNAME.bbmap.out",
        intable     => 0,
      },
      "ssu_coll_fasta",
      {
        description => "FASTA file of all full-length sequences and their closest database hits",
        discard     => 0,
        filename    => "$libraryNAME.SSU.collection.fasta",
        intable     => 1,
      },
      "all_final_fasta",
      {
        description => "FASTA file of full-length SSU sequences",
        discard     => 0,
        filename    => "$libraryNAME.all.final.fasta",
        intable     => 1,
      },
      "phyloFlash_log",
      {
        description => "phyloFlash log file",
        discard     => 0,
        filename    => "$libraryNAME.phyloFlash.log",
        intable     => 0,
      },
      "phyloFlash_archive",
      {
        description => "tar.gz archive of phyloFlash results",
        discard     => 0,
        filename    => "$libraryNAME.phyloFlash.tar.gz",
        intable     => 0,
      },
      "readsf_subsample",
      {
        description => "Subsample of the forward SSU reads for running nhmmer",
        discard     => 0, # Keep while testing
        filename    => "$libraryNAME.readsf.subsample.fasta",
        intable     => 0,
      },
      "nhmmer_tblout",
      {
        description => "Tabular output from aligning SSU HMM models with nhmmer on sample of reads",
        discard     => 0, # Keep while testing
        filename    => "$libraryNAME.nhmmer.tblout",
        intable     => 0,
      },
      "nhmmer_prok_histogram",
      {
        description => "Histogram of alignment position counts against prokaryotic HMM models",
        discard     => 0, # keep while testing
        filename    => "$libraryNAME.nhmmer.prok.histogram",
        intable     => 0,
      },
      "nhmmer_prok_histogram_svg",
      {
        description => "SVG graphic of histogram of alignment position counts against prokaryotic HMM model",
        discard     => 0,
        filename    => "$libraryNAME.nhmmer.prok.histogram.svg",
        intable     => 0,
      },
      "nhmmer_euk_histogram",
      {
        description => "Histogram of alignment position counts against eukaryotic HMM model",
        discard     => 0,
        filename    => "$libraryNAME.nhmmer.euk.histogram",
        intable     => 0,
      },
      "nhmmer_euk_histogram_svg",
      {
        description => "SVG graphic of histogram of alignment position counts against eukaryotic HMM model",
        discard     => 0,
        filename    => "$libraryNAME.nhmmer.euk.histogram.svg",
        intable     => 0,
      },
      "plotscript_out",
      {
        description => "Output stream from Plotscript",
        discard     => 1,
        filename    => "$libraryNAME.plotscript.out",
        intable     => 0,
      },
      "basecompositionhist",
      {
        description => "Base composition histogram from BBmap, from initial mapping of reads vs database",
        discard     => 1,
        filename    => "$libraryNAME.basecompositionhistogram",
        intable     => 0,
      },
      "gff_all",
      {
        description => "GFF of all gene models combined, for extracting SSU sequence from Fasta with fastaFromBed",
        discard     => 1,
        filename    => "$libraryNAME.scaffolds.final.gff",
        intable     => 0,
      },
      "fastaFromBed_out",
      {
        description => "Output stream from fastaFromBed",
        discard     => 1,
        filename    => "$libraryNAME.fastaFromBed.out",
        intable     => 0,
      },
      "reads_mapped_cat",
      {
        description => "Concatenation of both Fwd and Rev reads mapping to SSU rRNA database, for EMIRGE if read lengths are > 150 bp",
        discard     => 1,
        filename    => "$libraryNAME.SSU.all.fq",
        intable     => 0,
      },
      "reads_mapped_cat_rename",
      {
        description => "Concatenation of both Fwd and Rev reads mapping to SSU rRNA database, renamed, for EMIRGE if read lengths are > 150 bp",
        discard     => 1,
        filename    => "$libraryNAME.renamed.SSU.all.fq",
        intable     => 0,
      },
      "emirge_raw_fasta",
      {
        description => "Raw output Fasta file from EMIRGE iteration 40",
        discard     => 1,
        filename    => "$libraryNAME.emirge.result.fasta",
        intable     => 0,
      },
      "vsearch_out",
      {
        description => "Output stream from Vsearch",
        discard     => 1,
        filename    => "$libraryNAME.all.vsearch.out",
        intable     => 0,
      },
      "vsearch_clusterdb_out",
      {
        description => "Output stream from Vsearch cluster_fast",
        discard     => 1,
        filename    => "$libraryNAME.clusterdbhits.out",
        intable     => 0,
      },
      "reformat_out",
      {
        description => "Output stream from reformat.sh",
        discard     => 1,
        filename    => "$libraryNAME.reformat.out",
        intable     => 0,
      },
      "sam_remap_spades",
      {
        description => "SAM file of re-mapping extracted reads to SPAdes full-length sequences",
        discard     => 0,
        filename    => "$libraryNAME.$readsf.SSU_spades.sam",
        intable     => 0,
      },
      "bbmap_remap_log_spades",
      {
        description => "Log file from BBmap of re-mapping to SPAdes sequences",
        discard     => 0,
        filename    => "$libraryNAME.remap_spades.bbmap.out",
        intable     => 0,
      },
      "sam_remap_emirge",
      {
        description => "SAM file of re-mapping extracted reads to EMIRGE full-length sequences",
        discard     => 0,
        filename    => "$libraryNAME.$readsf.SSU_emirge.sam",
        intable     => 0,
      },
      "bbmap_remap_log_emirge",
      {
        description => "Log file from BBmap of re-mapping to EMIRGE sequences",
        discard     => 0,
        filename    => "$libraryNAME.remap_emirge.bbmap.out",
        intable     => 0,
      },
      "trusted_gff_bac",
      {
        description => "GFF file for Barrnap run of Bacteria model on trusted contigs",
        discard     => 1,
        filename    => "$libraryNAME.trusted.bac.gff",
        intable     => 0,
      },
      "trusted_gff_arch",
      {
        description => "GFF file for Barrnap run of Archaea model on trusted contigs",
        discard     => 1,
        filename    => "$libraryNAME.trusted.arch.gff",
        intable     => 0,
      },
      "trusted_gff_euk",
      {
        description => "GFF file for Barrnap run of Eukaryota model on trusted contigs",
        discard     => 1,
        filename    => "$libraryNAME.trusted.euk.gff",
        intable     => 0,
      },
      "trusted_gff_all",
      {
        description => "Concatenated GFF file of SSU rRNA sequences for all three models in trusted contigs",
        discard     => 0,
        filename    => "$libraryNAME.trusted.all.gff",
        intable     => 1,
      },
      "trusted_barrnap_log",
      {
        description => "Log file for Barrnap run of on trusted contigs",
        discard     => 1,
        filename    => "$libraryNAME.trusted.barrnap.log",
        intable     => 0,
      },
      "trusted_fastaFromBed",
      {
        description => "Fasta file of extracted SSU rRNA from trusted contigs from fastaFromBed",
        discard     => 1,
        filename    => "$libraryNAME.trusted_fastaFromBed.all.fasta",
        intable     => 0,
      },
      "trusted_fasta",
      {
        description => "Fasta file of extracted SSU rRNA from trusted contigs",
        discard     => 1,
        filename    => "$libraryNAME.trusted.all.fasta",
        intable     => 0,
      },
      "trusted_fastaFromBed_out",
      {
        description => "Output stream from fastaFromBed for trusted contigs",
        discard     => 1,
        filename    => "$libraryNAME.trusted.fastaFromBed.out",
        intable     => 0,
      },
      "sam_remap_trusted",
      {
        description => "SAM file of read mapping vs trusted contigs",
        discard     => 1,
        filename    => "$libraryNAME.trusted.bbmap.sam",
        intable     => 0,
      },
      "bbmap_remap_log_trusted",
      {
        description => "Log file from running bbmap on trusted contigs",
        discard     => 1,
        filename    => "$libraryNAME.trusted.bbmap.out",
        intable     => 0,
      },
      "reads_mapped_notrusted_f",
      {
        description => "Forward read file of reads not mapping to trusted contigs",
        discard     => 1,
        filename    => "$libraryNAME.trusted.bbmap.outu.fwd.fastq",
        intable     => 0,
      },
      "reads_mapped_notrusted_r",
      {
        description => "Reverse read file of reads not mapping to trusted contigs",
        discard     => 1,
        filename    => "$libraryNAME.trusted.bbmap.outu.rev.fastq",
        intable     => 0,
      },
      "reads_uncompressed",
      {
        description => "Uncompressed input reads for Sortmerna",
        discard     => 1,
        filename    => "$libraryNAME.$readsf.uncompressed.fastq",
        intable     => 0,
      },
      "sortmerna_sam",
      {
        description => "Original SAM file produced by Sortmerna",
        discard     => 1,
        filename    => "$libraryNAME.sortmerna.sam",
        intable     => 0,
      },
      "sortmerna_fastq",
      {
        description => "Original Fastq file produced by Sortmerna",
        discard     => 1,
        filename    => "$libraryNAME.sortmerna.fastq",
        intable     => 0,
      },
      "sortmerna_log",
      {
        description => "Log file produced by Sortmerna",
        discard     => 0,
        filename    => "$libraryNAME.sortmerna.log",
        intable     => 1,
      },
      #"",
      #{
      #  description => "",
      #  discard     => 1,
      #  filename    => "",
      #  intable     => 0,
      #},
    );
    return (\%hash);
}

=back

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2015 Elmar Pruesse <elmar.pruesse@ucdenver.edu>
                   Harald Gruber-Vodicka <hgruber@mpi-bremen.de>

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
along with this program.
If not, see L<http://www.gnu.org/licenses/>.

=cut

1; # keep require happy
