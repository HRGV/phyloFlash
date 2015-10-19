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

our $VERSION     = 2.00;
our @ISA         = qw(Exporter);
our @EXPORT      = qw(
  get_cpus
  msg
  err
  version_sort
  file_is_newer
  get_subdirs
  open_or_die
  csv_escape
  require_tools
  check_environment
  run_prog
  file_download
  fasta_copy_except
  fasta_copy_iupac_randomize
  cluster
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
sub msg {
    my $t = localtime;
    my $line = wrap ("[".$t->hms."] ", "           ",
                     join "\n", @_);
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
    return  map { $_->[0] }  
            sort { - ($a->[1] cmp $b->[1]) }
            map { [ $_, ($_ =~ s/.*\/[^\d]*//r) =~ s[(\d+)][pack "N", $1]ger ] } 
            @_; 
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
