#  PhyloFlash.pm
#
#  Copyright (C) 2015- Elmar Pruesse <elmar.pruesse@ucdenver.edu>
#                      Harald Gruber-Vodicka <hgruber@mpi-bremen.de>
#
#  LICENCE
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.
#  If not, see <http://www.gnu.org/licenses/>.
#
#
#  This file contains helper functions shared by the phyloFlash scripts
#


package PhyloFlash;
use strict;
use Exporter qw(import);
use Time::Piece;

our $VERSION     = 2.00;
our @ISA         = qw(Exporter);
our @EXPORT      = qw(
  get_cpus
  msg
  file_is_newer
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

sub get_cpus {
    my $cpus = `grep -c -P '^processor\\s+:' /proc/cpuinfo`;
    chomp($cpus);
    return $cpus;
}

# log a message to STDERR with time stamp prefix
sub msg {
  my $t = localtime;
  my $line = "[".$t->hms."] @_\n";
  print STDERR $line;
}

# check of arg1 was modified more recently than arg2
sub file_is_newer {
    my $file1 = shift;
    my $file2 = shift;

    return (stat($file1))[9] > (stat($file2))[9];
}

# open a filehandle
# dies with message on error
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
        or die "Failed to $msg '$fname': $!";
}

# escape array or scalar for CSV output
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

# add a hash of tools
sub require_tools {
    my %arg = @_;
    %progs = (%progs, %arg);
}

# verify that all required tools are available and locate their paths
sub check_environment {
    my $error = 0;
    foreach my $prog (keys %progs) {
        my $progname = $progs{$prog};
        if ($progs{$prog} = can_run($progname)) {
            msg("Using $prog found at \"".$progs{$prog}."\".");
        } else {
            $error = 1;
        }
    }
    if ($error == 1) {
        msg("Unable to find all required tools. These are missing:");
        foreach my $prog (keys %progs) {
            msg("  ".$prog) if (!defined($progs{$prog}));
        }
        die "Please make sure these are installed and in your PATH.\n\n";
    }
}

# run a tool
#   prog:         the name of the tool (see %progs at the top)
#   args:         command line arguments to be passed
#   redir_stdout: file or file descriptor to redirect stdout to
#   redir_stdin:  same for stderr
#
# does not return if the tool failed
sub run_prog {
    my ($prog, $args, $redir_stdout, $redir_stderr) = @_;

    if (not exists $progs{$prog}) {
        msg("trying to launch unknown tool \"$prog\". pls add to %progs");
        $progs{$prog} = can_run($prog) or die "Failed to find $prog";
    }
    my $cmd = $progs{$prog}." ".$args;
    $cmd .= " >".$redir_stdout if ($redir_stdout);
    $cmd .= " 2>".$redir_stderr if ($redir_stderr);

    msg("executing [$cmd]");
    system($cmd) == 0
        or die "Couldn't launch [$cmd]: $!/$?";

    # FIXME: print tail of stderr if redirected
}

# compute md5 sum on local file
# @parm   file: path to name of file
# @return       md5 sum as hex
sub file_md5 {
    my $file = shift;
    my $fh;
    open($fh, "<", $file)
        or die "Unable to open $file.";
    my $ctx = Digest::MD5->new;
    $ctx->addfile($fh);
    close($fh);
    return $ctx->hexdigest;
}

# fetch file from ftp into string
# @parm   ftp: Net:FTP object
# @parm   pat: string pattern in current path
# @return      content of file as string
sub ftp_read_var {
    my $ftp = shift;
    my $pat = shift;
    my $out = "";

    my $files = $ftp->ls($pat) or return "";
    return "" if (@$files == 0);

    my $fh;
    open($fh, '>', \$out);
    my $file = shift(@$files);
    $ftp->get($file, $fh)
        or die "Could not download $file", $ftp->message;
    close($fh);

    return $out;
}


# download a file from FTP
#
# Existing file is overwritten if
#  - byte sizes do not match
#  - *.md5 exists on server and md5 does not match
#
# Shows license promit if LICENSE.txt exists on server
#
# @param path: "ftp.somedomain.tld/some/dirs/pa.*ttern"
# @return      filename downloaded into CWD
#
sub file_download {
    my ($server,$path,$pat) = ($_[0] =~ /([^\/]*)(\/.*\/)([^\/]*)/);

    msg("  Connecting to $server");
    my $ftp = new Net::FTP($server,(
                               Passive => 1,
                               Debug => 0,
                               Timeout => 600
                           ))
        or die "Could not connect to $server";
    $ftp->login("anonymous", "-phyloFlash")
        or die "Could not login to $server:", $ftp->message;
    $ftp->binary();
    $ftp->pasv();

    msg("  Finding $path$pat");
    $ftp->cwd($path)
        or die "Could not enter path '\$path\': ", $ftp->message;
    my $files = $ftp->ls($pat)
        or die "Could not list files matching \'$pat'\ in \'$path\': ",
               $ftp->message;
    die "No files found?!"
        if (@$files == 0);
    msg("  Multiple files found?! Using first of ".join(@$files))
        if (@$files > 1);
    my $file = shift(@$files);
    my $file_size = $ftp->size($file)
        or die "Could not get file size:", $ftp->message;
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
        or die "Failed to download $file:", $ftp->message;
    print STDERR "\n";

    die "File size mismatch?!"
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

    msg("  MD5 sum mismatch: '$file_md5' != '$local_md5'");
    die;
}

# copies a FASTA file filtering out some sequences
# @param source: source file
# @param dest:   destination file
# @param accs:   accessions to be excluded
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

# create non-redundant DB using vsearch cluster_fast
# @param src: source FASTA
# @param dst: destination FASTA
# @param id:  identity threshold for clustering
sub cluster {
    my ($src, $dst, $id) = @_;
    msg("clustering database");
    run_prog("vsearch",
             "  --cluster_fast $src "
             . "--id $id "
             . "--centroids $dst "
             . "--notrunclabels ");
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

# normalize FASTA
# - remove alignment characters "." and ","
# - uppercase
# - DNAify (U->T)
# - replace IUPAC coded ambiguous bases with random matching base
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

sub run_stage {
    my %opts = @_;
}


{
    package Timer;
    use strict;
    use Time::Piece;
    use Time::Seconds;
    sub new {
        my ($class) = @_;
        my $self = bless {}, $class;
        $self->{time} = localtime;
        return $self;
    }
    sub minutes {
        my ($self) = @_;
        my $endtime = localtime;
        my $diff = $endtime - $self->{time};
        return sprintf "%.2f minutes", $diff->minutes;
    }
}


1; # keep require happy
