#!/usr/bin/perl -w
#   A simple script to prepare the phyloFlash databases
#   - Assumes the script is executed in a dir with write access and
#   - Dependencies: usearch v6, bbmap and Emirge
#   - This script generates approx. 5Gbytes of extracted databases

#    Copyright (C) 2014- by Harald Gruber-Vodicka mail:hgruber@mpi-bremen.de

#thanks goes to Torsten Seeman for the code to give timing informations. The code snippets for this were initially distributed in Prokka 1.7.1

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

#version 1.00 - last edit 27/08/2014

use strict;
use Time::Piece;
use Time::Seconds;
use Getopt::Long;
use Cwd;


# constants
my $syscpus = `cat /proc/cpuinfo | grep -c -P '^processor\\s+:'`;
chomp($syscpus);
my $dbsource;
my $silva_release;
my $cwd = getcwd;

my $starttime = localtime;


### MAIN ###
download_RefNR();
#download_univec();
getfilename();#always run!!!
#unpack_SSU();
#find_LSU();
#remove_LSU();
#mask_repeats();
#univec_trim();
#clustering();
bbmap_db();
#emirge_db();
finish();




### subroutines ###

sub msg {
  my $t = localtime;
  my $line = "[".$t->hms."] @_\n";
  print STDERR $line;
}

#downloads the current SSURef_NR99 file from the current silva release 
sub download_RefNR {
    msg("downloading latest SSU RefNR from www.arb-silva.de");
    my $wget_silva_cmd = "wget -nv -r -l1 -np -P ./download ftp://ftp.arb-silva.de/current/Exports/ -A '*_SSURef_Nr99_tax_silva_trunc.fasta.gz'";
    system($wget_silva_cmd) == 0
    or die "Couldn't launch [$wget_silva_cmd]: $! / $?\n Do you have wget and are you connected to the internet?";
    my $move_cmd="find ./download/ -name \"*.gz\" -exec cp '{}' ./download/ \\;";
    msg($move_cmd);
    system($move_cmd);
    system ("rm ./download/ftp.arb-silva.de -r");
}

sub download_univec{
  msg("downloading latest univec from ncbi");
    my $wget_ncbi_cmd = "wget -P ./download ftp://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec";
    system($wget_ncbi_cmd) == 0
    or die "Couldn't launch [$wget_ncbi_cmd]: $! / $?\n Do you have wget and are you connected to the internet?";
    my $move_cmd="find ./download/ -name \"UniVec\" -exec cp '{}' ./download/ \\;";
    msg($move_cmd);
    system($move_cmd);
}

#gets the current silva release number from the filename and creates output directory for the release
sub getfilename {
  my $fp = glob("./download/*SSURef_Nr99_tax_silva_trunc.fasta.gz");
  $dbsource = $fp if  $fp =~ m/SILVA.*/ ;
  msg ("you are working with $dbsource");
  if ($dbsource =~ m/SILVA_(\d+)_/) {
    $silva_release = $1;
    }
  msg ("This is silva release $silva_release");
  if (-e $silva_release) 
    {
        msg ("Silva release $silva_release directory exists");
    }
    else
    {
      msg ("creating directory for Silva release $silva_release");
      mkdir $silva_release or die "Error creating directory: $silva_release";
    }
}
  
#unpacking the downloaded SSURef_Nr file to the release folder
sub unpack_SSU {
  msg("unpacking SSU db");
  my $gunzip_cmd = "gunzip -c $dbsource > ./$silva_release/SILVA_SSU.fasta";
  system($gunzip_cmd) == 0
    or die "Couldn't launch [$gunzip_cmd]: $! / $?\n Do you have gunzip installed?";
}

#searching for LSU genes in the SSU RefNR using a modified barrnap version that only utilizes LSU hmm profiles
#barrnap was changed to use a different db folder that only contains the LSU profiles
sub find_LSU {
  msg("searching for LSU contamination in SSU RefNR");
  foreach ('bac', 'arch', 'euk') {
    my $barrnap_cmd = "./barrnap-HGV/bin/barrnap_HGV --kingdom $_ --threads $syscpus --evalue 1e-50 --reject 0.01 ./$silva_release/SILVA_SSU.fasta > tmp.barrnap_hits.$_.gff 2>tmp.barrnap_hits.$_.barrnap.out";
    system($barrnap_cmd) == 0
    or die "Couldn't launch [$barrnap_cmd]: $! / $?";
    my $egrep_cmd = "egrep '23S_rRNA\|28S_rRNA' tmp.barrnap_hits.$_.gff \>\> tmp.LSU_hits.collection.$silva_release.gff";
    system($egrep_cmd);
  }
}  

#all sequences with LSU predicitions are removed from the database
sub remove_LSU {
  msg("removing LSU contamination in SSU RefNR");
  my $cut_cmd = "cut -f -1 tmp.LSU_hits.collection.$silva_release.gff | sort -u > tmp.LSU_contamination_list.$silva_release.txt";
    system($cut_cmd) == 0
    or die "Couldn't launch [$cut_cmd]: $! / $?";
  my $fasomerec_cmd = "faSomeRecords ./$silva_release/SILVA_SSU.fasta tmp.LSU_contamination_list.$silva_release.txt ./$silva_release/SILVA_SSU.noLSU.fasta -exclude";
      system($fasomerec_cmd) == 0
    or die "Couldn't launch [$fasomerec_cmd]: $! / $?";
    unlink "./$silva_release/SILVA_SSU.fasta";
  
}


#bbmask masks low entropy regions and repeats in the fasta file
sub mask_repeats {
  msg("masking low entropy regions in SSU RefNR");
  my $bbmask_cmd = "bbmask.sh overwrite=t -Xmx4g in=./$silva_release/SILVA_SSU.noLSU.fasta out=./$silva_release/SILVA_SSU.noLSU.masked.fasta minkr=4 maxkr=8 mr=t minlen=20 minke=4 maxke=8 fastawrap=0";
    system($bbmask_cmd) == 0
    or die "Couldn't launch [$bbmask_cmd]: $! / $?";
    unlink "./$silva_release/SILVA_SSU.noLSU.fasta";
}

#db is screened against UniVec and hit bases are trimmed using bbduk
sub univec_trim {
  msg("removing UniVec contamination in SSU RefNR");
  my $bbduk_cmd = "bbduk.sh ref=./download/UniVec fastawrap=0 ktrim=r ow=t minlength=800 mink=11 hdist=1 in=./$silva_release/SILVA_SSU.noLSU.masked.fasta out=./$silva_release/SILVA_SSU.noLSU.masked.trimmed.fasta stats=tmp.UniVec_contamination_stats.$silva_release.txt";
    system($bbduk_cmd) == 0
    or die "Couldn't launch [$bbduk_cmd]: $! / $?";
    unlink "./$silva_release/SILVA_SSU.noLSU.masked.fasta";
}

# the cleaned, masked and trimmed databases are clustered
# 1) at 99id with full labels for bbmap
# 2) at 96id for emirge
sub clustering {
  msg("clustering database");
  my $cluster_99 = "vsearch --cluster_fast ./$silva_release/SILVA_SSU.noLSU.masked.trimmed.fasta --id 0.99 --centroids ./$silva_release/SILVA_SSU.noLSU.masked.trimmed.NR99.fasta --notrunclabels";
    system($cluster_99) == 0
    or die "Couldn't launch [$cluster_99]: $! / $?";
  my $cluster_96 = "vsearch --cluster_fast ./$silva_release/SILVA_SSU.noLSU.masked.trimmed.NR99.fasta --id 0.96 --centroids ./$silva_release/SILVA_SSU.noLSU.masked.trimmed.NR96.fasta";
    system($cluster_96) == 0
    or die "Couldn't launch [$cluster_96]: $! / $?";
}

#the NR99 file is fixed and a reference file is generated for bbmap
sub bbmap_db {
  msg("preparing bbmap database");
  my $fix_cmd = "python2 fix_nonstandard_chars_X.py \< ./$silva_release/SILVA_SSU.noLSU.masked.trimmed.NR99.fasta \> ./$silva_release/SILVA_SSU.noLSU.masked.trimmed.NR99.fixed.fasta";
    system($fix_cmd) == 0
    or die "Couldn't launch [$fix_cmd]: $! / $?\n Do you have python2 available?";
    msg("creating bbmap reference"); 
  my $bbmap_ref_cmd = "bbmap.sh -Xmx4g ref=./$silva_release/SILVA_SSU.noLSU.masked.trimmed.NR99.fixed.fasta path=./$silva_release/";
  system($bbmap_ref_cmd) == 0
    or die "Couldn't launch [$bbmap_ref_cmd]: $! / $?\n Do you have bbmap installed?";
}

#the NR96 file is fixed and a bowtie index is built    
sub emirge_db{
  msg("preparing emirge database");
  my $fix_cmd = "python2 fix_nonstandard_chars_X.py \< ./$silva_release/SILVA_SSU.noLSU.masked.trimmed.NR96.fasta \> ./$silva_release/SILVA_SSU.noLSU.masked.trimmed.NR96.fixed.fasta";
    system($fix_cmd) == 0
    or die "Couldn't launch [$fix_cmd]: $! / $?\n Do you have python2 available?";
 
  msg("creating emirge bowtie index");    
  my $bowtie_index_cmd = "bowtie-build ./$silva_release/SILVA_SSU.noLSU.masked.trimmed.NR96.fixed.fasta ./$silva_release/SILVA_SSU.noLSU.masked.trimmed.NR96.fixed.bt -q";
    system($bowtie_index_cmd) == 0
    or die "Couldn't launch [$bowtie_index_cmd]: $! / $?\n Do you have bowtie (without the 2) installed?";
  
}

#------------------------------------------------- final timing stats and goodbye
sub finish {
  my $endtime = localtime;
  my $walltime = $endtime - $starttime;
  my $prettytime = sprintf "%.2f minutes", $walltime->minutes;
  msg("Walltime used: $prettytime");;
  msg("The databases for Silva release $silva_release are ready for phyloFlash\n
    Please provide your location of
    the databases with -dbhome: /path/to/your/databases/
    or change script line XXX accordingly");#Fixme
}