#!/usr/bin/env python2
"""
USAGE: fix_nonstandard_chars.py < input.fasta > output.fasta

This script is designed to replace all bases not in {ACTG} in a fasta
file with a semi-random base within {ATCG}.  This is necessary for
bowtie, which will not map reads to bases with non-standard chars.
From the bowtie manual
(http://bowtie-bio.sourceforge.net/manual.shtml#the-bowtie-aligner):

"Alignments involving one or more ambiguous reference characters (N,
-, R, Y, etc.) are considered invalid by Bowtie. This is true only for
ambiguous characters in the reference; alignments involving ambiguous
characters in the read are legal, subject to the alignment
policy. Ambiguous characters in the read mismatch all other
characters. Alignments that"fall off" the reference sequence are not considered valid."

Where IUPAC codes specifying only some of the possible bases are
given, one from the appropriate subset of bases will replace the
ambiguous base.  Otherwise (N or unknown char), bases are chosen at
random from {ATCG}.



"""

import sys
from Bio import SeqIO, Seq
from random import choice


def main(inf = None, outf = None):
    replace_dict = {"-": "",
                    ".": "",
                    "X": "ACGT",		
                    "A": "A",	
                    "C": "C",
                    "G": "G",
                    "T": "T",
                    "U": "T",
                    "R": "AG", #	Purine (A or G)
                    "Y": "CT", #	Pyrimidine (C, T, or U)
                    "M": "CA", #	C or A
                    "K": "TG", #	T, U, or G
                    "W": "TA", #	T, U, or A
                    "S": "CG", #	C or G
                    "B": "CTG",#	C, T, U, or G (not A)
                    "D": "ATG",#	A, T, U, or G (not C)
                    "H": "ATC",#	A, T, U, or C (not G)
                    "V": "ACG",#	A, C, or G (not T, not U)
                    "N": "ACTG" }#	Any base (A, C, G, T, or U)
    if inf is None:
        inf = sys.stdin
    if outf is None:
        outf = sys.stdout

    joinstring = "".join
    _str = str
    _Seq = Seq.Seq
    outf_write = outf.write
    _choice = choice
    get = replace_dict.get

    for record in SeqIO.parse(inf, 'fasta'):
        newseqstring = joinstring([_choice(get(x, "ACTG")) for x in _str(record.seq.upper())])
        record.seq =_Seq(newseqstring)
        outf_write(record.format("fasta"))

if __name__ == "__main__":
    main()

