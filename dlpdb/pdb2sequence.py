#!/usr/bin/env python

"""
 Usage:

    pdb2sequence.py < PDB_FILE

 To select a limited
    pdb2sequence.py ChainID1 SeqNum1 ICode1 ChainID2 SeqNum2 ICode2 < PDB_FILE

 This program extracts the sequence of residue types from a PDB file, 
 and prints it the the standard out (all on one line).
 Each residue is converted to a 1-letter amino-acid code.
 The 20 standard amino-acids, and 5 nucleic acid types are supported.
 Unknown residue types are denoted 'x'.
 (To add additional residue types, edit the code and add entries to the
 "three2one" lookup table.)

 WARNING: Residues represented as HETATM records will be ignored (skipped)!

     Optional interval arguments:

 Users can limit the PDB file to residues which lie within an interval.

 Intervals are pairs of amino acids signifying the first and last
 residue in the interval.
 In order to specify any amino acid in a PDB file, you must provide 3
 identifiers:
   the ChainID a single letter specifying
   the SeqNum an integer indicating the location within that chain
   the ICode ("insert code" usually 0.  I don't know why this number is 
              necessary. ..For the record, I never loved the PDB file format.
              BUT THIS PROGRAM EXPECTS AN ICode, SO YOU MUST SUPPLY ONE.)
 All 3 identifiers are needed for both the starting and ending residues.
   ---->
 SO, this program expects either 0 or 6 arguments: 
    (to be read from the command line)
 ChainID_first SeqNum_first ICode_first ChainID_last SeqNum_last ICode_last
"""

import sys
from operator import attrgetter

try:
    from .resid import *
except ImportError:
    from resid import *


def main():
    if len(sys.argv) == 1:
        use_all_residues = True
    elif len(sys.argv) == 7:
        use_all_residues = False
        first = ResID(sys.argv[1], int(sys.argv[2]), sys.argv[3])
        last  = ResID(sys.argv[4], int(sys.argv[5]), sys.argv[6])
    else:
        sys.stderr.write("Error: This program requires either 0 or 6 arguments.\n"
                         "       By default, the the sequence is extracted from the entire PDB file.\n"
                         "       In that case, no arguments are required.\n"
                         "       Alternately, you can limit the selection to a single interval of\n"
                         "       residues from one of the chains in the PDB file.\n"
                         "       To specify an interval, you must passing 6 arguments to this program.\n"
                         "       This program requires a pair of residues to designate the first and\n"
                         "       last members of the interval.  Each residue requires 3 identifiers.\n"
                         "       Consequently the six arguments needed are:\n"
                         "ChainID_first SeqNum_first ICode_first ChainID_last SeqNum_last ICode_last\n")
        exit(-1)

    three2one = {\
        'GLY':'G',\
        'ALA':'A',\
        'SER':'S',\
        'CYS':'C',\
        'VAL':'V',\
            \
        'THR':'T',\
        'ILE':'I',\
        'PRO':'P',\
        'MET':'M',\
        'ASP':'D',\
            \
        'ASN':'N',\
        'LEU':'L',\
        'LYS':'K',\
        'GLU':'E',\
        'GLN':'Q',\
            \
        'ARG':'R',\
        'HIS':'H',\
        'PHE':'F',\
        'TYR':'Y',\
        'TRP':'W',
            \
        ' DA':'A',\
        ' DC':'C',\
        ' DG':'G',\
        ' DT':'T',\
            \
        '  A':'A',\
        '  U':'U', \
        '  G':'G',\
        '  C':'C'}



    resID2type3 = {}

    for line in sys.stdin:
        #if line[0:6] in ("ATOM  ", "HETATM"):
        if line[0:6] in "ATOM  ":
            #atomID    = int(line[6:11])
            #atomType  = line[12:16]
            #altLoc    = line[16:17]
            resType   = line[17:20]
            chainID   = line[21:22]
            resSeq    = line[22:26]
            iCode     = line[26:27]
            resID = ResID(chainID, int(resSeq), iCode)
            if (use_all_residues or ((first<=resID) and (resID<=last))):
                resID2type3[resID] = resType

    # Extract an (unordered) list of the resIDs of the residues in the sequence
    resIDs = [resID for resID in resID2type3]

    # Residues in PDB files are often not listed in order.
    # Consequently, we must sort the list by chainID, seqNum, and finnaly iCode:
    sequence_of_resIDs = sorted(resIDs, key=attrgetter('chainID','seqNum','iCode'))


    # Now loop through the sequence of resIDs, and 
    # lookup the residue type for each residue (ie. 'ALA', 'PRO, 'GLY', 'VAL' etc...)
    # and convert it to a 1-letter residue name ('A', 'P', 'G', 'V'...)
    # and store the sequence in a string 'APGV...'

    sequence = ''
    for resID in sequence_of_resIDs:
        three_letter_code = resID2type3[resID]
        if three_letter_code in three2one:
            c = three2one[three_letter_code]
        else:
            c = 'x' # default character used for unknown/non-standard residues
        sequence += c

    # Now, finally print out the sequence:
    sys.stdout.write(sequence+'\n')


if __name__ == "__main__":
    main()
