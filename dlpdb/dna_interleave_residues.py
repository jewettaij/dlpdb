#!/usr/bin/env python
"""
Explanation:

This program reads 2 pdb files (each of which should contain the same number
of residues).  It sorts them both according to the 3-part ResID identifier
(chainID, resSeq, iCode).  Then it prints out a new PDB file, alternately
printing 1 residue from 1 PDB file, with a residues from the other file,
(listed IN REVERSE ORDER).  The new chainID is set to 'X', and the resSeq
numbers will increase consecutively from 1 to 2N, where N is the number of
residues in either file.  This program was intended to be used to make it
easier to associate each nucleotide in double-stranded DNA with it's partner.

The user has to prepare the two PDB files to make sure that they contain two
different DNA strands which completely hybridize with each other (no overhangs),
AND the order must be reversed.  (Because the two strands in DNA are physically
oriented in opposite directions.)

"""

import sys
from operator import attrgetter
from collections import defaultdict

try:
    from .resid import *
except ImportError:
    from resid import *


g_program_name = __file__.split('/')[-1]
g_date_str     = '2017-4-18'
g_version_str  = '0.9.0'


usage_string = "Usage:\n\n  "+g_program_name+" file1.pdb file2.pdb > file_merged.pdb\n\n"

usage_string += __doc__



def ResID2Text_from_pdb_file(file_name):

    """
    This function reads a pdb file, and creates a dictionary containing
    the text for each residue, indexed by it's identifier.  Identifiers
    in PDB files (ResIDs) have 3 parts (chainID,resSeq,iCode).

    """

    pdb_file = open(file_name, 'r')

    resID2text = defaultdict(list)
    
    chain_text = {}

    for line in pdb_file:
        #if line[0:6] in ("ATOM  ", "HETATM"):
        if line[0:6] in "ATOM  ":
            #atomID    = int(line[6:11])
            #atomType  = line[12:16]
            #altLoc    = line[16:17]
            #resType   = line[17:20]
            chainID   = line[21:22]
            resSeq    = line[22:26]
            iCode     = line[26:27]
            resID = ResID(chainID, int(resSeq), iCode)
            ##xyz_str=line[30:54] <- not safe. spaces not always between numbers
            #x_str     = line[30:38].strip()
            #y_str     = line[38:46].strip()
            #z_str     = line[46:54].strip()

            resID2text[resID].append(line)

    pdb_file.close()
    return resID2text




def SortedResidueText_from_pdb_file(file_name):

    """
    This function reads a pdb file, and creates a list of lists of strings.
    Each sub-list of strings contains the text for a particular residue.
    The residues are sorted according to their 3-part ResID identifiers
    (see above).

    """

    resID2text = ResID2Text_from_pdb_file(file_name)

    # Extract an (unordered) list of the resIDs of the residues in the sequence
    resIDs = [resID for resID in resID2text]

    # Residues in PDB files are often not listed in order.
    # Consequently, we must sort the list by chainID, seqNum, and finnaly iCode:
    resIDs_sorted = sorted(resIDs, 
                           key=attrgetter('chainID','seqNum','iCode'))

    index2text = []
    for resID in resIDs_sorted:
        index2text.append(resID2text[resID])

    return index2text


def main():
    if len(sys.argv) != 3:
        sys.stderr.write("Error: Expected two arguments (two pdb files)\n")
        sys.stderr.write(usage_string)
        exit(1)

    sys.stderr.write(g_program_name+' v'+g_version_str+' '+g_date_str+'\n')

    sys.stderr.write('  Reading file \"'+sys.argv[1]+'\"\n')
    strand1 = SortedResidueText_from_pdb_file(sys.argv[1])
    sys.stderr.write('  Reading file \"'+sys.argv[2]+'\"\n')
    strand2 = SortedResidueText_from_pdb_file(sys.argv[2])
    if (len(strand1) != len(strand2)):
        sys.stderr.write("Error: The two PDB files must have the same length.\n")
        exit(1)


    N = len(strand1)
    new_chainID=' '
    i = 1;
    for n in range(0, N):
        for line in strand1[n]:
            new_resSeq = str(i).rjust(4)
            new_line = (line[0:21] + new_chainID + new_resSeq + line[26:])
            # (leave the iCode in line[26:27] alone)
            sys.stdout.write(new_line)
        i += 1
        for line in strand2[N-n-1]:
            new_resSeq = str(i).rjust(4)
            new_line = (line[0:21] + new_chainID + new_resSeq + line[26:])
            # (leave the iCode in line[26:27] alone)
            sys.stdout.write(new_line)
        i += 1




if __name__ == "__main__":

    main()
