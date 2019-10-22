#!/usr/bin/env python

"""
 This is a crude script which prints out the coordinates for the average 
 position of each residue in a PDB file (three numbers per line).
 (Note: This program does not calculate the center-of-mass.  All atoms 
        are given the same weight when calculating the average position.)

 Typical usage: 

 pdb2coords_ave.py < PDB_FILE

 Exceptions:
   Backbone atoms (other than alpha-carbon " CA " atoms) are ignored.
   Glycine residues are also ignored.

 These settings can not be changed without editing the code for this script.
 (unlike other scripts whose behavior can be changed by command-line arguments)

 Output:
   The number of lines of output should match the number of residues.
   For glycines (and residues lacking elligible atoms) a blank line is printed.
   For all other residues, the x, y, and z coordinates of the average
   position are printed (3 numbers), one line per residue.

"""

import sys
from operator import attrgetter

try:
    from .resid import *
except ImportError:
    from resid import *


# --- THE FOLLOWING FEATURES (interval restrictions) may be removed later:--
#
# The entire PDB file does not have to be used.  Instead, the user can pass 
# 6 arguments denoting the first and last residue of an interval in the chain.
# In order to specify any amino acid in a PDB file, you must provide 3
# identifiers:
#   the ChainID a single letter specifying
#   the SeqNum an integer indicating the location within that chain
#   the ICode (an insert code, usually " ")
#
# All 3 identifiers are needed for both the starting and ending residues.
# Consequently this program expects 6 arguments: 
#    (to be read from the command line)
# ChainID_first SeqNum_first ICode_first ChainID_last SeqNum_last ICode_last


# The following atoms are not included in the average:
ignore_these_atoms    = [' N  ', ' H  ', ' C  ', ' O  ', ' CA ']
# The following residues are completely excluded:
ignore_these_residues = set(['GLY', 'PRO'])



def main():
    use_all_residues = True

    if len(sys.argv) > 1:
        i = 1
        while i < len(sys.argv):

            if (len(sys.argv[i]) == 3):
                # Add the string to the list of amino acids we want to ignore
                resType = sys.argv[i]
                ignore_these_residues.add(resType)
                i += 1

            else: # if not is_digit(sys.argv[i][0])

                # If the next argument is a number, then interpret this number
                # as a residue sequence number
                if len(sys.argv) < i+6:
                    sys.stderr.write("Error: Not enough arguments or argument type error:\n"
                                     "       Offending arguemt #"+str(i)+": \""+sys.argv[i]+"\"\n")

                    if len(sys.argv[i]) == 1:
                        sys.stderr.write("\n"
                                         "       Is argument \""+sys.argv[i]+"\" a chain ID letter?\n"
                                         "       (Chain IDs are only passed as arguments when you want to limit the\n"
                                         "        residues considered to within an interval in the PDB file.\n"
                                         "        To specify an interval, you must provide 5 more aruments.  See below.)\n"
                                         "\n"
                                         "       Note: one-letter residue codes (eg A,V,L,...) are not considered\n"
                                         "             valid residue types by this program.\n"
                                         "             You must specify a 3-letter equivalent (eg ALA, VAL, LYS,...).\n"
                                         "             One-letter arguments are interpreted as chain-IDs. (See below.)\n"
                                         "       \n"
                                         "       By default, the the sequence is extracted from the entire PDB file.\n"
                                         "       In that case, no arguments are required.\n"
                                         "       Alternately, you can limit the selection to a single interval of\n"
                                         "       residues from one of the chains in the PDB file.\n"
                                         "       To specify an interval, you must passing 6 arguments to this program.\n"
                                         "       This program requires a pair of residues to designate the first and\n"
                                         "       last members of the interval.  Each residue requires 3 identifiers.\n"
                                         "       Consequently the six arguments needed are:\n"
                                         "\n"
                                         "ChainID_first SeqNum_first ICode_first ChainID_last SeqNum_last ICode_last\n"
                                         "\n")
                    else:
                        sys.stderr.write("\n"
                                         "Note: PDB files are not passed as arguments\n"
                                         "      but are read from the standard input, for example using the notation\n"
                                         "\n"
                                         "      "+sys.argv[0]+" < file.pdb\n")
                    exit(-1)
                else:
                    use_all_residues = False
                    first = ResID(sys.argv[i], int(sys.argv[i+1]), sys.argv[i+2])
                    last  = ResID(sys.argv[i+3], int(sys.argv[i+4]), sys.argv[i+5])
                    #sys.stderr.write('  Interval selected: (\"'+first.chainID+'\", '+str(first.seqNum)+', \"'+first.iCode+'\") ... (\"'+last.chainID+'\", '+str(last.seqNum)+', \"'+last.iCode+'\")\n')
                    i += 6

    resID2pos = {}

    for line in sys.stdin:
        if (line[0:6] == "ATOM  "):
            #atomID    = int(line[6:11])
            atomType  = line[12:16]
            #altLoc    = line[16:17]
            resType   = line[17:20]
            chainID   = line[21:22]
            resSeq    = line[22:26]
            iCode     = line[26:27]
            resID = ResID(chainID, int(resSeq), iCode)
            #xyz_str   = line[30:54]<-not safe. spaces not always in between numbers
            x_str     = line[30:38].strip()
            y_str     = line[38:46].strip()
            z_str     = line[46:54].strip()
            x = float(x_str)
            y = float(y_str)
            z = float(z_str)
            if (use_all_residues or ((first<=resID) and (resID<=last))):
                if resID not in resID2pos:
                    resID2pos[resID] = []

                # Ignore atoms on the backbone (other than CA), and also 
                # ignore all atoms which are not heavy atoms (hydrogen atoms).
                if ((resType not in ignore_these_residues) and
                    (atomType not in ignore_these_atoms)):
                    resID2pos[resID].append([x, y, z, atomType])

    # Extract an (unordered) list of the resIDs of the residues in the sequence
    resIDs = [resID for resID in resID2pos]

    # Residues in PDB files are often not listed in order.
    # Consequently, we must sort the list by chainID, seqNum, and finnaly iCode:
    sequence_of_resIDs = sorted(resIDs, key=attrgetter('chainID','seqNum','iCode'))

    # Now loop through the sequence of resIDs, and calculate 
    # the average position of the atoms in that residue.

    for resID in sequence_of_resIDs:
        coords = resID2pos[resID]
        if coords == []:
            sys.stdout.write('\n')
        else:
            xyz_tot = [0.0 for d in range(0,3)]
            num_atoms_in_res = len(coords)
            for i in range(0, num_atoms_in_res):

                for d in range(0,3):
                    xyz_tot[d] += coords[i][d]

            x_ave = xyz_tot[0] / num_atoms_in_res
            y_ave = xyz_tot[1] / num_atoms_in_res
            z_ave = xyz_tot[2] / num_atoms_in_res


            x_str = str(x_ave)
            y_str = str(y_ave)
            z_str = str(z_ave)

            # Alternately, I could have used:
            # x_str = "%5.3" % x_ave
            # y_str = "%5.3" % y_ave
            # z_str = "%5.3" % z_ave

            sys.stdout.write(x_str+' '+y_str+' '+z_str+'\n')


if __name__ == "__main__":
    main()
