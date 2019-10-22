#!/usr/bin/env python

"""
This program extracts the information from a PDB file
which lies within an interval in the PDB file.

Usage:

select_interval.py select_interval.py A 3 " " A 10 z < 309d.pdb > 309d_A3-10.pdb

(This selects residues 3-10 from chain A of the PDB file)

Intervals are pairs of amino acids signifying the first and last
residue in the interval.
In order to specify any amino acid in a PDB file, you must provide 3
identifiers:
  the ChainID a single letter specifying (typically " " or A-Z)
  the SeqNum an integer indicating the location within that chain
  the ICode (a single character code for residue insertions.
             Usually this charcter is the space character " ".
             The smallest possible value is " ".  The largest is "z".
             For the record, I hate this file format.)

All 3 identifiers are needed for both the starting and ending residues.
Consequently this program expects 6 arguments: 
   (to be read from the command line)
ChainID_first SeqNum_first ICode_first ChainID_last SeqNum_last ICode_last
"""

import sys
# Sometimes this program pipes its output to other programs which stops reading
# the PDB file prematurely (such as when multiple MODEL records are present).
# Below we silently suppress the ugly "Broken pipe" message this generates:
import signal
signal.signal(signal.SIGPIPE, signal.SIG_DFL)

try:
    from .resid import *
except ImportError:
    from resid import *


def main():
    if len(sys.argv) != 7:
        sys.stderr.write("Error: This program requires 6 arguments.\n"
                         "       This program requires a pair of residues to designate the first and\n"
                         "       last members of the interval.  Each residue requires 3 identifiers.\n"
                         "       Consequently the six arguments needed are:\n"
                         "ChainID_first SeqNum_first ICode_first ChainID_last SeqNum_last ICode_last\n")
        exit(-1)

    print('\"'+sys.argv[1] + '\" \"' + sys.argv[2] + '\" \"'+ sys.argv[3]+'\"')
    first = ResID(sys.argv[1], int(sys.argv[2]), sys.argv[3])
    last  = ResID(sys.argv[4], int(sys.argv[5]), sys.argv[6])

    for line in sys.stdin:
        line_type = line[0:6]

        if line_type in set(["ATOM  ", "HETATM", "ANISOU", "SIGATM", "SIGUIJ"]):
            #atomID    = int(line[6:11])
            #atomType  = line[12:16]
            #altLoc    = line[16:17]
            #resType   = line[17:20]
            chainID   = line[21:22]
            seqNum    = line[22:26]
            iCode     = line[26:27]
            resID = ResID(chainID, int(seqNum), iCode)
            if (first <= resID <= last):
                sys.stdout.write(line)

        elif (line_type == "HET   "):
            #hetID     = line[7:10]
            chainID   = line[12:13]
            seqNum    = line[13:17]
            iCode     = line[17:18]
            #numHETATMs = int(line[20:25])
            #descriptor = line[30:70]
            resID = ResID(chainID, int(seqNum), iCode)
            if (first <= resID <= last):
                sys.stdout.write(line)

        elif (line_type == "HELIX "):
            initChainID = line[19:20]
            initSeqNum  = int(line[21:25])
            initICode   = line[25:26]
            initID = ResID(initChainID, int(initSeqNum), initICode)
            endChainID  = line[31:32]
            endSeqNum   = int(line[33:37])
            endICode    = line[37:38]
            endID = ResID(endChainID, int(endSeqNum), endICode)
            if (first <= initID <= last) and (first <= endID <= last):
                sys.stdout.write(line)

        elif (line_type == "SHEET "):
            initChainID = line[21:22]
            initSeqNum  = int(line[22:26])
            initICode   = line[26:27]
            initID = ResID(initChainID, int(initSeqNum), initICode)
            endChainID  = line[32:33]
            endSeqNum   = int(line[33:37])
            endICode    = line[37:38]
            endID = ResID(endChainID, int(endSeqNum), endICode)
            if (first <= initID <= last) and (first <= endID <= last):
                sys.stdout.write(line)

        elif (line_type == "TURN  "):
            initChainID = line[19:20]
            initSeqNum  = int(line[20:24])
            initICode   = line[24:25]
            initID = ResID(initChainID, int(initSeqNum), initICode)
            endChainID  = line[30:31]
            endSeqNum   = int(line[31:35])
            endICode    = line[35:36]
            endID = ResID(endChainID, int(endSeqNum), endICode)
            if (first <= initID <= last) and (first <= endID <= last):
                sys.stdout.write(line)

        elif line_type == "SEQRES":
            chainID = line[11:12]
            if (first.chainID <= chainID) and (chainID <= last.chainID):
                sys.stdout.write(line)

        elif line_type == "TER   ":
            chainID = line[21:22]
            seqNum  = int(line[22:26])
            iCode   = line[26:27]
            terRes  = ResID(chainID, seqNum, iCode)
            if (first <= resID <= last):
                sys.stdout.write(line)

        else:
            sys.stdout.write(line)

if __name__ == "__main__":
    main()
