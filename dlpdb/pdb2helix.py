#!/usr/bin/env python

"""
This program finds the location of each HELIX in a PDB file.
(For each helix, it prints out the starting and ending residue
 using a format which is recognized by the "select_interval.py" script.
 That script is typically used for extracting excerpts of a PDB file.)
This script is not really intended to be used on it's own.
"""

import sys

def main():
    for line in sys.stdin:
        if (line[0:6] == "HELIX "):
            initChainID = line[19:20]
            initSeqNum  = int(line[21:25])
            initICode   = line[25:26]
            endChainID  = line[31:32]
            endSeqNum   = int(line[33:37])
            endICode    = line[37:38]
            sys.stdout.write("\""+initChainID+"\" "+str(initSeqNum)+" \""+initICode+
                             "\"  \""+
                             endChainID+"\" "+str(endSeqNum)+" \""+endICode+"\"\n")


if __name__ == "__main__":
    main()
