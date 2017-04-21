#!/usr/bin/env python

"""
This program finds the location of each short SHEET interval in a PDB file.
(For each sheet interval, it prints out the starting and ending residue
 using a format which is recognized by the "select_interval.py" script.
 That script is typically used for extracting excerpts of a PDB file.)
This script is not really intended to be used on it's own.
"""

import sys

def main():
    for line in sys.stdin:
        if (line[0:6] == "SHEET "):
            initChainID = line[21:22]
            initSeqNum  = int(line[22:26])
            initICode   = line[26:27]
            endChainID  = line[32:33]
            endSeqNum   = int(line[33:37])
            endICode    = line[37:38]
            sys.stdout.write("\""+initChainID+"\" "+str(initSeqNum)+" \""+initICode+
                             "\"  \""+
                             endChainID+"\" "+str(endSeqNum)+" \""+endICode+"\"\n")


if __name__ == "__main__":
    main()
