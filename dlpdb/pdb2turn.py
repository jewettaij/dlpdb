#!/usr/bin/env python

# This program finds the location of each TURN in a PDB file.
# (For each turn, it prints out the starting and ending residue
#  using a format which is recognized by the "select_interval.py" script.
#  That script is typically used for extracting excerpts of a PDB file.)
# This script is not really intended to be used on it's own.

import sys

def main():
    for line in sys.stdin:
        if (line[0:6] == "TURN  "):
            initChainID = line[19:20]
            initSeqNum  = int(line[20:24])
            initICode   = line[24:25]
            endChainID  = line[30:31]
            endSeqNum   = int(line[31:35])
            endICode    = line[35:36]
            sys.stdout.write("\""+initChainID+"\" "+str(initSeqNum)+" \""+initICode+
                             "\"  \""+
                             endChainID+"\" "+str(endSeqNum)+" \""+endICode+"\"\n")


if __name__ == "__main__":
    main()
