#!/usr/bin/env python

"""
This program currently does not do anything which grep could not do. 
This returns prints nothing, but it returns (using an exit code)
whether or not a PDB file contains any secondary structure 
information. If it doesn't it usually means that the original authors of 
the PDB file did not indicate where the they are.

This program was not meant to be invoked explicitly by the user.
(In the future, perhaps I will eliminate it.)
"""

import sys

def main():
    secondary_str_found = False;
    for line in sys.stdin:
        if ((line[0:6] == "HELIX ") or
            (line[0:6] == "SHEET ") or
            (line[0:5] == "TURN ")):
            secondary_str_found = True;

    if (secondary_str_found):
        exit(0) #normal termination indicates a helix/sheet/turn was found
    else:
        exit(1) #non-zero (abnormal termination) exit code indicates this
                #PDB file does not contain any secondary structure information

if __name__ == "__main__":
    main()
