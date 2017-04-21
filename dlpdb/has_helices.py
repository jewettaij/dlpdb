#!/usr/bin/env python

"""
This program currently does not do anything which grep could not do.
This returns prints nothing, but it returns (using an exit code)
whether or not a PDB file contains a "HELIX" record.
If it doesn't it could either be because the original authors of the 
PDB file neglected to indicate where the helices are, or because
there simply arne't any.

This program was not meant to be invoked explicitly by the user.
(In the future, perhaps I will eliminate it.)
"""

import sys

def main():
    helix_found = False;
    for line in sys.stdin:
        if (line[0:6] == "HELIX "):
            helix_found = True;

    if (helix_found):
        exit(0)  # normal termination indicates a helix was found
    else:
        exit(1)  # non-zero (abnormal termination) exit code indicates 
                 # this PDB file does not contain a helix.

if __name__ == "__main__":
    main()
