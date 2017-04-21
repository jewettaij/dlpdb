#!/usr/bin/env python

"""
This program currently does not do anything which grep could not do.
This returns prints nothing, but it returns (using an exit code)
whether or not a PDB file contains a "TURN" record.
If it doesn't it could either be because the original authors of the 
PDB file neglected to indicate where the turns are, or because
there simply arne't any.

This program was not meant to be invoked explicitly by the user.
(In the future, perhaps I will eliminate it.)
"""

import sys

def main():
    turn_found = False;
    for line in sys.stdin:
        if (line[0:5] == "TURN "):
            turn_found = True;

    if (turn_found):
        exit(0)  # normal termination indicates a turn was found
    else:
        exit(1)  # non-zero (abnormal termination) exit code indicates 
                 # this PDB file does not contain a turn.

if __name__ == "__main__":
    main()
