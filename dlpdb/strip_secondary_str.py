#!/usr/bin/env python

"""
This removes lines from pdb files which define 
where HELIX, SHEET, or TURNs are located.
"""

import sys

def main():
    for line in sys.stdin:
        if ((line[0:6] != "HELIX ") and 
            (line[0:6] != "SHEET ") and
            (line[0:5] != "TURN ")):
            sys.stdout.write(line)

if __name__ == "__main__":
    main()
