#!/usr/bin/env python

"""
This program reads a 6-column numeric text file and prints a list of distances.
The 6 numbers on each line represent the x,y,z coordinates of 2 atoms.
For each line this program computes the distance between these atom pairs,
and prints it to the standard out.
(When a line contains the wrong number of numbers, this program does not crash.
 Instead an impossible value, -1.0 is printed.)

"""


import sys
from math import sqrt
# Sometimes this program pipes its output to other programs which halt early.
# Below we silently suppress the ugly "Broken pipe" message this generates:
import signal
signal.signal(signal.SIGPIPE, signal.SIG_DFL)


def length_v(r):
    lsqd = 0.0
    for d in range(0,len(r)):
        lsqd += r[d]*r[d]
    return sqrt(lsqd)


def main():
    if (len(sys.argv) > 3):
        sys.stderr.write('Error (coords2distances): number of arguments should not exceed 2.\n'\
                          '    (The two arguments correspond to the number of lines of\n'\
                          '     text to omit from the beginning and end of the file, respectively.)\n'\
                          '     If one argument is passed, then both are assumed to be the same.\n'\
                          '     If no argument is passed, then by default, no data is ignored.\nExiting...\n\n')
        sys.exit(-1)

    # NOTE: The "truncate" arguments are not really supported any more.  Instead
    #       use other scripts to post-process the results printed by this program.
    elif (len(sys.argv) == 3):
        truncate_a = int(sys.argv[1])
        truncate_b = int(sys.argv[2])
    elif (len(sys.argv) == 2):
        truncate_a = int(sys.argv[1])
        truncate_b = truncate_a
    else:
        truncate_a = 0
        truncate_b = 0


    coords_list = []
    # Read the file
    for line in sys.stdin:
        line = line.strip()
        # Each line should contain a list of 6 numbers separated by whitespace.
        # If so, store the 6 numbers in a list variable (named xyz), and append
        # it to the list of coordinates.
        # However some lines might also be blank, in which case we append the
        # empty list [] to the list of coordinates.
        if line == '':
            coords = []
        else:
            # Each line should contain a list of 3 numbers separated by whitespace.
            coords = list(map(float, line.split()))
            if len(coords) != 2*3:
                sys.stderr.write('Error(coords2distances.py):\n'+'Each line should either contain 6 numbers or be blank.\n')
                sys.exit(-1)
        coords_list.append(coords)


    # Truncate the data we don't want.
    # (Why?  The residues at the beginning and ending of helices 
    #  are less trustworthy then the residues in the middle.)
    coords_list = coords_list[truncate_a:len(coords_list)-truncate_b]

    N = len(coords_list)
    for i in range(0,N):
        if len(coords_list[i]) == 2*3:
            r10 = [0.0, 0.0, 0.0]
            for d in range(0,3):
                r10[d] = coords_list[i][3*1+d] - coords_list[i][3*0+d]
            l10 = length_v(r10)

            sys.stdout.write(str(l10)+'\n')

        else:
            # Otherwise, we write out an impossible value (-1.0) to let the caller 
            # know that this particular distance could not be computed
            sys.stdout.write('-1.0\n')


if __name__ == "__main__":
    main()
