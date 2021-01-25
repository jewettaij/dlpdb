#!/usr/bin/env python

"""
This program reads in x,y,z coordinates from a file
and computes the "helixAngleOmega"
(defined in "helix_Omega_angle_derivation.png" in the "doc" subdirectory)
This has been used to infer the periodicity of alpha-helices within proteins.
"""


import sys
from math import sqrt, cos, sin, tan, acos, asin, atan, pi
try:
    from .helixAngleOmega import CalcOmegaFromThetaPhi, CalcOmega
except (ImportError, SystemError, ValueError):
    # not installed as a package
    from helixAngleOmega import CalcOmegaFromThetaPhi, CalcOmega





def main():
    if (len(sys.argv) > 3):
        sys.stderr.write('Error (coords2angles_Omega): number of arguments should not exceed 2.\n'\
                             '    (The two arguments correspond to the number of lines of\n'\
                             '     text to omit from the beginning and end of the file, respectively.)\n'\
                             '     If one argument is passed, then both are assumed to be the same.\n'\
                             '     If no argument is passed, then by default, no data is ignored.\nExiting...\n\n')
        sys.exit(-1)

    elif (len(sys.argv) == 3):
        truncate_a = int(sys.argv[1])
        truncate_b = int(sys.argv[2])
    elif (len(sys.argv) == 2):
        truncate_a = int(sys.argv[1])
        truncate_b = truncate_a
    else:
        truncate_a = 0
        truncate_b = 0


    r_i = []
    # Read the file
    for line in sys.stdin:
        line = line.strip()
        # Each line should contain a list of 3 numbers separated by whitespace.
        # If so, store the 3 numbers in a list variable (named xyz), and append
        # it to the list of coordinates.
        # However some lines might also be blank, in which case we append the
        # empty list [] to the list of coordinates.
        if line == '':
            xyz = []
        else:
            # Each line should contain a list of 3 numbers separated by whitespace.
            xyz = list(map(float, line.split()))
        r_i.append(xyz)


    # Truncate the data we don't want.
    # (Why?  The residues at the beginning and ending of helices 
    #  are less trustworthy then the residues in the middle.)
    r_i = r_i[truncate_a:len(r_i)-truncate_b]

    N = len(r_i)
    for i in range(0,N-4+1):

        if (len(r_i[i])==3)and(len(r_i[i+1])==3)and(len(r_i[i+2])==3)and(len(r_i[i+3])==3):

            Omega = CalcOmega(r_i[i], r_i[i+1], r_i[i+2], r_i[i+3])

            sys.stdout.write(str(Omega*180.0/pi))

        else:
            # Otherwise, we write out an impossible value (-720) to let the caller 
            # know that this particular angle could not be computed
            sys.stdout.write('-720')

        if i+1 < N-4+1:
            sys.stdout.write(' ')

    sys.stdout.write('\n')


if __name__ == "__main__":
    main()
