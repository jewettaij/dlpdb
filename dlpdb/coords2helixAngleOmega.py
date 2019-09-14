#!/usr/bin/env python

"""
This program reads in x,y,z coordinates from a file
and computes the "helixAngleOmega"
(defined in "helix_Omega_angle_derivation.png" in the "doc" subdirectory)
This has been used to infer the periodicity of alpha-helices within proteins.
"""


from math import sqrt, cos, sin, tan, acos, asin, atan, pi
import sys
from helixAngleOmega import CalcOmega



def length_v(r):
    lsqd = 0.0
    for d in range(0,len(r)):
        lsqd += r[d]*r[d]
    return sqrt(lsqd)


def inner_prod_v(r1,r2):
    result = 0.0
    for d in range(0,len(r1)):
        result += r1[d]*r2[d]
    return result


def cross_prod_v3(a,b):
    c = [0.0,0.0,0.0]
    c[0] = a[1]*b[2] - a[2]*b[1]
    c[1] = a[2]*b[0] - a[0]*b[2]
    c[2] = a[0]*b[1] - a[1]*b[0]
    return c


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

            r10 = [0.0, 0.0, 0.0]
            r21 = [0.0, 0.0, 0.0]
            r32 = [0.0, 0.0, 0.0]
            for d in range(0,3):
                r10[d] = r_i[i+1][d] - r_i[i][d]
                r21[d] = r_i[i+2][d] - r_i[i+1][d]
                r32[d] = r_i[i+3][d] - r_i[i+2][d]
            l10 = length_v(r10)
            l21 = length_v(r21)
            l32 = length_v(r32)

            n012 = cross_prod_v3(r10, r21)
            n123 = cross_prod_v3(r21, r32)


            # The torsion-angle or 4-body angle is named "angle0124"
            cos_phi = inner_prod_v(n012, n123) /(length_v(n012)*length_v(n123))

            # There is a problem whenever 4 consecutive atoms are coplanar:
            #
            #            *---*
            #                |      (all 4 atoms are coplanar, and phi = 0)
            #            *---*
            #
            # In this example, the torsion angle phi is well defined and =0.
            # The problem is that, due to floating point roundoff
            # "cos_phi" can sometimes slightly exceed 1.
            # This causes a NAN when you calculate acos(cos_phi).

            if (cos_phi > 1.0):
                cos_phi = 1.0
            elif (cos_phi < -1.0):
                cos_phi = -1.0

            phi = acos(cos_phi)

            # This formula does not distinguish positive and negative phi.
            #
            # Negative torsion angles:
            #
            # Check if  the position of atom i+3 is above the phi=0 plane
            # (in the region of positive phi), or below the phi=0 plane.
            # It is above the phi=0 plane if the bond from atom i+2 to i+3
            # points in the same direction as the negative-phi-tangent-vector 
            # for atom i (not i+3)  (...which points in the n012 direction)
            if inner_prod_v(n012, r32) < 0.0:
                phi = -phi

            # The two bond-angles or 3-body angles are named "angle012" and "angle123"
            angle012 = acos( -inner_prod_v(r10, r21) / (l10 * l21) )
            angle123 = acos( -inner_prod_v(r21, r32) / (l21 * l32) )
            # (The negative sign above comes from the fact that we are really 
            #  interested in the angle between r21 and r01 (which is -r10).)

            # Convert these angles to degrees, and print it out
            #sys.stdout.write(str(angle012*180.0/pi)+' ')
            #sys.stdout.write(str(angle123*180.0/pi)+' ')

            # Let theta = the average of the two 3-body angles
            theta = 0.5 * (angle012 + angle123) 

            # Omega (usually 360/3.6 ~= 100 degrees) is the helix rotation angle.
            Omega = CalcOmega(theta, phi)

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
