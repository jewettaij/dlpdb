#!/usr/bin/env python

"""
Typical usage:

    amino_acid_energy.py [candidate.txt] < all_sequences.txt

In this example, "all_sequences.txt" contains a large list of sequences
from which the probability of each character appearing is determined.
The natural logarithm of the probability of each type of character 
is returned to the user, along with the probability that the corresponding
amino acid would occur at random.  This is written to the standard-error.

An optional argument ("candidate.txt") contains a sequence.
If you provide an optional argument (in this example "candidate.txt"),
this program assumes the file contains a sequence of letters, and
that you wish to calculate the probability of that sequence 
appearing by chance (using the probabilities you calculated above).
This program calculates the natural logarithm of this probability
and writes it to the standard-out.
"""


import sys
import math



def main():

    type_count= {} #keep track of the number of times 
                   #each character in the sequence 
                   #appears in the entire file.

    #lines = sys.stdin.readlines()
    #for line in lines:

    for line in sys.stdin:
        line = line.strip()

        # Count the number of residues of each type
        for c in line:
            if (c in type_count):
                type_count[c] += 1
            else:
                type_count[c] = 1


    # What is the probability of finding each type of amino acid (P[c]) ?
    type_count_tot = 0
    for c in type_count:
        type_count_tot += type_count[c]

    P = {}
    deltaP = {}
    for c in type_count:
        P[c] = float(type_count[c])/type_count_tot
        # Suppose x = type_count[c] / N
        # is the probability of encountering a certain letter (c)
        # in a huge list of N letters (in this example type_count_tot = N).
        # It's easy to show that the uncertainty in the estimate for 
        # the probability is sqrt(x*(1-x)/N) for large N.
        deltaP[c] = math.sqrt((P[c]*(1.0-P[c])) / type_count_tot)





    # This information is somewhat meaningless because
    # some amino acids simply occur more often than others 
    # even in the absence of evolutionary pressure.
    #
    # I want to divide this probability (P[c]) by the probability that 
    # this amino acid (c) would occure due to a random mutation in the DNA.

    # Pevo[c] is the probability that an amino acid would be type c
    # due to random mutations (in the absence of evolutionary pressure).
    # I assume (perhaps incorrectly), that this is mostly determined by
    # the number of codons that code for that amino acid.

    Pevo = {}
    num_codons = {'F':2,
                  'L':6,
                  'I':3,
                  'M':1,
                  'V':4,
                  'R':6,
                  'K':2,
                  'D':2,
                  'E':2,
                  'N':2,
                  'Q':2,
                  'H':2,
                  'S':6,
                  'P':4,
                  'T':4,
                  'A':4,
                  'Y':2,
                  'W':1,
                  'C':2,
                  'G':4}
    tot_num_codons = 0 # should end up 4^3-3 = 61
    for c in num_codons:
        tot_num_codons += num_codons[c]

    for c in num_codons:
        Pevo[c] = float(num_codons[c]) / float(tot_num_codons)


    ##############################################################################
    kB = 0.0019872 #Boltzmann's constant for converting T in Kelvin into kcal/mole
    T  = 300.0


    sys.stderr.write('\nProbabilities for each type (sum = 1.0):\n')
    sys.stderr.write('AA  P(AA) deltaP Pevo(AA) ln(Pevo(AA)/P(AA)) delta(ln(P(AA)))\n')

    ener = {}
    delta_ener = {}
    for c in type_count:
        sys.stderr.write(''+c+' '+str(P[c])+' '+str(deltaP[c]))
        if c in Pevo:  #<-if c is one of the standard 20 amino acids)
            sys.stderr.write(' '+str(Pevo[c]))
            if ((P[c] > 0.0) and (Pevo[c] > 0.0)): #<-being needlessly cautions here
                # The next variable can be used to estimate energetic penalties 
                # as a result of substituting an amino acid of type c into 
                # the environment from which you collected your sequence data.
                # Such energetic penalties may explain why some amino acids
                # appear more or less often than you would expect by 
                # considering the number of codons that code for them.
                ener[c] = math.log(Pevo[c] / P[c])
                delta_ener[c] = deltaP[c] / P[c]
                sys.stderr.write(' '+str(ener[c]*kB*T)+' '+str(delta_ener[c]*kB*T))
        sys.stderr.write('\n')




    # If an argument to this program was given, then we assume it is
    # the name of a file containing a candidate sequence.
    # Assuming that the energy of every amino acid is additive,
    # we can estimate the energy of the candidate sequence.
    # Ignoring correlations, this is proportional to the the negative logarithm
    # of the probability of a sequence of that length evolving at random.

    if len(sys.argv) > 1:
        seqfilename = sys.argv[1]
        file = open(seqfilename,'r')

        tot_ener = 0.0
        variance_tot_ener = 0.0
        for line in file:
            line = line.strip()
            for c in line:
                if c in ener:
                    tot_ener += ener[c]
                    variance_tot_ener += (delta_ener[c]**2)
                else:
                    sys.stderr.write('\nError: sequence argument contains non-standard amino acid: \''+c+'\'\n')
                    sys.stdout.write('1.0e38\n')
                    exit(-1)

        sys.stdout.write(str(tot_ener*kB*T)+' '+str(math.sqrt(variance_tot_ener)*kB*T)+'\n')




if __name__ == "__main__":
    main()
