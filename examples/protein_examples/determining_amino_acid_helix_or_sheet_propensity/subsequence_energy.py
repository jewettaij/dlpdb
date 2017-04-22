#!/usr/bin/env python

"""
Typical usage:

    subsequence_energy.py N [candidate.txt] < all_sequences.txt

 
In this example, "sequences.txt" contains a large list of sequences
from which the probability of each subsequence appearing is determined.
The natural logarithm of the probability of each pattern
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

    if len(sys.argv) == 1:
        sys.stderr.write('Error: expected at least 1 argument.\n'\
                         '       SEE CODE FOR USAGE EXPLANATION.\n')

    # This next variable is the size of the subsequences whose
    # frequency we will be counting.
    subseq_size = int(sys.argv[1])

    subseq_count= {} #keep track of the number of times 
                   #each subsequence appears in the entire file.

    #lines = sys.stdin.readlines()
    #for line in lines:

    for line in sys.stdin:
        line = line.strip()

        # Count the number of residues of each type
        for i in range(0,1+len(line)-subseq_size):
            s = line[i:i+subseq_size]
            if (s in subseq_count):
                subseq_count[s] += 1
            else:
                subseq_count[s] = 1


    # What is the probability of finding each type of amino acid (P[c]) ?
    subseq_count_tot = 0
    for s in subseq_count:
        subseq_count_tot += subseq_count[s]

    P = {}
    deltaP = {}
    for s in subseq_count:
        P[s] = float(subseq_count[s])/subseq_count_tot
        # Suppose x = subseq_count[c] / N
        # is the probability of encountering a certain letter (c)
        # in a huge list of N letters (in this example subseq_count_tot = N).
        # It's easy to show that the uncertainty in the estimate for 
        # the probability is sqrt(x*(1-x)/N) for large N.
        deltaP[s] = math.sqrt((P[s]*(1.0-P[s])) / subseq_count_tot)





    # This information is somewhat meaningless because
    # some amino acids simply occur more often than others 
    # even in the absence of evolutionary pressure.
    #
    # I want to divide this probability (P[s]) by the probability that 
    # this subsequence (s) would occure due to a random mutation in the DNA.

    # Pevo[c] is the probability that an amino acid would be type c
    # due to random mutations (in the absence of evolutionary pressure).
    # I assume (perhaps incorrectly), that this is mostly determined by
    # the number of codons that code for that amino acid.
    #
    # And PevoS[s] is the probability that a subsequence, s, 
    # appears randomly    __
    # I assume PevoS[s] = ||  Pevo[s_i]  
    #                      i
    # where s_i is the ith amino acid in the sequence s

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

    # Now calculate PevoS[s] from Pevo[],
    # for all the sequences s in subseq_count[]
    PevoS = {}
    for s in subseq_count:
        all_standard_amino_acids = True
        Pevo_of_s = 1.0
        for c in s:
            if c in Pevo:
                Pevo_of_s *= Pevo[c]
            else:
                all_standard_amino_acids = False
                Pevo_of_s = 0.0
        if all_standard_amino_acids:
            PevoS[s] = Pevo_of_s


    ##############################################################################
    kB = 0.0019872 #Boltzmann's constant for converting T in Kelvin into kcal/mole
    T  = 300.0

    ##sys.stdout.write('\nProbabilities for each subsequence (sum = 1.0):\n')
    ##sys.stdout.write('seq  P(seq) deltaP Pevo(seq) ln(Pevo(seq)/P(seq)) delta(ln(P(seq)))\n')

    ener = {}
    delta_ener = {}
    for s in subseq_count:
        ##sys.stdout.write(''+s+' '+str(P[s])+' '+str(deltaP[s]))
        if s in PevoS:  #<-(if s consists of only the 20 standard amino acids)
            ##sys.stdout.write(' '+str(PevoS[s]))
            if ((P[s] > 0.0) and (PevoS[s] > 0.0)): #<-being needlessly cautions here
                # The next variable can be used to estimate energetic penalties 
                # as a result of substituting an amino acid of type c into 
                # the environment from which you collected your sequence data.
                # Such energetic penalties may explain why some amino acids
                # appear more or less often than you would expect by 
                # considering the number of codons that code for them.
                ener[s] = math.log(PevoS[s] / P[s])
                delta_ener[s] = deltaP[s] / P[s]
                ##sys.stdout.write(' '+str(ener[s])+' '+str(delta_ener[s]))
        ##sys.stdout.write('\n')




    # If an argument to this program was given, then we assume it is
    # the name of a file containing a candidate sequence.
    # Assuming that the energy of every amino acid is additive,
    # we can estimate the energy of the candidate sequence.
    # Ignoring correlations, this is proportional to the the negative logarithm
    # of the probability of a sequence of that length evolving at random.

    if len(sys.argv) > 2:
        seqfilename = sys.argv[2]
        file = open(seqfilename,'r')

        for line in file:
            line = line.strip()

            tot_ener = 0.0
            variance_tot_ener = 0.0
            sys.stderr.write('---------- calculating energy of sequence \''+line+'\' ----------\n')
            line = line.strip()
            unknown_subsequence = False
            for i in range(0,1+len(line)-subseq_size):
                s = line[i:i+subseq_size]
                if s in ener:
                    tot_ener += ener[s]
                    variance_tot_ener += (delta_ener[s]**2)
                    sys.stderr.write('U(\''+s+'\') = '+str(ener[s]*kB*T)+' +/- '+str(delta_ener[s]*kB*T)+' kcal\mole\n')
                else:
                    sys.stderr.write('\nError: A new subsequence was found: \''+s+'\'\n'\
                                     '         (consequently its probability/energy can not be estimated).\n')
                    unknown_subsequence = True

            sys.stderr.write('---------------------\n'\
                             'sum(\''+line+'\') = ')
            if unknown_subsequence:
                sys.stdout.write('1.0e38 +/- 1.0e38 kcal/mole\n')
                sys.stderr.write(' divided by '+str(len(line)-subseq_size)+':\n')
                sys.stderr.write('1.0e38 +/- 1.0e38 kcal/mole\n')
            else:
                sys.stdout.write(str(tot_ener*kB*T)+' +/- '+str(math.sqrt(variance_tot_ener)*kB*T)+' kcal/mole\n')
                sys.stderr.write(' divided by '+str(len(line)-subseq_size)+':\n')
                sys.stderr.write(str(tot_ener*kB*T/(len(line)-subseq_size))+' +/- '+str(math.sqrt(variance_tot_ener)*kB*T/(len(line)-subseq_size))+' kcal/mole\n')




if __name__ == "__main__":
    main()
