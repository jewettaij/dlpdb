#!/usr/bin/env python

"""
 Typical usage:
count_NPO_patterns.py N < sequences.txt
"""


import sys


# global variable:

D = 3 #reduced alphabet size (3 types of residues: N, P, O)



# Convert a number from 0,1,2 into a 1-character string: 'N','P', or 'O'
def Digit2Char(d):
    if (d == 0):
        return 'N'
    elif (d == 1):
        return 'P'
    else:
        return 'O'


# Each amino acid is categorized as 
#        'N' non-polar, 
#        'P' polar, or 
#        'O' "other"
# according to the criteria explained in:
#   West, MW; Hecht, MH. 
#   "Binary patterning of polar and nonpolar amino acids in the sequences and structures of native proteins"
#   Protein Sci. 1995, Oct, Vol 4, No 10, pp 2032-2039

translate_AAtoReduced = {'F':'N', # ^
                         'L':'N', # |
                         'I':'N', # non-polar amino acids
                         'M':'N', # |
                         'V':'N', # v

                         'R':'P', # ^
                         'K':'P', # |
                         'D':'P', # |
                         'E':'P', # polar amino acids
                         'N':'P', # |
                         'Q':'P', # |
                         'H':'P', # v

                         'S':'O', # ^
                         'P':'O', # |
                         'T':'O', # |
                         'A':'O', # other
                         'Y':'O', # |
                         'W':'O', # |
                         'C':'O', # |
                         'G':'O'} # V




#def TranslateAAtoReduced(c):
#    if (c in translate_AAtoReduced):
#        return translate_AAtoReduced[c]
#    else:
#        return 'O'





def SequenceUsesNormalAA(s):
    for c in s:
        if (not (c in translate_AAtoReduced)):
            sys.stderr.write('  sequence \"'+s+'\"\n')
            sys.stderr.write('  contains nonstandard char:\"'+c+'\"\n'
                             '---------\n')
            return False
    return True










def main():

    if (len(sys.argv) != 2):
        sys.stderr.write('Error: expected one argument: N\n'
                         '       N = the size of the patterns you want to search for.\n'
                         '       For example, N=5: all penta-peptide patterns are counted\n'
                         '       and ranked in popularity.\n')
        exit(-1)

    N = int(sys.argv[1])



    pattern_num  =[0 for n in range(0,N)]         #stores a sequence as a list of
                                                  #numbers in the range 0 to D-1.
                                                  #Example: [0,1,1,0,0,2]
    pattern_chars=[Digit2Char(0) for n in range(0,N)]#this list is a short sequence
                                                  #of characters, represented using
                                                  #the reduced (D-char) alphabet:
                                                  #eg: ['N','P','P','N','N','O']
                                                  #(it's not a string like 'NPPNNO')
    #print('pattern_chars=',pattern_chars)





    pattern_count = {}





    ## ---------------------------
    ## Optional: Fill pattern_count{} with every possible pattern.
    ## Here, we initialize every possible pattern to have an initial count of 0
    ##           (Probably not a good idea. I'll comment this out later)
    ##
    #
    ## The following function increments an N-digit number in base D.
    ## (For example, [1,0,2,2,2,2] -> [1,1,0,0,0,0] (when D=3)
    ## It returns whether or not overflow occurs.
    #def IncrPattern(pattern_num):
    #    N=len(pattern_num)
    #    n = 0
    #    while (n < N):
    #        if n == N:
    #            break
    #        pattern_num[n] += 1
    #        if (pattern_num[n] < D):
    #            break
    #        else:
    #            pattern_num[n] = 0
    #        n += 1
    #    return (n == N)
    #
    #while (True):
    #    for n in range(0,N):
    #        #sys.stderr.write('  num['+str(n)+']='+str(pattern_num[n])+'\n')
    #        #sys.stderr.write('chars['+str(n)+']='+Digit2Char(pattern_num[n])+'\n')
    #        pattern_chars[n] = Digit2Char(pattern_num[n])
    #
    #    pattern_str = ''.join(pattern_chars) #the dictionary key is a python string
    #    pattern_count[pattern_str] = 0
    #    #sys.stderr.write('possible pattern #'+str(len(pattern_count))+
    #    #                 ': \"'+pattern_str+'\"\n')
    #    if (IncrPattern(pattern_num)):
    #        break;
    ## ---------------------------



    num_patterns = 0


    type_count= {} #keep track of the number of times each 
                   #character in the reduced alphabet appears
                   #in the entire file.

    #lines = sys.stdin.readlines()
    #for line in lines:

    for line in sys.stdin:
        line = line.strip()

        # Count the number of residues of each type
        for c in line:
            if c in translate_AAtoReduced:
                c = translate_AAtoReduced[c]
                if (c in type_count):
                    type_count[c] += 1
                else:
                    type_count[c] = 1



        # Now, consider each contiguous susbtring of length N in the line, convert
        # each character in the string to a character in the reduced alphabet,
        # (which has fewer types of residues)
        # and count the number of times that (reduced) pattern appears.
        for offset in range(0,len(line)+1-N):
            pattern_orig = line[offset:offset+N]

            # If the pattern_orig contains any non-standard amino acids, ignore it
            if SequenceUsesNormalAA(pattern_orig):
                # Otherwise, translate the pattern_orig into the reduced alphabet
                # (I have to store this as a list because strings are immutable.)
                pattern_translated = []
                for c in pattern_orig:
                    pattern_translated.append(translate_AAtoReduced[c])

                #The "pattern_count" dictionary uses strings for keys, not lists,
                #so I must convert this list of characters back into a string:
                pattern = ''.join(pattern_translated)
                #sys.stderr.write('encountered pattern=\"'+pattern+'\"\n')

                # Have we encountered this pattern before?
                if (pattern in pattern_count):
                    pattern_count[pattern] += 1
                    num_patterns += 1
                else:
                    pattern_count[pattern] = 1


    # What is the probability of finding each type of amino acid?
    type_count_total = 0
    for c in type_count:
        type_count_total += type_count[c]

    type_prob = {}
    sys.stderr.write('\nProbabilities for each type (sum = 1.0):\n')
    for c in type_count:
        type_prob[c] = float(type_count[c])/type_count_total
        sys.stderr.write('\''+c+'\' '+str(type_prob[c])+'\n')


    # For reference:
    # ---
    #   compare these to the probabilities reported in
    #      West+Hecht, Protein Science, 1995:
    # type_prob = {'N':0.2687,
    #              'P':0.3295,
    #              'O':0.4018}
    # Note: The probabilities computed by this program may be more accurate
    #       than the ones reported in the paper since they were using a very small
    #       early snapshot of the PDB repository containing only ~150 proteins.
    # ---
    # More recently (2/06/2009), using the culledpdb list from the PISCES server
    # (sequence homology<=40%, resolution<=2.5A, R<=1.0),  I got:
    #
    #  for helices:
    # type_prob = {'P':0.370930557862,
    #              'O':0.344594099281,
    #              'N':0.284475342857}
    #
    #  for sheets:
    # type_prob = {'P':0.246765741807,
    #              'O':0.333002936995,
    #              'N':0.420231321198}



    sys.stderr.write('\nReweighting the probabilities for each pattern according to composition.\n')

    # Now compute the relative probability ("rp") that such a pattern would
    # be encountered.  The relative probability is defined as:
    #
    #                probability pattern encountered
    #   rp = ------------------------------------------------
    #         probability pattern encountered drawing randomly
    #          from a pool of letters with wieghted probability
    #
    # The numerator is just pattern_count[pattern] / num_patterns
    # The denominator is the product of 

    pattern_rp = {}
    for pattern in pattern_count:
        prob = float(pattern_count[pattern]) / num_patterns
        rand_prob = 1.0
        for c in pattern:
            rand_prob *= type_prob[c]
        pattern_rp[pattern] = prob / rand_prob


    # Now sort the list of patterns by relative probability (rp) and
    # write out the list of patterns (and probabilities) standard out:

    sys.stderr.write('\nSorting patterns by the weighted relative probability.\n')

    import operator
    for pc in sorted(pattern_rp.iteritems(),
                     key=operator.itemgetter(1),
                     reverse=True):

        pattern = pc[0]
        count = pattern_count[pattern]
        prob = float(count)/num_patterns

        sys.stdout.write(pattern
                         + ' ' + str(prob)
                         + ' ' + str(pattern_rp[pattern])
                         + '\n')



    # Old code: (sorts by probability, not by relative probability (rp))
    #import operator
    #for pc in sorted(pattern_count.iteritems(),
    #                 key=operator.itemgetter(1),
    #                 reverse=True):
    #    pattern = pc[0]
    #    rand_prob = 1.0
    #    for c in pattern:
    #        rand_prob *= type_prob[c]
    #
    #    count = pc[1]
    #    prob = float(pc[1])/num_patterns
    #    sys.stdout.write(pattern
    #                     + ' ' + str(prob)
    #                     + ' ' + str(prob/rand_prob)
    #                     + '\n')




if __name__ == "__main__":
    main()
