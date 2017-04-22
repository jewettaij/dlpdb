#!/usr/bin/env python

import math, cmath, sys

# author: Andrew Jewett
g_program_name = __file__.split('/')[-1]
g_date_str = '2017-4-21'
g_version_str = '0.5.1'



#    MOTIVATION:
#
# Given two lists (denoted x_i, and y_i),
# we want to calculate the Fourier transform of y(x)
# according to the convention:
#
#               1       / infty         -i*k*x    <-- ("i" = sqrt(-1) here)
# ft(k)  =   -------    |         y(x)  e
#           sqrt(2*pi)  /-infty
#
# where:
#            __N_
#            \
#   y(x)  ~= /___  y_i * delta(x - x_i)
#              i
#
# ...which leads to:
#                      __N_
#               1      \            -i*k*x_i
# ft(k)  =   -------   /___  y_i * e 
#           sqrt(2*pi)   i
#
# where:
# x_i and y_i are lists, passed as separate arguments to this function
# (If x_i is omitted, the non-negative integers are substituted.)
#
# 
#     PROBLEMS WITH THIS APPROACH:
#
# 1) Due to the fact that the average signal <y> is not typically zero,
#    there is a huge peak in the frequency distribution
#    near k=0.  This is partially a finite-size effect,
#    since width of this peak seems to decay as the length of each
#    sequence increases, but it is also a consequence of the fact that,
#    there's no reason to expect <y> to be zero.
#    You could ignore this, but the size of the peak is typically huge.
#           To get rid of this, we need to 
#    subtract off the average value of <y> from y_i...
#    ...and <y> should be calculated as the average value of 
#    considering all the sequences.  (<y> should not be calculated
#    separately for each sequence as this would give a peak at k=0 again.)
#
# 2) The magnitude of ft(k) as defined above, 
#    depends on the length of each sequence.
#    Because not all sequences are the same length, 
#    it makes no sense to take the average of them.
#    This causes the uncertainty in the spread of |ft(k)|^2 values
#    to appear artificially high.
#    We need to normalize each the value of ft(k) for each
#    sequence (by dividing by the sequence length, N), 
#    before we add them together to take the average.
#
# 3) It's more convenient to interpret the frequency distribution if
#    the frequencies are the same as the inverse of the period.
#    In other words, a frequency of 1/4 would correspond to a pattern
#    that repeats itself once every 4 letters.
#    In other words, we want to find the frequency 
#    distribution as a function of:
#    f = w / 2pi (which has units of 1/length along the sequence), 
#    not w (which has units of radians/length).
#    
#
#     We use this function instead:
#
#                __N_
#                \            -i*2*pi*f*x_i
#     ft(f)  =   /___  y_i * e 
#                  i


def ftnorm(f, y_i, x_i=[]):
    N = len(y_i)

    # If user neglected to specify the list of x_i values, 
    # then fill it with a list of integers from 0 to N-1
    if (len(x_i) == 0):
        x_i = [i for i in range(0,N)]

    sum = 0.0 + 0.0j
    k = math.pi * 2.0 * f
    for i in range(0,N):
        sum += y_i[i] * cmath.exp(-1j*(k*x_i[i]))

    return sum





def main():

    if (len(sys.argv) != 3):
        sys.stderr.write('Error: expected two arguments: min_seq_length, num_k_values\n'
                         '\n'
                         '  Sytax (typical usage):'
                         '\n'
                         '  '+g_program_name+' 5 100 < sequences.txt > ft.txt\n'
                         '\n'
                         'exiting...\n')
        exit(-1)

    min_seq_length = int(sys.argv[1])
    num_f_values =  int(sys.argv[2])


    # Global variables:

    # The next variable indicates the set of k values at which we will
    # calculate the Fourier Transform.

    f_values = [0.5*(i+1.0)/num_f_values for i in range(0,num_f_values)]

    # The next set of tables determines what we are taking the Fourier transform of.


    #  Amino acid hydrophobicity scale from
    #  the supplementary information for "Recognition of transmembrane 
    #  helices by the endoplasmic reticulum translocon," Hessa et al., 
    #  Nature 433:377 (2005).
    #
    #  More negative means more hydrophobic.  
    #  The met value is also assigned to selenomethionine.
    #  For more information, see:
    #http://www.cgl.ucsf.edu/chimera/docs/UsersGuide/midas/hydrophob.html


    hydrophobicity_hh_lookup = {'D':3.49,  #asp
                                'E':2.68,  #glu
                                'N':2.05,  #asn
                                'Q':2.36,  #gln
                                'K':2.71,  #lys
                                'R':2.58,  #arg
                                'H':2.06,  #his
                                'G':0.74,  #gly
                                'P':2.23,  #pro
                                'S':0.84,  #ser
                                'T':0.52,  #thr
                                'C':-0.13, #cys
                                'M':-0.10, #met
                                #'mse':-0.10,
                                'A':0.11,  #ala
                                'V':-0.31, #val
                                'I':-0.60, #ile
                                'L':-0.55, #leu
                                'F':-0.32, #phe
                                'W':0.30,  #trp
                                'Y':0.68}  #tyr

    #  Amino acid hydrophobicity scale from
    #  "Experimentally determined hydrophobicity scale for proteins at
    #  membrane interfaces," Wimley and White, Nat Struct Biol 3:842 (1996).
    #
    #  More positive means more hydrophobic.
    #  Values for the ionized forms of asp, glu, his are used here, and the
    #  met value is also assigned to selenomethionine.
    #  For more information, see:
    #http://www.cgl.ucsf.edu/chimera/docs/UsersGuide/midas/hydrophob.html

    hydrophobicity_ww_lookup = {'D':-1.23, #asp
                                'E':-2.02, #glu
                                'N':-0.42, #asn
                                'Q':-0.58, #gln
                                'K':-0.99, #lys
                                'R':-0.81, #arg
                                'H':-0.96, #his
                                'G':-0.01, #gly
                                'P':-0.45, #pro
                                'S':-0.13, #ser
                                'T':-0.14, #thr
                                'C':0.24,  #cys
                                'M':0.23,  #met
                                #'mse':0.23,
                                'A':-0.17, #ala
                                'V':-0.07, #val
                                'I':0.31,  #ile
                                'L':0.56,  #leu
                                'F':1.13,  #phe
                                'W':1.85,  #trp
                                'Y':0.94}  #tyr

    #  Amino acid hydrophobicity scale from
    #J. Kyte and R.F. Doolittle, "A Simple Method for Displaying the Hydropathic Character of a Protein" J Mol Biol 157:105 (1982). 
    #  For more information, see:
    #http://www.cgl.ucsf.edu/chimera/docs/UsersGuide/midas/hydrophob.html

    hydrophobicity_lookup = {'D':-3.5, #asp
                             'E':-3.5, #glu
                             'N':-3.5, #asn
                             'Q':-3.5, #gln
                             'K':-3.9, #lys
                             'R':-4.5, #arg
                             'H':-3.2, #his
                             'G':-0.4, #gly
                             'P':-1.6, #pro
                             'S':-0.8, #ser
                             'T':-0.7, #thr
                             'C':2.5,  #cys
                             'M':1.9,  #met
                             #'mse':0.23,
                             'A':1.8, #ala
                             'V':4.2, #val
                             'I':4.5,  #ile
                             'L':3.8,  #leu
                             'F':2.8,  #phe
                             'W':-0.9,  #trp
                             'Y':-1.3}  #tyr


    # The next variable is a copy of one of the lookup tables above.
    # The lookup table will be used in a couple places later on.
    # We want to use the same lookup table each time, so it's safer to choose
    # which one we want at the beginning and refer to it later.
    sequence2signal = hydrophobicity_hh_lookup


    # Main program:
    #
    # Read in all the sequences (lines of text)
    # For each sequence of characters (one sequence per line),
    # 1) convert it to a numerical signal based on the characters on that line
    # 2) Take the Fourier transform of that signal (evaluated at one of the 
    #    frequencies, f, normalized by length)
    #    and compute the power spectrum (the square modulus)
    #    Do this for a variety of differnet spacial frequencies in the list
    #    Call this I_f
    # 4) Keep track of the average, and standard deviation of I_f.
    #    To do this, we need to keep track of the sum of the I_f values and the
    #    sum of the squares of the I_f values.

    # I_f is defined as |y_f|^2, where y_f is the 
    # Fourier transform of y_x (evaluated at one of the frequencies in f_val[]).
    I_f_sum     = [0.0 for f in range(0, num_f_values)]
    # We keep track of the average I_f, as well as its variance:
    I_f_sqr_sum = [0.0 for f in range(0, num_f_values)]

    num_sequences = 0


    #Read the signal in:

    #for line in sys.stdin:  <-- old version (single pass)

    # Wait, just realized I have to subtract off the average signal before I
    # compute the Fourier transforms.  (Otherwise, I get a big
    # spike at k = 0.)  If I subtract off each sequence's signal average
    # before I take the Fourier transform, the graph of the Fourier transform
    # will be pulled down to zero as k approaches 0.
    # I don't want this either.
    # I think best thing to do (which gives you the flattest prettiest graphs)
    # is to compute the average signal over all of the data sets
    # and subtract that from each signal.


    #I have to do this in two passes:
    lines = sys.stdin.readlines()

    # In the first pass, I try to find the average signal
    # (over all the lines in the file)
    y_sum = 0.0
    num_data = 0
    sys.stderr.write('Pass 1: determining the average signal = ')
    for line in lines:
        line = line.strip()
        line_has_strange_data = False
        for c in line:
            if not (c in sequence2signal):
                line_has_strange_data = True;
                sys.stderr.write('WARNING: Unrecognized character in sequence: \''+
                                 c+'\'\n'+
                                 '  sequence = \"'+line+'\"\n'+
                                 '  This sequence will be discarded.\n')

        if ((len(line) >= min_seq_length) and
            (not line_has_strange_data)):

            for c in line:
                y_sum += sequence2signal[c]

            num_data += len(line)

    y_ave = y_sum / num_data
    sys.stderr.write(str(y_ave)+'\n\n')
    sys.stderr.write('Pass 2: calculating the Fourier transform for each sequence.\n')


    # In the second pass, I calculate the Fourier transform of each sequence's 
    # signal, after subtracting off the average signal.  (This was averaged over 
    # all sequences.  This is not the average for each sequence.)
    for line in lines:
        line = line.strip()
        line_has_strange_data = False
        for c in line:
            if not (c in sequence2signal):
                line_has_strange_data = True;

        if ((len(line) >= min_seq_length) and
            (not line_has_strange_data)):

            #sys.stderr.write('sequence = \"'+line+'\"\n')

            y_x = [0.0 for x in range(0,len(line))]
            for x in range(0,len(line)):
                y_x[x] = sequence2signal[line[x]]

            # subtract off the average value from y_x
            for x in range(0,len(line)):
                y_x[x] -= y_ave

            # this list stores the fourier transform of y_of_x

            for f in range(0,num_f_values):
                f_val = f_values[f]
                y_f = ftnorm(f_val, y_x)

                # I_f = |y_f|^2
                I_f = abs(y_f) 
                I_f *= I_f  #this way, I_f will be internally represented as a float
                I_f_sqr         = I_f * I_f
                I_f_sum[f]     += I_f
                I_f_sqr_sum[f] += I_f_sqr

            num_sequences += 1


    # Now print out the results to the standard output.
    # This is a three column file:
    #   column 1: f value (lies in the range (0, 0.5]
    #   column 2: <I_f>   (the average value of the power spectrum at f)
    #   column 3: the fluctuations in I_f (computed using sqrt(<I_f^2> - <I_f>^2))
    for f in range(0,num_f_values):
        I_f_ave     =     I_f_sum[f] / num_sequences
        I_f_sqr_ave = I_f_sqr_sum[f] / num_sequences

        I_f_var = I_f_sqr_ave - I_f_ave*I_f_ave
        # paranoid check:
        if (I_f_var < 0.0):#<--In principle, this could happen due to numerical 
            I_f_var = 0.0  #   roundoff error but it is extremely unlikely

        if (num_sequences > 1):
            I_f_var_unbiased = I_f_var * (num_sequences / (num_sequences - 1.0))
        else:
            I_f_var_unbiased = 1.0  #(any really huge value will do)
        I_f_unc = math.sqrt(I_f_var_unbiased)

        sys.stdout.write(str(f_values[f]) + ' ' + 
                         str(I_f_ave) + ' ' + 
                         str(I_f_unc) + '\n')



if __name__ == "__main__":
    main()
