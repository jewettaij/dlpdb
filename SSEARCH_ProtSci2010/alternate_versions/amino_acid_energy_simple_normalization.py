#!/usr/bin/env python

"""
Typical usage:

 amino_acid_energy.py sub_sequences.txt all_sequences.txt [candidate.txt]
 
In this example, "sub_sequences.txt" 
contains a large list of sequences from which the probability 
of each character appearing is determined. 
The natural logarithm of the probability of each type of amino acid
in "sub_sequences.txt" is returned to the user.
The negative-logarithm of this probability can be interpreted as the
chemical potential energy of inserting this type of amino acid into 
the context that these sub-sequences have in common (in units of kT).  
(For example, if the sub-sequences are all beta sheets, then this 
 logarithm may crudely reflect energy if inserting this kind of amino 
 acid into a beta-sheet.)

Problem:
   Even in the absence of evolutionary selection, without worrying
about how each amino acid effects stability, 
some amino acids are more common than others.
For one thing, some amino acids correspond to more DNA codons than others.
(Actually, there is another version of this program which considers
 the number of condons for each amino acid.)
However in this program, we make fewer assumptions.
In this program, the natural (unbiased) frequency of each amino acid 
in proteins is not estimated by the number of codons.
Instead the probability that each type of amino acid would be selected
at random in the absence of evolutionary considerations, is determined 
from their frequency in some other set.
This is the SECOND argument for this program.

Typically, the first argument (sub_sequences.txt in this example)
contains sequences which are subsequences of the second argument.
The second argument (all_sequences.txt) is typically list of sequences
corresponding to a representative sample of some large set of proteins
(for example the set of all known cytosolic proteins).
I myself don't know if this is better than looking at the number of codons
for each amino acid.

-->  In any case, the probability of belonging to the first set
--> is divided by the probability of belonging to the second set,
--> and then the (negative) logarithm is taken and returned to the user.
--> This is done separately for each amino acid, 
--> and printed to the standard error (along with uncertainties).
 
An optional third argument ("candidate.txt") contains a sequence.
If you provide an optional argument (in this example "candidate.txt"),
this program assumes the file contains a sequence of letters, and
that you wish to calculate the probability of that sequence 
appearing by chance (using the probabilities you calculated above).
This program calculates the natural logarithm of this probability
and writes it to the standard-out.
"""

import sys
import math



def CalcSequenceEnergy(s,energy):
    tot_ener = 0.0
    for c in s:
        if c in energy:
            tot_ener += energy[c]
        else:
            sys.stderr.write('\nError: sequence argument contains non-standard amino acid: \''+c+'\'\n')
            exit(-1)
    return tot_ener



def CalcSequenceEnergyErr(s,energy,delta_energy):
    tot_ener = 0.0
    var_tot_ener = 0.0
    for c in s:
        if c in energy:
            tot_ener += energy[c]
            var_tot_ener += (delta_energy[c]**2)
        else:
            sys.stderr.write('\nError: sequence argument contains non-standard amino acid: \''+c+'\'\n')
            exit(-1)
    return tot_ener,math.sqrt(var_tot_ener)



def main():
    if len(sys.argv) <= 2 :
        sys.stderr.write('Error: expected at least 2 arguments.\n'\
                         '       SEE CODE FOR USAGE EXPLANATION.\n')
        sys.exit(-1)

    filename = sys.argv[1]
    evofilename = sys.argv[2]


    type_count= {} #the number of times each character appears in the file
    type_count_tot = 0 #the number of characters in the file

    file = open(filename,'r')
    for line in file:
        line = line.strip()

        # Count the number of residues of each type
        for c in line:
            if (c in type_count):
                type_count[c] += 1
            else:
                type_count[c] = 1
            type_count_tot += 1
    file.close()


    ########## OLD SIMPLE VERSION ############
    #for c in type_count:
    #    P[c] = float(type_count[c])/type_count_tot
    #    # Suppose x = type_count[c] / N
    #    # is the probability of encountering a certain letter (c)
    #    # in a huge list of N letters (in this example type_count_tot = N).
    #    # It's easy to show that the uncertainty in the estimate for 
    #    # the probability is sqrt(x*(1-x)/N) for large N.
    #    deltaP[c] = math.sqrt((P[c]*(1.0-P[c])) / type_count_tot)


    ########## NEW VERSION ##############
    # The probabilities should sum up to the probability that
    # a randomly chosen amino acid happens to be in a helix
    #NAA_tot=2123360 # number of amino acids in the database
    #NAA_hel=942456  # number of amino acids in the database in helices
    #NAA_she=509097  # number of amino acids in the database in sheets
    # After throwing away the PDB files containing "MEMBRANE" OR " NMR", we get
    NAA_tot=1986157  # number of amino acids in the database
    NAA_hel=886428   # number of amino acids in the database in helices
    NAA_she=471407   # number of amino acids in the database in sheets
    prior__helix = (float(NAA_hel)/float(NAA_tot))
    prior__sheet = (float(NAA_she)/float(NAA_tot))
    P = {}
    deltaP = {}
    for c in type_count:
        conditional_probability__c = float(type_count[c])/type_count_tot

        #######################################
        # Giovanni's method:
        P[c] = conditional_probability__c
        deltaP[c] = math.sqrt((P[c]*(1.0-P[c])) / type_count_tot)

        #######################################
        # Andrew's method:
        #
        #joint_probability__c_helix = conditional_probability__c * prior__helix
        #joint_probability__c_sheet = conditional_probability__c * prior__sheet
        #P[c] = joint_probability__c_helix
        #P[c] = joint_probability__c_sheet

        # Uncertainty:
        # Suppose x = 
        # is the probability of encountering a certain letter (c)
        # in a huge list of N letters (in this example N=NAA_tot).
        # It's easy to show that the uncertainty in the estimate for 
        # the probability is sqrt(x*(1-x)/N) for large N.
        #deltaP[c] = math.sqrt((P[c]*(1.0-P[c])) / NAA_tot)
        #deltaP[c] = math.sqrt((P[c]*(1.0-P[c])) / type_count_tot)

        #######################################


    # This information is somewhat meaningless because
    # some amino acids simply occur more often than others 
    # even in the absence of evolutionary pressure.
    #
    # I want to divide this probability (P[c]) by the probability that 
    # this amino acid (c) would occure due to a random mutation in the DNA.
    # 
    # To estimate that, (in this version of the software, at least),
    # I consider the probability that this amino acid appears
    # in a larger set of amino acids which is representative of
    # proteins as a whole.
    type_count_evo = {}
    type_count_evo_tot = 0
    evofile = open(evofilename,'r')
    for line in evofile:
        line = line.strip()
        # Count the number of residues of each type
        for c in line:
            if (c in type_count_evo):
                type_count_evo[c] += 1
            else:
                type_count_evo[c] = 1
            type_count_evo_tot += 1
    evofile.close()


    Pevo = {}
    deltaPevo = {}
    for c in type_count:
        # Pevo[c] is the probability that an amino acid would be type c
        # due to random mutations (in the absence of evolutionary pressure).
        Pevo[c] = float(type_count_evo[c])/type_count_evo_tot
        # Suppose x = type_count[c] / N
        # is the probability of encountering a certain letter (c)
        # in a huge list of N letters (in this example type_count_evo_tot = N).
        # It's easy to show that the uncertainty in the estimate for 
        # the probability is sqrt(x*(1-x)/N) for large N.
        deltaPevo[c] = math.sqrt((Pevo[c]*(1.0-Pevo[c])) / type_count_evo_tot)



    ener = {}
    delta_ener = {}
    for c in P:
        assert(P[c] > 0.0)
        if c in Pevo:
            assert(Pevo[c] > 0.0)
            # The next variable can be used to estimate energetic penalties 
            # as a result of substituting an amino acid of type c into 
            # the environment from which you collected your sequence data.
            # Such energetic penalties may explain why some amino acids
            # appear more or less often than you would expect by 
            # considering the number of codons that code for them.
            ener[c] = math.log(Pevo[c] / P[c])
            delta_lnP = deltaP[c] / P[c]
            delta_lnPevo = deltaPevo[c] / Pevo[c]
            delta_ener[c] = math.sqrt(delta_lnP**2 + delta_lnPevo**2)



    # Amino acid 1-letter and 3-letter codes:
    OneToThree = {'E':'Glu',
                  'K':'Lys',
                  'D':'Asp',
                  'A':'Ala',
                  'Q':'Gln',
                  'F':'Phe',
                  'N':'Asn',
                  'M':'Met',
                  'L':'Leu',
                  'I':'Ile',
                  'Y':'Tyr',
                  'V':'Val',
                  'W':'Trp',
                  'G':'Gly',
                  'T':'Thr',
                  'H':'His',
                  'R':'Arg',
                  'S':'Ser',
                  'P':'Pro',
                  'C':'Cys'}

    ThreeToOne = dict([[v,k] for k,v in OneToThree.items()])

    # Experimental measures of helix propensity:

    # This one from (Horovitz,Matthews,&Fersht JMB 1992 227:560-568) is in kCal/mole
    helixPropensityHorovitz_3ltr = {'Ala':0.00,
                                   'Arg':0.14,
                                   'Lys':0.19,
                                   'Met':0.31,
                                   'Leu':0.35,
                                   'Ser':0.41,
                                   'Gln':0.48,
                                   'Glu':0.55,
                                   'Asn':0.66,
                                   'Phe':0.69,
                                   'Asp':0.71,
                                   'His':0.78,
                                   'Thr':0.79,
                                   'Ile':0.81,
                                   'Tyr':0.82,
                                   'Val':0.88,
                                   'Gly':0.91,
                                   'Trp':0.98,
                                   'Cys':1.00,
                                   'Pro':4.08}

    helixPropensityHorovitz = {}
    for AA in helixPropensityHorovitz_3ltr:
        c = ThreeToOne[AA]
        helixPropensityHorovitz[c] = helixPropensityHorovitz_3ltr[AA]


    # Alternate helix propensity metric:
    # (from O'Neil and Degrado, Science 1990, 250:646-651)
    # Again, units are in kCal/mole.
    #  The experiment substituted residues into the following sequence:
    #  Ac-EWEALEKKLAALE-X-KLQALEKKLEALEHG at position "X"
    # This is a designed protein, not a natural protein

    helixPropensityONeil_3ltr = {'Ala':-0.77,
                                'Arg':-0.68,
                                'Lys':-0.65,
                                'Met':-0.50,
                                'Leu':-0.62,
                                'Ser':-0.35,
                                'Gln':-0.33,
                                'Glu':-0.27,
                                'Asn':-0.07,
                                'Phe':-0.41,
                                'Asp':-0.15,
                                'His':-0.06,
                                'Thr':-0.11,
                                'Ile':-0.23,
                                'Tyr':-0.17,
                                'Val':-0.14,
                                'Gly':0.00,
                                'Trp':-0.45,
                                'Cys':-0.23,
                                'Pro': 3.00}

    helixPropensityONeil = {}
    for AA in helixPropensityONeil_3ltr:
        c = ThreeToOne[AA]
        helixPropensityONeil[c] = helixPropensityONeil_3ltr[AA]


    helixPropensityBlaber44_3ltr = {'Ala':-0.96,
                                    'Leu':-0.92,
                                    'Met':-0.86,
                                    'Ile':-0.84,
                                    'Gln':-0.80,
                                    'Arg':-0.77,
                                    'Lys':-0.73,
                                    'Tyr':-0.72,
                                    'Val':-0.63,
                                    'Phe':-0.59,
                                    'Trp':-0.58,
                                    'His':-0.57,
                                    'Thr':-0.54,
                                    'Glu':-0.53,
                                    'Ser':-0.53,
                                    'Asp':-0.42,
                                    'Cys':-0.42,
                                    'Asn':-0.39,
                                    'Gly':0.00,
                                    'Pro':2.50}


    helixPropensityBlaber44 = {}
    for AA in helixPropensityBlaber44_3ltr:
        c = ThreeToOne[AA]
        helixPropensityBlaber44[c] = helixPropensityBlaber44_3ltr[AA]


    # The next couple tables were taken from table 1 of 
    # Myers+Pace+Scholtz, Biochemistry 1997,36:10923-10929
    #
    # (Incidentally, what they called the "wild type" sequence
    #  based on residues 13-29 of RNase T1:
    #  SSDVSTAQAAGYKLHED
    #  However they muted residue 23 to A 
    #  to stabilize the helix.  This resulted in:
    #  SSDVSTAQAAAYKLHED
    #  They then replaced the residue at position 21
    #  with 19 different amino acids (all except proline)
    #
    #  SSDVSTAQxAAYKLHED
    #          ^
    #    mutated here

    helixPropensityMyersProteinPh2p5 = {'A':0.00,
                                        'C':0.74,
                                        'D':-0.33,
                                        'E':-0.05,
                                        'F':0.57,
                                        'G':0.90,
                                        'H':0.56,
                                        'I':0.44,
                                        'K':0.51,
                                        'L':0.13,
                                        'M':0.15,
                                        'N':-0.34,
                                        'Q':0.26,
                                        'R':0.41,
                                        'S':0.49,
                                        'T':0.57,
                                        'V':0.66,
                                        'W':0.30,
                                        'Y':0.39}

    # pH=7 data taken from table 1 of Myers+Pace+Scholtz, Biochemistry 1997
    #      (right-most column) when available.  If not available, this data
    #      was taken from the pH=2.5 column.  (Perhaps this is okay because, 
    #      if I understand correctly, only two of the amino acid's protonation 
    #      states, D,E was effected by the change in pH, and the others
    #      did not change that much. However I remain confused by this
    #      I think this is what these authors did when they plotted 
    #      figure 3.  They didn't explain which pH they used, but I assume that
    #      it is reflected in the protonation state of the amino acids in the
    #      Still, some amino acids, like H,Q,S have different propensities 
    #      at pH=2.5 and 7, and I don't know which value they used
    #      but I'm guessing they used the pH=7 value.  This was confusing.)

    helixPropensityMyersProteinPh7 = helixPropensityMyersProteinPh2p5

    helixPropensityMyersProteinPh7_corrections = {'A':0.0,
                                                  'D':0.71,
                                                  'E':0.69,
                                                  'H':0.17,
                                                  'Q':0.40,
                                                  'S':0.40}
    for c in helixPropensityMyersProteinPh7_corrections:
        helixPropensityMyersProteinPh7[c] = helixPropensityMyersProteinPh7_corrections[c]


    # Taken from table 2 of Myers+Pace+Scholtz, Biochemistry 1997,36:10923-10929
    helixPropensityMyersPeptidePh2p5 = {'A':0.0,
                                        'C':0.53,
                                        'D':0.66,
                                        'E':0.17,
                                        'F':0.61,
                                        'G':0.98,
                                        'H':1.2,
                                        'I':0.38,
                                        'K':0.45,
                                        'L':0.25,
                                        'M':0.18,
                                        'N':0.66,
                                        'P':1.1,
                                        'Q':0.31,
                                        'R':0.56,
                                        'S':0.51,
                                        'T':0.71,
                                        'V':0.66,
                                        'W':0.14,
                                        'Y':0.45}



    helixPropensityMyersPeptidePh7 = {'A':0.0,
                                      'C':0.51,
                                      'D':0.68,
                                      'E':0.31,
                                      'F':0.59,
                                      'G':0.95,
                                      'H':0.67,
                                      'I':0.29,
                                      'K':0.30,
                                      'L':0.19,
                                      'M':0.12,
                                      'N':0.58,
                                      'P':1.1,
                                      'Q':0.29,
                                      'R':0.38,
                                      'S':0.42,
                                      'T':0.59,
                                      'V':0.61,
                                      'W':0.02,
                                      'Y':0.31}


    # The help the reader, I also subtract off the average of each energy.
    # (First I must calculate these averages.)

    ave_ener = 0.0
    n_temp = 0
    for c in ener:
        if c != 'P': # According to Horovitz, Proline is an outlier. exclude it
            ave_ener += ener[c]
            n_temp += 1
    if n_temp > 0:
        ave_ener /= n_temp

    avePropHorovitz = 0.0
    for AA in helixPropensityHorovitz_3ltr:
        if AA != 'Pro': # According to Horovitz, Proline is an outlier. exclude it
            avePropHorovitz += helixPropensityHorovitz_3ltr[AA]
    avePropHorovitz /= (len(helixPropensityHorovitz_3ltr) - 1)

    avePropONeil = 0.0
    # The help the reader, subtract off the average. (First calculate the average)
    for AA in helixPropensityONeil_3ltr:
        if AA != 'Pro': # According to Horovitz, Proline is an outlier. exclude it
            avePropONeil += helixPropensityONeil_3ltr[AA]
    avePropONeil /= (len(helixPropensityONeil_3ltr) - 1)

    count = 0
    avePropBlaber44 = 0.0
    # The help the reader, subtract off the average. (First calculate the average)
    for c in helixPropensityBlaber44:
        if c != 'P': # According to Horovitz, Proline is an outlier. exclude it
            avePropBlaber44 += helixPropensityBlaber44[c]
            count += 1
    avePropBlaber44 /= count

    count = 0
    avePropMyersProtein= 0.0
    # The help the reader, subtract off the average. (First calculate the average)
    for c in helixPropensityMyersProteinPh7:
        if c != 'P': # According to Horovitz, Proline is an outlier. exclude it
            avePropMyersProtein += helixPropensityMyersProteinPh7[c]
            count += 1
    avePropMyersProtein /= count

    count = 0
    avePropMyersPeptide= 0.0
    # The help the reader, subtract off the average. (First calculate the average)
    for c in helixPropensityMyersPeptidePh7:
        if c != 'P': # According to Horovitz, Proline is an outlier. exclude it
            avePropMyersPeptide += helixPropensityMyersPeptidePh7[c]
            count += 1
    avePropMyersPeptide /= count


    ##############################################################################
    # Now sort the amino acids by energy, and write out the table
    kB = 0.0019872 #Boltzmann's constant for converting T in Kelvin into kcal/mole
    T  = 300.0
    sys.stdout.write('\nProbabilities for each type (sum = '+str(prior__helix)+'):\n')
    #sys.stdout.write('C AA Horovitz1992 ONeil1990 ln(Pevo(AA)/P(AA))*kB*T delta(ln(P(AA)))*kB*T P(AA) deltaP Pevo(AA)\n')
    sys.stdout.write('C AA P(AA) deltaP Pevo(AA) ln(Pevo(AA)/P(AA))*kB*T delta(ln(P(AA)))*kB*T ONeil1990 Horovitz1992 Blaber1993 Myers1997pep Myers1997pro\n')

    #for c in ener:
    #    sys.stdout.write(c)


    #sort by value (energy):
    #ener_sorted = sorted(ener.iteritems(), key=lambda(k,v):v)
    #sort by key (amino acid):
    #ener_sorted = sorted(ener.iteritems(), key=lambda(k,v):k)
    ## alternately, to sort by Horovitsz's energy values:
    ##ener_sorted = sorted(helixPropensityHorovitz.iteritems(), key=lambda(k,v):v)
    #for ener_pair in ener_sorted:
    #    c = ener_pair[0]

    AAsorted = sorted(ThreeToOne.iteritems(), key=lambda(k,v):k)
    for AA_c_pair in AAsorted:
        c = AA_c_pair[1]

        #sys.stdout.write(c)
        #sys.stdout.write(c+' ')
        if c in OneToThree:
            sys.stdout.write(OneToThree[c])
        else:
            sys.stdout.write('???')

        if c in P:
         #sys.stdout.write(' '+str(P[c])+' '+str(deltaP[c]))
            sys.stdout.write((' %5.3f' % (P[c]*100.0))+(' %5.3f'%(deltaP[c]*100.0)))
        else:
            sys.stdout.write((' %5.3f'%0.0)+(' %4.2f'%0.0))

        if c in Pevo:
            #sys.stdout.write(' '+str(Pevo[c]))
            sys.stdout.write((' %5.3f' % (Pevo[c]*100.0))+(' %5.3f' % (deltaPevo[c]*100.0)))

        if c in ener:
            assert(c in P)
            ave_ener=0.0
            #sys.stdout.write((' %5.s3f'%((ener[c]-ave_ener)*kB*T)) + (' %5.3f'%(delta_ener[c]*kB*T)))
            #sys.stdout.write(('& %5.s3f'%((ener[c]-ave_ener)*kB*T)) + (' %5.3f'%(delta_ener[c]*kB*T)))
            sys.stdout.write(' %5.3f'%((ener[c]-ave_ener)*kB*T))

        #if c in helixPropensityONeil:
        #    sys.stdout.write(' %4.2f' % (helixPropensityONeil[c]-avePropONeil))

        #if c in helixPropensityHorovitz:
        #    sys.stdout.write(' %4.2f' % (helixPropensityHorovitz[c]-avePropHorovitz))

        #if c in helixPropensityBlaber44:
        #    sys.stdout.write(' %4.2f' % (helixPropensityBlaber44[c]-avePropBlaber44))

        #if c in helixPropensityMyersPeptidePh7:
        #    sys.stdout.write(' %4.2f' % (helixPropensityMyersPeptidePh7[c]-avePropMyersPeptide))

        #if c in helixPropensityMyersProteinPh7:
        #    sys.stdout.write(' %4.2f' % (helixPropensityMyersProteinPh7[c]-avePropMyersProtein))

        sys.stdout.write('\n')


    # Now, calculate the Pearlson covariance (R) between the different metrics:
    cov11_Horovitz = 0.0
    cov22_Horovitz = 0.0
    cov12_Horovitz = 0.0
    cov11_ONeil = 0.0
    cov22_ONeil = 0.0
    cov12_ONeil = 0.0
    cov11_Blaber44 = 0.0
    cov22_Blaber44 = 0.0
    cov12_Blaber44 = 0.0
    cov11_MyersProtein = 0.0
    cov22_MyersProtein = 0.0
    cov12_MyersProtein = 0.0
    cov11_MyersPeptide = 0.0
    cov22_MyersPeptide = 0.0
    cov12_MyersPeptide = 0.0
    num_common_AA_Horovitz = 0
    num_common_AA_ONeil = 0
    num_common_AA_Blaber44 = 0
    num_common_AA_MyersProtein = 0
    num_common_AA_MyersPeptide = 0
    for c in ener:
        if (c in helixPropensityHorovitz) and (c != 'P'): #exclude proline
            # Note: I left in the kB*T factor even though it cancels out later
            num_common_AA_Horovitz += 1
            cov11_Horovitz += (kB*T*(ener[c]-ave_ener))**2
            cov22_Horovitz += (helixPropensityHorovitz[c]-avePropHorovitz)**2
            cov12_Horovitz += kB*T*(ener[c]-ave_ener)*(helixPropensityHorovitz[c]-avePropHorovitz)
        if (c in helixPropensityONeil) and (c != 'P'): #exclude proline:
            # Note: I left in the kB*T factor even though it cancels out later
            num_common_AA_ONeil += 1
            cov11_ONeil += (kB*T*(ener[c]-ave_ener))**2
            cov22_ONeil += (helixPropensityONeil[c]-avePropONeil)**2
            cov12_ONeil += kB*T*(ener[c]-ave_ener)*(helixPropensityONeil[c]-avePropONeil)
        if (c in helixPropensityBlaber44) and (c != 'P'): #exclude proline:
            # Note: I left in the kB*T factor even though it cancels out later
            num_common_AA_Blaber44 += 1
            cov11_Blaber44 += (kB*T*(ener[c]-ave_ener))**2
            cov22_Blaber44 += (helixPropensityBlaber44[c]-avePropBlaber44)**2
            cov12_Blaber44 += kB*T*(ener[c]-ave_ener)*(helixPropensityBlaber44[c]-avePropBlaber44)
        if (c in helixPropensityMyersPeptidePh7) and (c != 'P'): #exclude proline:
            num_common_AA_MyersPeptide += 1
            cov11_MyersPeptide += (kB*T*(ener[c]-ave_ener))**2
            cov22_MyersPeptide += (helixPropensityMyersPeptidePh7[c]-avePropMyersPeptide)**2
            cov12_MyersPeptide += kB*T*(ener[c]-ave_ener)*(helixPropensityMyersPeptidePh7[c]-avePropMyersPeptide)
        if (c in helixPropensityMyersProteinPh7) and (c != 'P'): #exclude proline:
            num_common_AA_MyersProtein += 1
            cov11_MyersProtein += (kB*T*(ener[c]-ave_ener))**2
            cov22_MyersProtein += (helixPropensityMyersProteinPh7[c]-avePropMyersProtein)**2
            cov12_MyersProtein += kB*T*(ener[c]-ave_ener)*(helixPropensityMyersProteinPh7[c]-avePropMyersProtein)


    sys.stdout.write('\n-------------------------------------\n')
    sys.stdout.write('         sigma = '+str(math.sqrt(cov11_Horovitz/num_common_AA_Horovitz))+'\n')
    sys.stdout.write(' sigmaHorovitz = '+str(math.sqrt(cov22_Horovitz/num_common_AA_Horovitz))+'\n')
    sys.stdout.write('             R = '+str(math.sqrt((cov12_Horovitz*cov12_Horovitz)/(cov11_Horovitz*cov22_Horovitz)))+'\n')
    sys.stdout.write('-------------------------------------\n')
    sys.stdout.write('         sigma = '+str(math.sqrt(cov11_ONeil/num_common_AA_ONeil))+'\n')
    sys.stdout.write('    sigmaONeil = '+str(math.sqrt(cov22_ONeil/num_common_AA_ONeil))+'\n')
    sys.stdout.write('             R = '+str(math.sqrt((cov12_ONeil*cov12_ONeil)/(cov11_ONeil*cov22_ONeil)))+'\n')
    sys.stdout.write('-------------------------------------\n')
    sys.stdout.write('         sigma = '+str(math.sqrt(cov11_Blaber44/num_common_AA_Blaber44))+'\n')
    sys.stdout.write('    sigmaBlaber44 = '+str(math.sqrt(cov22_Blaber44/num_common_AA_Blaber44))+'\n')
    sys.stdout.write('             R = '+str(math.sqrt((cov12_Blaber44*cov12_Blaber44)/(cov11_Blaber44*cov22_Blaber44)))+'\n')
    sys.stdout.write('-------------------------------------\n')
    sys.stdout.write('         sigma = '+str(math.sqrt(cov11_MyersProtein/num_common_AA_MyersProtein))+'\n')
    sys.stdout.write('    sigmaMyersProtein = '+str(math.sqrt(cov22_MyersProtein/num_common_AA_MyersProtein))+'\n')
    sys.stdout.write('             R = '+str(math.sqrt((cov12_MyersProtein*cov12_MyersProtein)/(cov11_MyersProtein*cov22_MyersProtein)))+'\n')
    sys.stdout.write('-------------------------------------\n')
    sys.stdout.write('         sigma = '+str(math.sqrt(cov11_MyersPeptide/num_common_AA_MyersPeptide))+'\n')
    sys.stdout.write('    sigmaMyersPeptide = '+str(math.sqrt(cov22_MyersPeptide/num_common_AA_MyersPeptide))+'\n')
    sys.stdout.write('             R = '+str(math.sqrt((cov12_MyersPeptide*cov12_MyersPeptide)/(cov11_MyersPeptide*cov22_MyersPeptide)))+'\n')
    sys.stdout.write('-------------------------------------\n')



    # If an argument to this program was given, then we assume it is
    # the name of a file containing a candidate sequence.
    # Assuming that the energy of every amino acid is additive,
    # we can estimate the energy of the candidate sequence.
    # Ignoring correlations, this is proportional to the the negative logarithm
    # of the probability of a sequence of that length evolving at random.

    if len(sys.argv) > 3:
        seqfilename = sys.argv[3]
        file = open(seqfilename,'r')

        for line in file:
            line = line.strip()

            sys.stdout.write('Using energies inferred from the database:\n')
            U, deltaU = CalcSequenceEnergyErr(line,ener,delta_ener)
            sys.stdout.write(str(U*kB*T)+' '+str(deltaU*kB*T)+'\n')

            sys.stdout.write('Using Horovitz1992:\n')
            U = CalcSequenceEnergy(line,helixPropensityHorovitz)
            sys.stdout.write(str(U)+'\n')

            sys.stdout.write('Using ONeil1990:\n')
            U = CalcSequenceEnergy(line,helixPropensityONeil)
            sys.stdout.write(str(U)+'\n')

            sys.stdout.write('Using Blaber44_1993:\n')
            U = CalcSequenceEnergy(line,helixPropensityBlaber44)
            sys.stdout.write(str(U)+'\n')

            sys.stdout.write('Using Myers1997_peptide_pH7:\n')
            U = CalcSequenceEnergy(line,helixPropensityMyersPeptidePh7)
            sys.stdout.write(str(U)+'\n')

            sys.stdout.write('Using Myers1997_protein_pH7:\n')
            U = CalcSequenceEnergy(line,helixPropensityMyersProteinPh7)
            sys.stdout.write(str(U)+'\n')




if __name__ == "__main__":
    main()
