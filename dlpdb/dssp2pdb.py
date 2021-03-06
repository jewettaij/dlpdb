#!/usr/bin/env python

"""
 Usage:

   dssp2pdb.py PDB_FILE < DSSP_FILE > helix_sheet_records.pdb

 This program reads a DSSP file, and generates lines of text which can
 be inserted into a PDB file.  These lines of text are
 "HELIX" and "SHEET" records, which can be interpreted by other programs.

 Information on the DSSP format was obtained from:
 http://swift.cmbi.ru.nl/gv/dssp/DSSP_3.html

 DSSP files for proteins with 4-letter PDB/RCSB identifiers
 can be downloaded from
 ftp://ftp.cmbi.kun.nl/pub/molbio/data/dssp/????.dssp
 For example:
 ftp://ftp.cmbi.kun.nl/pub/molbio/data/dssp/1est.dssp

 Information on the PDB HELIX/SHEET format was obtained from:
 http://www.wwpdb.org/documentation/format32/sect5.html

"""

import sys

# author: Andrew Jewett
g_program_name = __file__.split('/')[-1]
g_date_str = '2012-12-10'
g_version_str = '0.8.0'



def Int2Digits(i, base):
    """
    The next function converts the number into a list of digits
    For example:
    if base=10, then Int2Digits(126, 10) = [1, 2, 6]
    if base=2,  then Int2Digits(126, 2)  = [1, 1, 1, 1, 1, 1, 0]
    (This function is called by the Number2String() function.)
    """

    if (i < base):
        return [i]
    else:
        return Int2Digits(i//base, base) + [i % base]



# The next function converts 1-letter amino acid names into 3-letter names.
# If the user doesn not specify a PDB file (with 3-letter res names), use this.
def ResNamesFrom1Char(c):
    lookup = {'A':'ALA',
              'R':'ARG',
              'N':'ASN',
              'D':'ASP',
              'C':'CYS',
              'E':'GLU',
              'Q':'GLN',
              'G':'GLY',
              'H':'HIS',
              'I':'ILE',
              'L':'LEU',
              'K':'LYS',
              'M':'MET',
              'F':'PHE',
              'P':'PRO',
              'S':'SER',
              'T':'THR',
              'W':'TRP',
              'Y':'TYR',
              'V':'VAL'}
    if (lookup.has_key(c)):
        return lookup[c]
    else:
        return c+c+c



# For some reason, I also have to give each helix/strand a 3-character name.
def Number2String(i):
    str = ''
    for d in Int2Digits(i, 26):
        c = chr(ord('A') + d)
        str += c
    return str




def IsHelix(secondary_type):
    return (secondary_type in 'HIG')

def IsStrand(secondary_type):
    return (secondary_type in 'E')

def SameType(secondary_type1, secondary_type2):
    return (secondary_type1 == secondary_type2)

#def SameType(secondary_type1, secondary_type2):
#    return ((IsHelix(secondary_type1) and IsHelix(secondary_type2))
#            or
#            (IsStrand(secondary_type1) and IsStrand(secondary_type2)))


def PrintPDB(secondary_type,
             start_chainID,
             start_resSeq,
             start_iCode,
             start_resName,
             stop_chainID,
             stop_resSeq,
             stop_iCode,
             stop_resName,
             average_chirality,
             length,
             helix_counter,
             sheet_counter):

    if (IsHelix(secondary_type)):
        comment = ' generated by DSSP & dssp2pdb '
        class_code = ' 1' #PDB code corresponding to a right-handed alpha helix
        if (average_chirality < 0.0):
            class_code = ' 7' # left handed alpha helix
        if (secondary_type == 'I'):
            class_code = ' 3' # right-handed pi helix
        elif (secondary_type == 'G'):
            class_code = ' 5' # right-handed 3-10 helix
        sys.stdout.write("HELIX  "+
                         str(helix_counter).rjust(3)+" "+
                         str(helix_counter).rjust(3)+" "+
                         #Number2String(helix_counter-1).rjust(3)+" "+
                         start_resName+" "+start_chainID+" "+
                         str(start_resSeq).rjust(4)+start_iCode+" "+
                         stop_resName+" "+stop_chainID+" "+
                         str(stop_resSeq).rjust(4)+stop_iCode+
                         class_code + " " + 
                         comment + str(length).rjust(5)+"\n")

    if (IsStrand(secondary_type)):
        num_strands_in_sheet = 1
        sense = ' 0'
        sys.stdout.write("SHEET  "+
                         str(sheet_counter).rjust(3)+" "+
                         #str(sheet_counter).rjust(3)+" "+
                         Number2String(sheet_counter-1).rjust(3)+
                         str(num_strands_in_sheet).rjust(2)+" "+
                         start_resName+" "+start_chainID+
                         str(start_resSeq).rjust(4)+start_iCode+" "+
                         stop_resName+" "+stop_chainID+
                         str(stop_resSeq).rjust(4)+stop_iCode+
                         sense+
                         "                                   "+"\n")
    



def main():

    min_length_strand = 0
    min_length_helix  = 0

    resNamesFromPDBid = {}

    if (len(sys.argv) == 2):
        pdb_file_name = sys.argv[1]
        # Then open the PDB file and load the names of each residue 
        for line in open(pdb_file_name,'r'):
            if (line[0:6] == 'ATOM  '):
                # the next 3 data identify the residue to which the atom belongs
                chainID = line[21]
                resSeq  = int(line[22:26])
                iCode   = line[26]
                # this next short string specifies the type of this residue
                resName = line[17:20]

                resNamesFromPDBid[(chainID, resSeq, iCode)] = resName


    elif (len(sys.argv) != 1):
        sys.stderr("Error: wrong number of arguments.\n")
        exit(-1)



    # Now, read the DSSP file:

    dssp_lines = sys.stdin.readlines()
    prev_secondary_type = ' '



    # First, skip over the description at the begining of the DSSP file
    skip = 0
    for i in range(0,len(dssp_lines)):
        line = dssp_lines[i]
        if (line[0:25] == '  #  RESIDUE AA STRUCTURE'):
            skip = i+1

    dssp_lines = dssp_lines[skip:]


    helix_counter  = 0
    sheet_counter  = 0


    # Finally, scan over the (remaining lines of the) DSSP file:

    for i in range(0,len(dssp_lines)):

        line = dssp_lines[i]
        #sys.stderr.write("line = \""+line.rstrip()+"\"\n")

        line_contains_junk = False
        line_contains_digits = False
        for c in line[6:10]:
            if (c in '0123456789'):
                line_contains_digits = True
            elif (c != ' '):
                line_contains_junk = True
        if (not line_contains_digits):
            line_contains_junk = True


        if (not line_contains_junk):

            # Lookup the 3-part identifier for this residue (used in PDB files)
            chainID = line[11]       #To which chain does this residue belong?
            resSeq  = int(line[6:10])#This integer gives a rough idea of the 
                                     #location of a residue in the chain's sequence.
            iCode   = line[10]       #A 1-character code for inserting residues

            secondary_type = line[16]
            chirality_code = line[22]

            # What is the name of this residue? (This should be a 3-character string)
            if (len(resNamesFromPDBid) > 0):
                resName = resNamesFromPDBid[(chainID,resSeq,iCode)]#lookup from PDB file
            else:
                resName = ResNamesFrom1Char(line[13]) #convert 1-lttr to 3-lttr code

            # Check for my stupidity
            if (len(resName) != 3):
                sys.stderr.write("Internal Error: Illegal residue name: \""+name+"\"\n"\
                                     "  for residue (chain\""+chainID+"\", res="+str(resSeq)+" iCode=\""+iCode+"\")\n"\
                                     "  Please report this bug (along with input files).\n")
                exit(-1)



            # Does this new residue mark the end of a helix or a strand?

            # If the previous residue belonged to a helix or strand...
            if (IsHelix(prev_secondary_type) or
                IsStrand(prev_secondary_type)):

                # ...and if the current residue does -not- also share the same
                # kind of secondary structure (either helix or strand),
                # ...(or if it belongs to a different chain)
                # ...THEN this new residue marks the end of that helix/strand
                if ((not SameType(secondary_type, 
                                  prev_secondary_type)) or
                    (chainID != prev_chainID)):

                    # ... Check for one more thing:
                    #     if the helix/sheet is too short, then ignore it
                    if (IsHelix(prev_secondary_type) and
                        (length < min_length_helix)):
                        helix_counter -= 1 #-1 because we should not have counted it
                    elif (IsStrand(prev_secondary_type) and
                        (length < min_length_strand)):
                        sheet_counter -= 1 #-1 because we should not have counted it
                    else:
                        # If there are no problems, print out that helix/strand
                        PrintPDB(prev_secondary_type,
                                 start_chainID,
                                 start_resSeq,
                                 start_iCode,
                                 start_resName,
                                 prev_chainID,
                                 prev_resSeq,
                                 prev_iCode,
                                 prev_resName,
                                 total_chirality / length,
                                 length,
                                 helix_counter,
                                 sheet_counter)
                else:
                    # Keep track of the length of each helix/strand
                    length += 1
                    # Was it right handed '+' or left handed '-' ?
                    if (chirality_code != '-'):
                        total_chirality += 1.0
                    else:
                        total_chirality -= 1.0




            # Is this the beginning of a new helix/strand?
            if ((IsHelix(secondary_type) or 
                 IsStrand(secondary_type)) and
                (not SameType(secondary_type, prev_secondary_type))):

                start_resSeq  = resSeq
                start_iCode   = iCode
                start_chainID = chainID
                start_resName = resName
                length = 1
                if (IsHelix(secondary_type)):
                    helix_counter += 1

                if (IsStrand(secondary_type)):
                    sheet_counter += 1

                total_chirality = 0.0


            # keep track of the properties of this residue
            prev_secondary_type = secondary_type
            prev_resSeq    = resSeq
            prev_iCode     = iCode
            prev_chainID   = chainID
            prev_resName   = resName




    # If we reached the end of the file / molecule, then we have obviously
    # reached the end of whatever helix or strand which was located at the end of
    # that molecule.  We should get around to printing out that helix/strand now:

    if (IsHelix(secondary_type) or IsStrand(secondary_type)):
        PrintPDB(prev_secondary_type,
                 start_resSeq,
                 start_iCode,
                 start_chainID,
                 start_resName,
                 prev_resSeq,
                 prev_iCode,
                 prev_chainID,
                 prev_resName,
                 total_chirality / length,
                 length,
                 helix_counter,
                 sheet_counter)

if __name__ == "__main__":
    main()
