#!/usr/bin/env python

"""
This program extracts a series of x,y,z coordinates from a PDB file (stdin).
The user must supply a list of 4-character atom names. (Usually these names 
include spaces and thus must be enclosed in quotes when passed to the shell.)

For every residue in the PDB file (which falls within the user's selection),
the coordinates of these atoms are printed to the standard out (stdout),
one line per residue (multiple atoms per line).  The order of the atoms
printed on each line should match the order they were passed as arguments
to this program (not the order they were read from the pdb file).
If any of the requested atom names are not present in the PDB file, then an
"? ? ?" is printed instead at the appropriate location on that line.
(...unless the "-blank" argument is passed, in which case a blank line is 
printed.  Either way the # of output lines equals the # of residues.)

Results can be limited to a specific residue interval in the PDB file.
Intervals are pairs of amino acids signifying the first and last
residue in the interval.  Specifying an interval is optional.
In order to specify any amino acid in a PDB file, you must provide 3
identifiers:
  the ChainID a single letter specifying
  the SeqNum an integer indicating the location within that chain
  the ICode ("insert code" usually 0.  I don't know why this number is 
             necessary. ..For the record, I never loved the PDB file format.)

All 3 identifiers are needed for both the starting and ending residues.
Consequently this program expects 6 arguments: 
   (to be read from the command line)
ChainID_first SeqNum_first ICode_first ChainID_last SeqNum_last ICode_last

 ---  Mixing atoms from different residues together ---
Sometimes bonded interactions occur between atoms in in different 
residues from the same chain.  (For example, the " C  " atom
of one residue binds to the " N  " of the next residue.)
It's convenient to print the coordinates of the atoms involved in these 
interactions together on the same line.  To do that, insert the "i+1"
argument before the " N  " argument, as in:
 " C  "  i+1  " N  "
After the i+1 argument, all of the atoms following it are assumed to come
from the next residue.  (The " N  " atom in this example.)
(Note: If for some reason you want the third atom to come from the original
 residue, you can switch back to it by adding "i+0" after the " N  " atom:
 " C  "  i+1  " N  " i+0 " CA  "
 More generally, you can also use "i+2", or "i+3", ..., and "i-1", "i-2",...
 to mix atoms from more distant residues together.)

Sometimes this program pipes its output to other programs which halt early.
Below we silently suppress the ugly "Broken pipe" message this generates:
"""

# author: Andrew Jewett
g_program_name = __file__.split('/')[-1]
g_date_str = '2017-4-21'
g_version_str = '0.5.1'



import sys
from operator import attrgetter
from math import sqrt
# Sometimes this program pipes its output to other programs which stops reading
# the PDB file prematurely (such as when multiple MODEL records are present).
# Below we silently suppress the ugly "Broken pipe" message this generates:
import signal
signal.signal(signal.SIGPIPE, signal.SIG_DFL)

try:
    from .resid import *
except ImportError:
    from resid import *


# Ignore atoms on the backbone (other than CA), 
RAVE_exclude_atoms = [' N  ', ' H  ', ' C  ', ' O  ', ' CA ']
RAVE_exclude_residues = set(['GLY', 'PRO']) #Glycine residues do not have CB atoms


def main():
    atoms_needed = []
    atoms_res_offsets = []
    atoms_res_offset = 0
    omit_incomplete = False
    use_all_residues = True
    firstR = None
    lastR = None

    final_range_a = None
    final_range_b = None
    final_slice_incr = 1



    if len(sys.argv) > 1:
        i = 1
        while i < len(sys.argv):

            if ((sys.argv[i] == '-blank') or (sys.argv[i] == 'blank')):
                omit_incomplete = True
                i += 1

            elif ((sys.argv[i] == 'i') or
                  ((sys.argv[i][:1] == 'i') and (len(sys.argv[i]) >= 3) and
                   (sys.argv[i][1:2] in ['+','-','=']))):
                if sys.argv[i] == 'i':
                    atoms_res_offset = 0
                else:
                    atoms_res_offset = int(sys.argv[i][1:].lstrip('='))
                i += 1

            elif ((len(sys.argv[i]) >= 2) and
                  (sys.argv[i][0] == '[') and (sys.argv[i][-1] == ']')):

                tokens=(sys.argv[i][1:-1]).split(':')
                if len(tokens[0].strip()) > 0:
                    final_range_a = int(tokens[0])
                if len(tokens[1].strip()) > 0:
                    final_range_b = int(tokens[1])
                if len(tokens[2].strip()) > 0:
                    final_slice_incr = int(tokens[2])
                i += 1

            elif (len(sys.argv[i]) == 4):
                # Keep track of the atom type name and the order it appeared
                atoms_needed.append(sys.argv[i])
                # Keep track of from which residue it comes from
                atoms_res_offsets.append(atoms_res_offset)
                i += 1

            else: # if not is_digit(sys.argv[i][0])

                # If the next argument is a number, then interpret this number
                # as a residue sequence number
                if len(sys.argv) < i+6:
                    sys.stderr.write("Error: Not enough arguments or argument type error:\n"
                                     "       Offending arguemt #"+str(i)+": \""+sys.argv[i]+"\"\n")

                    if len(sys.argv[i]) == 1:
                        sys.stderr.write("\n"
                                         "       Is argument \""+sys.argv[i]+"\" a chain ID letter?\n"
                                         "       (Chain IDs are only passed as arguments when you want to limit the\n"
                                         "        residues considered to within an interval in the PDB file.\n"
                                         "        To specify an interval, you must provide 5 more aruments.  See below.)\n"
                                         "\n"
                                         "       Note: one-letter residue codes (eg A,V,L,...) are not considered\n"
                                         "             valid residue types by this program.\n"
                                         "             You must specify a 3-letter equivalent (eg ALA, VAL, LYS,...).\n"
                                         "             One-letter arguments are interpreted as chain-IDs. (See below.)\n"
                                         "       \n"
                                         "       By default, the the sequence is extracted from the entire PDB file.\n"
                                         "       In that case, no arguments are required.\n"
                                         "       Alternately, you can limit the selection to a single interval of\n"
                                         "       residues from one of the chains in the PDB file.\n"
                                         "       To specify an interval, you must passing 6 arguments to this program.\n"
                                         "       This program requires a pair of residues to designate the first and\n"
                                         "       last members of the interval.  Each residue requires 3 identifiers.\n"
                                         "       Consequently the six arguments needed are:\n"
                                         "\n"
                                         "ChainID_first SeqNum_first ICode_first ChainID_last SeqNum_last ICode_last\n"
                                         "\n")
                    else:
                        sys.stderr.write("\n"
                                         "Note: Atom-type names should be exactly 4 characters long, and typically\n"
                                         "      contain spaces.  For example: \" C  \", \" O2 \", or \" C1\'\"\n"
                                         "\n"
                                         "Note: PDB files are not passed as arguments\n"
                                         "      but are read from the standard input, for example using the notation\n"
                                         "\n"
                                         "      "+sys.argv[0]+" < file.pdb\n")
                    exit(-1)
                else:
                    use_all_residues = False
                    firstR = ResID(sys.argv[i], int(sys.argv[i+1]), sys.argv[i+2])
                    lastR  = ResID(sys.argv[i+3], int(sys.argv[i+4]), sys.argv[i+5])
                    #sys.stderr.write('  Interval selected: (\"'+firstR.chainID+'\", '+str(firstR.seqNum)+', \"'+firstR.iCode+'\") ... (\"'+lastR.chainID+'\", '+str(lastR.seqNum)+', \"'+lastR.iCode+'\")\n')
                    i += 6

    resID2coords   = {}
    resID2CrdTot   = {}
    resID2CrdSqTot = {}
    resID2CrdNum   = {}
    model_ID = None

    for line in sys.stdin:
        #print("\""+line.rstrip()+"\"")
        if line[0:6] == "MODEL ":
            if model_ID == None:
                model_ID = line[10:14]
            else:
                # Otherwise, if it's not the first model quit
                sys.stderr.write('  Warning(pdb2coords.py): Omitted alternate models from pdb file.\n')
                break
        elif line[0:6] == "ATOM  ":
            #atomID    = int(line[6:11])
            atomType  = line[12:16]
            altLoc    = line[16:17]
            resType   = line[17:20]
            chainID   = line[21:22]
            resSeq    = line[22:26]
            iCode     = line[26:27]
            resID = ResID(chainID, int(resSeq), iCode)
            #xyz_str   = line[30:54]<-not safe. spaces not always in between numbers
            x_str     = line[30:38].strip()
            y_str     = line[38:46].strip()
            z_str     = line[46:54].strip()
            if (use_all_residues or ((firstR<=resID) and (resID<=lastR))):
                if resID not in resID2coords:
                    resID2coords[resID]   = [None for i in range(0,len(atoms_needed))]
                    resID2CrdTot[resID]   = [0.0, 0.0, 0.0]
                    resID2CrdSqTot[resID] = [0.0, 0.0, 0.0]
                    resID2CrdNum[resID]   = 0

                if altLoc == ' ':  # (Ignore "alternate" atom records)
                    for i in range(0, len(atoms_needed)):

                        if atomType == atoms_needed[i]:
                            resID2coords[resID][i] = [x_str, y_str, z_str]

                        if ((not (atomType in RAVE_exclude_atoms)) and
                            (not (resType in RAVE_exclude_residues))):
                            resID2CrdNum[resID] += 1
                            crd_tot = resID2CrdTot[resID]
                            crd_tot[0] += float(x_str)
                            crd_tot[1] += float(y_str)
                            crd_tot[2] += float(z_str)
                            crd_sq_tot = resID2CrdSqTot[resID]
                            crd_sq_tot[0] += float(x_str)**2
                            crd_sq_tot[1] += float(y_str)**2
                            crd_sq_tot[2] += float(z_str)**2



    # Extract an (unordered) list of the resIDs of the residues in the sequence
    resIDs = [resID for resID in resID2coords]

    # Residues in PDB files are often not listed in order.
    # Consequently, we must sort the list by chainID, seqNum, and finnaly iCode:
    sequence_of_resIDs = sorted(resIDs, key=attrgetter('chainID','seqNum','iCode'))

    positions        = [ [None for j in range(0,len(atoms_needed))] 
                         for i in range(0,len(sequence_of_resIDs)) ]

    assert(len(atoms_needed) == len(atoms_res_offsets))

    if len(atoms_needed) > 0:
        min_offset = None
        max_offset = None
        for offset in atoms_res_offsets:
            if ((min_offset == None) or (min_offset > offset)):
                min_offset = offset
            if ((max_offset == None) or (max_offset < offset)):
                max_offset = offset
    else:
        sys.stderr.write('  Warning(pdb2coords.py): NO ATOM TYPES SELECTED.\n')
        min_offset = 0
        max_offset = 0

    # If the user requested any offsets (so that some atoms on the same line 
    # of output belong to different residues), then this means that by
    # default (if omit_incomplete is False), additional lines of coordinate
    # data will be printed which contain references to missing coordinates
    # (because some of the needed atoms are out of index range).
    # The "num_extra" variable equals the number the number of additional lines
    # of output.
    num_extra = max_offset - min_offset
    N = len(sequence_of_resIDs)

    for i in range(0, len(atoms_res_offsets)):
        atoms_res_offsets[i] -= min_offset

    sorted_positions = [ [None for j in range(0,len(atoms_needed))] 
                         for i in range(0, N+num_extra) ]

    for i in range(0, N):
        resID = sequence_of_resIDs[i]
        positions[i] = resID2coords[resID]

        # now deal with atoms of type "RAVE"
        for j in range(0, len(atoms_needed)):
            if atoms_needed[j] == 'RAVE':
                RAVE_num_atoms = resID2CrdNum[resID]
                if RAVE_num_atoms > 0:
                    positions[i][j] = [str(resID2CrdTot[resID][0] / RAVE_num_atoms),
                                       str(resID2CrdTot[resID][1] / RAVE_num_atoms),
                                       str(resID2CrdTot[resID][2] / RAVE_num_atoms)]
                else:
                    assert(positions[i][j] == None)
            elif atoms_needed[j] == 'RGYR':
                RAVE_num_atoms = resID2CrdNum[resID]
                if RAVE_num_atoms > 0:
                    rave   = ((resID2CrdTot[resID][0] / RAVE_num_atoms),
                              (resID2CrdTot[resID][1] / RAVE_num_atoms),
                              (resID2CrdTot[resID][2] / RAVE_num_atoms))
                    rsqave = ((resID2CrdSqTot[resID][0] / RAVE_num_atoms),
                              (resID2CrdSqTot[resID][1] / RAVE_num_atoms),
                              (resID2CrdSqTot[resID][2] / RAVE_num_atoms))
                    rgyrsq = 0.0
                    for d in range(0,3):
                        rgyrsq += rsqave[d] - (rave[d]*rave[d])
                    rgyr = sqrt(rgyrsq)
                    positions[i][j] = [str(rgyr), '', '']
                    if rgyr == 0.0:
                        print('resID=\"'+str(resID)+'\"')
                        exit(0)
                else:
                    assert(positions[i][j] == None)


    for i in range(-num_extra, N):
        for j in range(0,len(atoms_needed)):
            I = i + atoms_res_offsets[j]
            assert(I >= i)
            if (0 <= I) and (I < N):
                sorted_positions[i+num_extra][j] = positions[I][j]

    if final_range_a == None:
        final_range_a = 0
    else:
        final_range_a += max_offset

    if final_range_b == None:
        final_range_b = N+num_extra
    else:
        final_range_b += max_offset

    for i in range(final_range_a, final_range_b, final_slice_incr):
        coords_str_list = []
        for j in range(0,len(sorted_positions[i])):
            if sorted_positions[i][j] == None:
                if omit_incomplete:
                    coords_str_list = []
                    break
                else:
                    coords_str_list += ['?','?','?']
            else:
                coords_str_list += sorted_positions[i][j]

        # Finally, write out the coordinates.
        sys.stdout.write(' '.join(coords_str_list) + '\n')


if __name__ == "__main__":
    main()
