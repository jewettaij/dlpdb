#!/usr/bin/env python

"""
This returns whether or not a PDB file contains all of the required backbone
heavy atoms: " C1'",...," C5'"," O3'"," O4'"," O5'"," P  "," OP1"," OP2"
"""

import sys

res_types = (['  A',
              '  G',
              '  C',
              '  U'])


atoms_found = {' C1\'':False,
               ' C2\'':False,
               ' C3\'':False,
               ' C4\'':False,
               ' C5\'':False,
               ' O3\'':False,
               ' O4\'':False,
               ' O5\'':False,
               ' P  ':False,
               ' OP1':False,
               ' OP2':False}

def main():
    for line in sys.stdin:
        if (line[0:5] == 'ATOM '):
            atom_type=line[12:16]
            res_type=line[17:20]
            if (res_type in res_types):
                #if ((atom_type in atoms_found) and (not atoms_found[atom_type])):
                #    print('found atom_type=\"'+atom_type+'\", in res_type=\"'+res_type+'\"')
                atoms_found[atom_type] = True

    search_criteria_satisfied = True
    for atom_type in atoms_found:
        if (not atoms_found[atom_type]):
            search_criteria_satisfied = False

    if (search_criteria_satisfied):
        exit(0)  # normal termination indicates all atoms were found
    else:
        exit(1)  # non-zero (abnormal termination) exit code indicates this PDB
                 # file does not contain very many heavy atoms in amino acids

if __name__ == "__main__":
    main()
