#!/usr/bin/env python

"""
This returns whether or not a PDB file contains all of the require
backbone heavy atoms: ' CA ', ' CB ', ' C  ', ' N  ', ' O  '
"""

import sys

res_types = (['ALA',
              'ARG',
              'ASN',
              'ASP',
              'CYS',
              'GLU',
              'GLN',
              'GLY',
              'HIS',
              'ILE',
              'LEU',
              'LYS',
              'MET',
              'PHE',
              'PRO',
              'SER',
              'THR',
              'TRP',
              'TYR',
              'VAL'])


atoms_found = {' CA ':False,
               ' CB ':False,
               ' C  ':False,
               ' N  ':False,
               ' O  ':False}

def main():
    for line in sys.stdin:
        if (line[0:5] == 'ATOM '):
            atom_type=line[12:16]
            res_type=line[17:20]
            #print('atom_type=\"'+atom_type+'\", res_type=\"'+res_type+'\"')
            if (res_type in res_types):
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
