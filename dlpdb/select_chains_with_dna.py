#!/usr/bin/env python

"""
This program finds each chain in a PDB file containing DNA residues (with 
backbone atoms: " C1'",...," C5'"," O3'"," O4'"," O5'"," P  "," OP1"," OP2")
... and create a new file containing only that chain.
PDB-files lacking DNA residues will be ignored.
The files which are created will match the name of the original file, with
the .pdb or .PDB extension stripped off (if present), and replaced with 
_?.pdb,  Where "?" is the chain letter.
"""

import sys

g_program_name  = __file__.split('/')[-1]
g_date_str     = '2015-8-17'
g_version_str  = '0.5.0'


usage_string = "Usage:\n\n  "+g_program_name+" file1.pdb [file2.pdb file3.pdb ...]\n\n"

usage_string += \
"""
Explanation:

"""
usage_string += __doc__


res_types = ([' DA',
              ' DG',
              ' DC',
              ' DT'])


heavy_atoms = (' C1\'',
               ' C2\'',
               ' C3\'',
               ' C4\'',
               ' C5\'',
               ' O3\'',
               ' O4\'',
               ' O5\'',
               ' P  ',
               ' OP1',
               ' OP2')


def ExtractChainText(file_name):
    """
    ExtractChainText() finds all of the chains in a PDB file which contain DNA
    nucleotides and creates a new PDB file for that chain with a similar name.

    """

    pdb_file = open(file_name, 'r')

    chain_atoms_found = {}
    chain_text = {}

    for line in pdb_file:
        if line[0:6] == 'ATOM  ':
            chainID = line[21:22]
            if (not chainID in chain_atoms_found):
                chain_atoms_found[chainID] = {}
                chain_text[chainID] = []
                for atom_name in heavy_atoms:
                    chain_atoms_found[chainID][atom_name] = False

            chain_text[chainID].append(line)

            atom_type = line[12:16]
            res_type = line[17:20]
            if res_type in res_types:
                #if ((atom_type in atoms_found) and (not atoms_found[atom_type])):
                #    print('found atom_type=\"'+atom_type+'\", in res_type=\"'+res_type+'\"')
                chain_atoms_found[chainID][atom_type] = True


    for chainID in chain_atoms_found:
        search_criteria_satisfied = True
        for atom_type in chain_atoms_found[chainID]:
            if (not chain_atoms_found[chainID][atom_type]):
                search_criteria_satisfied = False
        if search_criteria_satisfied:
            sys.stderr.write("  Chain \""+chainID+"\" contains DNA.\n")
            # Then create a new PDB file with a name similar to the original:
            pdb_file_chain_name = file_name
            i = pdb_file_chain_name.lower().rfind('.pdb')
            if i != -1:
                pdb_file_chain_name = (pdb_file_chain_name[:i] +
                                       '_' + chainID +
                                       pdb_file_chain_name[i:])
            else:
                pdb_file_chain_name = file_name + '_' + chainID
            sys.stderr.write('    Creating file \"'+pdb_file_chain_name+'\"\n')
            pdb_file_chain = open(pdb_file_chain_name, 'w')
            pdb_file_chain.write(''.join(chain_text[chainID]))
            pdb_file_chain.close()

    pdb_file.close()



def main():
    if len(sys.argv) < 2:
        sys.stderr.write("Error: Expected at least one argument.\n")
        sys.stderr.write(usage_string)
        exit(1)

    sys.stderr.write(g_program_name+' v'+g_version_str+' '+g_date_str+'\n')

    for file_name in sys.argv[1:]:
        sys.stderr.write('  Reading file \"'+file_name+'\"\n')
        ExtractChainText(file_name)

if __name__ == "__main__":
    main()
