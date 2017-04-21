The main script in this directory is "pdb2coords.py".
  It extracts coordinates from every residue according to atom name,
  (and sorts them according to chain, residue-number, and insert-code).

The documentation for this script is the "README_pdb2coords.txt" file.

Most of the other scripts in this directory are simple scripts which 
calculate distances and angles BETWEEN ATOMS WHICH ARE NEARBY IN THE CHAIN
using the coordinates extracted by "pdb2coords.py".
(As of 2012-10-24, I have not yet added the capability to calculate distances 
 or angles between sequence-distant atoms, tertiary contacts.)

-- Most of these scripts allow you to specify the atoms you want to use in
-- the calculation.  You must refer to them by their full 4-character 
-- PDB atom name.  These atom names typically contain spaces.

Common atoms are shown below
    " CA "   alpha-carbon atom
    " CB "   beta-carbon atom
    " N  "   nitrogen backbone atom
    " O  "   oxygen backbone atom
    " H  "   amide hydrogen
    " C  "   carbon backbone atom
The names for other common PDB atom types are listed in 
JL Markley et. al, Pure & Appl. Chem., 70(1):117-142 (1998)
(Available at http://icnm.cerm.unifi.it/iupac.pdf)
See the documentation for README_pdb2coords.txt for details.
