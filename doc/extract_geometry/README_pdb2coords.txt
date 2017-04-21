pdb2coords.py ATOM_LIST [-blank] < PDB_FILE

This script extracts the coordinates of the atoms named in ATOM_LIST and prints them to the standard out.
For every residue in the PDB_FILE, a line of text containing the coordinates of the requested atoms in ATOM_LIST are printed to the standard out.
(If N atom names appear in the list, then 3 x N numbers will be printed on every line.)
Each line printed by pdb2coords.py contains a list of numbers which is typically passed to other programs which calculate distances, angles, or other geometric quantities, one line at a time.

 ---  DISCARDED OR ABSENT ATOMS: ---

If one of the requested atom types is not present in a residue, then "? ? ?" is printed at the appropriate place in the list.
The "? ? ?" blank symbols are also printed at the beginning and ending of every PDB file whenever the "i+1" symbol is used. (See below.)
Atoms are also discarded if there are multiple alternate locations in the PDB file for the same atom, as indicated by a non-blank alternate location indicator on that line.  These atom coordinates are not used.  See http://deposit.rcsb.org/adit/docs/pdb_atom_format.html
If the optional "-blank" argument is included, then lines which would containing "? ? ?" symbols are replaced with blank lines.

 ---  ATOM_LIST arguments ---

ATOM_LIST is a list of arguments, each one is a quoted 4-character 
atom-type-name, included whitespace. YOU MUST INCLUDE THE QUOTES AND SPACES. 
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
In addition, there is a special atom type:
    "RAVE"   Average position of all non-backbone atoms (ecludes, C CA N O H)
             (For GLY and PRO residues, RAVE coords are undefined & left-blank.)
              

 ---  Extracting coordinates for atoms in subsequent residues  ---    

The coordinates on each line do not all have to come from the same residue.
ATOM_LIST can be separated by the "i+1" symbol.  This indicate that the next atom in the list comes from the next residue.  General separator symbols such as "i-1", "i+2", "i-2" are also supported.  In this way, the same atom's coordinates may appear on multiple lines in different places.  For example:

pdb2coords.py " CA " i+1 " CA " < PDB_FILE

In this example, each line printed by pdb2coords.py contains 6 numbers (or ? blank symbols), the alpha carbon for this residue, and the alpha carbon for the next residue.  (In this example, for a protein which has n residues, n+1 lines will be printed. The first and last lines will contain blanks.)
