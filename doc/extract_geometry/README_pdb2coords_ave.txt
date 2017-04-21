 This is a crude script which prints out the coordinates for the average 
 position of each residue in a PDB file (three numbers per line).
 (Note: This program does not calculate the center-of-mass.  All atoms 
        are given the same weight when calculating the average position.)

 Typical usage: 

 pdb2coords_ave.py < PDB_FILE

 Exceptions:
   Backbone atoms (other than alpha-carbon " CA " atoms) are ignored.
   Glycine residues are also ignored.

 These settings can not be changed without editing the code for this script.
 (unlike other scripts whose behavior can be changed by command-line arguments)

 Output:
   The number of lines of output should match the number of residues.
   For glycines (and residues lacking eligible atoms) a blank line is printed.
   For all other residues, the x, y, and z coordinates of the average
   position are printed (3 numbers), one line per residue.
