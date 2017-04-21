This program extracts the information from a PDB file
which lies within an interval in the PDB file.

Usage:

select_interval.py select_interval.py A 3 " " A 10 z < 309d.pdb > 309d_A3-10.pdb

(This selects residues 3-10 from chain A of the PDB file)

Intervals are pairs of amino acids signifying the first and last
residue in the interval.
In order to specify any amino acid in a PDB file, you must provide 3
identifiers:
  the ChainID a single letter specifying (typically " " or A-Z)
  the SeqNum an integer indicating the location within that chain
  the ICode (a single character code for residue insertions.
             Usually this charcter is the space character " ".
             The smallest possible value is " ".  The largest is "z".
             For the record, I hate this file format.)

All 3 identifiers are needed for both the starting and ending residues.
Consequently this program expects 6 arguments: 
   (to be read from the command line)
ChainID_first SeqNum_first ICode_first ChainID_last SeqNum_last ICode_last
