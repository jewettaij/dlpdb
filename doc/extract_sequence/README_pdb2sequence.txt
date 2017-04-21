    pdb2sequence.py

 Typical usage:

    pdb2sequence.py < file.pdb

    (..where "file.pdb" is the name of a PDB file, eg: 3FG9.pdb.  Note: You 
     must download the PDB file yourself.  This program will not retrieve it.)

 This program extracts the sequence of residue types from a PDB file, 
 and prints it the the standard out (all on one line).
 Each residue is converted to a 1-letter amino-acid code.
 The 20 standard amino-acids, and 5 nucleic acid types are supported.
 Unknown residue types are denoted 'x'.
 (To add additional residue types, edit the code and add entries to the
 "three2one" lookup table.)

 WARNING: Residues represented as HETATM records will be ignored (skipped)!

     Optional interval arguments:

 Users can limit the PDB file to residues which lie within an interval.

 Intervals are pairs of amino acids signifying the first and last
 residue in the interval.
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

