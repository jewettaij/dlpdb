    dssp2pdb.py

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

