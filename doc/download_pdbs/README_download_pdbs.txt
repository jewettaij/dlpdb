---
   This file contains more details about how to use "download_pdbs.py".

--- Overview ---
When attempting to ask statistical questions about proteins, 
it is preferable to use an unbiased sample of PDB files.

    The dlppdb.py script is a very simple script which downloads a list of
PDB files from the official rscb.org repository.

 Typical usage (example):

    download_pdbs.py < list_of_pdb_codes.txt

 Where the input file contains a list of pdb-codes, one per line:
100d
101d
102d
103d
104d
107d
108d
109d
 Every line in this file contains a 4-letter pdb identifier.
Each 4-letter pdb code is compared with a list of pdb codes downloaded so far
(which is initially read from the file "pdbs_most_recent.txt", if present).
If it is new 4-letter pdb code, then the corresponding pdb file is downloaded,
and appended to the end of the "pdbs_most_recent.txt" file.
The corresponding DSSP file is also downloaded, if available.

For every line in the input file, this script
a) downloads a the relevant file in gzipped format, for example
   7odc.pdb.gz
b) unzips the file
c) This script also attempts to download the corresponding dssp file from 
   ftp.cmbi.kun.nl/pub/molbio/data/dssp
   (These DSSP files can be useful for older PDB files which lack secondary
    structure information.  However in most cases, they are not needed.)

This script keeps track of which files have been downloaded so far
so that you can run it again without having to delete your old PDB files.
The following files keep track of this information:

pdbs_most_recent.txt This is a list of PDB files (4-letter codes) which were
                     successfully downloaded the most recent time you ran
                     dlppdb.py.
pdbs_old.txt        A list of PDB files which were downloaded in the current
                    attempt using download_pdbs.py, as well as all earlier
                    attempts.

If this script gets interrupted before completing, you can safely run it again.
