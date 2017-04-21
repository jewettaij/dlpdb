  --- dlpisces Overview ---
When attempting to ask statistical questions about proteins, 
it is preferable to use an unbiased sample of PDB files.  Unfortunately,
the raw PDB library contains many redundant and similar PDB structures,
which can weight the statistics towards proteins which are easy to crystallize
(or heavily studied proteins with many artificial mutant variants).

    The dlpisces.py script downloads a "representative sample" of the entire 
PDB library (using the filtering method provided by the PISCES server).
The PISCES server attempts to remove most of this bias.
  -------------------------

1) Download a list of the names of culled pdb files from the pisces server

http://dunbrack.fccc.edu/PISCES.php

I usually choose from a list of previously-generated files.
The names of the files indicate the desired resolution, allowable sequence
overlap, R-value limits, and date.

2) Unzip the file (if necessary)

3) ---- run dlpisces.py ----

  Download the PDB files on this list using "dlpisces.py".  Typical usage:

  dlpisces.py < cullpdb_pc30_res1.6_R0.25_d120723_chains2924

  (The name of your "cullpdb_" file will be different.)
  This file contains a list of chains from PDB files with less than 30% sequence
  identity, a resolution better than 1.6Angstroms, and an R-value below 0.25.

  For every line in the "cullpdb_" file, this script
  a) downloads a the pdb relevant pdb file, for example
     7odc.pdb.gz
  b) unzips the file
  c) extracts the relevant chain from the pdb file and creates a new pdb file:
     7odc_chainA.pdb
  d) This script also downloads the corresponding dssp file from 
     ftp.cmbi.kun.nl/pub/molbio/data/dssp
     These DSSP files are useful for older PDB files which 
     lack secondary structure information.

  If this script gets interrupted before completing, you can safely 
  run it again, and it should continue from where it left off.
  dlpisces.py creates the following file:

  pdb_files.txt     This stores a list of names of PDB files created by
                    dlpisces.py which contain only the chains referred to
                    in the "cullpdb" file downloaded from the pisces server.
                    THESE ARE THE PDB FILES YOU WOULD PASS TO OTHER SCRIPTS
                    WHICH ANALYZE STATISTICAL PATTERNS IN THE PDB DATABASE.
                    THESE FILES CONTAIN A "REPRESENTATIVE SAMPLE" OF THE 
                    CHAINS IN THE PDB LIBRARY.  
                    (You may need to use use other criteria to further
                    refine the list of pdb files you want, but this is
                    a good starting point.)

  See "README_dlpisces.txt" for details.


---------------------------------------------------------------------------
---------        END OF MAIN DOCUMENTATION
--------- At this point, you are "done" downloading the PDB files.
--------- It's not usually necessary to worry about steps 4 and 5 below.
---------------------------------------------------------------------------


4) -- Optional: Eliminate the PDB files you don't want, --
   --           and clean up the existing PDB files     --

 a) Do the PDB-files lack heavy atoms?
    Scan all the pdb files to see if they lack heavy atoms.
    If so, move the ones that do to a separate directory.
    To do this, run:

    move_missing_heavy_atoms.sh pdbs_missing_heavy_atoms < pdb_files.txt

    (Note: This step may not be necessary.  I think the PISCES server allows 
           you to filter out PDB files lacking heavy atoms.)

 b) Do the PDB-files contain secondary structure information (helices & sheets)?
    (If not, throw them away or use the DSSP files to generate the information.)
    The following script moves the files which lack HELIX or SHEET records
    to a subdirectory named "pdbs_no_helices_or_sheets".  To do this, run

  move_missing_secondary_str.sh pdbs_no_helices_or_sheets < pdb_files.txt

    If you want to, you can use the dssp files to add HELIX and SHEET records.

  cd pdbs_no_helices_or_sheets
  mv ../*.dssp .
  ls -f1 *.pdb | replace_missing_secondary_str.sh
  cd ../

    (Personally, I don't bother doing this.  There are usually very few 
     PDB files with this problem.  I just throw them away, along with the
     .dssp files.  Furthermore, recently (2015) I've noticed that sometimes
     the old dssp server fails to respond, and the .dssp files fail to 
     download.  I suspect that this server may be taken offline soon.)


 c) Do the PDB-files contain membrane proteins?
    (I usually throw them away if the contain the words MEMBRANE or MICELLE.)
    To move them to a separate directory ("pdbs_membrane") run:

  move_membrane_proteins.sh pdbs_membrane < pdb_files.txt


 d) Were the PDB-files obtained by NMR?
    Usually, the PISCES server excludes structures obtained by NMR.
    To be sure, you can check which structures were obtained by XRAY, NMR, or
    some other method by reading the cullpdb file (the third column).
    If instead you only want structures obtained by NMR, you can 
    try something like this:

    awk '{if ($3=="NMR") print $0}' < cullpdb_orig > cullpdb_nmr_only
    dlpisces.py < cullpdb_nmr_only

    Make sure you manually instruct the PISCES server to include 
    PDB files which were obtained by other methods besides XRAY
    (by default NMR structures are excluded).

    Alternately, you can use the "move_nmr_structures.sh" script, this way:
    ls -f1 *.pdb | move_nmr_structures.sh pdbs_nmr


   (Note: single-chain PDB files ending in "_chain?.pdb" also contain all 
          descriptive remarks, from the original PDB file, so you can either 
          scan for text in the original PDB files, or the "_chain?.pdb" files, 
          or both.)

5) ---- Optional: Make a final list of the PDB files you want to use. -----

If this is the first time you have run "dlpisces.py", you're done. Stop here.

If you have run dlpisces.py before in the same directory using a different list,
then it is possible that there may be other PDB files/chains present in that 
directory as well, (which aren't referred to in the current "cullpdb" file).
If so, you can create a new file ("pdb_files_final.txt") 
which excludes these PDB files this way:




rm -f pdb_files_final.txt
for file_pdb in pdb_files.txt; do
    if [ -s "${file_pdb}" ]; then
        echo file_pdb >> pdb_files_final.txt
    fi
done



(You can use a similar for-loop to move these files to a different directory.
 Note: These commands assume you are using the bourne-shell.
 If the loop above does not work, then try typing:
   echo $SHELL
 If it returns:
   tcsh
 then type /bin/sh or /bin/bash first and try again.)

