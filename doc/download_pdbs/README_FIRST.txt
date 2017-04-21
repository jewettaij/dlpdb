This directory contains the following scripts for 
downloading classifying, and cleaning up PDB files.

dlpisces.py      Download a subset of the PDB library using the PISCES server
                 (See README_dlpisces.txt)

dssp2pdb.py      Convert helix/sheet/turn records from dssp file into PDB format
                 (README_dssp2pdb.txt)

move_missing_heavy_atoms.sh    Move PDB files missing heavy atoms somewhere else

move_nmr_structures.sh         Move PDB files containing the string " NMR"

move_membrane_proteins.sh      Move PDB files with "MEMBRANE" or "MICELLE"

move_missing_secondary_str.sh  Move PDB files missing helix/sheet records

------- Typical usage for the "move" scripts: -----------
You pass these scripts a list of PDB files (through the standard-in)
and supply an argument which is the name of the directory where you
want these PDB files to go.  (The script will attempt to create it
if it does not exist.)  For example:

ls -f1 *.pdb | move_membrane_proteins.sh pdbs_membrane 
ls -f1 *.pdb | move_missing_heavy_atoms.sh pdbs_missing_heavy_atoms
ls -f1 *.pdb | move_missing_secondary_str.sh pdbs_no_helices_or_sheets
ls -f1 *.pdb | move_nmr_structures.sh pdbs_nmr

The effected PDB files will be moved to these directories.
