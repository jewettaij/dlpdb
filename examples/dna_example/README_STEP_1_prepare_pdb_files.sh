# This README file contains a combination of instructions and bash commands to
# download and select PDB files containing nothing but duplex DNA, and convert
# them into a format that other tools in the "dlpdb/extract_geometry" directory
# can understand (such as "extract_distances.sh", "extract_angles.sh", ...)
# This README file also contains instructions to discard PDB files derived from
# NMR, PDB files containing non-standard nucleotides, PDB files with redundant
# nucleotide sequences, as well as PDB files containing other molecules (which
# might disturb the DNA conformation)
#
# WARNING: For DNA PDB files, there is no way to automate this process entirely.
#          At the end, you must verify each structure visually and discard
#          the files which are obviously not naked duplex B-form DNA.
#          (The goal here was to just reduce the number of PDB files you have
#            to check manually as much as possible.)
#
#   -----------------------
#
# Note: These commands assume you are using the bourne-shell.
# If any of the for-loops fail to work, then try typing this into the terminal:
#   echo $SHELL
# If it returns:
#   tcsh
# enter this command into the terminal:
#   bash
#
#   -----------------------
#
# First, make sure the "dlpdb" tools are installed.  We will use them later:


cd dlpdb/
source README_install.sh
cd ../

#   -----------------------
#
# Figure out what PDB files you want to download.
#
# For a long list of PDB files containing duplex DNA (2-DNA chains),
# go to:
#https://www.cbi.cnptia.embrapa.br/cgi-bin/SMS/pdb_metrics/FinalSearch.cgi?nr_DNA_chains=2
# Click on "Show the list of all PBDs names" to get a list of PDB codes.
# (Incidentally, a more general list of PDB files can be found here
#  http://www.cbi.cnptia.embrapa.br/SMS/pdb_metrics)
#
# Note: This list contains some PDB files which contain RNA not DNA.
#       The number of PDB files which actually have DNA chains in them is
#       lower than the number of PDB files containing the word "DNA", which
#       you can calculate using:
# grep DNA *.pdb | awk -F':' '{ if (NR==1) {print $1; x=$1} else { if ($1!=x) {print $1; x=$1}}}' | wc
#
#
# I saved this list of 4-letter PDB-codes in a file named:
#
# pdbs_with_2_DNA_chains.txt
#
#
#
#   -----------------------
#
# Then download all of the corresponding PDB files using this command:


dlpdb.py  < pdbs_with_2_DNA_chains.txt

# dlpdb.py downloads many PDB files in compressed format and uncompresses them.
# Now delete all of the compressed files.

rm -f ????.pdb.gz

#(Ignore the warnings about DSSP files.  Don't worry if you get interrupted.
# If you repeat the command above, it will remember where it left off earlier.)
#
#   -----------------------
#
# Now we start to discard PDB files with undesirable properties.
# (IE PDB files containing RNA or proteins.  NMR structures, etc...)
# 
# For example, we only want to consider iles containing exactly 2 DNA chains.
# First throw away PDB files with proteins or RNA:

mkdir pdbs_with_non-DNA
ls -f1 *.pdb | move_non-dna.sh pdbs_with_non-DNA

#   -----------------------

mkdir pdbs_missing_dna_heavy_atoms
ls -f1 *.pdb | move_missing_dna_heavy_atoms.sh pdbs_missing_dna_heavy_atoms

#   -----------------------

mkdir pdbs_nmr_structures
ls -f1 *.pdb | move_nmr_structures.sh pdbs_nmr_structures

#   -----------------------

# Many of these sequences are redundant
# Now, try to make some effort to select a non-redundant subset of them.
# To do this, I am using the "trim" service provided by
# http://web.expasy.org/decrease_redundancy/
# I use it to generate a subset of sequences with no pair of sequences
# with more than (X) % sequence similarity.
#
# First, we need to create a FASTA file with the DNA sequence
# of ever PDB file contained on a different line.

./convert_dna_pdbs_to_fasta.sh

#   -----------------------

# Copy and paste the contents of the "sequences_some_redundancy.txt" file
# into the web page at:
# http://web.expasy.org/decrease_redundancy/
# Use it to generate a subset of sequences with no pair of sequences
# containing more than 90% sequence similarity.

# Save the resulting FASTA file in the file: "sequences_no_redundancy.txt"

#   -----------------------

# Discard the redundant PDB files.

# First, extract the strings following the ">" character in the FASTA
# file created above to obtain the PDB codes for the PDB files we want:
# Any of the 4-letter PDB codes following the ">" character would work,
# but I arbitrarily chose the first one in that list.
# (It tends to be the highest one, lexicographically, ie. "9xi4" > "113d")

awk '{if ((NF>0) && (substr($1,1,1)==">")) print substr($1,2)}' \
    < sequences_no_redundancy.txt \
    > pdbs_no_redundancy.txt

# Discard the PDB files which are not in this list
# (One way to do that is move the ones we want to keep into a directory...)
mkdir pdbs_no_redundancy
while read PDB_CODE; do
    mv $PDB_CODE.pdb pdbs_no_redundancy
done < pdbs_no_redundancy.txt

# (...and move the remaining PDB files into some other directory...)
mkdir pdbs_some_redundancy
mv *.pdb pdbs_some_redundancy/

# (...and move the PDB files we want back into the main directory)
mv pdbs_no_redundancy/*.pdb .
rmdir pdbs_no_redundancy

#   -----------------------
#
# CHECK FOR THE PRESENCE OF OTHER MOLECULES BINDING TO THE DNA:
# Any molecules which bind to the DNA are likely to modify its structure.
# We don't want to consider these.
# There should be protein or RNA content present in these PDB files,
# because we filtered these out earlier.  However there could be small
# organic molecules.  There is probably way to automatically detect these
# molecules, but there were so few PDB files left by now, I just looked
# through them visually.
#
#   -----------------------
# 
# CHECK FOR EXOTIC CONFORMATIONS OF DNA:
# (IE junctions, A-form DNA, Z-form DNA, non-Watson-Crick pairing, etc...)
# I moved PDB files with these problems into this directory:
# pdbs_exotic_DNA_VISUAL_INSPECTION
#
#   -----------------------
#
# REMOVE OVERHANGS:
# (IE. bits of single-stranded DNA hanging off of either strand.
# Make sure each nucleotide hybridizes with a partner on the opposite strand.
# If not, then the unpaired nucleotides this will shift the numbering
# of the residues in the PDB file.  This will confuse the scripts below,
# and will mess up all of the distance and angle calculations which follow.
# We need to eliminate the overhangs (or discard PDB files which have them).
# I could write a script to check that each nucleotide in the PDB file
# is hydrogen bonded with another nucleotide in the opposite strand.
# But there were so few PDB files left at this point, it was easier to look
# through them visually (and this was a very wise choice for other reaons...)
# I moved these PDB files to pdbs_with_overhangs/
# Optionally, we can extract out the duplex DNA from each of these files

mkdir pdbs_with_overhangs
mkdir pdbs_FINAL_individual_chains
# First, if you checked for overhangs and created any individual-chain
# PDB files, move them into the pdbs_FINAL_individual_chains/ directory.

cd pdbs_FINAL_individual_chains/
  mkdir pdbs_FINAL_USE_THESE_FOR_ANALYSIS
  cd pdbs_FINAL_USE_THESE_FOR_ANALYSIS/

    # Hopefully we should have many pairs of PDB files, in this directory.
    # (Each pair corresponding to different strands from the same duplex DNA.)
    # Now merge these two chains back into a single PDB file
    # but with the residues alternating from either strand
    # (This is necessary.  In order to extract distances and angles
    #  from each base-pair, the atoms in each base-pair must appear
    #  together, not on different residues in different chains.)
    select_interval.py A 3 " " A 10 z < ../../309d.pdb > ../309d_A3-10.pdb
    select_interval.py B 13 " " B 20 z < ../../309d.pdb > ../309d_B13-20.pdb
    dna_interleave_residues.py ../309d_A3-10.pdb ../309d_B13-20.pdb \
              > 309d_interleaved_12.pdb
    # Swap chain order to increase sampling:
    dna_interleave_residues.py ../309d_B13-20.pdb ../309d_A3-10.pdb \
              > 309d_interleaved_21.pdb

    select_interval.py A 3 " " A 10 z < ../../3l1q.pdb > ../3l1q_A3-10.pdb
    select_interval.py B 23 " " B 30 z < ../../3l1q.pdb > ../3l1q_B23-30.pdb
    mv -f 3l1q_*.pdb ../
    dna_interleave_residues.py ../3l1q_A3-10.pdb ../3l1q_B23-30.pdb \
              > 3l1q_interleaved_12.pdb
    # Swap chain order to increase sampling:
    dna_interleave_residues.py ../3l1q_B23-30.pdb ../3l1q_A3-10.pdb \
              > 3l1q_interleaved_21.pdb

  cd ../
  mv pdbs_FINAL_USE_THESE_FOR_ANALYSIS ../
cd ..


#   Now do the same thing with the rest of the PDB files.
#   (Because those filenames match a pattern, we can do this automatically.)

./interleave_nucleotides_in_duplex_dna.sh


#   -----------------------
#
# Optional: Remove the temporary files created by dlpdb.py

rm -f pdbs_old.txt pdbs_most_recent.txt


#   -----------------------
