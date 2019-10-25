#!/bin/sh

#   This is the global script reads a list of pdb files in 
# the current directory (*.pdb), and calculates lists of
# all the angles and distances which are relevant to 1-bead
# protein coarse-grained models.
#    These lists are calculated independently for sheets, 
# and for the entire protein chain.  (It's easy to calculate these
# distances and angles for helices and, turns as well.)

# (Histograms of these number lists show the distribution of
#  distances and angles in this set of pdb files.)
# The average and standard-deviation of these lists 
# of angles and distances is calculated and saved in
# separate files.
# 
#   These lists are generated separately for the whole chain,
# and for individual sheets in each PDB file.
#
#   --- output files ---
#
# "_raw.dat" files:
#
#    Glycine residues are omitted when calculating angles and distances
# between average side-chain atoms, and/or CB atoms.
# Angles and distances involving the positions of glycine residues are either
# deliberately set to impossible values (-1.0 for distances, -720 for angles),
# as is the case for "_raw.dat" files...
# ...or they are simply omitted (in all the other .dat files), 
# The data in "_raw.dat" files preserves the original order and 
# location of every residue in each PDB file (that is, after 
# being sorted by PDB-residue-ID number and PDB-insert-code).
#    These "_raw.dat" files also contain multiple numbers on each line,
# corresponding to the distances and angles that come from the
# the same chain (all the other "..._raw.dat" files)
#
# regular ".dat" files
#
#    Other .dat files (ie. not ending in "_raw.dat") have data corresponding to
# separate residues on separate lines, usually resulting in 1 number per line.
# The data in these files omits the angles and distances from distarded
# residues (for example glycine residues, or residues near the beginning
# and end of a helix or sheet, see below).


# ---------------------------------------------------------------
# ---- Calculate the distances and angles for entire chains -----
# ---------------------------------------------------------------

#    The DNA helix pdb files were altered so that each pair
# of nucleotides (from alternating strands) appear consecutively
# in the PDB file.
#    (To do this, I selected residues from the PDB files
# which belonged to double-stranded DNA, and discarded any 
# single-stranded overhangs.  Then I extracted the two DNA strands
# into separate PDB files (running in opposite directions), 
# and used the "pdb_interleave_residues.py" to merge them into
# a single PDB file.  These merged .pdb files were used below:

# The distance between C3' pairs in a base-pair can be calculated by
# looking at the distance between successive C3' atoms in the PDB file
ls -f1 *_12.pdb | extract_distances.sh "\"[0::2]\" \" C3'\" i+1 \" C3'\"" > distances_basepairs_C3p-C3p_raw.dat
# Note: This also works
# ls -f1 *_12.pdb | extract_distances.sh '[0::2] " C3'"'"'" i+1 " C3'"'"'"' > distances_basepairs_raw.dat

# To find the distance between C3' atoms on the same strand,
# you must alternately skip over the C3' atoms on the opposite strand.
# (So I use i+2 instead of i+1)
ls -f1 *_12.pdb | extract_distances.sh "\" C3'\" i+6 \" C3'\"" > distances_backbone_C3p-C3p_raw.dat


# Optional new addition (2016-3-18): "ZIGZAG" interactions between atoms
# which are not directly bonded together in real DNA. These interactions 
# are between distant atoms which can strongly effect bending rigidity.
# Models that use these distance constraints might be able to avoid using
# dihedral angle forces (and possibly even angle forces).
ls -f1 *_12.pdb | extract_distances.sh "\"[0::2]\" \" C3'\" i+13 \" C3'\"" > distances_zigzag_majorgroove_C3p-C3p_raw.dat
ls -f1 *_12.pdb | extract_distances.sh "\"[1::2]\" \" C3'\" i+11 \" C3'\"" > distances_zigzag_minorgroove_C3p-C3p_raw.dat



# Backbone Angles:
ls -f1 *_12.pdb | extract_angles.sh "\" C3'\" i+6 \" C3'\" i+12 \" C3'\"" > angles_backbone_C3p-C3p-C3p_raw.dat

# Each atom participates in 2 base-pair angles, but the relative numbering
# of the other two atoms participating in those interactions
# also depends on whether or not the atom is on an even-or-odd nucleotide.
# "[0::2]" are for odd nucleotides    "[1::2]" are for odd nucleotides
# (This is python slice notation for printing the odd and even entries of a list.)
# This means we have to have 4 rules (2 for even, 2 for odd).
ls -f1 *_12.pdb | extract_angles.sh "\"[0::2]\" \" C3'\" i+1 \" C3'\" i+7 \" C3'\"" > angles_basepairs_obtuse_C3p-C3p-C3p_raw.dat
ls -f1 *_12.pdb | extract_angles.sh "\"[0::2]\" \" C3'\" i+6 \" C3'\" i+7 \" C3'\"" >> angles_basepairs_obtuse_C3p-C3p-C3p_raw.dat
ls -f1 *_12.pdb | extract_angles.sh "\"[1::2]\" \" C3'\" i-1 \" C3'\" i+5 \" C3'\"" > angles_basepairs_acute_C3p-C3p-C3p_raw.dat
ls -f1 *_12.pdb | extract_angles.sh "\"[1::2]\" \" C3'\" i+6 \" C3'\" i+5 \" C3'\"" >> angles_basepairs_acute_C3p-C3p-C3p_raw.dat


# Backbone Dihedrals:
# COMMENTING OUT (backbone dihedrals become numerically unstable as 
# bond-angles approach 180, which they do in the "every_base_pair" model)
#ls -f1 *_12.pdb | extract_dihedrals.sh 180 "\" C3'\" i+6 \" C3'\" i+12 \" C3'\" i+18 \" C3'\"" > dihedrals_backbone_C3p-C3p-C3p-C3p_raw.dat


# Each atom participates in 1 base-pair dihedral, but the relative numbering
# of the other three atoms participating in those interactions
# depends on whether or not the atom is on an even-or-odd nucleotide.
# "[0::2]" are for odd nucleotides    "[1::2]" are for odd nucleotides
# (This is python notation for printing the odd and even entries of a list.)
# This means we have to have 2 rules:
ls -f1 *_12.pdb | extract_dihedrals.sh 180 "\"[0::2]\" \" C3'\" i+1 \" C3'\" i+7 \" C3'\" i+6 \" C3'\"" > dihedrals_basepairs_C3p-C3p-C3p-C3p_raw.dat
ls -f1 *_12.pdb | extract_dihedrals.sh 180 "\"[1::2]\" \" C3'\" i-1 \" C3'\" i+5 \" C3'\" i+6 \" C3'\"" >> dihedrals_basepairs_C3p-C3p-C3p-C3p_raw.dat


# Optional new addition (2016-3-18): "ZIGZAG" interactions between atoms
# which are not directly bonded together in real DNA. These interactions 
# are between distant atoms which can strongly effect bending rigidity.
ls -f1 *_12.pdb | extract_dihedrals.sh 180 "\"[0::2]\" \" C3'\" i+6 \" C3'\" i+7 \" C3'\" i+13 \" C3'\"" > dihedrals_zigzag_majorgroove_C3p-C3p-C3p-C3p_raw.dat
ls -f1 *_12.pdb | extract_dihedrals.sh 180 "\"[1::2]\" \" C3'\" i+6 \" C3'\" i+5 \" C3'\" i+11 \" C3'\"" > dihedrals_zigzag_minorgroove_C3p-C3p-C3p-C3p_raw.dat

ls -f1 *_12.pdb | extract_dihedrals.sh 180 "\"[0::2]\" \" C3'\" i-6 \" C3'\" i-5 \" C3'\" i+1 \" C3'\"" > dihedrals_zigzag_torsion_C3p-C3p-C3p-C3p_raw.dat
ls -f1 *_12.pdb | extract_dihedrals.sh 180 "\"[0::2]\" \" C3'\" i+6 \" C3'\" i+7 \" C3'\" i+1 \" C3'\"" >> dihedrals_zigzag_torsion_C3p-C3p-C3p-C3p_raw.dat


#   Now put all the numbers on separate lines, and throw away
#   impossible values (large or negative distances or angles)
awk '{for (i=1;i<=NF;i++){if ($i>=0.0) print $i}}' < distances_basepairs_C3p-C3p_raw.dat > distances_basepairs_C3p-C3p.dat
awk '{for (i=1;i<=NF;i++){if ($i>=0.0) print $i}}' < distances_backbone_C3p-C3p_raw.dat > distances_backbone_C3p-C3p.dat
awk '{for (i=1;i<=NF;i++){if ($i>=0.0) print $i}}' < distances_zigzag_majorgroove_C3p-C3p_raw.dat > distances_zigzag_majorgroove_C3p-C3p.dat
awk '{for (i=1;i<=NF;i++){if ($i>=0.0) print $i}}' < distances_zigzag_minorgroove_C3p-C3p_raw.dat > distances_zigzag_minorgroove_C3p-C3p.dat

awk '{for (i=1;i<=NF;i++){if ($i>=-180) print $i}}' < angles_backbone_C3p-C3p-C3p_raw.dat > angles_backbone_C3p-C3p-C3p.dat
awk '{for (i=1;i<=NF;i++){if ($i>=-180) print $i}}' < angles_basepairs_obtuse_C3p-C3p-C3p_raw.dat > angles_basepairs_obtuse_C3p-C3p-C3p.dat
awk '{for (i=1;i<=NF;i++){if ($i>=-180) print $i}}' < angles_basepairs_acute_C3p-C3p-C3p_raw.dat > angles_basepairs_acute_C3p-C3p-C3p.dat

# COMMENTING OUT (backbone dihedrals become numerically unstable as 
# bond-angles approach 180, which they do in the "every_base_pair" model)
#awk '{for (i=1;i<=NF;i++){if ($i>=-360) print $i}}' < dihedrals_backbone_C3p-C3p-C3p-C3p_raw.dat > dihedrals_backbone_C3p-C3p-C3p-C3p.dat
awk '{for (i=1;i<=NF;i++){if ($i>=-360) print $i}}' < dihedrals_basepairs_C3p-C3p-C3p-C3p_raw.dat > dihedrals_basepairs_C3p-C3p-C3p-C3p.dat
awk '{for (i=1;i<=NF;i++){if ($i>=-360) print $i}}' < dihedrals_zigzag_majorgroove_C3p-C3p-C3p-C3p_raw.dat > dihedrals_zigzag_majorgroove_C3p-C3p-C3p-C3p.dat
awk '{for (i=1;i<=NF;i++){if ($i>=-360) print $i}}' < dihedrals_zigzag_minorgroove_C3p-C3p-C3p-C3p_raw.dat > dihedrals_zigzag_minorgroove_C3p-C3p-C3p-C3p.dat
awk '{for (i=1;i<=NF;i++){if ($i>=-360) print $i}}' < dihedrals_zigzag_torsion_C3p-C3p-C3p-C3p_raw.dat > dihedrals_zigzag_torsion_C3p-C3p-C3p-C3p.dat


# Some nucleotides at the very start and the end of a helix are
# arranged differently than the nucleotides in the middle, 
# therefore, perhaps we should discard them.
# Create new lists, truncating the first and 
# last few nucleotides from every helix.
#
# This is also a good opportunity to strip out the 
# "impossible" angles,
# (whose values are less than -180 or -360 degrees).
# These correspond to angles which could not be calculated
# because at least one of the residues was invalid.
# (For example, if only the beta-carbon (CB) atoms are used 
# to determine the position of the beads, AND one of the 
# residues in this group of 3 or 4 residues lacks this atom.)



awk '{for (i=1;i<=NF;i++) {if ($i>=0.0) printf "%g ",$i} printf "\n"}' < distances_basepairs_C3p-C3p_raw.dat | truncate_tokens.py 1 1 | awk '{ for (i=1;i<=NF;i++){print $i}}' > distances_basepairs_C3p-C3p_trunc1-1.dat

awk '{for (i=1;i<=NF;i++) {if ($i>=0.0) printf "%g ",$i} printf "\n"}' < distances_backbone_C3p-C3p_raw.dat | truncate_tokens.py 1 1 | awk '{ for (i=1;i<=NF;i++){print $i}}' > distances_backbone_C3p-C3p_trunc1-1.dat

awk '{ for (i=1;i<=NF;i++) {if ($i>=0.0) printf "%g ",$i} printf "\n"}' < distances_zigzag_majorgroove_C3p-C3p_raw.dat | truncate_tokens.py 1 1 | awk '{ for (i=1;i<=NF;i++){print $i}}' > distances_zigzag_majorgroove_C3p-C3p_trunc1-1.dat

awk '{for (i=1;i<=NF;i++) {if ($i>=0.0) printf "%g ",$i} printf "\n"}' < distances_zigzag_minorgroove_C3p-C3p_raw.dat | truncate_tokens.py 1 1 | awk '{ for (i=1;i<=NF;i++){print $i}}' > distances_zigzag_minorgroove_C3p-C3p_trunc1-1.dat


awk '{for (i=1;i<=NF;i++){if ($i >=-180) printf "%g ",$i} printf "\n"}' < angles_backbone_C3p-C3p-C3p_raw.dat | truncate_tokens.py 1 1 | awk '{ for (i=1;i<=NF;i++){print $i}}' > angles_backbone_C3p-C3p-C3p_trunc1-1.dat

awk '{for (i=1;i<=NF;i++){if ($i >=-180) printf "%g ",$i} printf "\n"}' < angles_basepairs_obtuse_C3p-C3p-C3p_raw.dat | truncate_tokens.py 1 1 | awk '{ for (i=1;i<=NF;i++){print $i}}' > angles_basepairs_obtuse_C3p-C3p-C3p_trunc1-1.dat

awk '{for (i=1;i<=NF;i++){if ($i >=-180) printf "%g ",$i} printf "\n"}' < angles_basepairs_acute_C3p-C3p-C3p_raw.dat | truncate_tokens.py 1 1 | awk '{ for (i=1;i<=NF;i++){print $i}}' > angles_basepairs_acute_C3p-C3p-C3p_trunc1-1.dat



# COMMENTING OUT (backbone dihedrals become numerically unstable as 
# bond-angles approach 180, which they do in the "every_base_pair" model)
#awk '{for (i=1;i<=NF;i++){if ($i >=-360) printf "%g ",$i} printf "\n"}' < dihedrals_backbone_C3p-C3p-C3p-C3p_raw.dat | truncate_tokens.py 1 1 | awk '{ for (i=1;i<=NF;i++){print $i}}' > dihedrals_backbone_C3p-C3p-C3p-C3p_trunc1-1.dat

awk '{for (i=1;i<=NF;i++){if ($i >=-360) printf "%g ",$i} printf "\n"}' < dihedrals_basepairs_C3p-C3p-C3p-C3p_raw.dat | truncate_tokens.py 1 1 | awk '{ for (i=1;i<=NF;i++){print $i}}' > dihedrals_basepairs_C3p-C3p-C3p-C3p_trunc1-1.dat

awk '{for (i=1;i<=NF;i++){if ($i >=-360) printf "%g ",$i} printf "\n"}' < dihedrals_zigzag_majorgroove_C3p-C3p-C3p-C3p_raw.dat | truncate_tokens.py 1 1 | awk '{ for (i=1;i<=NF;i++){print $i}}' > dihedrals_zigzag_majorgroove_C3p-C3p-C3p-C3p_trunc1-1.dat

awk '{for (i=1;i<=NF;i++){if ($i >=-360) printf "%g ",$i} printf "\n"}' < dihedrals_zigzag_minorgroove_C3p-C3p-C3p-C3p_raw.dat | truncate_tokens.py 1 1 | awk '{ for (i=1;i<=NF;i++){print $i}}' > dihedrals_zigzag_minorgroove_C3p-C3p-C3p-C3p_trunc1-1.dat

awk '{for (i=1;i<=NF;i++){if ($i >=-360) printf "%g ",$i} printf "\n"}' < dihedrals_zigzag_torsion_C3p-C3p-C3p-C3p_raw.dat | truncate_tokens.py 1 1 | awk '{ for (i=1;i<=NF;i++){print $i}}' > dihedrals_zigzag_torsion_C3p-C3p-C3p-C3p_trunc1-1.dat




# ---------------------------------------------------------------
# -- Finally, calculate their averages and standard-deviations --
# ---------------------------------------------------------------

awk '{for(i=1;i<=NF;++i){if ($i>=0.0){sum+=$i;sumsq+=$i*$i;n++}}} END{print sum/n " " sqrt((sumsq/n - (sum/n)*(sum/n))*(n/(n-1))) " " n}' \
 < distances_basepairs_C3p-C3p.dat \
 > distances_basepairs_C3p-C3p_ave_dev_n.dat
awk '{for(i=1;i<=NF;++i){if ($i>=0.0){sum+=$i;sumsq+=$i*$i;n++}}} END{print sum/n " " sqrt((sumsq/n - (sum/n)*(sum/n))*(n/(n-1))) " " n}' \
 < distances_basepairs_C3p-C3p_trunc1-1.dat \
 > distances_basepairs_C3p-C3p_trunc1-1_ave_dev_n.dat


awk '{for(i=1;i<=NF;++i){if ($i>=0.0){sum+=$i;sumsq+=$i*$i;n++}}} END{print sum/n " " sqrt((sumsq/n - (sum/n)*(sum/n))*(n/(n-1))) " " n}' \
 < distances_backbone_C3p-C3p.dat \
 > distances_backbone_C3p-C3p_ave_dev_n.dat
awk '{for(i=1;i<=NF;++i){if ($i>=0.0){sum+=$i;sumsq+=$i*$i;n++}}} END{print sum/n " " sqrt((sumsq/n - (sum/n)*(sum/n))*(n/(n-1))) " " n}' \
 < distances_backbone_C3p-C3p_trunc1-1.dat \
 > distances_backbone_C3p-C3p_trunc1-1_ave_dev_n.dat


awk '{for(i=1;i<=NF;++i){if ($i>=0.0){sum+=$i;sumsq+=$i*$i;n++}}} END{print sum/n " " sqrt((sumsq/n - (sum/n)*(sum/n))*(n/(n-1))) " " n}' \
 < distances_zigzag_majorgroove_C3p-C3p.dat \
 > distances_zigzag_majorgroove_C3p-C3p_ave_dev_n.dat
awk '{for(i=1;i<=NF;++i){if ($i>=0.0){sum+=$i;sumsq+=$i*$i;n++}}} END{print sum/n " " sqrt((sumsq/n - (sum/n)*(sum/n))*(n/(n-1))) " " n}' \
 < distances_zigzag_majorgroove_C3p-C3p_trunc1-1.dat \
 > distances_zigzag_majorgroove_C3p-C3p_trunc1-1_ave_dev_n.dat


awk '{for(i=1;i<=NF;++i){if ($i>=0.0){sum+=$i;sumsq+=$i*$i;n++}}} END{print sum/n " " sqrt((sumsq/n - (sum/n)*(sum/n))*(n/(n-1))) " " n}' \
 < distances_zigzag_minorgroove_C3p-C3p.dat \
 > distances_zigzag_minorgroove_C3p-C3p_ave_dev_n.dat
awk '{for(i=1;i<=NF;++i){if ($i>=0.0){sum+=$i;sumsq+=$i*$i;n++}}} END{print sum/n " " sqrt((sumsq/n - (sum/n)*(sum/n))*(n/(n-1))) " " n}' \
 < distances_zigzag_minorgroove_C3p-C3p_trunc1-1.dat \
 > distances_zigzag_minorgroove_C3p-C3p_trunc1-1_ave_dev_n.dat

awk '{for(i=1;i<=NF;++i){if ($i<=0.0){sum+=$i;sumsq+=$i*$i;n++}}} END{print sum/n " " sqrt((sumsq/n - (sum/n)*(sum/n))*(n/(n-1))) " " n}' \
 < dihedrals_zigzag_torsion_C3p-C3p-C3p-C3p.dat \
 > dihedrals_zigzag_torsion_C3p-C3p-C3p-C3p_ave_dev.dat
awk '{for(i=1;i<=NF;++i){if ($i<=0.0){sum+=$i;sumsq+=$i*$i;n++}}} END{print sum/n " " sqrt((sumsq/n - (sum/n)*(sum/n))*(n/(n-1))) " " n}' \
 < dihedrals_zigzag_torsion_C3p-C3p-C3p-C3p_trunc1-1.dat \
 > dihedrals_zigzag_torsion_C3p-C3p-C3p-C3p_trunc1-1_ave_dev.dat



awk '{for(i=1;i<=NF;++i){if ($i>=0.0){sum+=$i;sumsq+=$i*$i;n++}}} END{print sum/n " " sqrt((sumsq/n - (sum/n)*(sum/n))*(n/(n-1))) " " n}' \
 < angles_backbone_C3p-C3p-C3p.dat \
 > angles_backbone_C3p-C3p-C3p_ave_dev_n.dat
awk '{for(i=1;i<=NF;++i){if ($i>=0.0){sum+=$i;sumsq+=$i*$i;n++}}} END{print sum/n " " sqrt((sumsq/n - (sum/n)*(sum/n))*(n/(n-1))) " " n}' \
 < angles_backbone_C3p-C3p-C3p_trunc1-1.dat \
 > angles_backbone_C3p-C3p-C3p_trunc1-1_ave_dev_n.dat


awk '{for(i=1;i<=NF;++i){if ($i>=0.0){sum+=$i;sumsq+=$i*$i;n++}}} END{print sum/n " " sqrt((sumsq/n - (sum/n)*(sum/n))*(n/(n-1))) " " n}' \
 < angles_basepairs_obtuse_C3p-C3p-C3p.dat \
 > angles_basepairs_obtuse_C3p-C3p-C3p_ave_dev_n.dat
awk '{for(i=1;i<=NF;++i){if ($i>=0.0){sum+=$i;sumsq+=$i*$i;n++}}} END{print sum/n " " sqrt((sumsq/n - (sum/n)*(sum/n))*(n/(n-1))) " " n}' \
 < angles_basepairs_obtuse_C3p-C3p-C3p_trunc1-1.dat \
 > angles_basepairs_obtuse_C3p-C3p-C3p_trunc1-1_ave_dev_n.dat


awk '{for(i=1;i<=NF;++i){if ($i>=0.0){sum+=$i;sumsq+=$i*$i;n++}}} END{print sum/n " " sqrt((sumsq/n - (sum/n)*(sum/n))*(n/(n-1))) " " n}' \
 < angles_basepairs_acute_C3p-C3p-C3p.dat \
 > angles_basepairs_acute_C3p-C3p-C3p_ave_dev_n.dat
awk '{for(i=1;i<=NF;++i){if ($i>=0.0){sum+=$i;sumsq+=$i*$i;n++}}} END{print sum/n " " sqrt((sumsq/n - (sum/n)*(sum/n))*(n/(n-1))) " " n}' \
 < angles_basepairs_acute_C3p-C3p-C3p_trunc1-1.dat \
 > angles_basepairs_acute_C3p-C3p-C3p_trunc1-1_ave_dev_n.dat





# COMMENTING OUT (backbone dihedrals become numerically unstable as 
# bond-angles approach 180, which they do in the "every_base_pair" model)
#awk '{for(i=1;i<=NF;++i){if ($i>=0.0){sum+=$i;sumsq+=$i*$i;n++}}} END{print sum/n " " sqrt((sumsq/n - (sum/n)*(sum/n))*(n/(n-1))) " " n}' \
# < dihedrals_backbone_C3p-C3p-C3p-C3p.dat \
# > dihedrals_backbone_C3p-C3p-C3p-C3p_ave_dev_n.dat
#awk '{for(i=1;i<=NF;++i){if ($i>=0.0){sum+=$i;sumsq+=$i*$i;n++}}} END{print sum/n " " sqrt((sumsq/n - (sum/n)*(sum/n))*(n/(n-1))) " " n}'
# < dihedrals_backbone_C3p-C3p-C3p-C3p_trunc1-1.dat \
# > dihedrals_backbone_C3p-C3p-C3p-C3p_trunc1-1_ave_dev_n.dat


awk '{for(i=1;i<=NF;++i){if ($i>=0.0){sum+=$i;sumsq+=$i*$i;n++}}} END{print sum/n " " sqrt((sumsq/n - (sum/n)*(sum/n))*(n/(n-1))) " " n}' \
 < dihedrals_basepairs_C3p-C3p-C3p-C3p.dat \
 > dihedrals_basepairs_C3p-C3p-C3p-C3p_ave_dev_n.dat
awk '{for(i=1;i<=NF;++i){if ($i>=0.0){sum+=$i;sumsq+=$i*$i;n++}}} END{print sum/n " " sqrt((sumsq/n - (sum/n)*(sum/n))*(n/(n-1))) " " n}' \
 < dihedrals_basepairs_C3p-C3p-C3p-C3p_trunc1-1.dat \
 > dihedrals_basepairs_C3p-C3p-C3p-C3p_trunc1-1_ave_dev_n.dat


awk '{for(i=1;i<=NF;++i){if ($i>=0.0){sum+=$i;sumsq+=$i*$i;n++}}} END{print sum/n " " sqrt((sumsq/n - (sum/n)*(sum/n))*(n/(n-1))) " " n}' \
 < dihedrals_zigzag_majorgroove_C3p-C3p-C3p-C3p.dat \
 > dihedrals_zigzag_majorgroove_C3p-C3p-C3p-C3p_ave_dev_n.dat
awk '{for(i=1;i<=NF;++i){if ($i>=0.0){sum+=$i;sumsq+=$i*$i;n++}}} END{print sum/n " " sqrt((sumsq/n - (sum/n)*(sum/n))*(n/(n-1))) " " n}' \
 < dihedrals_zigzag_majorgroove_C3p-C3p-C3p-C3p_trunc1-1.dat \
 > dihedrals_zigzag_majorgroove_C3p-C3p-C3p-C3p_trunc1-1_ave_dev_n.dat


awk '{for(i=1;i<=NF;++i){if ($i>=0.0){sum+=$i;sumsq+=$i*$i;n++}}} END{print sum/n " " sqrt((sumsq/n - (sum/n)*(sum/n))*(n/(n-1))) " " n}' \
 < dihedrals_zigzag_minorgroove_C3p-C3p-C3p-C3p.dat \
 > dihedrals_zigzag_minorgroove_C3p-C3p-C3p-C3p_ave_dev_n.dat
awk '{for(i=1;i<=NF;++i){if ($i>=0.0){sum+=$i;sumsq+=$i*$i;n++}}} END{print sum/n " " sqrt((sumsq/n - (sum/n)*(sum/n))*(n/(n-1))) " " n}' \
 < dihedrals_zigzag_minorgroove_C3p-C3p-C3p-C3p_trunc1-1.dat \
 > dihedrals_zigzag_minorgroove_C3p-C3p-C3p-C3p_trunc1-1_ave_dev_n.dat


awk '{for(i=1;i<=NF;++i){if ($i<=0.0){sum+=$i;sumsq+=$i*$i;n++}}} END{print sum/n " " sqrt((sumsq/n - (sum/n)*(sum/n))*(n/(n-1))) " " n}' \
 < dihedrals_zigzag_torsion_C3p-C3p-C3p-C3p.dat \
 > dihedrals_zigzag_torsion_C3p-C3p-C3p-C3p_ave_dev_n.dat
awk '{for(i=1;i<=NF;++i){if ($i<=0.0){sum+=$i;sumsq+=$i*$i;n++}}} END{print sum/n " " sqrt((sumsq/n - (sum/n)*(sum/n))*(n/(n-1))) " " n}' \
 < dihedrals_zigzag_torsion_C3p-C3p-C3p-C3p_trunc1-1.dat \
 > dihedrals_zigzag_torsion_C3p-C3p-C3p-C3p_trunc1-1_ave_dev_n.dat

