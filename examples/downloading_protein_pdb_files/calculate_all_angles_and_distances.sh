#!/bin/sh

#   This is the global script reads a list of pdb files in 
# the current directory (*.pdb), and calculates lists of
# all the angles and distances which are relevant to 1-bead
# protein coarse-grained models.
#    These lists are calculated independently for sheets, 
# and for the entire protein chain.  (It's easy to calculate these
# distances and angles for helices and, turns as well.)
#
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
# the same helix (in files named "helix_ ..._raw.dat"), or
# the same sheet (in files named "sheet_ ..._raw.dat"), or
# the same turn  (in files named "turn_ ... _raw.dat"), or
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

ls -f1 *_chain?.pdb | extract_resAveDistances.sh > resAveDistances_raw.dat
ls -f1 *_chain?.pdb | extract_distances.sh '" N  " " CB "' > distances_N-CB_raw.dat
ls -f1 *_chain?.pdb | extract_distances.sh '" CB " " O  "' > distances_CB-O_raw.dat
ls -f1 *_chain?.pdb | extract_distances.sh '" O  " i+1 " N  "' > distances_O-N_raw.dat
ls -f1 *_chain?.pdb | extract_distances.sh '" N  " " H  "' > distances_N-H_raw.dat
ls -f1 *_chain?.pdb | extract_angles.sh '" N  " " CB " " O  "' > angles_N-CB-O_raw.dat
ls -f1 *_chain?.pdb | extract_angles.sh '" CB " " O  " i+1 " N  "' > angles_CB-O-N_raw.dat
ls -f1 *_chain?.pdb | extract_angles.sh '" O  " i+1 " N  " " CB "' > angles_O-N-CB_raw.dat
ls -f1 *_chain?.pdb | extract_angles.sh '" H  " " N  " i-1 " O  "' > angles_H-N-O_raw.dat
ls -f1 *_chain?.pdb | extract_angles.sh '" H  " " N  " " CB "' > angles_H-N-CB_raw.dat
ls -f1 *_chain?.pdb | extract_dihedrals.sh '" N  " " CB " " O  " i+1 " N  "' > dihedrals_N-CB-O-N_raw.dat
ls -f1 *_chain?.pdb | extract_dihedrals.sh '" CB " " O  " i+1 " N  " " CB "' > dihedrals_CB-O-N-CB_raw.dat
ls -f1 *_chain?.pdb | extract_dihedrals.sh '" O  " i+1 " N  " " CB " " O  "' > dihedrals_O-N-CB-O_raw.dat
ls -f1 *_chain?.pdb | extract_dihedrals.sh 180 '" H  " " N  " " CB " " O  "'     > dihedrals_H-N-CB-O_raw.dat
ls -f1 *_chain?.pdb | extract_dihedrals.sh 180 '" H  " " N  " i-1 " O  " " CB "' > dihedrals_H-N-O-CB_raw.dat

# ---------------------------------------------------------------
# ---- Calculate the distances and angles for sheets        -----
# ---------------------------------------------------------------
ls -f1 *_chain?.pdb | extract_sheet_resAveDistances.sh > sheet_resAveDistances_raw.dat
ls -f1 *_chain?.pdb | extract_sheet_distances.sh '" N  " " CB "' > sheet_distances_N-CB_raw.dat
ls -f1 *_chain?.pdb | extract_sheet_distances.sh '" CB " " O  "' > sheet_distances_CB-O_raw.dat
ls -f1 *_chain?.pdb | extract_sheet_distances.sh '" O  " i+1 " N  "' > sheet_distances_O-N_raw.dat
ls -f1 *_chain?.pdb | extract_sheet_distances.sh '" N  " " H  "' > sheet_distances_N-H_raw.dat
ls -f1 *_chain?.pdb | extract_sheet_angles.sh '" N  " " CB " " O  "' > sheet_angles_N-CB-O_raw.dat
ls -f1 *_chain?.pdb | extract_sheet_angles.sh '" CB " " O  " i+1 " N  "' > sheet_angles_CB-O-N_raw.dat
ls -f1 *_chain?.pdb | extract_sheet_angles.sh '" O  " i+1 " N  " " CB "' > sheet_angles_O-N-CB_raw.dat
ls -f1 *_chain?.pdb | extract_sheet_angles.sh '" H  " " N  " i-1 " O  "' > sheet_angles_H-N-O_raw.dat
ls -f1 *_chain?.pdb | extract_sheet_angles.sh '" H  " " N  " " CB "' > sheet_angles_H-N-CB_raw.dat

ls -f1 *_chain?.pdb | extract_sheet_angles.sh '" N  " i+1 " N  " i+2 " N  "' > sheet_angles_N-N-N_raw.dat
ls -f1 *_chain?.pdb | extract_sheet_angles.sh '" H  " " N  " i+1 " N  "' > sheet_angles_H-N-N_raw.dat
ls -f1 *_chain?.pdb | extract_sheet_angles.sh '" N  " i+1 " N  " " H  "' > sheet_angles_N-N-H_raw.dat

ls -f1 *_chain?.pdb | extract_sheet_dihedrals.sh '" N  " " CB " " O  " i+1 " N  "' > sheet_dihedrals_N-CB-O-N_raw.dat
ls -f1 *_chain?.pdb | extract_sheet_dihedrals.sh '" CB " " O  " i+1 " N  " " CB "' > sheet_dihedrals_CB-O-N-CB_raw.dat
ls -f1 *_chain?.pdb | extract_sheet_dihedrals.sh '" O  " i+1 " N  " " CB " " O  "' > sheet_dihedrals_O-N-CB-O_raw.dat
ls -f1 *_chain?.pdb | extract_sheet_dihedrals.sh 180 '" H  " " N  " " CB " " O  "'     > sheet_dihedrals_H-N-CB-O_raw.dat
ls -f1 *_chain?.pdb | extract_sheet_dihedrals.sh 180 '" H  " " N  " i-1 " O  " " CB "' > sheet_dihedrals_H-N-O-CB_raw.dat



#   Now put all the numbers on separate lines, and throw away
#   impossible values (large or negative distances or angles)
awk '{for (i=1;i<=NF;i++){if ($i>=0.0) print $i}}' < resAveDistances_raw.dat > resAveDistances.dat
awk '{for (i=1;i<=NF;i++){if ($i>=0.0) print $i}}' < distances_N-CB_raw.dat > distances_N-CB.dat
awk '{for (i=1;i<=NF;i++){if ($i>=0.0) print $i}}' < distances_CB-O_raw.dat > distances_CB-O.dat
awk '{for (i=1;i<=NF;i++){if ($i>=0.0) print $i}}' < distances_O-N_raw.dat  > distances_O-N.dat
awk '{for (i=1;i<=NF;i++){if ($i>=0.0) print $i}}' < distances_N-H_raw.dat  > distances_N-H.dat
awk '{for (i=1;i<=NF;i++){if ($i>=-360) print $i}}' < angles_N-CB-O_raw.dat > angles_N-CB-O.dat
awk '{for (i=1;i<=NF;i++){if ($i>=-360) print $i}}' < angles_CB-O-N_raw.dat > angles_CB-O-N.dat
awk '{for (i=1;i<=NF;i++){if ($i>=-360) print $i}}' < angles_O-N-CB_raw.dat > angles_O-N-CB.dat
awk '{for (i=1;i<=NF;i++){if ($i>=-360) print $i}}' < angles_H-N-O_raw.dat > angles_H-N-O.dat
awk '{for (i=1;i<=NF;i++){if ($i>=-360) print $i}}' < angles_H-N-CB_raw.dat > angles_H-N-CB.dat
awk '{for (i=1;i<=NF;i++){if ($i>=-360) print $i}}' < dihedrals_N-CB-O-N_raw.dat > dihedrals_N-CB-O-N.dat
awk '{for (i=1;i<=NF;i++){if ($i>=-360) print $i}}' < dihedrals_CB-O-N-CB_raw.dat > dihedrals_CB-O-N-CB.dat
awk '{for (i=1;i<=NF;i++){if ($i>=-360) print $i}}' < dihedrals_O-N-CB-O_raw.dat > dihedrals_O-N-CB-O.dat
awk '{for (i=1;i<=NF;i++){if ($i>=-360) print $i}}' < dihedrals_H-N-CB-O_raw.dat > dihedrals_H-N-CB-O.dat
awk '{for (i=1;i<=NF;i++){if ($i>=-360) print $i}}' < dihedrals_H-N-O-CB_raw.dat > dihedrals_H-N-O-CB.dat

awk '{for (i=1;i<=NF;i++){if ($i>=0.0) print $i}}' < sheet_resAveDistances_raw.dat > sheet_resAveDistances.dat
awk '{for (i=1;i<=NF;i++){if ($i>=0.0) print $i}}' < sheet_distances_N-CB_raw.dat > sheet_distances_N-CB.dat
awk '{for (i=1;i<=NF;i++){if ($i>=0.0) print $i}}' < sheet_distances_CB-O_raw.dat > sheet_distances_CB-O.dat
awk '{for (i=1;i<=NF;i++){if ($i>=0.0) print $i}}' < sheet_distances_O-N_raw.dat  > sheet_distances_O-N.dat
awk '{for (i=1;i<=NF;i++){if ($i>=0.0) print $i}}' < sheet_distances_N-H_raw.dat  > sheet_distances_N-H.dat
awk '{for (i=1;i<=NF;i++){if ($i>=-360) print $i}}' < sheet_angles_N-CB-O_raw.dat > sheet_angles_N-CB-O.dat
awk '{for (i=1;i<=NF;i++){if ($i>=-360) print $i}}' < sheet_angles_CB-O-N_raw.dat > sheet_angles_CB-O-N.dat
awk '{for (i=1;i<=NF;i++){if ($i>=-360) print $i}}' < sheet_angles_O-N-CB_raw.dat > sheet_angles_O-N-CB.dat
awk '{for (i=1;i<=NF;i++){if ($i>=-360) print $i}}' < sheet_angles_H-N-O_raw.dat > sheet_angles_H-N-O.dat
awk '{for (i=1;i<=NF;i++){if ($i>=-360) print $i}}' < sheet_angles_H-N-CB_raw.dat > sheet_angles_H-N-CB.dat

awk '{for (i=1;i<=NF;i++){if ($i>=-360) print $i}}' < sheet_angles_N-N-N_raw.dat > sheet_angles_N-N-N.dat
awk '{for (i=1;i<=NF;i++){if ($i>=-360) print $i}}' < sheet_angles_H-N-N_raw.dat > sheet_angles_H-N-N.dat
awk '{for (i=1;i<=NF;i++){if ($i>=-360) print $i}}' < sheet_angles_N-N-H_raw.dat > sheet_angles_N-N-H.dat

awk '{for (i=1;i<=NF;i++){if ($i>=-360) print $i}}' < sheet_dihedrals_N-CB-O-N_raw.dat > sheet_dihedrals_N-CB-O-N.dat
awk '{for (i=1;i<=NF;i++){if ($i>=-360) print $i}}' < sheet_dihedrals_CB-O-N-CB_raw.dat > sheet_dihedrals_CB-O-N-CB.dat
awk '{for (i=1;i<=NF;i++){if ($i>=-360) print $i}}' < sheet_dihedrals_O-N-CB-O_raw.dat > sheet_dihedrals_O-N-CB-O.dat
awk '{for (i=1;i<=NF;i++){if ($i>=-360) print $i}}' < sheet_dihedrals_H-N-CB-O_raw.dat > sheet_dihedrals_H-N-CB-O.dat
awk '{for (i=1;i<=NF;i++){if ($i>=-360) print $i}}' < sheet_dihedrals_H-N-O-CB_raw.dat > sheet_dihedrals_H-N-O-CB.dat


# Some residues at the very start and the end of a helix or sheet are 
# arranged differently than the residues in the middle, 
# therefore, perhaps we should discard them.
# Create new lists, truncating the first and 
# last few residues from every helix.
#
# This is also a good opportunity to strip out the 
# "impossible" angles,
# (whose values are less than -360 degrees).
# These correspond to angles which could not be calculated
# because at least one of the residues was invalid.
# (For example, if only the beta-carbon (CB) atoms are used 
# to determine the position of the beads, AND one of the 
# residues in this group of 3 or 4 residues lacks this atom.)
truncate_tokens.py 1 1 < sheet_resAveDistances_raw.dat | awk '{ for (i=1;i<=NF;i++) {if ($i>=0.0) print $i}}' > sheet_resAveDistances_trunc1_1.dat
truncate_tokens.py 1 1 < sheet_distances_N-CB_raw.dat | awk '{ for (i=1;i<=NF;i++) {if ($i>=0.0) print $i}}' > sheet_distances_N-CB_trunc1_1.dat
truncate_tokens.py 1 1 < sheet_distances_CB-O_raw.dat | awk '{ for (i=1;i<=NF;i++) {if ($i>=0.0) print $i}}' > sheet_distances_CB-O_trunc1_1.dat
truncate_tokens.py 1 1 < sheet_distances_O-N_raw.dat | awk '{ for (i=1;i<=NF;i++) {if ($i>=0.0) print $i}}' > sheet_distances_O-N_trunc1_1.dat
truncate_tokens.py 1 1 < sheet_distances_N-H_raw.dat | awk '{ for (i=1;i<=NF;i++) {if ($i>=0.0) print $i}}' > sheet_distances_N-H_trunc1_1.dat
truncate_tokens.py 1 1 < sheet_angles_N-CB-O_raw.dat | awk '{ for (i=1;i<=NF;i++){if ($i >=-360) print $i}}' > sheet_angles_N-CB-O_trunc1_1.dat
truncate_tokens.py 1 1 < sheet_angles_CB-O-N_raw.dat | awk '{ for (i=1;i<=NF;i++){if ($i >=-360) print $i}}' > sheet_angles_CB-O-N_trunc1_1.dat
truncate_tokens.py 1 1 < sheet_angles_O-N-CB_raw.dat | awk '{ for (i=1;i<=NF;i++){if ($i >=-360) print $i}}' > sheet_angles_O-N-CB_trunc1_1.dat
truncate_tokens.py 1 1 < sheet_angles_H-N-O_raw.dat | awk '{ for (i=1;i<=NF;i++){if ($i >=-360) print $i}}' > sheet_angles_H-N-O_trunc1_1.dat
truncate_tokens.py 1 1 < sheet_angles_H-N-CB_raw.dat | awk '{ for (i=1;i<=NF;i++){if ($i >=-360) print $i}}' > sheet_angles_H-N-CB_trunc1_1.dat

truncate_tokens.py 1 1 < sheet_angles_N-N-N_raw.dat | awk '{ for (i=1;i<=NF;i++){if ($i >=-360) print $i}}' > sheet_angles_N-N-N_trunc1_1.dat
truncate_tokens.py 1 1 < sheet_angles_H-N-N_raw.dat | awk '{ for (i=1;i<=NF;i++){if ($i >=-360) print $i}}' > sheet_angles_H-N-N_trunc1_1.dat
truncate_tokens.py 1 1 < sheet_angles_N-N-H_raw.dat | awk '{ for (i=1;i<=NF;i++){if ($i >=-360) print $i}}' > sheet_angles_N-N-H_trunc1_1.dat

truncate_tokens.py 1 1 < sheet_dihedrals_N-CB-O-N_raw.dat | awk '{ for (i=1;i<=NF;i++){if ($i >=-360) print $i}}' > sheet_dihedrals_N-CB-O-N_trunc1_1.dat
truncate_tokens.py 1 1 < sheet_dihedrals_CB-O-N-CB_raw.dat | awk '{ for (i=1;i<=NF;i++){if ($i >=-360) print $i}}' > sheet_dihedrals_CB-O-N-CB_trunc1_1.dat
truncate_tokens.py 1 1 < sheet_dihedrals_O-N-CB-O_raw.dat | awk '{ for (i=1;i<=NF;i++){if ($i >=-360) print $i}}' > sheet_dihedrals_O-N-CB-O_trunc1_1.dat
#truncate_tokens.py 1 1 < sheet_dihedrals_H-N-CB-O_raw.dat | awk '{ for (i=1;i<=NF;i++){if ($i >=-360) print $i}}' > sheet_dihedrals_H-N-CB-O_trunc1_1.dat
#truncate_tokens.py 1 1 < sheet_dihedrals_H-N-O-CB_raw.dat | awk '{ for (i=1;i<=NF;i++){if ($i >=-360) print $i}}' > sheet_dihedrals_H-N-O-CB_trunc1_1.dat
truncate_tokens.py 1 1 < sheet_dihedrals_H-N-CB-O_raw.dat | awk '{ for (i=1;i<=NF;i++){if ($i >=-360) { if ($i < 180.0) {print $i} else {print $i-360}} }}' > sheet_dihedrals_H-N-CB-O_trunc1_1.dat
truncate_tokens.py 1 1 < sheet_dihedrals_H-N-O-CB_raw.dat | awk '{ for (i=1;i<=NF;i++){if ($i >=-360) { if ($i < 180.0) {print $i} else {print $i-360}} }}' > sheet_dihedrals_H-N-O-CB_trunc1_1.dat

# ---------------------------------------------------------------
# -- Finally, calculate their averages and standard-deviations --
# ---------------------------------------------------------------

awk '{for(i=1;i<=NF;++i){if ($i>=0.0){sum+=$i;sumsq+=$i*$i;n++}}} END{print sum/n " " sqrt((sumsq/n - (sum/n)*(sum/n))*(n/(n-1)))}' < resAveDistances.dat > resAveDistances_ave_dev.dat
awk '{for(i=1;i<=NF;++i){if ($i>=0.0){sum+=$i;sumsq+=$i*$i;n++}}} END{print sum/n " " sqrt((sumsq/n - (sum/n)*(sum/n))*(n/(n-1)))}' < distances_N-CB.dat > distances_N-CB_ave_dev.dat
awk '{for(i=1;i<=NF;++i){if ($i>=0.0){sum+=$i;sumsq+=$i*$i;n++}}} END{print sum/n " " sqrt((sumsq/n - (sum/n)*(sum/n))*(n/(n-1)))}' < distances_CB-O.dat > distances_CB-O_ave_dev.dat
awk '{for(i=1;i<=NF;++i){if ($i>=0.0){sum+=$i;sumsq+=$i*$i;n++}}} END{print sum/n " " sqrt((sumsq/n - (sum/n)*(sum/n))*(n/(n-1)))}' < distances_O-N.dat > distances_O-N_ave_dev.dat
awk '{for(i=1;i<=NF;++i){if ($i>=0.0){sum+=$i;sumsq+=$i*$i;n++}}} END{print sum/n " " sqrt((sumsq/n - (sum/n)*(sum/n))*(n/(n-1)))}' < distances_N-H.dat > distances_N-H_ave_dev.dat
awk '{if ((-360.0<=$1) && ($1<=360.0)) {n++; sum+=$1; sumsqr+=$1*$1}} END{print sum/n " " sqrt((sumsqr/n - (sum/n)*(sum/n))*(n/(n-1)))}' < angles_N-CB-O.dat > angles_N-CB-O_ave_dev.dat
awk '{if ((-360.0<=$1) && ($1<=360.0)) {n++; sum+=$1; sumsqr+=$1*$1}} END{print sum/n " " sqrt((sumsqr/n - (sum/n)*(sum/n))*(n/(n-1)))}' < angles_CB-O-N.dat > angles_CB-O-N_ave_dev.dat
awk '{if ((-360.0<=$1) && ($1<=360.0)) {n++; sum+=$1; sumsqr+=$1*$1}} END{print sum/n " " sqrt((sumsqr/n - (sum/n)*(sum/n))*(n/(n-1)))}' < angles_O-N-CB.dat > angles_O-N-CB_ave_dev.dat
awk '{if ((-360.0<=$1) && ($1<=360.0)) {n++; sum+=$1; sumsqr+=$1*$1}} END{print sum/n " " sqrt((sumsqr/n - (sum/n)*(sum/n))*(n/(n-1)))}' < angles_H-N-O.dat > angles_H-N-O_ave_dev.dat
awk '{if ((-360.0<=$1) && ($1<=360.0)) {n++; sum+=$1; sumsqr+=$1*$1}} END{print sum/n " " sqrt((sumsqr/n - (sum/n)*(sum/n))*(n/(n-1)))}' < angles_H-N-CB.dat > angles_H-N-CB_ave_dev.dat
awk '{if ((-360.0<=$1) && ($1<=360.0)) {n++; sum+=$1; sumsqr+=$1*$1}} END{print sum/n " " sqrt((sumsqr/n - (sum/n)*(sum/n))*(n/(n-1)))}' < dihedrals_N-CB-O-N.dat > dihedrals_N-CB-O-N_ave_dev.dat
awk '{if ((-360.0<=$1) && ($1<=360.0)) {n++; sum+=$1; sumsqr+=$1*$1}} END{print sum/n " " sqrt((sumsqr/n - (sum/n)*(sum/n))*(n/(n-1)))}' < dihedrals_CB-O-N-CB.dat > dihedrals_CB-O-N-CB_ave_dev.dat
awk '{if ((-360.0<=$1) && ($1<=360.0)) {n++; sum+=$1; sumsqr+=$1*$1}} END{print sum/n " " sqrt((sumsqr/n - (sum/n)*(sum/n))*(n/(n-1)))}' < dihedrals_O-N-CB-O.dat > dihedrals_O-N-CB-O_ave_dev.dat
awk '{if ((-360.0<=$1) && ($1<=360.0)) {n++; sum+=$1; sumsqr+=$1*$1}} END{print sum/n " " sqrt((sumsqr/n - (sum/n)*(sum/n))*(n/(n-1)))}' < dihedrals_H-N-CB-O.dat > dihedrals_H-N-CB-O_ave_dev.dat
awk '{if ((-360.0<=$1) && ($1<=360.0)) {n++; sum+=$1; sumsqr+=$1*$1}} END{print sum/n " " sqrt((sumsqr/n - (sum/n)*(sum/n))*(n/(n-1)))}' < dihedrals_H-N-O-CB.dat > dihedrals_H-N-O-CB_ave_dev.dat


awk '{for(i=1;i<=NF;++i){if ($i>=0.0){sum+=$i;sumsq+=$i*$i;n++}}} END{print sum/n " " sqrt((sumsq/n - (sum/n)*(sum/n))*(n/(n-1)))}' < sheet_resAveDistances_trunc1_1.dat > sheet_resAveDistances_trunc1_1_ave_dev.dat
awk '{for(i=1;i<=NF;++i){if ($i>=0.0){sum+=$i;sumsq+=$i*$i;n++}}} END{print sum/n " " sqrt((sumsq/n - (sum/n)*(sum/n))*(n/(n-1)))}' < sheet_distances_N-CB_trunc1_1.dat > sheet_distances_N-CB_trunc1_1_ave_dev.dat
awk '{for(i=1;i<=NF;++i){if ($i>=0.0){sum+=$i;sumsq+=$i*$i;n++}}} END{print sum/n " " sqrt((sumsq/n - (sum/n)*(sum/n))*(n/(n-1)))}' < sheet_distances_CB-O_trunc1_1.dat > sheet_distances_CB-O_trunc1_1_ave_dev.dat
awk '{for(i=1;i<=NF;++i){if ($i>=0.0){sum+=$i;sumsq+=$i*$i;n++}}} END{print sum/n " " sqrt((sumsq/n - (sum/n)*(sum/n))*(n/(n-1)))}' < sheet_distances_O-N_trunc1_1.dat > sheet_distances_O-N_trunc1_1_ave_dev.dat
awk '{for(i=1;i<=NF;++i){if ($i>=0.0){sum+=$i;sumsq+=$i*$i;n++}}} END{print sum/n " " sqrt((sumsq/n - (sum/n)*(sum/n))*(n/(n-1)))}' < sheet_distances_N-H_trunc1_1.dat > sheet_distances_N-H_trunc1_1_ave_dev.dat
awk '{if ((-360.0<=$1) && ($1<=360.0)) {n++; sum+=$1; sumsqr+=$1*$1}} END{print sum/n " " sqrt((sumsqr/n - (sum/n)*(sum/n))*(n/(n-1)))}' < sheet_angles_N-CB-O_trunc1_1.dat > sheet_angles_N-CB-O_trunc1_1_ave_dev.dat
awk '{if ((-360.0<=$1) && ($1<=360.0)) {n++; sum+=$1; sumsqr+=$1*$1}} END{print sum/n " " sqrt((sumsqr/n - (sum/n)*(sum/n))*(n/(n-1)))}' < sheet_angles_CB-O-N_trunc1_1.dat > sheet_angles_CB-O-N_trunc1_1_ave_dev.dat
awk '{if ((-360.0<=$1) && ($1<=360.0)) {n++; sum+=$1; sumsqr+=$1*$1}} END{print sum/n " " sqrt((sumsqr/n - (sum/n)*(sum/n))*(n/(n-1)))}' < sheet_angles_O-N-CB_trunc1_1.dat > sheet_angles_O-N-CB_trunc1_1_ave_dev.dat
awk '{if ((-360.0<=$1) && ($1<=360.0)) {n++; sum+=$1; sumsqr+=$1*$1}} END{print sum/n " " sqrt((sumsqr/n - (sum/n)*(sum/n))*(n/(n-1)))}' < sheet_angles_H-N-O_trunc1_1.dat > sheet_angles_H-N-O_trunc1_1_ave_dev.dat
awk '{if ((-360.0<=$1) && ($1<=360.0)) {n++; sum+=$1; sumsqr+=$1*$1}} END{print sum/n " " sqrt((sumsqr/n - (sum/n)*(sum/n))*(n/(n-1)))}' < sheet_angles_H-N-CB_trunc1_1.dat > sheet_angles_H-N-CB_trunc1_1_ave_dev.dat

awk '{if ((-360.0<=$1) && ($1<=360.0)) {n++; sum+=$1; sumsqr+=$1*$1}} END{print sum/n " " sqrt((sumsqr/n - (sum/n)*(sum/n))*(n/(n-1)))}' < sheet_angles_N-N-N_trunc1_1.dat > sheet_angles_N-N-N_trunc1_1_ave_dev.dat
awk '{if ((-360.0<=$1) && ($1<=360.0)) {n++; sum+=$1; sumsqr+=$1*$1}} END{print sum/n " " sqrt((sumsqr/n - (sum/n)*(sum/n))*(n/(n-1)))}' < sheet_angles_H-N-N_trunc1_1.dat > sheet_angles_H-N-N_trunc1_1_ave_dev.dat
awk '{if ((-360.0<=$1) && ($1<=360.0)) {n++; sum+=$1; sumsqr+=$1*$1}} END{print sum/n " " sqrt((sumsqr/n - (sum/n)*(sum/n))*(n/(n-1)))}' < sheet_angles_N-N-H_trunc1_1.dat > sheet_angles_N-N-H_trunc1_1_ave_dev.dat

awk '{if ((-360.0<=$1) && ($1<=360.0)) {n++; sum+=$1; sumsqr+=$1*$1}} END{print sum/n " " sqrt((sumsqr/n - (sum/n)*(sum/n))*(n/(n-1)))}' < sheet_dihedrals_N-CB-O-N_trunc1_1.dat > sheet_dihedrals_N-CB-O-N_trunc1_1_ave_dev.dat
awk '{if ((-360.0<=$1) && ($1<=360.0)) {n++; sum+=$1; sumsqr+=$1*$1}} END{print sum/n " " sqrt((sumsqr/n - (sum/n)*(sum/n))*(n/(n-1)))}' < sheet_dihedrals_CB-O-N-CB_trunc1_1.dat > sheet_dihedrals_CB-O-N-CB_trunc1_1_ave_dev.dat
awk '{if ((-360.0<=$1) && ($1<=360.0)) {n++; sum+=$1; sumsqr+=$1*$1}} END{print sum/n " " sqrt((sumsqr/n - (sum/n)*(sum/n))*(n/(n-1)))}' < sheet_dihedrals_O-N-CB-O_trunc1_1.dat > sheet_dihedrals_O-N-CB-O_trunc1_1_ave_dev.dat
awk '{if ((-360.0<=$1) && ($1<=360.0)) {n++; sum+=$1; sumsqr+=$1*$1}} END{print sum/n " " sqrt((sumsqr/n - (sum/n)*(sum/n))*(n/(n-1)))}' < sheet_dihedrals_H-N-CB-O_trunc1_1.dat > sheet_dihedrals_H-N-CB-O_trunc1_1_ave_dev.dat
awk '{if ((-360.0<=$1) && ($1<=360.0)) {n++; sum+=$1; sumsqr+=$1*$1}} END{print sum/n " " sqrt((sumsqr/n - (sum/n)*(sum/n))*(n/(n-1)))}' < sheet_dihedrals_H-N-O-CB_trunc1_1.dat > sheet_dihedrals_H-N-O-CB_trunc1_1_ave_dev.dat

