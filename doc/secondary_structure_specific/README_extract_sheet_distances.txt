#    extract_sheet_distances.sh
#
# simple usage:
#
# ls -f1 *.pdb | extract_sheet_distances.sh '" CB " i+1 " CB "' \
#              > sheet_distances.dat
#
# Note: Be sure to enclose the atom type arguments with quotes, as above.
#
# This script extracts coordinates of three specified atom types from every
# residue in every sheet of every PDB file.  It then invokes coords2distances.py
# to calculate the distance between these 2 atoms.  This is done 
# repetitively for every residue in every PDB file. 
# (In the example above, we calculate the distance between " CB " atoms in 
#  successive residues.  The i+1 flag tells the script to use the next residue.)
#
# Distances corresponding to the same sheet are printed on
# the same line.  Different sheets on different lines.
# (Short sheets sometimes produce blank lines.)
#
# Later, you might want to post-process the results of this script 
# (sheet_distances.dat) using the "truncate_tokens.py" script:
#  truncate_tokens.py 1 1 < sheet_distances.dat \
#              | awk '{ for (i=1;i<=NF;i++) print $i }' \
#              | awk '{if ($1<0) print 360+$1; else print $1}' \
#              > sheet_distances_trunc1_1.dat
#  This will effectively truncate the beginning and ending of each sheet 
#  both by 1, as well as converting the list of numbers to a single-column 
#  ascii file. (Residues residues near the ends of a sheet are not trustworthy.)
