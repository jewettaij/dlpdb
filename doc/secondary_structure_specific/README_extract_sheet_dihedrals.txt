#    extract_sheet_dihedrals.txt
#
# simple usage:
#
# ls -f1 *.pdb|extract_sheet_dihedrals.sh 180 '" N  " " CB " " O  " i+1 " N  "'\
#             > sheet_dihedrals.dat
#
# Note: Be sure to enclose atom type arguments with quotes, exactly as above.
# 
# Note: The "180" is an optional argument. It causes dihedral angles 
#       to range from -180 to 180. (Any angle can be used. Default: 360)
# 
# Explanation:
# This script extracts coordinates of four specified atom types from every
# residue in every sheet of every PDB file.  It then invokes coords2dihedrals.py
# to calculate the 4-body dihedral angle between these 4 atoms.  This is done 
# repetitively for every residue in every PDB file. 
# (In the example above, we calculate the dihedral angle between the " N  ", 
#  " CB ", and " O  " atoms, and the " CB " atom in the next residue.)
#
# Dihedral angles corresponding to the same sheet are printed on
# the same line.  Different sheets on different lines.
# (Short sheets sometimes produce blank lines.)
#
# Later, you might want to post-process the results of this script 
# (sheet_dihedrals.dat) using the "truncate_tokens.py" script:
#  truncate_tokens.py 1 1 < sheet_dihedrals.dat \
#              | awk '{ for (i=1;i<=NF;i++) print $i }' \
#              | awk '{if ($1<0) print 360+$1; else print $1}' \
#              > sheet_dihedrals_trunc1_1.dat
#  This will effectively truncate the beginning and ending of each sheet 
#  both by 1, as well as converting the list of numbers to a single-column 
#  ascii file. (Residues residues near the ends of a sheet are not trustworthy.)
