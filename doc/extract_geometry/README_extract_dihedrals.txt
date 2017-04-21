#    extract_dihedrals.sh
#
# This script extracts coordinates of four specified atom types from every
# residue in a list of PDB files.  It then invokes coords2dihedrals.py to
# calculate the 4-body dihedral angle between these 4 atoms.  This is done 
# repetitively for every residue in every PDB file. (Angles between atoms in 
# the same PDB file are printed on the same line.  There is one line of output
# per PDB file.  To list these angles on separate lines, see the example below.)
# When a needed atom type is missing, an impossible value (-720) is printed.
# Angles between atoms in different (successive) residues can be calculated
# using the "i+1" (or "i+2", "i+3", ...) flags as shown below:
# 
# simple usage:
# ls -f1 *.pdb | extract_dihedrals.sh 180 '" N  " " CB " " O  " i+1 " N  "' \
#              > dihedrals.dat
#
# Note: Be sure to enclose atom type arguments with quotes, exactly as above.
# 
# Note: The "180" is an optional argument. It causes dihedral angles 
#       to range from -180 to 180. (Any angle can be used. Default: 360)
# 
# Explanation: This will calculate the dihedral angle between these atoms:
#            the " N  " atom
#            the " CB " atom, 
#            the " O  " atom, and
#            the " N  " atom (from the next residue).
#    for every residue in all of the pdb files in this directory.
# 
#   If you want all of the numbers to appear on different lines
#   (irrespective of which PDB file they come from), then pipe the result to
#          awk '{ for (i=1;i<=NF;i++) print $i }'    (or   sed -e 's/\s\+/\n/g')
#   as shown below:
# 
# ls -f1 *.pdb | extract_dihedrals.sh 180 '" N  " " CB " " O  " i+1 " N  "' \
#              | awk '{ for (i=1;i<=NF;i++) print $i }' > dihedrals.dat
