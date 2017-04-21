#    extract_angles.sh
#
# This script extracts coordinates of three specified atom types from every
# residue in a list of PDB files.  It then invokes coords2angles.py to
# calculate the 3-body angle between these 3 atoms.  This is done 
# repetitively for every residue in every PDB file. (Angles between atoms in 
# the same PDB file are printed on the same line.  There is one line of output
# per PDB file.  To list these angles on separate lines, see the example below.)
# When a needed atom type is missing, an impossible (negative) value is printed.
# Angles between atoms in different (successive) residues can be calculated
# using the "i+1" (or "i+2", "i+3", ...) flags as shown below:
#
# simple usage:
# ls -f1 *.pdb | extract_angles.sh '" CB " " O  " i+1 " N  "' > angles.dat
#
# Note: Be sure to enclose atom type arguments with quotes, exactly as above.
#
#    This will calculate the angle between these atoms:
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
# ls -f1 *.pdb | extract_angles.sh '" CB " " O  " i+1 " N  "' \
#              | awk '{ for (i=1;i<=NF;i++) print $i }' > angles.dat
