#    extract_resAveDistances.sh
#
# This script extracts the coordinates for all of the atoms 
# in all of the residues in a list of PDB files.
# -For each residue in a PDB file, it calculates the average position 
#  for the atoms in that residue.  
# -Then it calculates a list of distances between consecutive residue centers
#  and prints these distances to the standard out (all on the same line).
#  There is one line of output per PDB file. (To list these distances on
#  separate lines, regardless of PDB file, see the example below.)
#
#   simple usage:
#
# ls -f1 *.pdb | extract_resAveDistances.sh  > resAveDistances.dat
#
#   If you want all of the numbers to appear on different lines
#   (irrespective of which PDB file they come from), then pipe the result to
#          awk '{ for (i=1;i<=NF;i++) print $i }'    (or   sed -e 's/\s\+/\n/g')
#   as shown below:
#
# ls -f1 *.pdb | extract_resAveDistances.sh \ 
#              | awk '{ for (i=1;i<=NF;i++) print $i }' > resAveDistances.dat 
