#    extract_sheet_resAveDistances.sh
#
# simple usage:
# ls -f1 *.pdb | extract_sheet_resAveDistances.sh  > sheet_resAveDistances.dat
#
# Distances corresponding to the same sheet are printed on
# the same line.  Different sheets on different lines.
# (Short sheets sometimes produce blank lines.)
#
# Residues near the ends of a sheet are not trustworthy.
# You may want to discard them later by processing the output of this script.
# You can do this using:
# 
#  truncate_tokens.py 1 1 < sheet_resAveDistances.dat | awk '{ for (i=1;i<=NF;i++) print $i }' > sheet_resAveDistances_trunc1_1.dat
#  This will truncate the beginning and ending of each sheet both by 1,
#  as well as converting the list of numbers to a single-column ascii file.
