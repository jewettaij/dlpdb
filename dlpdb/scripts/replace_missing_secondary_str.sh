#!/bin/sh

# syntax (example):
#    This script scans through all the pdb files in the current directory.
#    For each file it:
#      i) checks to see if there is any helix and sheet information
#     ii) If not, it creates a new pdb file named "oldname+helixsheet.pdb"
#         containing helix and sheet records generated by DSSP
#         (and converted to pdb format using dssp2pdb).

# Typical usage:
#  
#  ls -f1 *.pdb | replace_missing_secondary_str.sh 

while read file_pdb; do
    filenoext="${file_pdb%.pdb}"
    pdb_code=`echo "$file_pdb" | awk '{print substr($0,0,4)}'`
    file_dssp="${pdb_code}.dssp"

    #echo "dssp_file = $file_dssp"

    if has_secondary_str.py < $file_pdb; then
        echo "${file_pdb} has helix/sheet info" >&2
    else
        echo "${file_pdb} is missing helix/sheet info" >&2
        dssp2pdb.py "${file_pdb}" < "${file_dssp}" > helix_sheet.pdb.tmp
        cat "${file_pdb}" helix_sheet.pdb.tmp > "${filenoext}+helixsheet.pdb"
        rm -f helix_sheet.pdb.tmp
    fi
done
