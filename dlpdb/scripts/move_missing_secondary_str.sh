#!/bin/sh

# syntax (example):
#    This script scans through all the pdb files in the current directory.
#    For each file it:
#      i) checks to see if there is any helix and sheet information
#     ii) If not, moves it to the directory indicated by argument $1

# Typical usage:
#  
#  ls -f1 *.pdb | move_missing_secondary_str.sh pdbs_no_helices_or_sheets

if [ "$#" -lt 1 ]; then
    echo "Error: move_missing_secondary_str.sh requires a directory argument" >&2
fi

if [ ! -d "$1" ]; then
    echo "creating directory $1" >&2
    mkdir "$1"
fi

while read file_pdb; do
    if [ -s "$file_pdb" ]; then
        if has_secondary_str.py < "$file_pdb"; then
            echo "$file_pdb has helix/sheet info" >&2
        else
            echo "$file_pdb lacks helices/sheets. Moving to $1" >&2
            mv -f "$file_pdb" "$1"
        fi
    fi
done
