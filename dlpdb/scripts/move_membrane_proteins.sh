#!/bin/sh

# Run this script in the directory you store your PDB 
# files right after downloading them.

# Typical usage:
#  
#  ls -f1 *.pdb | move_membrane_proteins.sh pdbs_membrane

if [ "$#" -lt 1 ]; then
    echo "move_membrane_proteins.sh requires a directory argument" >&2
fi

if [ ! -d "$1" ]; then
    echo "creating directory $1" >&2
    mkdir "$1"
fi

while read file_pdb; do
    if [ -s "$file_pdb" ]; then
        if grep -i "MEMBRANE\|MICELLE" "$file_pdb"; then
            echo "$file_pdb appears to be a membrane protein" >&2
            mv -f "$file_pdb" "$1"
        else
            echo "$file_pdb does not appear to be a membrane protein" >&2
        fi
    fi
done

