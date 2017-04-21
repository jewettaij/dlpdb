#!/bin/sh

# Run this script in the directory you store your PDB 
# files right after downloading them.

# Typical usage:
#  
#  ls -f1 *.pdb | move_missing_dna_heavy_atoms.sh pdbs_missing_heavy_atoms

if [ "$#" -lt 1 ]; then
    echo "Error: move_missing_dna_heavy_atoms.sh requires a directory argument" >&2
fi

if [ ! -d "$1" ]; then
    echo "creating directory $1" >&2
    mkdir "$1"
fi

while read file_pdb; do
    if [ -s "$file_pdb" ]; then
        if has_dna_heavy_atoms.py < "$file_pdb"; then
            echo "$file_pdb has heavy atoms" >&2
	else
            echo "$file_pdb lacks heavy atoms. Moving to $1" >&2
	    mv -f "$file_pdb" "$1"
	fi
    fi
done
