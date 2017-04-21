#!/usr/bin/env sh

# Run this script in the directory you store your PDB 
# files right after downloading them.

# Typical usage:
#  
# ls -f1 *.pdb | move_non-dna.sh pdbs_with_non-DNA

if [ "$#" -lt 1 ]; then
    echo "Error: move_non-dna.sh requires a directory argument" >&2
fi

if [ ! -d "$1" ]; then
    echo "creating directory $1" >&2
    mkdir "$1"
fi

while read file_pdb; do
    #echo "processing $file_pdb" >&2
    if [ -s "$file_pdb" ]; then
	if has_rna_heavy_atoms.py < "$file_pdb"; then
            echo "$file_pdb appears to contain RNA. Moving to $1" >&2
	    mv -f "$file_pdb" "$1"
	elif has_protein_heavy_atoms.py < "$file_pdb"; then
            echo "$file_pdb appears to contain proteins. Moving to $1" >&2
	    mv -f "$file_pdb" "$1"
	# One final check: Does the first line of the HEADER
	# only contain the word "DNA" ?
	elif `awk '{if ((NR==1) && (NF==4) && ($1=="HEADER") && ($2=="DNA")) {exit 1}}' < $file_pdb`; then
            echo "$file_pdb might contain things that are not DNA. Moving to $1" # >&2
	    mv -f "$file_pdb" "$1"
	else
            echo "$file_pdb appears to contain nothing but DNA. Keeping." >&2
	fi
    fi
done
