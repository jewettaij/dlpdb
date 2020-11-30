#!/usr/bin/env sh

mkdir pdbs_FINAL_individual_chains
mkdir pdbs_FINAL_USE_THESE_FOR_ANALYSIS

for f in *.pdb; do
    echo "extracting the chains from $f" >&2
    # Split each PDB file into two chains:
    # This creates 2 pdb files whose suffix is _X.pdb, where X
    # is the Chain-ID letter for the corresponding strand
    select_chains_with_dna.py "$f"
    PDB_CODE=`basename "$f" .pdb`
    # Who knows what these chain-ID letters are?  A? B?  Depends on the file.
    # To remove ambiguity, rename them to _1.pdb and _2.pdb
    i=1
    for g in ${PDB_CODE}_?.pdb; do
	mv "$g" ${PDB_CODE}_${i}.pdb
	i=$((i+1))
    done
    # (Note: You can also use: select_chains_with_dna.py *.pdb)
    dna_interleave_residues.py ${PDB_CODE}_1.pdb  ${PDB_CODE}_2.pdb \
       > pdbs_FINAL_USE_THESE_FOR_ANALYSIS/${PDB_CODE}_interleaved_12.pdb
    # Optional: Create a second PDB file with the order of the chains swapped.
    # (Why? If you swap the order of the chains before interleaving the residues
    #  then later, when you extract angles from atom positions in the resulting
    #  interleaved double stranded DNA structure, the new angles will be
    #  different.  ...But also valid.  Doing this effectively doubles the amount
    #  of data you can collect from the orignal PDB files you started with.)
    dna_interleave_residues.py ${PDB_CODE}_2.pdb  ${PDB_CODE}_1.pdb \
       > pdbs_FINAL_USE_THESE_FOR_ANALYSIS/${PDB_CODE}_interleaved_21.pdb
    mv ${PDB_CODE}_?.pdb pdbs_FINAL_individual_chains/
done
    
mkdir pdbs_FINAL
mv *.pdb pdbs_FINAL/
mv pdbs_with_overhangs pdbs_FINAL/
