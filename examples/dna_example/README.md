Extracting distances and angles from double-stranded DNA
=======
This directory contains the scripts I used to extract distances and angles
from PDB files containing double-stranded DNA.
I was using these distances and angles to derive force-field parameters
for the simulated coarse-grained DNA model explained here
[here](https://github.com/jewettaij/moltemplate/tree/master/examples/coarse_grained/DNA_models/dsDNA_only/2strands/3bp_2particles/simple_dna_example).

The distances and angles I measured are shown here:
![distances_and_angles_from_DNA](./statistics_keeping_every_3rd_base_pair/dsDNA_3to1_C3p.png)

### "Interleaving" the DNA monomers

The scripts that were included with *dlpdb* were originally intended to be
used to extract distances and angles from PDB files containing proteins.
As you walk along the length of a protein, starting from the N terminal,
the atoms that you encounter will belong to monomers that are numbered
sequentially in the PDB file. This is also true for single-stranded DNA and RNA.

However it is not true for double-stranded DNA.
PDB files containing double-stranded DNA usually consist of two chains,
and each chain represents one strand of DNA.
In addition, the two strands are pointing in opposite directions,
so as you walk from the 5' end to the 3' end from one DNA strand,
you travel in the opposite direction on the other DNA strand.
This makes it impossible to extract distances or angles between atoms
using the scripts included with *dlpdb* if these atoms are in opposite strands.

So before we can use the *dlpdb* scripts to extract geometry from these
PDB files, we need to convert the PDB file containing two different chains
into a PDB containing only one chain.  In this new PDB file, the single chain
contains one "polymer".  Each "monomer" in this polymer contains a base-pair,
and they are numbered sequentially from one end of the double-stranded DNA
molecule to the other.  This is accomplished using this script:
[interleave_nucleotides_in_duplex_dna.sh](interleave_nucleotides_in_duplex_dna.sh)

*Note: Before this can be accomplished, any overhangs or bubbles that exist in
the DNA molecule must be removed.  (Currently this must be done manually.)*
