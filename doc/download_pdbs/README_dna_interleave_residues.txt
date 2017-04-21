dna_interleave_residues.py

This program reads two PDB files (containing a single chain each)
and creates a new PDB file with a single chain which "interleaves"
the residues of the two PDB files.  The two PDB files should contain
chains of the same length, and in opposite directions.
This program was designed to be used to manipulate PDB files containing
duplex DNA.  It is the user's responsibility to insure that there
are no overhangs.  (See example below using "select_interval.py")


## Typical usage:

# First select out chains A and B from "309.pdb", selecting the residues
# we want carefully (A3-10 and B13-20) to avoid overhangs::

   select_interval.py A  3 " " A 10 z < 309d.pdb > 309d_A3-10.pdb
   select_interval.py B 13 " " B 20 z < 309d.pdb > 309d_B13-20.pdb

# Then combine them into a single PDB file:

   dna_interleave_residues.py ../309d_A3-10.pdb ../309d_B13-20.pdb \
              > 309d_interleaved.pdb



## Motivation:

Most (but not all) PDB files containing DNA structures with the sequence numbers listed in this order:

Chain  chain
  A      B


i=1 --- i=10
 |       |
i=2 --- i=9
 |       |
i=3 --- i=8
 |       |
i=4 --- i=7
 |       |
i=5 --- i=6

In other words, the residues in one of the DNA strands (in either chain) are listed in order, followed by residues in the other strand (as if it were a hairpin).

The residue ordering above makes it somewhat tricky to extract distances and angles between atoms if they don't belong to the same strand, because those atoms are likely stored in very different places in the PDB file.  This script renumbers and reorders the residues so that it is easier to figure out which pairs of residues are neighbors from opposite strands.

After running this script, the residues in the PDB file will be (hopefully)
be numbered this way.  (Again, be sure to check.)


  chain X
  

i=1 --- i=2
 |       |
i=3 --- i=4
 |       |
i=5 --- i=6
 |       |
i=7 --- i=8
 |       |
i=9 --- i=10
     :

This makes it easier for other scripts to find pairs of atoms in opposite strands, because they will now be located nearby in the PDB file.  You can then use other scripts like pdb2coords.py to extract the coordinate in each base pair (and then send these coordinates to other programs like coords2distances.py).

WARNING: Not all PDB structures with DNA are numbered in this way.  Be sure to check before running this program that the order of the residues matches this pattern, and be sure to check the results are reasonable after running this probram.)

-----------------
(NOTE2: Many other software tools written by other people that assume that consecutively numbered residues are connected by a backbone bond.  Some visualization software will attempt to draw backbone bonds connecting successively numbered residues.  i=1,2,3,4,5,6...  However after renumbering, this is no longer true, and the resulting structure will be ugly to look at (although the coordinates should still be accurate.)  It is up to you to use these renumbered PDB files in an appropriate way.)
-----------------
