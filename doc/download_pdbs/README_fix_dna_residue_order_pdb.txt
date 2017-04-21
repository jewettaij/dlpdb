I am writing this documentation several months after I wrote the original program named "fix_dna_residue_order_pdb.py".  If my memory is correct, then the goal of that program was to overcome a difficulty found in the order that residues in double-stranded DNA appear in a PDB file.

Many (but not all) PDB files containing DNA structures with the sequence numbers listed in this order:


i=1 --- i=10
 |       |
i=2 --- i=9
 |       |
i=3 --- i=8
 |       |
i=4 --- i=7
 |       |
i=5 --- i=6

In other words, the residues in one of the two DNA strands chain are listed in order, followed by residues in the other strand (as if it were a hairpin).


-----------------
(NOTE1: Sequence numbers need not begin at "1", as they did in this example.  If any overhangs exist in either strand, it is the user's responsibility to manually delete them beforehand.  Not all PDB structures with DNA are numbered in this way.  Be sure to check before running this program that the order of the residues matches this pattern, and be sure to check the results are reasonable after running this probram.)

-----------------

The residue ordering above makes it somewhat tricky to extract distances and angles between atoms if they don't belong to the same strand, because those atoms are likely stored in very different places in the PDB file.  This script renumbers and reorders the residues so that it is easier to figure out which pairs of residues are neighbors from opposite strands.

After running this script, the residues in the PDB file will be (hopefully)
be numbered this way.  (Again, be sure to check.)


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

-----------------
(NOTE2: Many other software tools written by other people that assume that consecutively numbered residues are connected by a backbone bond.  Some visualization software will attempt to draw backbone bonds connecting successively numbered residues.  i=1,2,3,4,5,6...  However after renumbering, this is no longer true, and the resulting structure will be ugly to look at (although the coordinates should still be accurate.)  It is up to you to use these renumbered PDB files in an appropriate way.)
-----------------
