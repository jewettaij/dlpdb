This program is not used for extracting coordinates from PDB files.  Instead, it is useful for extracting distances and angles from polymer coordinate data in trajectory files (in ".RAW" format).  I include it in this directory because eventually, you may want to compare simulated conformations with PDB structures, and this program is compatible with the other scripts in this directory used for measuring angles and distances.

This script will print coordinates from atoms which participate in the distance or angle you are interested all on the same line.  (This mimics the way coordinates are extracted using the "pdb2coords.py" script.)  Then these coordinates can be processed using other simple scripts like "coords2distances.py", "coords2angles.py", which reads one line at a time and calculates distances and angles, respectively.  This script was designed to work with polymers whose monomers contain the same number of atoms in every monomer (this is often true of many simple coarse-grained models.)

As an example, I will demonstrate what happens when I apply this program to the simple trajectory file ("coordinates.raw")

Suppose we have a short polymer (4 monomers), with two atoms in each monomer ("-p 2").  Suppose we want to measure the angle between the backbone atoms in 3 successive monomers.  First we extract their coordinates using:

merge_lines_periodic.py -2 0 2 -p 2 < coordinates.raw

For reference, "coordinates.raw" is a 3-column numeric text file containing coordinate data, with blank lines delimiting data recorded at different times during the simulation.
----- coordinates.raw -----
#snapshot #1 from trajectory:
-0.301 -3.622 2.722
-0.793 -3.543 2.273
-0.466 -3.883 2.040
-0.007 -3.824 2.566
-0.425 -3.548 2.719
-0.579 -3.734 2.352
0.0557 -3.812 2.580
-0.281 -3.404 2.870

#snapshot #2 from trajectory:
-1.883 -2.285 0.371
-1.770 -2.263 0.579
-1.957 -2.287 -0.003
-1.906 -2.387 0.167
-1.867 -2.163 0.470
-1.963 -1.763 -0.050
-1.912 -2.141 -0.383
-1.830 -2.541 0.131
----------------------
(Note: Comment-only lines are ignored.  They are not treated as blank lines.)

When applied to this file, this command:

merge_lines_periodic.py -2 0 2 -p 2 < coordinates.raw

...divides each snapshot of the trajectory into "blocks" of size 2 lines each (corresponding to the 2 atoms in each monomer as indicated by "-p 2").  In each "block", these atoms will be selected: "-2 0 2".  This means:

-2  <==>  the first atom in the previous block
0   <==>  the first atom from this block  (indexing begins at "0" not "1")
2   <==>  the first atom in the next block

Then it will print coordinates for these 3 atoms on the same line,
once for each successive monomer in the chain, repetitively for every
snapshot in the trajectory.

The polymer has length 4 monomers (8 atoms), and this corresponds to 2 backbone angles (between monomers 1-2-3, and 2-3-4).  Consequently, it will print 2 lists of numbers for each snapshot:

--------------------------
-0.301 -3.622 2.722 -0.466 -3.883 2.040 -0.425 -3.548 2.719
-0.466 -3.883 2.040 -0.425 -3.548 2.719 0.0557 -3.812 2.580

-1.883 -2.285 0.371 -1.957 -2.287 -0.003 -1.867 -2.163 0.470
-1.957 -2.287 -0.003 -1.867 -2.163 0.470 -1.912 -2.141 -0.383
--------------------------
(The first two lines are from the first snapshot.  The next two lines come from the second snapshot.  A blank line is inserted between data from different snapshots.)

This text can then be passed to "coords2angles.py" to calculate the angle betwen the 3 atom on each line.  When reading this data "coords2angles.py" will be confused by the blank line:

merge_lines_periodic.py -2 0 2 -p 2 < coordinates.raw | ./coords2angles.py
This prints out:

10.9760620712
67.6525732686
-720
14.1473978473
17.6722430771
"-720" is an impossible value meant to indicate something went wrong (the blank line).  To avoid that, pipe this text through awk '{if (NF!=0) print $0}' beforehand:

merge_lines_periodic.py -2 0 2 -p 2 < coordinates.raw \
       | awk '{if (NF!=0) print $0}' \
       | ./coords2angles.py

This prints out:

10.9760620712
67.6525732686
14.1473978473
17.6722430771


