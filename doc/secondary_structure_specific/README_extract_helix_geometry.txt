The documentation for

extract_helix_angles.sh
extract_helix_dihedrals.sh
extract_helix_distances.sh
extract_helix_resAveDistances.sh

is identical to:

README_extract_sheet_angles.txt
README_extract_sheet_dihedrals.txt
README_extract_sheet_distances.txt
README_extract_sheet_resAveDistances.txt

(Just replace the word "sheet" by "helix" everywhere.)

   There is one minor difference in extract_helix_dihedrals.sh
Dihedral angles returned by extract_helix_dihedrals.sh
lie in the range from -180.0 to 180.0 by default
(instead of 0 to 360 as they are extract_sheet_dihedrals.sh by default).
However this can be overridden by supplying a numeric argument before
the atom types. (See the 180 argument in README_extract_sheet_dihedrals.txt.)

