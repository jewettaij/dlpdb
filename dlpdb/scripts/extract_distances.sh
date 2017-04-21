#!/bin/sh

ATOM_SELECTION="$@"

if [ -z "${EXTRACTCOORDS}" ]
then
   EXTRACTCOORDS="pdb2coords.py -blank"
fi

while read pdb_file_name; do
    echo "${0##*/} processing $pdb_file_name" >&2
    print_coords_command="${EXTRACTCOORDS} ${ATOM_SELECTION} < $pdb_file_name"
    eval "$print_coords_command" | coords2distances.py | tr "\n" " "
    # You can pipe the results to  sed -e 's/\s\+/\n/g' to put on separate lines
    echo ""
done

