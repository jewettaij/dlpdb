#!/bin/sh

ATOM_SELECTION="$@"

if [ -z "${EXTRACTCOORDS}" ]
then
   EXTRACTCOORDS="pdb2coords.py ${ATOM_SELECTION} -blank"
fi


while read pdb_file_name; do
    #echo "pdb_file = \"${pdb_file_name}\""
    echo "${0##*/} processing $pdb_file_name" >&2
    pdb2turn.py < "$pdb_file_name" > turn_intervals_angles.tmp
    while read interval; do
        command="select_interval.py $interval < \"$pdb_file_name\" | ${EXTRACTCOORDS} | coords2angles.py | tr \"\n\" \" \""
        #echo "${command}"
        eval "$command"
	echo ""
    done < turn_intervals_angles.tmp
done

rm -f turn_intervals_angles.tmp

