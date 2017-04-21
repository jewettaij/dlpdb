#!/bin/sh

if [ -z "${EXTRACTCOORDS}" ]
then
   EXTRACTCOORDS="pdb2coords_ave.py"
fi


while read pdb_file_name; do
    echo "${0##*/} processing $pdb_file_name" >&2
    "${EXTRACTCOORDS}" < "$pdb_file_name" | awk 'BEGIN{pNF=0} {if (NF==3) {if (pNF==3) {print px" "py" "pz"  "$1" "$2" "$3} px=$1; py=$2; pz=$3} pNF=NF}' | coords2distances.py | tr "\n" " "
    # You can pipe the results to  sed -e 's/\s\+/\n/g' to put on separate lines
    echo ""
done


