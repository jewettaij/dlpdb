#!/bin/sh

if [ -z "${EXTRACTCOORDS}" ]
then
   EXTRACTCOORDS="pdb2coords_ave.py"
fi

while read pdb_file_name; do
    #echo "pdb_file = \"${pdb_file_name}\""
    echo "${0##*/} processing $pdb_file_name" >&2
    pdb2turn.py < "$pdb_file_name" > turn_intervals_resAveDistances.tmp
    while read interval; do
	command="select_interval.py $interval < $pdb_file_name | ${EXTRACTCOORDS} | awk 'BEGIN{pNF=0} {if (NF==3) {if (pNF==3) {print px\" \"py\" \"pz\"  \"\$1\" \"\$2\" \"\$3} px=\$1; py=\$2; pz=\$3} pNF=NF}' | coords2distances.py | tr \"\n\" \" \""
        #echo "$command"
        eval $command
	echo ""
    done < turn_intervals_resAveDistances.tmp
done

rm -f turn_intervals_resAveDistances.tmp

