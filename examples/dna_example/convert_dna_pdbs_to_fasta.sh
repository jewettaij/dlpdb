#!/usr/bin/env sh

mkdir pdbs_too_short
mkdir pdbs_GC_content_too_high
mkdir pdbs_AT_content_too_high
mkdir pdbs_with_non-standard_bases
rm -f sequences_some_redundant.txt

for f in *.pdb; do
    echo "extracting the sequence from $f" >&2
    PDB_CODE=`basename "$f" .pdb`
    SEQ=`pdb2sequence.py < "$f"`
    # (Note: When run this way, pdb2sequence.py will output the
    #        sequence of both strands, back-to-back.
    #        You CAN select individual strands, but I did not want
    #        to be bothered figuring out the chain-ID letters for
    #        each strand, which differs from PDB file to PDB file.)
    LEN_SUM=`expr length "$SEQ"` #sum of lengths of BOTH strands
    if [ $LEN_SUM -lt 20 ]; then
        echo "$PDB_CODE sequence too short ($LEN_SUM).  Discarding." >&2
        mv $f pdbs_too_short
    else
        # Discard if more than 70% GC content
	NUM_X=$(grep -o "x" <<< "$SEQ" | wc -l)
        NUM_G=$(grep -o "G" <<< "$SEQ" | wc -l)
        NUM_C=$(grep -o "C" <<< "$SEQ" | wc -l)
        GC_CONTENT=`echo "scale = 10; ($NUM_G+$NUM_C)/$LEN_SUM" | bc`
        NUM_A=$(grep -o "A" <<< "$SEQ" | wc -l)
        NUM_T=$(grep -o "T" <<< "$SEQ" | wc -l)
        AT_CONTENT=`echo "scale = 10; ($NUM_A+$NUM_T)/$LEN_SUM" | bc`
	if [ $NUM_X -gt 0 ]; then
            echo "$PDB_CODE contsians non-standard bases.  Discarding." >&2
            mv $f pdbs_with_non-standard_bases
        elif (( $(echo "$GC_CONTENT > 0.6" | bc) )); then
            echo "$PDB_CODE GC conntent too high ($GC_CONTENT).  Discarding." >&2
            mv $f pdbs_GC_content_too_high
        elif (( $(echo "$AT_CONTENT > 0.6" | bc) )); then
            echo "$PDB_CODE AT conntent too high ($GC_CONTENT).  Discarding." >&2
            mv $f pdbs_AT_content_too_high
        else
            echo ">$PDB_CODE" >> sequences_some_redundant.txt
                 echo "$SEQ" >> sequences_some_redundant.txt
        fi
    fi
done
