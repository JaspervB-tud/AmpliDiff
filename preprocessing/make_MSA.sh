SEQFOLDER="$1"
FILE="$SEQFOLDER"/sequences_aligned.fasta

##Check if alignment already exists in same folder as sequences, otherwise perform MSA with MAFFT and output:
##  - sequences_aligned.fasta (MSA)
##  - runtime_alignment.txt (Runtime for MSA)

if [ -f "$FILE" ]; then
	echo "$FILE exists."
else
	SECONDS=0
	echo "$FILE does not exist."
	mafft --auto --leavegappyregion "$SEQFOLDER"/sequences.fasta > "$SEQFOLDER"/sequences_aligned.fasta
	echo "Elapsed: $(($SECONDS / 3600)) hours, $((($SECONDS / 60) % 60)) minutes, $((SECONDS % 60)) seconds" > "$outfolder"/runtime_alignment.txt
fi