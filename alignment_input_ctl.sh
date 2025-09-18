#!/bin bash

# Directories
GENOME="Bio_Project/bowtie2_index/mouse.mm9.genome"
FASTQ_DIR="Bio_Project"

#Samples list
SAMPLES=("SRR5409170.fastq.gz:WT"
"SRR5409172.fastq.gz:Set8KO"
"SRR5409178.fastq.gz:SDS1 WT"
"SRR5409180.fastq.gz:SDS1 Set8KO")

# Align reads, convert, sort, and index for each sample
for SAMPLE in "${SAMPLES[@]}"; do

	# Extract filename and sample name
	FILE=$(echo "$SAMPLE" | cut -d':' -f1)
	NAME=$(echo "$SAMPLE" | cut -d':' -f2)
	
	echo "Processing sample: $NAME ($FILE)..."

	# Align reads
	bowtie2 -x "$GENOME" -U "$FASTQ_DIR/$FILE" -S "$FASTQ_DIR/$NAME.sam" --very-
	sensitive -k 1

	# Convert SAM to BAM
	samtools view -bS "$FASTQ_DIR/$NAME.sam" > "$FASTQ_DIR/$NAME.bam"

	# Sort BAM file
	samtools sort -o "$FASTQ_DIR/$NAME.sorted.bam" "$FASTQ_DIR/$NAME.bam"

	# Index sorted BAM file
	samtools index "$FASTQ_DIR/$NAME.sorted.bam"
	
	echo "Finished processing $NAME."

done

echo "All samples processed successfully."
