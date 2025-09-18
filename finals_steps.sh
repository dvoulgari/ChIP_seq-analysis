#!/bin/bash

# Bowtie2 index and genome files
BOWTIE2_INDEX="Bio_Project/bowtie2_index/mm9"
GENOME="Bio_Project/bowtie2_index/mouse.mm9.genome"
GENOME_FA="Bio_Project/bowtie2_index/mm9.fa"

# List of samples
SAMPLES=("sample22"
"sample24"
"sample26"
"sample28"
"sample30"
"sample32"
"sample38"
"sample40")

# Function to process each sample step by step
process_step() {
	local step="$1"
	for sample in "${SAMPLES[@]}"; do

	sample_name="$sample"
	echo "Processing $step for $sample_name..."

	case "$step" in
		"summits")

		# Step 9: Process summits
		echo "Processing summits for $sample_name..."
		slopBed -i "${sample_name}_summits.bed" -g "$GENOME" -b 20 > "$
{sample_name}_summits-b20.bed"
		fastaFromBed -fi "$GENOME_FA" -bed "${sample_name}_summits-
b20.bed" > "${sample_name}_summits-b20.fa"
		;;

	"meme")
		# Step 10: Run MEME on summits
		echo "Running MEME on summits for $sample_name..."
		meme "${sample_name}_summits-b20.fa" -o "${sample_name}_meme" -
dna
		;;
	esac
done
}

# Process all steps for all samples
echo "Processing summits for all samples..."
process_step "summits"

echo "Running MEME for all samples..."
process_step "meme"

echo "All samples processed."
