#!/bin/bash
 
# Directory containing the FASTQ files
FASTQ_DIR="Bio_Project"
# Bowtie2 index and genome files
BOWTIE2_INDEX="Bio_Project/bowtie2_index/mm9"
GENOME="Bio_Project/bowtie2_index/mouse.mm9.genome"
GENOME_FA="Bio_Project/bowtie2_index/mm9.fa"
 
# List of samples
SAMPLES=("SRR5409158.fastq.gz:sample22"
    "SRR5409160.fastq.gz:sample24"
    "SRR5409162.fastq.gz:sample26"
    "SRR5409164.fastq.gz:sample28"
    "SRR5409166.fastq.gz:sample30"
    "SRR5409168.fastq.gz:sample32"
    "SRR5409170.fastq.gz:sample34"
    "SRR5409172.fastq.gz:sample36"
    "SRR5409174.fastq.gz:sample38"
    "SRR5409176.fastq.gz:sample40"
    "SRR5409178.fastq.gz:sample42"
    "SRR5409180.fastq.gz:sample42")
 
# Function to process each sample step by step
process_step() {
    local step="$1"
    local file_suffix="$2"
   
    for sample in "${SAMPLES[@]}"; do
        fastq_file="${sample%%:*}"
        sample_name="${sample##*:}"
       
        echo "Processing $step for $sample_name..."
       
        case "$step" in
            "fastqc")
                # Step 0: Run FASTQC on the FASTQ file
                echo "Running FASTQC on $fastq_file for $sample_name..."
                fastqc "$FASTQ_DIR/$fastq_file" -o "$FASTQ_DIR/fastqc_reports"
                ;;
            "align")
                # Step 1: Align with Bowtie2 and create SAM file
                echo "Aligning $fastq_file with Bowtie2 for $sample_name..."
                bowtie2 -p 4 -x "$BOWTIE2_INDEX" -U "$FASTQ_DIR/$fastq_file" -S "$sample_name.sam"
                ;;
            "convert_bam")
                # Step 2: Convert SAM to BAM
                echo "Converting SAM to BAM for $sample_name..."
                samtools view -bSo "$sample_name.bam" "$sample_name.sam"
                # Remove SAM file after conversion
                rm -f "$sample_name.sam"
                ;;
            "sort_bam")
                # Step 3: Sort BAM file
                echo "Sorting BAM file for $sample_name..."
                samtools sort -o "$sample_name.sorted.bam" "$sample_name.bam"
                rm -f "$sample_name.bam"  # Remove original BAM after sorting
                ;;
           "index_bam")
                # Step 4: Index sorted BAM file
                echo "Indexing sorted BAM file for $sample_name..."
                samtools index "$sample_name.sorted.bam"
                ;;
            "flagstat")
                # Step 5: Generate flagstat statistics
                echo "Generating flagstat for $sample_name..."
                samtools flagstat "$sample_name.sorted.bam" > "$sample_name.flagstat.txt"
                ;;
            "bedgraph")
                # Step 6: Generate BedGraph file
                echo "Generating BedGraph for $sample_name..."
                genomeCoverageBed -bg -ibam "$sample_name.sorted.bam" -g "$GENOME" > "$sample_name.bedgraph"
                ;;
            "bigwig")
                # Step 7: Convert BedGraph to BigWig
                echo "Converting BedGraph to BigWig for $sample_name..."
                bedGraphToBigWig "$sample_name.bedgraph" "$GENOME" "$sample_name.bw"
                ;;
           "macs")
                # Step 8: Run MACS2 for peak calling
                echo "Running MACS2 for peak calling on $sample_name..."
                macs2 callpeak -t "$sample_name.sorted.bam" -c "gfp.sorted.bam" --format BAM --name "$sample_name" --gsize 138000000 --tsize 26 --diag --wig
                ;;
            "summits")
                # Step 9: Process summits
                echo "Processing summits for $sample_name..."
                slopBed -i "${sample_name}_summits.bed" -g "$GENOME" -b 20 > "${sample_name}_summits-b20.bed"
                fastaFromBed -fi "$GENOME_FA" -bed "${sample_name}_summits-b20.bed" > "${sample_name}_summits-b20.fa"
                ;;
            "meme")
                # Step 10: Run MEME on summits
                echo "Running MEME on summits for $sample_name..."
                meme "${sample_name}_summits-b20.fa" -o "${sample_name}_meme" -dna
                ;;
            "cleanup")
                # Step 11: Clean up after processing
                echo "Cleaning up intermediate files for $sample_name..."
                rm -f "$sample_name.bedgraph" "$sample_name.flagstat.txt"
                ;;
        esac
    done
}
 
# Step 1: Build Bowtie2 index for mm9 genome (run this only once)
echo "Building Bowtie2 index for mm9 genome..."
bowtie2-build "$GENOME_FA" "$BOWTIE2_INDEX"
echo "Bowtie2 index for mm9 genome built."
 
# Step 2: Build Bowtie2 index for GFP control (run this only once)
echo "Building Bowtie2 index for GFP control..."
bowtie2-build "$FASTQ_DIR/gfp.fastq" "$FASTQ_DIR/gfp_index"
echo "Bowtie2 index for GFP control built."
 
# Step 3: Create a directory for FASTQC reports
mkdir -p "$FASTQ_DIR/fastqc_reports"
 
# Step 4: Run FASTQC for all samples
echo "Running FASTQC for all samples..."
process_step "fastqc"
 
# Step 5: Process all steps for all samples
echo "Starting alignment for all samples..."
process_step "align"
 
echo "Converting SAM to BAM for all samples..."
process_step "convert_bam"
 
echo "Sorting BAM files for all samples..."
process_step "sort_bam"
 
echo "Indexing sorted BAM files for all samples..."
process_step "index_bam"
 
echo "Generating flagstat for all samples..."
process_step "flagstat"
 
echo "Generating BedGraph for all samples..."
process_step "bedgraph"
 
echo "Converting BedGraph to BigWig for all samples..."
process_step "bigwig"
 
echo "Running MACS2 peak calling for all samples..."
process_step "macs"
 
echo "Processing summits for all samples..."
process_step "summits"
 
echo "Running MEME for all samples..."
process_step "meme"
 
echo "Cleaning up intermediate files for all samples..."
process_step "cleanup"
 
echo "All samples processed."
