#!/bin/bash
# Run MACS for peak calling
echo "Running MACS for peak calling..."

# Loop through the samples
for sample_name in sample22 sample26 sample30 sample24 sample28 sample32
sample38 sample40
do
	# Determine the control file based on sample name
	case "$sample_name" in
		sample22|sample26|sample30)
			control="WT.sorted.bam"
			;;
		sample24|sample28|sample32)
			control="Set8KO.sorted.bam"
			;;
		sample38)
			control="SDS1_WT.sorted.bam"
			;;
		sample40)
			control="SDS1_Set8KO.sorted.bam"
			;;
		*)
			echo "No control file defined for $sample_name!"
			continue
			;;
esac

# Run MACS with the determined control
echo "Running MACS for $sample_name with control $control..."
macs -t "${sample_name}.sorted.bam" -c "$control" --format BAM --name
"$sample_name" --gsize 138000000 --tsize 26 --diag --wig

done
