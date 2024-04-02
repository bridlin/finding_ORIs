#!/bin/bash
#
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err
#SBATCH --mail-type END
#SBATCH --mail-user b-barckmann@chu-montpellier.fr
#
#
#SBATCH --partition fast
#SBATCH --cpus-per-task 6
#SBATCH --mem 64GB


module load samtools/1.13
module load  macs2/2.2.7.1
module load  bedtools/2.30.0




source scripts/finding_ORIs/config_finding-ORIs.txt

# seperation of the minus and plus strand read pairs:

mkdir $output_dir

for sample in "${input_list[@]}"; do
	samtools view \
		-b \
		-f 128 \
		-F 16 \
		$read_directory/$sample\_$file_prefix\.bam > $read_directory/$sample\_F2.bam &&
	samtools view \
		-b \
		-f 80 \
		$read_directory/$sample\_$file_prefix\.bam > $read_directory/$sample\_R1.bam &&
	samtools merge \
		-f $read_directory/$sample\_F2R1_$file_prefix\.bam \
		$read_directory/$sample\_F2.bam \
		$read_directory/$sample\_R1.bam &&
	samtools index $read_directory/$sample\_F2R1_$file_prefix\.bam &&

	samtools view \
		-b \
		-f 144 \
		$read_directory/$sample\_$file_prefix\.bam > $read_directory/$sample\_R2.bam &&
	samtools view \
		-b \
		-f 64 \
		-F 16 \
		$read_directory/$sample\_$file_prefix\.bam > $read_directory/$sample\_F1.bam &&
	samtools merge \
		-f $read_directory/$sample\_F1R2_$file_prefix\.bam \
		$read_directory/$sample\_R2.bam \
		$read_directory/$sample\_F1.bam &&
	samtools index $read_directory/$sample\_F1R2_$file_prefix\.bam \
; done

#PEAK CALLING ALONE:

for sample in "${input_list[@]}"; do
	macs2 callpeak  \
		--bdg  -t $read_directory/$sample\_F2R1_$file_prefix\.bam   \
		-f BAMPE \
		-n $sample\-alone_Minus_bowtie2_trimmed_uniq_dupsre_narrow_p005   \
		--outdir $output_dir/ \
		-p 5e-2 \
		-s 130 \
		-m 10 30 \
		--gsize 2.5e7 &&
	macs2 callpeak  \
		--bdg  -t $read_directory/$sample\_F1R2_$file_prefix\.bam  \
		-f BAMPE \
		-n $sample\-alone_Plus_bowtie2_trimmed_uniq_dupsre_narrow_p005  \
		--outdir $output_dir/ \  
		-p 5e-2 \
		-s 130 \
		-m 10 30 \
		--gsize 2.5e7 \
; done


#PEAK SELECTION :
#1) Sorting out overlapping peaks:

for overlap in "${overlap_list[@]}" ; do for sample in "${input_list[@]}" ; do
	bedtools intersect \
		-wa  \
		-e  \
		-v \
		-f 0.$overlap  \
		-F 0.$overlap  \
		-a $output_dir/$sample\-alone_Minus_bowtie2_trimmed_uniq_dupsre_narrow_p005_peaks.narrowPeak \
		-b $output_dir/$sample\-alone_Plus_bowtie2_trimmed_uniq_dupsre_narrow_p005_peaks.narrowPeak \
	> $output_dir/$sample\-alone_Minus_nonoverlap$overlap\_narrow_p005_peaks.narrowPeak &&
	bedtools intersect \
		-wa  \
		-e  \
		-v \
		-f 0.$overlap  \
		-F 0.$overlap  \
		-a $output_dir/$sample\-alone_Plus_bowtie2_trimmed_uniq_dupsre_narrow_p005_peaks.narrowPeak \
		-b $output_dir/$sample\-alone_Minus_bowtie2_trimmed_uniq_dupsre_narrow_p005_peaks.narrowPeak \
		> $output_dir/$sample\-alone_Plus_nonoverlap$overlap\_narrow_p005_peaks.narrowPeak \
;done ;done


#2) Selecting peaks in a window:

for window in "${window_list[@]}"; do for sample in "${input_list[@]}" ; do for overlap in "${overlap_list[@]}" ; do
	bedtools window \
		-w $window \
		-a $output_dir/$sample\-alone_Minus_nonoverlap$overlap\_narrow_p005_peaks.narrowPeak \
		-b $output_dir/$sample\-alone_Plus_nonoverlap$overlap\_narrow_p005_peaks.narrowPeak \
		> $output_dir/union$window\_$sample\-alone_Minus-Plus-nonoverlap$overlap\_narrow_p005_peaks.narrowPeak  &&
	cat $output_dir/union$window\_$sample\-alone_Minus-Plus-nonoverlap$overlap\_narrow_p005_peaks.narrowPeak \
	| cut -f 1,2,3,4,5,6,7,8,9,10 \
	| sort -k1,1 -k2,2n  \
	| uniq \
	| awk 'BEGIN { OFS="\t" } {if ($4 ~ /_Minus_/) print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" "-" "\t" $7 "\t"$8 "\t" $9 "\t" $10 ;else print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" "+" "\t" $7 "\t"$8 "\t" $9 "\t" $10}' \
	> $output_dir/Minus-union$window\_$sample\-alone_nonoverlap$overlap\_narrow_p005_peaks.narrowPeak &&
	cat $output_dir/union$window\_$sample\-alone_Minus-Plus-nonoverlap$overlap\_narrow_p005_peaks.narrowPeak \
	| cut -f 11,12,13,14,15,16,17,18,19,20 \
	| sort -k1,1 -k2,2n  \
	| uniq \
	| awk 'BEGIN { OFS="\t" } {if ($4 ~ /_Minus_/) print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" "-" "\t" $7 "\t"$8 "\t" $9 "\t" $10 ;else print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" "+" "\t" $7 "\t"$8 "\t" $9 "\t" $10}'  \
	> $output_dir/Plus-union$window\_$sample\-alone_nonoverlap$overlap\_narrow_p005_peaks.narrowPeak &&
	cat \
	$output_dir/Minus-union$window\_$sample\-alone_nonoverlap$overlap\_narrow_p005_peaks.narrowPeak  \
	$output_dir/Plus-union$window\_$sample\-alone_nonoverlap$overlap\_narrow_p005_peaks.narrowPeak \
	| sort -k1,1 -k2,2n \
	> $output_dir/catMinus-Plus_union$window\_$sample\-alone_nonoverlap$overlap\_narrow_p005_peaks.narrowPeak  \
;done ;done ;done

#ORI SELECTION
#Ori-selection strand-seperated Minus-Plus direction 

for window in "${window_list[@]}"; do for sample in "${input_list[@]}" ; do for overlap in "${overlap_list[@]}" ; do
	bedtools closest \
		-iu \
		-D ref \
		-t all \
		-a $output_dir/Minus-union$window\_$sample\-alone_nonoverlap$overlap\_narrow_p005_peaks.narrowPeak  \
		-b $output_dir/Plus-union$window\_$sample\-alone_nonoverlap$overlap\_narrow_p005_peaks.narrowPeak \
	| awk -v dist="$window"  'BEGIN { OFS="\t" } {if ($12 != -1 && $21 < dist ) {print $0} }'\
	> $output_dir/closest-iu_Minus-Plus_union$window\_$sample\-alone_nonoverlap$overlap_narrow_p005_peaks.narrowPeak &&
	bedtools closest \
		-id \
		-D ref \
		-t all \
		-a $output_dir/Plus-union$window\_$sample\-alone_nonoverlap$overlap\_narrow_p005_peaks.narrowPeak \
		-b $output_dir/Minus-union$window\_$sample\-alone_nonoverlap$overlap\_narrow_p005_peaks.narrowPeak  \
	| awk -v dist="$window"  'BEGIN { OFS="\t" } {if ($12 != -1 && $21 > -dist ) {print $0} }'  \
	> $output_dir/closest-id_Plus-Minus_union$window\_$sample\-alone_nonoverlap$overlap\_narrow_p005_peaks.narrowPeak &&
	cat $output_dir/closest-iu_Minus-Plus_union$window\_$sample\-alone_nonoverlap$overlap\_narrow_p005_peaks.narrowPeak \
	| awk 'BEGIN { OFS="\t" } {if ($21 > 0 ) {print $1 "\t" $3 "\t" $12 "\t" $4 "\t" "." "\t" $6 "\t" $21}   else if ($21 == 0 && $2 < $12 ) {print $1 "\t" $12 "\t" $3 "\t" $4 "\t" "." "\t" $6 "\t" $21} }' \
	> $output_dir/Minus_oris_union$window\_$sample\-alone_nonoverlap$overlap\_narrow_p005_peaks.bed &&
	cat $output_dir/closest-id_Plus-Minus_union$window\_$sample\-alone_nonoverlap$overlap\_narrow_p005_peaks.narrowPeak \
	| awk 'BEGIN { OFS="\t" } {if ($21 < 0 ) {print $1 "\t" $13 "\t" $2 "\t" $4 "\t" "." "\t" "-" "\t" $21}  else if ($21 == 0 && $2 > $12 ) {print $1 "\t" $2 "\t" $13 "\t" $4 "\t" "." "\t" "-" "\t" $21} }' \
	> $output_dir/Plus_oris_union$window\_$sample\-alone_nonoverlap$overlap\_narrow_p005_peaks.bed &&
	cat \
	$output_dir/Plus_oris_union$window\_$sample\-alone_nonoverlap$overlap\_narrow_p005_peaks.bed  \
	$output_dir/Minus_oris_union$window\_$sample\-alone_nonoverlap$overlap\_narrow_p005_peaks.bed \
	| sort -k1,1 -k2,2n \
	| awk '!a[$1 $2 $3]++' | sort -rk2 \
	| awk '!seen[$1 $3]++' | sort -k3  \
	| awk '!seen[$1 $2]++' | sort -k1,1 -k2,2n \
	> $output_dir/ORI_$sample\-alone_union$window\_nonoverlap$overlap\_narrow_p005_peaks_withStrand.bed   \
	awk ' {print  $1 "\t" $2  "\t" $3  "\t" $4  "\t" $5  "\t""." "\t" $7}' $output_dir/ORI_$sample\-alone_union$window\_nonoverlap$overlap\_narrow_p005_peaks_withStrand.bed  > $output_dir/ORI_$sample\-alone_union$window\_nonoverlap$overlap\_narrow_p005_peaks.bed
;done ;done ;done
