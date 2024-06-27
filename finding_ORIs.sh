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

### seperation of the minus and plus strand read pairs: The read pairs are seperated by the flag information in the bam file. F2R1 for minus read pairs and F1R2 for plus read pairs. The bam files are indexed for further processing.

mkdir $output_dir
mkdir $output_dir/peak_calling
mkdir $output_dir/peak_filtering
mkdir $output_dir/oris

# for sample in "${input_list[@]}"; do
# 	samtools view \
# 		-b \
# 		-f 128 \
# 		-F 16 \
# 		$read_directory/$sample\_$file_prefix\.bam > $read_directory/$sample\_F2.bam &&
# 	samtools view \
# 		-b \
# 		-f 80 \
# 		$read_directory/$sample\_$file_prefix\.bam > $read_directory/$sample\_R1.bam &&
# 	samtools merge \
# 		-f $read_directory/$sample\_F2R1_$file_prefix\.bam \
# 		$read_directory/$sample\_F2.bam \
# 		$read_directory/$sample\_R1.bam &&
# 	samtools index $read_directory/$sample\_F2R1_$file_prefix\.bam &&

# 	samtools view \
# 		-b \
# 		-f 144 \
# 		$read_directory/$sample\_$file_prefix\.bam > $read_directory/$sample\_R2.bam &&
# 	samtools view \
# 		-b \
# 		-f 64 \
# 		-F 16 \
# 		$read_directory/$sample\_$file_prefix\.bam > $read_directory/$sample\_F1.bam &&
# 	samtools merge \
# 		-f $read_directory/$sample\_F1R2_$file_prefix\.bam \
# 		$read_directory/$sample\_R2.bam \
# 		$read_directory/$sample\_F1.bam &&
# 	samtools index $read_directory/$sample\_F1R2_$file_prefix\.bam \
# ; done

###PEAK CALLING without control (alone): done seperatly for minus and plus strand originating read pairs. Narrow peaks are called with a p-value of 5e-2. The effective genome size is set to 2.5e7 bp for T.brucei.



for sample in "${input_list[@]}"; do
	macs2 callpeak  \
		--bdg  \
		-t $read_directory/$sample\_F2R1_$file_prefix\.bam   \
		-c $read_directory/$sample\_control_F2R1_$file_prefix\.bam  \
		-f BAMPE \
		-n $sample\-alone_Minus_bowtie2_trimmed_uniq_dupsre_narrow   \
		--outdir $output_dir/peak_calling/ \
		-s 130 \
		-p 5e-5 \
		-m 10 30 \
		--slocal 10000\
		--llocal 100000\
		--gsize 2.5e7 &&
	macs2 callpeak  \
		--bdg  \
		-t $read_directory/$sample\_F1R2_$file_prefix\.bam  \
		-c $read_directory/$sample\_control_F1R2_$file_prefix\.bam  \
		-f BAMPE \
		-n $sample\-alone_Plus_bowtie2_trimmed_uniq_dupsre_narrow  \
		--outdir $output_dir/peak_calling/ \
		-s 130 \
		-p 5e-5 \
		-m 10 30 \
		--slocal 10000\
		--llocal 100000\
		--gsize 2.5e7 \
;done

#--nolambda \
#--nolambda \
#-s size of the reads
###PEAK FILTERING: 
###1) Sorting out overlapping peaks: we used 50% of maximal overlap for at least one peak for the selection of non-overlapping peaks. 

printf "peak filtering\n"

for overlap in "${overlap_list[@]}" ; do for sample in "${input_list[@]}" ; do
	bedtools intersect \
		-wa  \
		-e  \
		-v \
		-f 0.$overlap  \
		-F 0.$overlap  \
		-a $output_dir/peak_calling/$sample\-alone_Minus_bowtie2_trimmed_uniq_dupsre_narrow_peaks.narrowPeak \
		-b $output_dir/peak_calling/$sample\-alone_Plus_bowtie2_trimmed_uniq_dupsre_narrow_peaks.narrowPeak \
	> $output_dir/peak_filtering/$sample\-alone_Minus_nonoverlap$overlap\_narrow_peaks.narrowPeak &&
	bedtools intersect \
		-wa  \
		-e  \
		-v \
		-f 0.$overlap  \
		-F 0.$overlap  \
		-a $output_dir/peak_calling/$sample\-alone_Plus_bowtie2_trimmed_uniq_dupsre_narrow_peaks.narrowPeak \
		-b $output_dir/peak_calling/$sample\-alone_Minus_bowtie2_trimmed_uniq_dupsre_narrow_peaks.narrowPeak \
	> $output_dir/peak_filtering/$sample\-alone_Plus_nonoverlap$overlap\_narrow_peaks.narrowPeak \
;done ;done

printf "overlapping peak filtered\n"


###2) Selecting peaks in a selected window: we used 500 bp max distance. 

for window in "${window_list[@]}"; do for sample in "${input_list[@]}" ; do for overlap in "${overlap_list[@]}" ; do
	bedtools window \
		-w $window \
		-a $output_dir/peak_filtering/$sample\-alone_Minus_nonoverlap$overlap\_narrow_peaks.narrowPeak \
		-b $output_dir/peak_filtering/$sample\-alone_Plus_nonoverlap$overlap\_narrow_peaks.narrowPeak \
		> $output_dir/peak_filtering/union$window\_$sample\-alone_Minus-Plus-nonoverlap$overlap\_narrow_peaks.narrowPeak  &&
	# seperating the bedtools window output in Minus and Plus peaks tables
	cat $output_dir/peak_filtering/union$window\_$sample\-alone_Minus-Plus-nonoverlap$overlap\_narrow_peaks.narrowPeak \
	| cut -f 1,2,3,4,5,6,7,8,9,10 \
	| sort -k1,1 -k2,2n  \
	| uniq \
	| awk 'BEGIN { OFS="\t" } {if ($4 ~ /_Minus_/) print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" "-" "\t" $7 "\t"$8 "\t" $9 "\t" $10 ;else print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" "+" "\t" $7 "\t"$8 "\t" $9 "\t" $10}' \
	> $output_dir/peak_filtering/Minus-union$window\_$sample\-alone_nonoverlap$overlap\_narrow_peaks.narrowPeak &&
	cat $output_dir/peak_filtering/union$window\_$sample\-alone_Minus-Plus-nonoverlap$overlap\_narrow_peaks.narrowPeak \
	| cut -f 11,12,13,14,15,16,17,18,19,20 \
	| sort -k1,1 -k2,2n  \
	| uniq \
	| awk 'BEGIN { OFS="\t" } {if ($4 ~ /_Minus_/) print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" "-" "\t" $7 "\t"$8 "\t" $9 "\t" $10 ;else print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" "+" "\t" $7 "\t"$8 "\t" $9 "\t" $10}'  \
	> $output_dir/peak_filtering/Plus-union$window\_$sample\-alone_nonoverlap$overlap\_narrow_peaks.narrowPeak \
;done ;done ;done

printf "peak pairs selected\n"

###3)last filtering step : selecting peak pairs in minus upstream plus downstram direction
###ORI SELECTION  defining origina as the region between filtered peak pairs

###  (a) closest iu: closest peak downstream  to the minus strand peaks, awk line to exclude all distances bigger than the window size

for window in "${window_list[@]}"; do for sample in "${input_list[@]}" ; do for overlap in "${overlap_list[@]}" ; do
	bedtools closest \
		-iu \
		-D ref \
		-t all \
		-a $output_dir/peak_filtering/Minus-union$window\_$sample\-alone_nonoverlap$overlap\_narrow_peaks.narrowPeak  \
		-b $output_dir/peak_filtering/Plus-union$window\_$sample\-alone_nonoverlap$overlap\_narrow_peaks.narrowPeak | \
		awk -v dist="$window"  'BEGIN { OFS="\t" } {if ($12 != -1 && $21 < dist ) {print $0} }'\
	> $output_dir/peak_filtering/closest-iu_Minus-Plus_union$window\_$sample\-alone_nonoverlap$overlap\_narrow_peaks.narrowPeak &&
### (b) closest id: closest peak upstream to the plus strand peaks, awk line to exclude all distances bigger than the window size	
	bedtools closest \
		-id \
		-D ref \
		-t all \
		-a $output_dir/peak_filtering/Plus-union$window\_$sample\-alone_nonoverlap$overlap\_narrow_peaks.narrowPeak \
		-b $output_dir/peak_filtering/Minus-union$window\_$sample\-alone_nonoverlap$overlap\_narrow_peaks.narrowPeak |  \
		awk -v dist="$window"  'BEGIN { OFS="\t" } {if ($12 != -1 && $21 > -dist ) {print $0} }'  \
	> $output_dir/peak_filtering/closest-id_Plus-Minus_union$window\_$sample\-alone_nonoverlap$overlap\_narrow_peaks.narrowPeak &&
### (c) selection of the origins as the region between a minus downstream and a following plus peaks from the closest output table:
### to be considered: the closest iu option is not ignoring upstream peaks that are overlapping. Those are filtered out by distance ($21) smaller equal 0 and start coordinate of a ($2) bigger than start coordinate of b ($12).
	cat $output_dir/peak_filtering/closest-iu_Minus-Plus_union$window\_$sample\-alone_nonoverlap$overlap\_narrow_peaks.narrowPeak |\
	awk 'BEGIN { OFS="\t" } {if ($21 > 0 ) {print $1 "\t" $3 "\t" $12 "\t" $4 "\t" "." "\t" $6 "\t" $21}   else if ($21 == 0 && $2 < $12 ) {print $1 "\t" $12 "\t" $3 "\t" $4 "\t" "." "\t" $6 "\t" $21} }' \
	> $output_dir/peak_filtering/Minus_oris_union$window\_$sample\-alone_nonoverlap$overlap\_narrow_peaks.bed &&
### (d) selection of the Origins as the region between a plus upstream and a following Minus peaks from the closest output table:
### to be considered: the closest id option is not ignoring downstream peaks that are overlapping. Those are filtered out by distance ($21) equal  0 and start coordinate of a ($2) smaller than start coordinate of b ($12).	
	cat $output_dir/peak_filtering/closest-id_Plus-Minus_union$window\_$sample\-alone_nonoverlap$overlap\_narrow_peaks.narrowPeak |\
	awk 'BEGIN { OFS="\t" } {if ($21 < 0 ) {print $1 "\t" $13 "\t" $2 "\t" $4 "\t" "." "\t" "-" "\t" $21}   else if ($21 == 0 && $2 > $12 ) {print $1 "\t" $2 "\t" $13 "\t" $4 "\t" "." "\t" "-" "\t" $21} }' \
	> $output_dir/peak_filtering/Plus_oris_union$window\_$sample\-alone_nonoverlap$overlap\_narrow_peaks.bed &&
### (e) union of the two Ori tables, most of the Ori regions are identical, so the union is done by sorting the regions by start coordinate and removing duplicates and region with tha same sart or end coordinate keeping the smaller one
	echo $output_dir/peak_filtering/Plus_oris_union$window\_$sample\-alone_nonoverlap$overlap\_narrow_peaks.bed && \
	echo $output_dir/peak_filtering/Minus_oris_union$window\_$sample\-alone_nonoverlap$overlap\_narrow_peaks.bed && \
	cat $output_dir/peak_filtering/Plus_oris_union$window\_$sample\-alone_nonoverlap$overlap\_narrow_peaks.bed $output_dir/peak_filtering/Minus_oris_union$window\_$sample\-alone_nonoverlap$overlap\_narrow_peaks.bed| sort -k1,1 -k2,2n | awk '!a[$1 $2 $3]++' | sort -rk2 | awk '!seen[$1 $3]++' | sort -k3 | awk '!seen[$1 $2]++' | sort -k1,1 -k2,2n> $output_dir/oris/ORI_$sample\-alone_union$window\_nonoverlap$overlap\_narrow_peaks_withStrand.bed &&  
	cat $output_dir/oris/ORI_$sample\-alone_union$window\_nonoverlap$overlap\_narrow_peaks_withStrand.bed | \
	awk ' {print  $1 "\t" $2  "\t" $3  "\t" $4  "\t" $5  "\t" "." "\t" $7}'  \
	> $output_dir/oris/ORI_$sample\-alone_union$window\_nonoverlap$overlap\_narrow_peaks.bed \
;done ;done ;done

printf "ori selection completed\n"