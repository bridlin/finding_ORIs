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




source config_finding-ORIs.txt

# seperation of the minus and plus strand read pairs:

for x in "${input_list[@]}"; do
	samtools view -b -f 128 -F 16 $directory/$x\_$file_prefix\.bam > $directory/$x\_F2.bam &&
	samtools view -b -f 80 $directory/$x\_$file_prefix\.bam > $directory/$x\_R1.bam &&
	samtools merge -f $directory/$x\_F2R1_$file_prefix\.bam $directory/$x\_F2.bam $directory/$x\_R1.bam &&
	samtools index $directory/$x\_F2R1_$file_prefix\.bam &&

	samtools view -b -f 144 $directory/$x\_$file_prefix\.bam > $directory/$x\_R2.bam &&
	samtools view -b -f 64 -F 16 $directory/$x\_$file_prefix\.bam > $directory/$x\_F1.bam &&
	samtools merge -f $directory/$x\_F1R2_$file_prefix\.bam $directory/$x\_R2.bam $directory/$x\_F1.bam &&
	samtools index $directory/$x\_F1R2_$file_prefix\.bam \
; done

#PEAK CALLING ALONE:

for x in "${input_list[@]}"; do
	macs2 callpeak  --bdg  -t $directory/$x\_F2R1_$file_prefix\.bam   -f BAMPE -n $x\-alone_Minus_bowtie2_trimmed_uniq_dupsre_narrow_p005   --outdir $output_dir/ -p 5e-2 -s 130 -m 10 30 --gsize 2.5e7 &&
	macs2 callpeak  --bdg  -t $directory/$x\_F1R2_$file_prefix\.bam  -f BAMPE -n $x\-alone_Plus_bowtie2_trimmed_uniq_dupsre_narrow_p005  --outdir $output_dir/ -p 5e-2 -s 130 -m 10 30 --gsize 2.5e7 \
; done


#PEAK SELECTION :
#1) Sorting out overlapping peaks:

for i in 4 5 6 9 ; do for x in "${input_list[@]}" ; do
	bedtools intersect -wa  -e  -v -f 0.$i  -F 0.$i  -a $output_dir/$x\-alone_Minus_bowtie2_trimmed_uniq_dupsre_narrow_p005_peaks.narrowPeak -b $output_dir/$x\-alone_Plus_bowtie2_trimmed_uniq_dupsre_narrow_p005_peaks.narrowPeak > $output_dir/$x\-alone_Minus_nonoverlap$i\0_narrow_p005_peaks.narrowPeak &&
	bedtools intersect -wa  -e  -v -f 0.$i  -F 0.$i  -a $output_dir/$x\-alone_Plus_bowtie2_trimmed_uniq_dupsre_narrow_p005_peaks.narrowPeak -b $output_dir/$x\-alone_Minus_bowtie2_trimmed_uniq_dupsre_narrow_p005_peaks.narrowPeak > $output_dir/$x\-alone_Plus_nonoverlap$i\0_narrow_p005_peaks.narrowPeak \
;done ;done


#2) Selecting peaks in a window:

for i in 500 1000 1500 2000 2500; do for x in "${input_list[@]}" ; do for y in 4 5 6 9 ; do
	bedtools window -w $i -a $output_dir/$x\-alone_Minus_nonoverlap$y\0_narrow_p005_peaks.narrowPeak -b $output_dir/$x\-alone_Plus_nonoverlap$y\0_narrow_p005_peaks.narrowPeak > $output_dir/union$i\_$x\-alone_Minus-Plus-nonoverlap$y\0_narrow_p005_peaks.narrowPeak  &&
	cat $output_dir/union$i\_$x\-alone_Minus-Plus-nonoverlap$y\0_narrow_p005_peaks.narrowPeak | cut -f 1,2,3,4,5,6,7,8,9,10 | sort -k1,1 -k2,2n  | uniq | awk 'BEGIN { OFS="\t" } {if ($4 ~ /_Minus_/) print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" "-" "\t" $7 "\t"$8 "\t" $9 "\t" $10 ;else print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" "+" "\t" $7 "\t"$8 "\t" $9 "\t" $10}' > $output_dir/Minus-union$i\_$x\-alone_nonoverlap$y\0_narrow_p005_peaks.narrowPeak &&
	cat $output_dir/union$i\_$x\-alone_Minus-Plus-nonoverlap$y\0_narrow_p005_peaks.narrowPeak | cut -f 11,12,13,14,15,16,17,18,19,20 | sort -k1,1 -k2,2n  | uniq | awk 'BEGIN { OFS="\t" } {if ($4 ~ /_Minus_/) print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" "-" "\t" $7 "\t"$8 "\t" $9 "\t" $10 ;else print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" "+" "\t" $7 "\t"$8 "\t" $9 "\t" $10}'  > $output_dir/Plus-union$i\_$x\-alone_nonoverlap$y\0_narrow_p005_peaks.narrowPeak &&
	cat $output_dir/Minus-union$i\_$x\-alone_nonoverlap$y\0_narrow_p005_peaks.narrowPeak  $output_dir/Plus-union$i\_$x\-alone_nonoverlap$y\0_narrow_p005_peaks.narrowPeak | sort -k1,1 -k2,2n > $output_dir/catMinus-Plus_union$i\_$x\-alone_nonoverlap$y\0_narrow_p005_peaks.narrowPeak  \
;done ;done ;done

#ORI SELECTION
#Ori-selection strand-seperated Minus-Plus direction 

for i in 500 1000 1500 2000 2500 ; do for x in "${input_list[@]}"; do for y in 4 5 6 9 ; do
	bedtools closest -iu -D ref -t all -a $output_dir/Minus-union$i\_$x\-alone_nonoverlap$y\0_narrow_p005_peaks.narrowPeak  -b $output_dir/Plus-union$i\_$x\-alone_nonoverlap$y\0_narrow_p005_peaks.narrowPeak | awk -v dist="$i"  'BEGIN { OFS="\t" } {if ($12 != -1 && $21 < dist ) {print $0} }'> $output_dir/closest-iu_Minus-Plus_union$i\_$x\-alone_nonoverlap$y\0_narrow_p005_peaks.narrowPeak &&
	bedtools closest -id -D ref -t all -a $output_dir/Plus-union$i\_$x\-alone_nonoverlap$y\0_narrow_p005_peaks.narrowPeak -b $output_dir/Minus-union$i\_$x\-alone_nonoverlap$y\0_narrow_p005_peaks.narrowPeak  | awk -v dist="$i"  'BEGIN { OFS="\t" } {if ($12 != -1 && $21 > -dist ) {print $0} }'  > $output_dir/closest-id_Plus-Minus_union$i\_$x\-alone_nonoverlap$y\0_narrow_p005_peaks.narrowPeak &&
	cat $output_dir/closest-iu_Minus-Plus_union$i\_$x\-alone_nonoverlap$y\0_narrow_p005_peaks.narrowPeak | awk 'BEGIN { OFS="\t" } {if ($21 > 0 ) {print $1 "\t" $3 "\t" $12 "\t" $4 "\t" "." "\t" $6 "\t" $21}   else if ($21 == 0 && $2 < $12 ) {print $1 "\t" $12 "\t" $3 "\t" $4 "\t" "." "\t" $6 "\t" $21} }' > $output_dir/Minus_oris_union$i\_$x\-alone_nonoverlap$y\0_narrow_p005_peaks.bed &&
	cat $output_dir/closest-id_Plus-Minus_union$i\_$x\-alone_nonoverlap$y\0_narrow_p005_peaks.narrowPeak | awk 'BEGIN { OFS="\t" } {if ($21 < 0 ) {print $1 "\t" $13 "\t" $2 "\t" $4 "\t" "." "\t" "-" "\t" $21}  else if ($21 == 0 && $2 > $12 ) {print $1 "\t" $2 "\t" $13 "\t" $4 "\t" "." "\t" "-" "\t" $21} }' > $output_dir/Plus_oris_union$i\_$x\-alone_nonoverlap$y\0_narrow_p005_peaks.bed &&
	cat $output_dir/Plus_oris_union$i\_$x\-alone_nonoverlap$y\0_narrow_p005_peaks.bed  $output_dir/Minus_oris_union$i\_$x\-alone_nonoverlap$y\0_narrow_p005_peaks.bed | sort -k1,1 -k2,2n |  awk '!a[$1 $2 $3]++' | sort -rk2 | awk '!seen[$1 $3]++'|  sort -k3  | awk '!seen[$1 $2]++' |   sort -k1,1 -k2,2n >   $output_dir/ORI_$x\-alone_union$i\_nonoverlap$y\0_narrow_p005_peaks.bed   \
;done ;done ;done
