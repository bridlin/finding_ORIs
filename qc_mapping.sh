#!/bin/bash
#
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err
#SBATCH --mail-type END
#SBATCH --mail-user b-barckmann@chu-montpellier.fr
#
#
#SBATCH --partition fast
#SBATCH --cpus-per-task 4
#SBATCH --mem  128GB


module load cutadapt/4.0
module load trimmomatic/0.39
module load fastqc/0.11.9
module load bowtie2/2.4.1
module load samtools/1.13
module load picard/2.23.5
module load multiqc/1.13


source scripts/finding_ORIs/config_mapping.txt

mkdir $output_dir
mkdir $output_dir/fastqc
mkdir $output_dir/reports


for sample in "${input_list[@]}"; do
fastqc $fastq_directory/$sample$readname_postfix_R1 --outdir $output_dir/fastqc &&
fastqc $fastq_directory/$sample$readname_postfix_R2 --outdir $output_dir/fastqc &&
cutadapt -u -10 -u 10   -U -10 -U 10  \
    -o $fastq_directory/$sample\R1_5tailtrimmed.fastq.gz  \
    -p $fastq_directory/$sample\R2_5tailtrimmed.fastq.gz  \
    $fastq_directory/$sample$readname_postfix_R1 $fastq_directory/$sample$readname_postfix_R2 &&
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA  -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT  \
    -o $fastq_directory/$sample\R1_5-3trimmed.fastq.gz \
    -p $fastq_directory/$sample\R2_5-3trimmed.fastq.gz  \
    $fastq_directory/$sample\R1_5tailtrimmed.fastq.gz  $fastq_directory/$sample\R2_5tailtrimmed.fastq.gz \
    --minimum-length 30 > $output_dir//reports/$sample\_cutadapt_report.txt &&
trimmomatic PE \
    -threads 4 \
    -trimlog $output_dir/$sample\trim \
    $fastq_directory/$sample\R1_5-3trimmed.fastq.gz $fastq_directory/$sample\R2_5-3trimmed.fastq.gz \
    $fastq_directory/$sample\R1_5-3trimmed_q20.fastq.gz   $fastq_directory/$sample\R1_5-3trimmed_q20_un.fastq.gz \
    $fastq_directory/$sample\R2_5-3trimmed_q20.fastq.gz   $fastq_directory/$sample\R2_5-3trimmed_q20_un.fastq.gz \
    SLIDINGWINDOW:4:20 MINLEN:40 &&
fastqc $fastq_directory/$sample\R1_5-3trimmed_q20.fastq.gz --outdir $output_dir/fastqc &&
fastqc $fastq_directory/$sample\R2_5-3trimmed_q20.fastq.gz --outdir $output_dir/fastqc && 
bowtie2 \
    -k1 \
    -x $genome \
    -1 $fastq_directory/$sample\R1_5-3trimmed_q20.fastq.gz \
    -2 $fastq_directory/$sample\R2_5-3trimmed_q20.fastq.gz   \
    -S $output_dir/$sample\aln-pe_$genome_prefix\.sam \
    2> $output_dir/reports/$sample\_bowtie.log &&
samtools view -S -b $output_dir/$sample\aln-pe_$genome_prefix\.sam > $output_dir/$sample\aln-pe_$genome_prefix\.sam.bam &&
samtools sort $output_dir/$sample\aln-pe_$genome_prefix\.sam.bam -o $output_dir/$sample\aln-pe_$genome_prefix\_sorted.bam &&
samtools reheader -c 'grep -v ^@PG' $output_dir/$sample\aln-pe_$genome_prefix\_sorted.bam  > $output_dir/$sample\aln-pe_$genome_prefix\_sorted_reheadered.bam &&
picard CollectInsertSizeMetrics \
    -I $output_dir/$sample\aln-pe_$genome_prefix\_sorted_reheadered.bam \
    -O $output_dir/reports/$sample\aln-pe_$genome_prefix\_sorted_reheadered_insert_size_metrics.txt \
    -H $output_dir/reports/$sample\aln-pe_$genome_prefix\_sorted_reheadered_insert_size_histogram.pdf \
    -M 0.5  &&
picard  MarkDuplicates \
    --REMOVE_DUPLICATES true \
    -I $output_dir/$sample\aln-pe_$genome_prefix\_sorted_reheadered.bam \
    -O $output_dir/$sample\aln-pe_$genome_prefix\_sorted_reheadered_dups-removed.bam \
    -M $output_dir/reports/$sample\aln-pe_$genome_prefix\_marked_dup_metrics.txt &&
samtools index $output_dir/$sample\aln-pe_$genome_prefix\_sorted_reheadered_dups-removed.bam &&
rm -f  $output_dir/$sample\aln-pe_$genome_prefix\_sorted_reheadered.bam &&
rm -f  $output_dir/$sample\aln-pe_$genome_prefix\.sam &&
rm -f  $output_dir/$sample\aln-pe_$genome_prefix\.sam.bam ;done

multiqc   \
    $output_dir \
    $output_dir/reports \
    $output_dir/fastqc  \
    --outdir $output_dir 