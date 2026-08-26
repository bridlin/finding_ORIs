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


source scripts/finding_ORIs/config.txt

mkdir $alined_reads_dir
mkdir ${alined_reads_dir}/fastqc
mkdir ${alined_reads_dir}/reports

printf "read qc and alignment\n"
printf "samples are ${input_list[@]}\n"
printf "read directory is ${fastq_directory}\n"
printf "aligned read directory is ${alined_reads_dir}\n"
printf "readname postfix is  ${readname_postfix_R1} and ${readname_postfix_R2} \n"
printf "genome used is ${genome} with prefix ${genome_prefix}\n"

for sample in "${input_list[@]}"; do
fastqc ${fastq_directory}/${sample}${readname_postfix_R1} --outdir ${alined_reads_dir}/fastqc &&
fastqc ${fastq_directory}/${sample}${readname_postfix_R2} --outdir ${alined_reads_dir}/fastqc &&
cutadapt -u -10 -u 10   -U -10 -U 10  \
    -o ${fastq_directory}/${sample}R1_5tailtrimmed.fastq.gz  \
    -p ${fastq_directory}/${sample}R2_5tailtrimmed.fastq.gz  \
    ${fastq_directory}/${sample}${readname_postfix_R1} ${fastq_directory}/${sample}${readname_postfix_R2} &&
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA  -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT  \
    -o ${fastq_directory}/${sample}R1_5-3trimmed.fastq.gz \
    -p ${fastq_directory}/${sample}R2_5-3trimmed.fastq.gz  \
    ${fastq_directory}/${sample}R1_5tailtrimmed.fastq.gz  ${fastq_directory}/${sample}R2_5tailtrimmed.fastq.gz \
    --minimum-length 30 > ${alined_reads_dir}//reports/${sample}\_cutadapt_report.txt &&
trimmomatic PE \
    -threads 4 \
    -trimlog ${alined_reads_dir}/${sample}trim \
    ${fastq_directory}/${sample}R1_5-3trimmed.fastq.gz ${fastq_directory}/${sample}R2_5-3trimmed.fastq.gz \
    ${fastq_directory}/${sample}R1_5-3trimmed_q20.fastq.gz   ${fastq_directory}/${sample}R1_5-3trimmed_q20_un.fastq.gz \
    ${fastq_directory}/${sample}R2_5-3trimmed_q20.fastq.gz   ${fastq_directory}/${sample}R2_5-3trimmed_q20_un.fastq.gz \
    SLIDINGWINDOW:4:20 MINLEN:40 &&
fastqc ${fastq_directory}/${sample}R1_5-3trimmed_q20.fastq.gz --outdir ${alined_reads_dir}/fastqc &&
fastqc ${fastq_directory}/${sample}R2_5-3trimmed_q20.fastq.gz --outdir ${alined_reads_dir}/fastqc && 
bowtie2 \
    -k1 \
    -x ${genome} \
    -1 ${fastq_directory}/${sample}R1_5-3trimmed_q20.fastq.gz \
    -2 ${fastq_directory}/${sample}R2_5-3trimmed_q20.fastq.gz   \
    -S ${alined_reads_dir}/${sample}aln-pe_${genome_prefix}.sam \
    2> ${alined_reads_dir}/reports/${sample}_bowtie.log &&
samtools view -S -b ${alined_reads_dir}/${sample}aln-pe_${genome_prefix}.sam > ${alined_reads_dir}/${sample}aln-pe_${genome_prefix}.sam.bam &&
samtools sort ${alined_reads_dir}/${sample}aln-pe_${genome_prefix}.sam.bam -o ${alined_reads_dir}/${sample}aln-pe_${genome_prefix}_sorted.bam &&
samtools reheader -c 'grep -v ^@PG' ${alined_reads_dir}/${sample}aln-pe_${genome_prefix}_sorted.bam  > ${alined_reads_dir}/${sample}aln-pe_${genome_prefix}_sorted_reheadered.bam &&
picard CollectInsertSizeMetrics \
    -I ${alined_reads_dir}/${sample}aln-pe_${genome_prefix}_sorted_reheadered.bam \
    -O ${alined_reads_dir}/reports/${sample}aln-pe_${genome_prefix}_sorted_reheadered_insert_size_metrics.txt \
    -H ${alined_reads_dir}/reports/${sample}aln-pe_${genome_prefix}_sorted_reheadered_insert_size_histogram.pdf \
    -M 0.5  &&
picard  MarkDuplicates \
    --REMOVE_DUPLICATES true \
    -I ${alined_reads_dir}/${sample}aln-pe_${genome_prefix}_sorted_reheadered.bam \
    -O ${alined_reads_dir}/${sample}aln-pe_${genome_prefix}_sorted_reheadered_dups-removed.bam \
    -M ${alined_reads_dir}/reports/${sample}aln-pe_${genome_prefix}_marked_dup_metrics.txt &&
samtools index ${alined_reads_dir}/${sample}aln-pe_${genome_prefix}_sorted_reheadered_dups-removed.bam &&
rm -f  ${alined_reads_dir}/${sample}aln-pe_${genome_prefix}_sorted_reheadered.bam &&
rm -f  ${alined_reads_dir}/${sample}aln-pe_${genome_prefix}.sam &&
rm -f  ${alined_reads_dir}/${sample}aln-pe_${genome_prefix}.sam.bam ;done

multiqc   \
    $alined_reads_dir \
    ${alined_reads_dir}/reports \
    ${alined_reads_dir}/fastqc  \
    --outdir ${alined_reads_dir} 