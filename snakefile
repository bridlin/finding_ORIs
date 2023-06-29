configfile: "config.yaml"

rule bam_strand-seperation :
    input:
        "directory/{sample}aln-pe_TB_trimmed_unique_sorted_reheadered_dups-removed.bam"
    output: 
        F2_output="directory/{sample}_F2.bam"
        R1_output="directory/{sample}_R1.bam"
        R2_output="directory/{sample}_R2.bam"
        F1_output="directory/{sample}_F1.bam"
    params:
        F2="-b --require-flags 128 --exclude-flags 16"
        R1="-b --require-flags 80"
        R2="-b --require-flags 144"
        F1="-b --require-flags 64 --exclude-flags 16"
    shell:
        """
        samtools view {params.F2}  {input} > {output.F2_output} &&
        samtools view {params.R1}  {input} > {output.R1_output} &&
        samtools view {params.R2}  {input} > {output.R2_output} &&
        samtools view {params.F1}  {input} > {output.F1_output} 
        """

rule F2R1_merge :
    input:
        F2_input="directory/{sample}_F2.bam"
        R1_input="directory/{sample}_R1.bam"
        R2_input="directory/{sample}_R2.bam"
        F1_input="directory/{sample}_F1.bam"
    output:
        F2R1_output="directory/F2R1_{sample}.bam"
        F1R2_output="directory/F1R2_{sample}.bam"
    params: "-f"
    shell:
        """
        samtools merge {params} {output.F2R1_output} {input.F2_input} {input.R1_input}
        samtools merge {params} {output.F1R2_output} {input.F1_input} {input.R2_input}
        """
