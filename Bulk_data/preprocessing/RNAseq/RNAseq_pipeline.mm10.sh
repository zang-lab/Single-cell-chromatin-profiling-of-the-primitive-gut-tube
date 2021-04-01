#!/bin/bash

if [ $# -lt 1 ];then
    echo "Need 1 parameters! <name> "
    exit
fi


R=$1

#create path
main=/scratch/hz9fq/toronto_project/YUN13891_and_YUN5339_RNA/
cd ${main}
mkdir ${R}
cd ${R}

##fastq-dump
#fastq-dump -I --split-files /scratch/hz9fq/sra/${R}.sra -O /scratch/hz9fq/RNAseq/${R}

### trim
module load gcc/7.1.0
module load cutadapt
#trim_galore --phred33 --clip_R1 30 --clip_R2 30 --three_prime_clip_R1 30 --three_prime_clip_R2 30 --fastqc --illumina --paired ${R}_1.fastq ${R}_2.fastq

### hisat2
module load hisat2/2.1.0
module load samtools
hisat2 -p 8 --dta -x /nv/vol190/zanglab/sh8tv/Data/Mapping_index/Hisat2/mm10 -1 ${main}fastq/${R}_1.fastq -2 ${main}fastq/${R}_2.fastq -S ${R}_raw_stringtie.sam  2>&1 >>/dev/null |tee -a ${R}_hisat2PE.out
#transfer sam to bam
samtools view -bS ${R}_raw_stringtie.sam > ${R}_raw_stringtie.bam

#get sam header
samtools view -H ${R}_raw_stringtie.bam > ${R}_filtered_stringtie.sam 
#get filtered sam content
samtools view -f 0x2 ${R}_raw_stringtie.bam | awk 'NR % 2 == 1{mapq=$5;forward=$0} NR % 2 == 0{if($5>=30 && mapq>=30 && substr($3,1,3)=="chr" && $7=="=" ) print forward"\n"$0}' >> ${R}_filtered_stringtie.sam 
samtools view -bS ${R}_filtered_stringtie.sam > ${R}_stringtie.bam
#samtools sort -@ 8 -o ${R}_stringtie.bam ${R}_filtered_stringtie.sam 

rm ${R}_raw_stringtie.sam
rm ${R}_filtered_stringtie.sam


samtools sort -@ 8 -o ${R}_sort1.bam ${R}_stringtie.bam
stringtie -e -B -p 8 -A ${R}_exp.txt -G /nv/vol190/zanglab/hz9fq/annotation/Mouse/mm10_refgenes.gtf -o ${R}_stringtie.gtf ${R}_sort1.bam

samtools sort -@ 8 -n -o ${R}_sort2.bam ${R}_stringtie.bam
htseq-count -f bam -r name ${R}_sort2.bam -i gene_id /nv/vol190/zanglab/hz9fq/annotation/Mouse/mm10_refgenes.gtf > ${R}_htseq_gene.txt


