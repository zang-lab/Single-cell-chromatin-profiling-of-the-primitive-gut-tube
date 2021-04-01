#!/bin/bash

if [ $# -lt 1 ];then
    echo "Need 1 parameters! <name> "
    exit
fi


R=$1
S=mm10

#create path
main=/scratch/hz9fq/toronto_project/public_ChIPseq2/
cd ${main}
mkdir ${R}
cd ${R}
mkdir macs14_p1e5
cd macs14_p1e5

##fastq-dump
#fastq-dump -I --split-files /scratch/hz9fq/sra/${R}.sra  -O /scratch/hz9fq/ChIPseq/${R}

### trim
#module load gcc/7.1.0
#module load cutadapt
#module load fastqc
#trim_galore --phred33 --fastqc --illumina ${R}.fastq

#module load gcc/7.1.0
#module load bowtie2/2.2.9
#bowtie2 -p 8 -x /nv/vol190/zanglab/hz9fq/index/mm10 -U ${main}fastq/${R}.fastq -S ${R}.sam 2>&1 >>/dev/null | tee ${R}_bowtie2SE.out
#bowtie2 -p 5 -x /nv/vol190/zanglab/sh8tv/Data/Mapping_index/bowtie2/${S}_clean -U /nv/vol190/zanglab/sh8tv/Project/tmp/MouseDev/Data/fastq/H3K4me3/${R}_p1.fastq -S ${R}.sam   2>&1 >>/dev/null |tee -a ${R}_bowtie2SE.out
#module load samtools
#samtools view -bS ${R}.sam > ${R}_raw.bam
#samtools view -q 30 -bS ${R}.sam > ${R}.bam

#rm ${R}.sam

#macs2 callpeak --SPMR -g mm -B --tsize 50 -q 0.1  --keep-dup 1 -f BAM  --nomodel --extsize 146 -t ${R}.bam -n ${R}_m2 --outdir macs2_q01
export PATH=$PATH:/home/sh8tv/bin/
export PYTHONPATH=$PYTHONPATH:/home/sh8tv/lib64/python2.7/site-packages/
macs14 -t ../${R}.bam -n ${R}_m14 -f BAM -g mm -p 1e-5 -B --single-profile --nomodel --shiftsize 73 --keep-dup 1 

#mv ${R}_m14_treat_pileup.bdg ${R}.bdg
#module load bedtools/2.26.0
#bdg2bw ${R}.bdg /scratch/hz9fq/data/mm10.chrom.sizes

