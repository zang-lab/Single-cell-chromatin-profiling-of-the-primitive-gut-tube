#!/bash/bin

if [ $# -lt 1 ];then
    echo "Need 1 parameters! <name> "
    exit
fi

N=$1
cd /scratch/hz9fq/toronto_project/atac_lung_panc/
mkdir ${N}
cd ${N}
####### mapping, formating
bowtie2 -p 8 -X 2000 -3 76 -x /nv/vol190/zanglab/sh8tv/Data/Mapping_index/bowtie2/mm10_clean -1 ../fastq/${N}_1.fastq -2 ../fastq/${N}_2.fastq -S ${N}.sam 2>&1 >>/dev/null | tee ${N}_bowtie2PE.out
module load samtools
#transfer sam to bam
samtools view -bS ${N}.sam > ${N}_raw.bam
#get sam header
samtools view -H ${N}_raw.bam > ${N}_filtered.sam
#get high quality reads
samtools view -f 0x2 ${N}_raw.bam | awk 'NR % 2 == 1{mapq=$5;forward=$0} NR % 2 == 0{if($5>=30 && mapq>=30 && substr($3,1,3)=="chr" && $7=="=" ) print forward"\n"$0}' >> ${N}_filtered.sam
samtools view -bS ${N}_filtered.sam > ${N}.bam
wait;
rm ${N}.sam ${N}_filtered.sam

module load gcc/7.1.0 bedtools/2.26.0
bamToBed -i ${N}.bam -bedpe | grep -v chrM | awk '{OFS="\t"; print $1,$2,$6}'|sort -k 1,1 -k 2,2g -k 3,3g |uniq  > ${N}_PEuniq.bed
