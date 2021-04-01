#!/bash/bin

if [ $# -lt 1 ];then
    echo "Need 1 parameters! <name> "
    exit
fi

N=$1
cd /scratch/hz9fq/toronto_project/E13bulkATACpeak
mkdir ${N}
cd ${N}

module load macs2
macs2 callpeak -f BED --SPMR -g mm -B -q 0.01 --keep-dup all --nomodel --extsize 100 -t ../E13bedfiles/readfiles/${N}_PEuniq.bed -n ${N}_m2 --outdir .
mv ${N}_m2_treat_pileup.bdg ${N}_m2.bdg

module load bedtools/2.26.0
sort -k1,1 -k2,2n ${N}_m2.bdg > ${N}.bdg
bdg2bw ${N}.bdg  /scratch/hz9fq/data_archive/annotation/mm10.chrom.sizes

