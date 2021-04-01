#!/bash/bin

if [ $# -lt 1 ];then
    echo "Need 1 parameters! <name> "
    exit
fi

N=$1
main=/scratch/hz9fq/toronto_project/E16bulkATACpeak
cd ${main}
mkdir ${N}
cd ${N}

module load macs2
#cmd <- sprintf("callpeak -g %s --name %s --treatment %s --outdir %s --format BED --call-summits --keep-dup all %s", 
#        genomeSize, basename(bedName), bedFile, dirname(bedName), additionalParams)
macs2 callpeak -f BED --SPMR -g mm -B -q 0.01 --keep-dup all --nomodel --extsize 100 -t ../readfiles/${N}_plus9_PEuniq.bed -n ${N}_m2 --outdir .
mv ${N}_m2_treat_pileup.bdg ${N}_m2.bdg

module load bedtools/2.26.0
sort -k1,1 -k2,2n ${N}_m2.bdg > ${N}.bdg
bdg2bw ${N}.bdg  /scratch/hz9fq/data_archive/annotation/mm10.chrom.sizes

