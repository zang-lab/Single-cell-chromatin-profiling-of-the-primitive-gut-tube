#1 get Sox2 motif center
awk -v OFS='\t' '{center=int(($2+$3)/2);print $1,center,center+1,$4,$5,$6}' /nv/vol190/zanglab/shared/Motif/sites/hg38_fimo_hocomoco/mappable_results/SOX2_map.bed > SOX2_map_center.bed

#2 get Sox2 motifs overlapped with pan-cancer peak set
#both hg38
module load bedtools/2.26.0
peak1=COAD_peakswithFC_bed6_top20k.bed
peak2=SOX2_map_center.bed

bedtools intersect -u -wa -a $peak2 -b $peak1 > Sox2motif_with_20k_peak.bed

#3 convert bed to bigwiggle
bedtools genomecov -bga -i Sox2motif_with_20k_peak.bed -g /nv/vol190/zanglab/hz9fq/annotation/Human/hg38.chrom.sizes > Sox2motif_with_20k_peak.bedgraph

sort -k1,1 -k2,2n Sox2motif_with_20k_peak.bedgraph > Sox2motif_with_20k_peak_sorted.bedgraph
bedGraphToBigWig Sox2motif_with_20k_peak_sorted.bedgraph /nv/vol190/zanglab/hz9fq/annotation/Human/hg38.chrom.sizes Sox2motif_with_20k_peak.bw

#4 rp
python RPCalc.py -b /sfs/lustre/bahamut/scratch/hz9fq/toronto_project/TCGA/ATAC-seq/RP/Sox2motif_with_20k_peak.bw -n Sox2motif_with_20k_RP.bed -g hg38 --cs /nv/vol190/zanglab/hz9fq/annotation/Human/hg38.chrom.sizes --tss /nv/vol190/zanglab/hz9fq/annotation/Human/hg38_gene_annotation_geneID_LenOrder_TSS_forRP.bed

