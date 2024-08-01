#!/bin/bash
# Requires-
# deeptools: https://deeptools.readthedocs.io/en/develop/index.html
# 10x genomics subset-bam tool: https://github.com/10XGenomics/subset-bam
#    - Barcodes files must be 1 barcode per line no other characters
# macs3: https://github.com/macs3-project/MACS/tree/master
# Sinto: https://github.com/timoast/sinto

## Extract reference genes from ensembl
python scr/gene_ref.py

## Reformat barcode files:
python format_barcodes.py

samples=("P1_T" \
		 "P1_N" \
		 "P2_T" \
		 "P2_N" \
		 "P3_T" \
		 "P4_T")

echo 'Starting'
for sample in ${samples[@]}
do 
	echo ${sample}
## Filter and reindex Bams using cell barcodes of corticotrophs identified by GEX analysis
	subset-bam --bam data/SN_MO/${sample}_outs/${sample}_outs_atac_possorted_bam.bam -c results/${sample}/${sample}_cort_barcodes_formatted.csv -o results/ATAC/bams/${sample}_atac_corts.bam
	samtools index results/ATAC/bams/${sample}_atac_corts.bam

## Create Coverage Bigwig for heatmaps : Paired-end: Reads with mates are always extended to match the fragment size defined by the two read mates. Unmated reads, mate reads that map too far apart (>4x fragment length) or even map to different chromosomes are treated like single-end reads. The input of a fragment length value is optional. If no value is specified, it is estimated from the data (mean of the fragment size of all mate reads).
	bamCoverage -b results/ATAC/bams/${sample}_atac_corts.bam -o results/ATAC/bigwigs/${sample}_corts_frags.bw -bs=1 --normalizeUsing RPKM -p=max/2 -e

## Compute Matrix for Peaks
	computeMatrix reference-point --referencePoint center -S results/ATAC/bigwigs/${sample}_corts_frags.bw -R data/SN_MO/${sample}_outs/${sample}_outs_atac_peaks.bed --missingDataAsZero -bs 1 -a 1000 -b 1000 -out results/ATAC/intermediates/${sample}_peak_matrix.gz 

## Plot Heatmap for Peaks
	plotHeatmap --refPointLabel peak  --regionsLabel peaks --xAxisLabel 'peak distance(bp)' -m results/ATAC/intermediates/${sample}_peak_matrix.gz  -out results/ATAC/plots/${sample}_peaks_WF.png --legendLocation none

## Compute matrix for genes
	computeMatrix scale-regions -S results/ATAC/bigwigs/${sample}_corts_frags.bw -R data/genes.bed --missingDataAsZero -bs 1 -a 1000 -b 1000 -out results/ATAC/intermediates/${sample}_gene_matrix.gz

## Plot Heatmap for Genes
	plotHeatmap -m results/ATAC/intermediates/${sample}_gene_matrix.gz -out results/ATAC/plots/${sample}_gene.png --legendLocation none

done

echo 'merging'
## Merge all Bams
samtools merge results/ATAC/bams/P1_T_atac_corts.bam results/ATAC/bams/P1_N_atac_corts.bam results/ATAC/bams/P2_T_atac_corts.bam results/ATAC/bams/P2_N_atac_corts.bam results/ATAC/bams/P3_T_atac_corts.bam results/ATAC/bams/P4_T_atac_corts.bam -o results/ATAC/bams/all_corts.bam
samtools index results/ATAC/bams/all_corts.bam

## Call peaks on merged psuedobulk for differential accessibility analsis
macs3 callpeak -t results/ATAC/bams/all_corts.bam -f BAM -n cort --outdir results/ATAC/macs/

## Annotate peaks with HOMER
annotatePeaks.pl results/ATAC/macs/cort_peaks.narrowPeak hg38 > results/ATAC/macs/cort_peaks_annotated.narrowpeak

## Merge tumor Bams
samtools merge results/ATAC/bams/P1_T_atac_corts.bam results/ATAC/bams/P2_T_atac_corts.bam results/ATAC/bams/P3_T_atac_corts.bam results/ATAC/bams/P4_T_atac_corts.bam -o results/ATAC/bams/tumor_corts.bam
samtools index results/ATAC/bams/tumor_corts.bam

## Merge normal Bams
samtools merge results/ATAC/bams/P1_N_atac_corts.bam results/ATAC/bams/P2_N_atac_corts.bam -o results/ATAC/bams/normal_corts.bam
samtools index results/ATAC/bams/normal_corts.bam

## Create Coverage Bigwig for Circos
bamCoverage -b results/ATAC/bams/tumor_corts.bam -o results/ATAC/bigwigs/tumor_corts_1mb_SF.bw -bs=1000000 --normalizeUsing None -p=max/2 -e --scaleFactor 6.965243e-5
bamCoverage -b results/ATAC/bams/normal_corts.bam -o results/ATAC/bigwigs/normal_corts_1mb_SF.bw -bs=1000000 --normalizeUsing None -p=max/2 -e --scaleFactor 6.1500615e-4


## Run DAC analyis
# python scr/DAC_analysis.py

## TF Motif enrichment with HOMER
findMotifsGenome.pl results/ATAC/homer/peaks_up.homer hg38 results/ATAC/homer/up/
findMotifsGenome.pl results/ATAC/homer/peaks_down.homer hg38 results/ATAC/homer/down/

## Create Heatmaps for differentially accessible peaks and genes
conditions=("tumor" "normal")
for condition in ${conditions[@]}
do
	echo ${condition}
## Create Coverage Bigwig for heatmaps
	bamCoverage -b results/ATAC/bams/${condition}_corts.bam -o results/ATAC/bigwigs/${condition}_corts.bw -bs=1 --normalizeUsing RPKM -p=max/2 -e

## Create Coverage Bigwig for Circos
	bamCoverage -b results/ATAC/bams/${condition}_corts.bam -o results/ATAC/bigwigs/${condition}_corts_1mb.bw -bs=1000000 --normalizeUsing RPKM -p=max/2 -e

## Compute Matrix for Peaks
	computeMatrix reference-point --referencePoint center -S results/ATAC/bigwigs/${condition}_corts.bw -R results/ATAC/macs/cort_peaks.narrowPeak --missingDataAsZero -bs 1 -a 1000 -b 1000 -out results/ATAC/intermediates/${condition}_peak_matrix.gz
	computeMatrix reference-point --referencePoint center -S results/ATAC/intermediates/${condition}_corts.bw -R results/ATAC/differential/cort_peaks_up.bed --missingDataAsZero -bs 1 -a 1000 -b 1000 -out results/ATAC/intermediates/${condition}_up_peak_matrix.gz 
	computeMatrix reference-point --referencePoint center -S results/ATAC/intermediates/${condition}_corts.bw -R results/ATAC/differential/cort_peaks_down.bed --missingDataAsZero -bs 1 -a 1000 -b 1000 -out results/ATAC/intermediates/${condition}_down_peak_matrix.gz
	
## Plot Heatmap for Peaks
	plotHeatmap --refPointLabel peak  --regionsLabel peaks --xAxisLabel 'peak distance(bp)' -m results/ATAC/intermediates/${condition}_peak_matrix.gz -out results/ATAC/plots/${condition}_peaks.png --legendLocation none --yMax 900 --zMax 2500
	plotHeatmap --refPointLabel peak  --regionsLabel peaks --xAxisLabel 'peak distance(bp)' -m results/ATAC/intermediates/${condition}_up_peak_matrix.gz  -out results/ATAC/plots/${condition}_peaks_up.png --legendLocation none
	plotHeatmap --refPointLabel peak  --regionsLabel peaks --xAxisLabel 'peak distance(bp)' -m results/ATAC/intermediates/${condition}_down_peak_matrix.gz -out results/ATAC/plots/${condition}_peaks_down.png --legendLocation none

## Compute matrix for genes
	computeMatrix scale-regions -S results/ATAC/bigwigs/${condition}_corts.bw -R data/genes.bed --missingDataAsZero -bs 1 -a 1000 -b 1000 -out results/ATAC/intermediates/${condition}_gene_matrix.gz

## Plot Heatmap for Genes
	plotHeatmap -m results/ATAC/intermediates/${condition}_gene_matrix.gz -out results/ATAC/plots/${condition}_gene.png --legendLocation none  --yMax 1200 --zMax 3000

done



