# FOOTPRINTING ###################
################################################

#hacer lo de diferenciar la orientacion del motivo
neurod2_peaks$dinucleotide_oriented<-substr(neurod2_peaks$ebox_seq, 3, 8) 
for (motif in unique(neurod2_peaks$dinucleotide_oriented)) {
  motifs<-neurod2_peaks %>% dplyr::filter(dinucleotide_oriented==motif)
  motifs$ebox_end<-as.numeric(motifs$ebox_start)+5
  motifs$ebox_start<-motifs$ebox_end-1
  motifs<-motifs %>% dplyr::select(chr,ebox_start,ebox_end)
  write.table(motifs,paste0("~/openchr/noack/footprinting/",motif,".bed"),sep="\t",col.names = F,row.names = F, quote = F, )
}

neurog2_peaks_motifs<-readRDS("~/neurog2/chipseq/neurog2_peaks_motifs.Rds")
neurog2_peaks_motifs$dinucleotide_oriented<-substr(neurog2_peaks_motifs$ebox_seq, 3, 8) 

library(dplyr)
for (motif in unique(neurog2_peaks_motifs$dinucleotide_oriented)) {
  motifs<-neurog2_peaks_motifs %>% dplyr::filter(dinucleotide_oriented==motif)
  motifs$ebox_end<-as.numeric(motifs$ebox_start)+5
  motifs$ebox_start<-motifs$ebox_end-1
  motifs<-motifs %>% dplyr::select(chr,ebox_start,ebox_end)
  write.table(motifs,paste0("~/openchr/noack/footprinting/neurog2_",motif,".bed"),sep="\t",col.names = F,row.names = F, quote = F, )
}

system ("
        cd ~/openchr/noack/footprinting
        for file in neurog2*CA**TG.bed
        do
        bedtools intersect -wa -a $file -b ~/openchr/noack/neural_atac_peaks_neurod2_intersect.bed > $file.filtered.tmp.bed
        done
        ")





#  heatmaps de los footprints
###
system ("
        cd ~/p/openchr/noack/footprinting
        sbatch  --wrap="bedtools intersect -wa -a $file -b ~/openchr/noack/neural_atac_peaks_neurod2_intersect.bed > $file.filtered.tmp.bed; "
        module load Miniconda3/4.12.0
        sbatch  --partition=short -n 50 --wrap="computeMatrix reference-point -p 50 --missingDataAsZero --binSize 1 -S ~/p/openchr/noack/cluster_bams/merged_clusters.bam.hint.corrected.bw \
        -R CA**.filtered.tmp.bed -out neurod2foot.matrix.tmp -a 20 -b 20 ; plotHeatmap --matrixFile neurod2foot.matrix.tmp --outFileName neurod2foot.matrix.tmp.pdf \
        --plotFileFormat pdf --perGroup --heatmapHeight 5  --colorMap viridis"
       #para neurog2
       sbatch  --partition=short -n 50 --wrap="computeMatrix reference-point -p 50 --missingDataAsZero --binSize 1 -S ~/p/openchr/noack/cluster_bams/merged_clusters.bam.hint.corrected.bw \
        -R neurog2*CA**.filtered.tmp.bed -out neurog2foot.matrix.tmp -a 20 -b 20 ; plotHeatmap --matrixFile neurog2foot.matrix.tmp --outFileName neurog2foot.matrix.tmp.pdf \
        --plotFileFormat pdf --perGroup --heatmapHeight 5  --colorMap viridis"
")


system("
  import_shiva p/openchr/noack/footprinting/\*pdf .
")




#hacer lo mismo pero con el genomicfeatures
system("
       import_berlin projects/openchr/noack/cluster_bams/merged_clusters.bam.hint.corrected.bw openchr/noack/cluster_bams
       import_berlin projects/openchr/noack/footprinting/CA**TG.bed.filtered.tmp.bed openchr/noack/footprinting
")

library(GenomicFeatures)
library(seqplots)
tracks <-'~/openchr/noack/cluster_bams/merged_clusters.bam.hint.corrected.bw'


plotmeth<-function(meth_matrix,main) {
  print(plotAverage(plotset=meth_matrix, xlim = c(-20,20),
                    ylim = c(0,120), main = main, xlab = "", ylab = '',
                    plotScale = "linear", type = "full", error.estimates = T, yaxs='i',xaxs='i',
                    legend = TRUE, legend_ext = FALSE, legend_pos = 'bottomright', labels = "",
                    legend_ext_pos = "topleft", cex.axis = 14, cex.lab = 16,xaxt='n',cex.legend = 12, pointsize = 12) )
}


## neurod2 motifs
for (motif in c("CAAATG","CAATTG","CACGTG","CAGCTG","CATATG","CATTTG","CAACTG","CACATG","CACTTG","CAGGTG",
                "CATCTG","CAAGTG","CACCTG","CAGATG","CAGTTG","CATGTG") ) {
  print(paste0("~/openchr/noack/footprinting/",motif,".bed.filtered.tmp.bed"))
  jj <- getPlotSetArray(tracks=tracks,features=paste0("~/openchr/noack/footprinting/",motif,".bed.filtered.tmp.bed"),
                        refgenome='mm10',type = 'mf',add_heatmap=F,xmin=20,xmax=20,bin = 1)
  pdf( paste0(motif,"_foot.pdf"), useDingbats = F, width = 10, height = 10)
  plotmeth(jj,main=motif)
  dev.off()
}

## neurog2 motifs
for (motif in c("CAAATG","CAATTG","CACGTG","CAGCTG","CATATG","CATTTG","CAACTG","CACATG","CACTTG","CAGGTG",
                "CATCTG","CAAGTG","CACCTG","CAGATG","CAGTTG","CATGTG") ) {
  print(paste0("~/openchr/noack/footprinting/neurog2",motif,".bed.filtered.tmp.bed"))
  jj <- getPlotSetArray(tracks=tracks,features=paste0("~/openchr/noack/footprinting/neurog2_",motif,".bed.filtered.tmp.bed"),
                        refgenome='mm10',type = 'mf',add_heatmap=F,xmin=20,xmax=20,bin = 1)
  pdf( paste0("neurog2",motif,"_foot.pdf"), useDingbats = F, width = 10, height = 10)
  plotmeth(jj,main=motif)
  dev.off()
}



