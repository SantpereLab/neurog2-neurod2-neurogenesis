# SPACING ######################
###############################################

neurod2_peaks<-readRDS("~/neurod2/neurod2_peaks.Rds")
neurog2_peaks_motifs<-readRDS("~/neurog2/chipseq/neurog2_peaks_motifs.Rds")

library(dplyr)

# Compute which are the proportions of spacing between all E-boxes.
proportion_spacing<-function(motifs) {
  for (chromosome in unique(motifs$chr)) {
    neurod2_peak_chr<-motifs %>% filter(chr==chromosome)
    motif1_locations<-as.numeric(neurod2_peak_chr$ebox_start)+2
    motif2_locations<-as.numeric(neurod2_peak_chr$ebox_start)+2
    #distancia de todos los motivos 1 con los 2.
    distance<-as.vector(outer(motif2_locations, motif1_locations, "-"))
    #restar 6, la longitud del ebox
    distance<-distance[-which(distance==0)]
    distance<-abs(distance)-6
    distance<-distance[which(distance<=20)] #que sea menor o igual a 20 para acotar
    distance<-distance[which(distance>=0)] #quitar los que son concatenados
    if (chromosome=="chr1") {distance_allchr<-distance }
    else {distance_allchr<-c(distance_allchr,distance) }
  }
  proportions_distances<-table(distance_allchr)
  proportions_distances<-proportions_distances/sum(proportions_distances)
  
  #añadir proporciones de 0 por si faltan datos
  proportions_distances_blank<-rep(0,21)
  names(proportions_distances_blank)<-0:20
  missing_distances<-proportions_distances_blank[-which(names(proportions_distances_blank) %in% names(proportions_distances) ) ]
  proportions_distances<-c(proportions_distances,missing_distances)
  proportions_distances<-proportions_distances[as.character(sort( as.numeric( names(proportions_distances) ) ) ) ]
  
  return(proportions_distances)
}
proportion_spacing(neurod2_peaks)


for (i in 1:10) {
  decil<-read.table(paste0("~/openchr/noack/neurod2_correlation_decil_",i,"_summits.txt"))$V5
  neurod2_peaks_decil<-neurod2_peaks %>% filter(peak_id %in% decil)
  proportion_spacing_quantil<-proportion_spacing(neurod2_peaks_decil)
  if (i==1) {
    proportion_spacing_matrix<-proportion_spacing_quantil
  }
  else {
    proportion_spacing_matrix<-rbind(proportion_spacing_matrix,proportion_spacing_quantil)
  }
}

rownames(proportion_spacing_matrix)<-c(1:10)
pdf("~/neurod2peaks_neurod2cor_spacing.pdf", useDingbats = F, height = 3)
pheatmap::pheatmap(proportion_spacing_matrix,
                   cluster_rows = F, cluster_cols = F,  color = rev(hcl.colors(50, "LaJolla")) )
dev.off()

# Proportion of spacings in the ATAC-seq peaks that do not intersect with proneural.
atac_motifs<-readRDS("~/openchr/noack/atac_peaks_motifs.Rds")
atac_overlap_proneural<-read.table("~/openchr/noack/proneural_table.bed")
library(bedr)
motifs<-bedr.join.region(bedr.sort.region(atac_motifs),bedr.sort.region(atac_overlap_proneural) )
atac_motifs_proneuralnt<-motifs %>% filter(V4=="1.proneuralnt")


atac_motifs_proneuralnt<-readRDS("~/p/atac_motifs_proneuralnt.RDS")
library(dplyr)
proportion_spacing_proneuralnt<-proportion_spacing(atac_motifs_proneuralnt)
saveRDS(proportion_spacing_proneuralnt,"~/p/proportion_spacing_proneuralnt.RDS")


proportion_spacing_proneuralnt<-readRDS("~/proportion_spacing_proneuralnt.RDS")
barplot(proportion_spacing_proneuralnt)



# The same but for the background genomic regions as outputted by HOMER.
background_peaks<-read.table("~/neurod2/background_peaks_homer/background.bed")
colnames(background_peaks)<-c("chr","start","end","peak_id")
background_peaks$summit<-"jii"
background_motifs<-find_motifs(background_peaks, peaks_file = "~/neurod2/background_peaks_homer/background.bed")

colnames(background_motifs)<-c("chr","ebox_start","ebox_end","motif_id")
proportion_spacing_homer<-proportion_spacing(background_motifs)
barplot(proportion_spacing_homer)

# A heatmap with everything together.
join_heatmap<-rbind(proportion_spacing_homer,proportion_spacing_proneuralnt,proportion_spacing_matrix)

pdf("~/neurod2peaks_neurod2cor_spacing_withbackground.pdf", useDingbats = F, height = 3)
pheatmap::pheatmap(join_heatmap,
                   cluster_rows = F, cluster_cols = F,  color = rev(hcl.colors(50, "LaJolla")) )
dev.off()






