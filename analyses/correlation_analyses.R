#### CORRELATION ################
##################################


setwd("~/openchr/noack")
cortex.atac <- readRDS('~/openchr/noack/cortex_atac_imputed_expression.RDS')
atac_counts<- readRDS('~/openchr/noack/mypeaks_featurematrix.RDS')
metadata<-cortex.atac@meta.data
barcodes<-rownames(metadata)
sct_expression<-cortex.atac@assays$RNA@data
rm(cortex.atac)
gc()

# Function that classifies peaks according to their correlation with the expression of the desired gene
correlation_expression_peaks<-function(gene, pseudobulk_size)  {
  require(dplyr)
  #sortear las celulas por la expresion del gen y hacer pseudobulks
  atac_counts<-atac_counts[,barcodes]
  expression_neurod2<-sct_expression[gene,]
  expression_neurod2<-sort(expression_neurod2)
  atac_counts<-atac_counts[,names(expression_neurod2)]
  gc()
  i<-1
  #loop que te hace los seudobulks
  while (i<length(colnames(atac_counts))) {
    print(i)
    if (i==1) { 
      pseudobulks<-rowSums(as.matrix(atac_counts[,i:(i+pseudobulk_size)]))
      neurod2_expression_pseudos<-mean(expression_neurod2[i:(i+pseudobulk_size)])
    }
    else if ( (i+pseudobulk_size) > length(colnames(atac_counts)) ) { 
      pseudobulk<-rowSums(as.matrix(atac_counts[ ,i:length(colnames(atac_counts)) ]))
      pseudobulks<-cbind(pseudobulks,pseudobulk)
      neurod2_expression_pseudo<-mean(expression_neurod2[ i:length(colnames(atac_counts)) ])
      neurod2_expression_pseudos<-c(neurod2_expression_pseudos,neurod2_expression_pseudo)
    }
    else {
      pseudobulk<-rowSums(as.matrix(atac_counts[,i:(i+pseudobulk_size)]))
      pseudobulks<-cbind(pseudobulks,pseudobulk)
      neurod2_expression_pseudo<-mean(expression_neurod2[i:(i+pseudobulk_size)])
      neurod2_expression_pseudos<-c(neurod2_expression_pseudos,neurod2_expression_pseudo)
    }
    i<-i+pseudobulk_size 
  }
  gc()
  rownames(pseudobulks)<-rownames(atac_counts)
  #cpms
  pseudobulks<-t(t(pseudobulks)/colSums(pseudobulks)*1000000)
  require(HybridMTest)
  sperma<-row.spearman(pseudobulks, neurod2_expression_pseudos)
  gc()
  sperma<-sperma %>% arrange(stat)
  rownames(sperma)<-paste0(rownames(sperma),"-",sperma$stat)
  sperma<-as.data.frame(stringr::str_split_fixed(rownames(sperma),"-",4) )
  colnames(sperma)<-c("chr","start","end","stat")
  sperma$chr<-as.character(sperma$chr)
  sperma$start<-as.numeric(sperma$start)
  sperma$end<-as.numeric(sperma$end)
  sperma$stat<-as.numeric(sperma$stat)
  return(sperma)
}

neurod2_correlation_atac<-correlation_expression_peaks(gene="Neurod2", pseudobulk_size=50)
neurog2_correlation_atac<-correlation_expression_peaks(gene="Neurog2", pseudobulk_size=50)

# Aggregate the expression of Neurod2 and Neurog2
Proneural<-colSums(rbind(sct_expression["Neurod2",],sct_expression["Neurog2",] ) )
sct_expression<-rbind(sct_expression,Proneural)
proneural_correlation_atac<-correlation_expression_peaks(gene="Proneural", pseudobulk_size=50)


# Classify the ATAC-seq peaks according to the intersection with the proneural factors.
system("
  cd ~/openchr/noack
  cat <(awk -v r='2.private_neurog2' '{print $1,$2,$3,r}' atac_private_neurog2.bed) <(awk -v r='4.private_neurod2' '{print $1,$2,$3,r}' atac_private_neurod2.bed) \
  <(awk -v r='3.shared' '{print $1,$2,$3,r}' atac_shared.bed) <(awk -v r='1.proneuralnt' '{print $1,$2,$3,r}' proneuralnt.bed) | tr ' ' '\t' > proneural_table.bed
")


#funcion que te saca para los picos ordenados por correlacion un stacked de los cuatro sets de picos
stalked<-function(correlation, nquantiles, factor) {
  require(bedr)
  require(ggplot2)
  proneural_table<-read.table("~/openchr/noack/proneural_table.bed")
  #sacar solo el midpoint para que no haya problemas con el itnersect
  proneural_table$V2<-proneural_table$V2+250
  proneural_table$V3<-proneural_table$V3-250
  intersect<-bedr.join.region(bedr.sort.region(correlation),bedr.sort.region(proneural_table) )
  intersect$stat<-as.numeric(intersect$stat)
  colnames(intersect)[8]<-"factor"
  intersect<-intersect %>% filter(factor!=".")
  intersect<-arrange(intersect,as.numeric(stat) )
  #poner los ids
  quantil<-as.integer(length(intersect$chr)/nquantiles)
  quantiles<-rep(1:nquantiles, each=quantil)
  quantiles<-append(quantiles,rep(nquantiles,length(intersect$chr)-length(quantiles)))
  intersect$quantiles<-quantiles
  summarised<-intersect %>% group_by(quantiles,factor) %>% summarise(n=n())
  colors<-rev(c("#CCFF00","#99FFFF","#FF6666","gray"))
  p<-ggplot(summarised, aes(fill=factor, y=n, x=quantiles)) +
    geom_bar(position="fill", stat="identity", colour="black", size=0.08) + theme_classic() + 
    scale_fill_manual(values=colors) + ggtitle(paste0("Correlation with ",factor))
  print(p)
  intersect_midpoints<-intersect %>% dplyr::select(chr,start, end)
  write.table(intersect_midpoints,file=paste0("~/openchr/noack/",factor,"_correlation_sorted_peaks.bed"),sep = "\t", col.names = F, row.names = F, quote = F)
  intersect_midpoints$start<-intersect_midpoints$start+(intersect_midpoints$end-intersect_midpoints$start)/2
  intersect_midpoints$end<-intersect_midpoints$start+1
  write.table(intersect_midpoints,file=paste0("~/openchr/noack/",factor,"_correlation_sorted_summits.bed"),sep = "\t", col.names = F, row.names = F, quote = F)
}


pdf("correlation_stacked.pdf", useDingbats = F, width = 8, height = 4)
stalked(correlation = neurog2_correlation_atac, nquantiles = 50, factor="Neurog2")
stalked(correlation = neurod2_correlation_atac, nquantiles = 50, factor="Neurod2")
stalked(correlation = proneural_correlation_atac, nquantiles = 50, factor="Proneural")
dev.off()

stalked_fc(correlation = neurog2_correlation_atac, nquantiles = 50, factor="Neurog2")
stalked_fc(correlation = neurod2_correlation_atac, nquantiles = 50, factor="Neurod2")
stalked_fc(correlation = proneural_correlation_atac, nquantiles = 50, factor="Proneural")


#funcion que te saca para los picos ordenados por correlacion un stacked de los promotores vs los enhancers
system("
  cd ~/openchr/noack
  cat promoters.bed enhancers.bed | sort -k1,1 -k2,2n | tr "_" "\t" | cut -f1,2,3,4 > promenh.bed
")

stalked_promenh<-function(correlation, nquantiles, factor) {
  require(bedr)
  require(ggplot2)
  proneural_table<-read.table("~/openchr/noack/promenh.bed")
  #sacar solo el midpoint para que no haya problemas con el itnersect
  intersect<-bedr.join.region(bedr.sort.region(correlation),bedr.sort.region(proneural_table) )
  intersect$stat<-as.numeric(intersect$stat)
  colnames(intersect)[8]<-"factor"
  intersect<-intersect %>% filter(factor!=".")
  intersect<-arrange(intersect,as.numeric(stat) )
  #poner los ids
  quantil<-as.integer(length(intersect$chr)/nquantiles)
  quantiles<-rep(1:nquantiles, each=quantil)
  quantiles<-append(quantiles,rep(nquantiles,length(intersect$chr)-length(quantiles)))
  intersect$quantiles<-quantiles
  summarised<-intersect %>% group_by(quantiles,factor) %>% summarise(n=n())
  colors<-c("#9999FF","#CC0000")
  p<-ggplot(summarised, aes(fill=factor, y=n, x=quantiles)) +
    geom_bar(position="fill", stat="identity", colour="black", size=0.08) + theme_classic() + 
    scale_fill_manual(values=colors) + ggtitle(paste0("Correlation with ",factor))
  print(p)
  intersect_midpoints<-intersect %>% dplyr::select(chr,start, end)
  write.table(intersect_midpoints,file=paste0("~/openchr/noack/",factor,"_correlation_sorted_peaks.bed"),sep = "\t", col.names = F, row.names = F, quote = F)
  intersect_midpoints$start<-intersect_midpoints$start+(intersect_midpoints$end-intersect_midpoints$start)/2
  intersect_midpoints$end<-intersect_midpoints$start+1
  write.table(intersect_midpoints,file=paste0("~/openchr/noack/",factor,"_correlation_sorted_summits.bed"),sep = "\t", col.names = F, row.names = F, quote = F)
}

stalked_promenh_intersect<-function(correlation, nquantiles, factor, intersect_factor, factor2) {
  require(bedr)
  require(ggplot2)
  proneural_table<-read.table("~/openchr/noack/promenh.bed")
  colnames(proneural_table)<-c("chr","start","end","factor")
  #sacar solo el midpoint para que no haya problemas con el itnersect
  intersect<-bedr.join.region(bedr.sort.region(correlation),bedr.sort.region(proneural_table) )
  intersect$stat<-as.numeric(intersect$stat)
  intersect<-intersect %>% filter(factor!=".")
  intersect<-arrange(intersect,as.numeric(stat) )
  #poner los ids
  quantil<-as.integer(length(intersect$chr)/nquantiles)
  quantiles<-rep(1:nquantiles, each=quantil)
  quantiles<-append(quantiles,rep(nquantiles,length(intersect$chr)-length(quantiles)))
  intersect$quantiles<-quantiles
  intersect<-bedr.join.region(bedr.sort.region(intersect),bedr.sort.region(intersect_factor) ) %>% filter(V5!=".")
  intersect$quantiles<-as.numeric(intersect$quantiles)
  summarised<-intersect %>% group_by(quantiles,factor) %>% summarise(n=n())
  colors<-c("#9999FF","#CC0000")
  p<-ggplot(summarised, aes(fill=factor, y=n, x=quantiles)) +
    geom_bar(position="fill", stat="identity", colour="black", size=0.08) + theme_classic() + 
    scale_fill_manual(values=colors) + ggtitle(paste0("Correlation with ",factor, ", peaks that intersect with", factor2))
  print(p)
  intersect_midpoints<-intersect %>% dplyr::select(chr,start, end)
  write.table(intersect_midpoints,file=paste0("~/openchr/noack/",factor,"_correlation_sorted_peaks.bed"),sep = "\t", col.names = F, row.names = F, quote = F)
  intersect_midpoints$start<-intersect_midpoints$start+(intersect_midpoints$end-intersect_midpoints$start)/2
  intersect_midpoints$end<-intersect_midpoints$start+1
  write.table(intersect_midpoints,file=paste0("~/openchr/noack/",factor,"_correlation_sorted_summits.bed"),sep = "\t", col.names = F, row.names = F, quote = F)
}


pdf("correlation_stacked_promoters_enhancers.pdf", useDingbats = F, width = 8, height = 2)
stalked_promenh(correlation = neurog2_correlation_atac, nquantiles = 50, factor="Neurog2")
stalked_promenh(correlation = neurod2_correlation_atac, nquantiles = 50, factor="Neurod2")
stalked_promenh(correlation = proneural_correlation_atac, nquantiles = 50, factor="Proneural")
dev.off()


neurod2_peaks<-read.table("~/neurod2/chipseq/e14.5.neurod2_highconfsubpeaks_summit.bed")
neurog2_peaks<-read.table("~/neurog2/chipseq/e14.5.neurog2_peaks.subpeaks.bed")
system("
  bedtools intersect -a ~/neurod2/chipseq/e14.5.neurod2_highconfsubpeaks_summit.bed -b ~/neurog2/chipseq/e14.5.neurog2_peaks.subpeaks.bed > ~/shared_peaks.bed
  bedtools intersect -v -a ~/neurod2/chipseq/e14.5.neurod2_highconfsubpeaks_summit.bed -b ~/neurog2/chipseq/e14.5.neurog2_peaks.subpeaks.bed > ~/neurod2_private.bed
  bedtools intersect -v -b ~/neurod2/chipseq/e14.5.neurod2_highconfsubpeaks_summit.bed -a ~/neurog2/chipseq/e14.5.neurog2_peaks.subpeaks.bed > ~/neurog2_shared.bed
")
shared_peaks<-read.table("~/shared_peaks.bed")


pdf("correlation_stacked_promoters_enhancers_factor_intersect.pdf", useDingbats = F, width = 8, height = 2)
stalked_promenh_intersect(correlation = neurog2_correlation_atac, nquantiles = 25, factor="Neurog2", intersect_factor = neurog2_peaks, factor2 = "Neurog2 private")
stalked_promenh_intersect(correlation = neurod2_correlation_atac, nquantiles = 25, factor="Neurod2", intersect_factor = neurod2_peaks, factor2 = "Neurod2 private")
stalked_promenh_intersect(correlation = proneural_correlation_atac, nquantiles = 25, factor="Proneural", intersect_factor = shared_peaks, factor2 = "shared peaks")
dev.off()






system("
  tunnel_shiva
  export_shiva openchr/noack/\*_correlation_sorted_summits.bed p/openchr/noack
")


system('
  cd ~/p/openchr/noack
  module load Miniconda3/4.12.0
  for file in Neurod2_correlation_sorted_summits.bed Neurog2_correlation_sorted_summits.bed Proneural_correlation_sorted_summits.bed
  
  #los pongo repetidos para luego quitar en el inkscape y asi tener una separacion
  computeMatrix reference-point -p 250 -S ~/p/openchr/noack/imputed_cells_pseudotime_pseudobulks/clusterpseudotime_pseudobulk_0_NSC.bam.bw ~/p/openchr/noack/imputed_cells_pseudotime_pseudobulks/clusterpseudotime_pseudobulk_1_NSC.bam.bw ~/p/openchr/noack/imputed_cells_pseudotime_pseudobulks/clusterpseudotime_pseudobulk_2_NSC.bam.bw ~/p/openchr/noack/imputed_cells_pseudotime_pseudobulks/clusterpseudotime_pseudobulk_3_NSC.bam.bw ~/p/openchr/noack/imputed_cells_pseudotime_pseudobulks/clusterpseudotime_pseudobulk_4_NSC.bam.bw ~/p/openchr/noack/imputed_cells_pseudotime_pseudobulks/clusterpseudotime_pseudobulk_5_NSC.bam.bw ~/p/openchr/noack/imputed_cells_pseudotime_pseudobulks/clusterpseudotime_pseudobulk_6_NSC.bam.bw ~/p/openchr/noack/imputed_cells_pseudotime_pseudobulks/clusterpseudotime_pseudobulk_7_NSC.bam.bw \
   ~/p/openchr/noack/imputed_cells_pseudotime_pseudobulks/clusterpseudotime_pseudobulk_0_IPC.bam.bw ~/p/openchr/noack/imputed_cells_pseudotime_pseudobulks/clusterpseudotime_pseudobulk_1_IPC.bam.bw ~/p/openchr/noack/imputed_cells_pseudotime_pseudobulks/clusterpseudotime_pseudobulk_2_IPC.bam.bw ~/p/openchr/noack/imputed_cells_pseudotime_pseudobulks/clusterpseudotime_pseudobulk_3_IPC.bam.bw ~/p/openchr/noack/imputed_cells_pseudotime_pseudobulks/clusterpseudotime_pseudobulk_4_IPC.bam.bw ~/p/openchr/noack/imputed_cells_pseudotime_pseudobulks/clusterpseudotime_pseudobulk_5_IPC.bam.bw \
   ~/p/openchr/noack/imputed_cells_pseudotime_pseudobulks/clusterpseudotime_pseudobulk_0_PN1.bam.bw ~/p/openchr/noack/imputed_cells_pseudotime_pseudobulks/clusterpseudotime_pseudobulk_1_PN1.bam.bw ~/p/openchr/noack/imputed_cells_pseudotime_pseudobulks/clusterpseudotime_pseudobulk_2_PN1.bam.bw ~/p/openchr/noack/imputed_cells_pseudotime_pseudobulks/clusterpseudotime_pseudobulk_3_PN1.bam.bw ~/p/openchr/noack/imputed_cells_pseudotime_pseudobulks/clusterpseudotime_pseudobulk_4_PN1.bam.bw ~/p/openchr/noack/imputed_cells_pseudotime_pseudobulks/clusterpseudotime_pseudobulk_5_PN1.bam.bw ~/p/openchr/noack/imputed_cells_pseudotime_pseudobulks/clusterpseudotime_pseudobulk_6_PN1.bam.bw \
   ~/p/openchr/noack/imputed_cells_pseudotime_pseudobulks/clusterpseudotime_pseudobulk_0_PN2.bam.bw ~/p/openchr/noack/imputed_cells_pseudotime_pseudobulks/clusterpseudotime_pseudobulk_1_PN2.bam.bw ~/p/openchr/noack/imputed_cells_pseudotime_pseudobulks/clusterpseudotime_pseudobulk_2_PN2.bam.bw ~/p/openchr/noack/imputed_cells_pseudotime_pseudobulks/clusterpseudotime_pseudobulk_3_PN2.bam.bw ~/p/openchr/noack/imputed_cells_pseudotime_pseudobulks/clusterpseudotime_pseudobulk_4_PN2.bam.bw ~/p/openchr/noack/imputed_cells_pseudotime_pseudobulks/clusterpseudotime_pseudobulk_5_PN2.bam.bw ~/p/openchr/noack/imputed_cells_pseudotime_pseudobulks/clusterpseudotime_pseudobulk_6_PN2.bam.bw ~/p/openchr/noack/imputed_cells_pseudotime_pseudobulks/clusterpseudotime_pseudobulk_7_PN2.bam.bw \
   ~/p/openchr/noack/imputed_cells_pseudotime_pseudobulks/clusterpseudotime_pseudobulk_0_PN3.bam.bw ~/p/openchr/noack/imputed_cells_pseudotime_pseudobulks/clusterpseudotime_pseudobulk_1_PN3.bam.bw \
   -R $file -o $file.tmp --sortRegions keep -a 750 -b 750

  for file in Neurod2_correlation_sorted_summits.bed Neurog2_correlation_sorted_summits.bed Proneural_correlation_sorted_summits.bed
  do
    plotHeatmap --matrixFile $file.tmp --outFileName $file.pdf --heatmapHeight 60 -T $file \
    --colorMap viridis --sortRegions keep -max 30 --whatToShow 'heatmap only' --plotFileFormat pdf  --sortRegions keep"
  done
')


#findear los motivos de los pios de atac
#las trayectorias las hice con subpicos longitud natural. y los heatmaps sobre los midpoints de esos picos
# luego el centrimo lo hice sobre los summits de esos mismos picos. que no se como saque esos summits pero bueno
#la file de los summits:
subpeaks_summits_centered<-read.table("~/openchr/noack/centrimo/all_summits.bed")
#extenderle un poco mas a los lados
subpeaks_summits_centered$V2<-subpeaks_summits_centered$V2-150
subpeaks_summits_centered$V3<-subpeaks_summits_centered$V3+150

subpeaks_summits_centered$summit<-subpeaks_summits_centered$V2+400
subpeaks_summits_centered$peak_id<-paste0("peak",c(1:length(subpeaks_summits_centered$V1)))
colnames(subpeaks_summits_centered)<-c("chr","start","end","summit","peak_id")

j<-1
while ((j+999)< length(subpeaks_summits_centered$chr)) {
  
  if (j==1) {
    atac_peaks_motif<-find_motifs(peaks=subpeaks_summits_centered[j:(j+999),],peaks_file="~/openchr/noack/all_summits_expanded.bed")
  }
  else {
    jj<-find_motifs(subpeaks_summits_centered[j:(j+999),],peaks_file="~/openchr/noack/all_summits_expanded.bed")
    atac_peaks_motif<-rbind(atac_peaks_motif,jj)
  }
  j<-j+1000
  print("");print("");print(j);print("");print("")
}
summit_centered_atac_motifs<-atac_peaks_motif
saveRDS(summit_centered_atac_motifs,"~/openchr/noack/atac_summits_expanded_motifs.Rds")
summit_centered_atac_motifs<-readRDS("~/openchr/noack/atac_summits_expanded_motifs.Rds")
# funcion que te saca los motivos del factor que quieras en los picos de atac sorteados que quieras
correlation_atac_motifs<-function(correlation,factor,nquantiles,boundary_1,boundary_2) {
  setwd("~/openchr/noack")
  cor<-read.table(correlation)
  colnames(cor)<-c("chr","start","end")
  #poner los ids de los quantiles
  quantil<-as.integer(length(cor$chr)/nquantiles)
  quantiles<-rep(1:nquantiles, each=quantil)
  quantiles<-append(quantiles,rep(nquantiles,length(cor$chr)-length(quantiles)))
  cor$quantil<-quantiles
  master<-read.table("proneural_table.bed")
  inter<-bedr.join.region(bedr.sort.region(cor),bedr.sort.region(master) )
  inter_factor<-inter[grep(factor,inter$V4),]
  #intersecar los picos de los quantiles con los motivos
  summit_centered_atac_motifs$end<-as.numeric(summit_centered_atac_motifs$summit)
  summit_centered_atac_motifs$start<-summit_centered_atac_motifs$end-1
  summit_centered_atac_motifs$chr<-as.character(summit_centered_atac_motifs$chr)
  interr<-bedr.join.region(bedr.sort.region(summit_centered_atac_motifs),bedr.sort.region(inter_factor) ) %>% filter(quantil!=".")
  #sacar por cada decil cuantos picos hay
  npeaks<-interr %>% select(quantil,peak_id) %>% distinct() %>% group_by(quantil) %>% summarise(npeaks=n())
  #agregar una extra data frame por si no hay picos
  npeaks_extra<-data.frame("quantil"=c(1:nquantiles), "npeaks"=rep(0,nquantiles) )
  npeaks<-rbind(npeaks,npeaks_extra) %>% group_by(quantil) %>% summarise(npeaks=sum(npeaks))
  npeaks<-arrange(npeaks, as.numeric(quantil))
  
  #filtrar los picos que no tienen ebox
  interr<-interr %>% filter(ebox_seq != ".")
  #dividir la tabla en los cercanos al summit y los lejanos
  interr$dist_summit<-abs(as.numeric(interr$summit)-as.numeric(interr$ebox_start)+5)
  interr_close<-interr %>% filter(dist_summit<boundary_1)
  interr_far<-interr %>% filter(dist_summit>boundary_2)
  
  interr_close<-interr_close %>% group_by(quantil,central_dinucleotide) %>% summarise(nmotifs=n())
  interr_far<-interr_far %>% group_by(quantil,central_dinucleotide) %>% summarise(nmotifs=n())
  
  #crear un dataframe con 0s para añadir los motivos que faltan
  indices_motivos<-as.character(rep(1:nquantiles, each=10))
  indices_motifs<-rep( c("CAT-CAT","CAT-CAG","CAG-CAG","CAG-CAC","CAT-CAC","CAC-CAC","CAA-CAA","CAA-CAG","CAA-CAT","CAA-CAC"), nquantiles )
  nmotivos<-rep(0,nquantiles)
  extraa<-data.frame("quantil"=indices_motivos,"central_dinucleotide"=indices_motifs,"nmotifs"=nmotivos)
  #añadir los 0s
  interr_close<-rbind(interr_close,extraa) %>% group_by(quantil,central_dinucleotide) %>% summarise(nmotifs=sum(nmotifs))
  interr_far<-rbind(interr_far,extraa) %>% group_by(quantil,central_dinucleotide) %>% summarise(nmotifs=sum(nmotifs))
  #añadir el background y el npeaks
  interr_close$nmotifs_background<-interr_far$nmotifs
  interr_close<-arrange(interr_close,as.numeric(quantil))
  interr_close$npeaks<-rep(npeaks$npeaks, each=10)
  
  interr_close$centrality<-(interr_close$nmotifs/(boundary_1*2))/(interr_close$nmotifs_background/((400-boundary_2)*2))
  interr_close$nmotifs_peak<-interr_close$nmotifs/interr_close$npeaks
  interr_close$quantil<-as.numeric(interr_close$quantil)
  
  require(ggplot2)
  p<-ggplot(interr_close, aes(x=quantil, y=nmotifs_peak, group=central_dinucleotide)) +
    geom_line(aes(color=central_dinucleotide)) + scale_color_manual(values=colors) + 
    geom_point(aes(color=central_dinucleotide,size=centrality)) + scale_x_continuous(n.breaks=10)
  print(p)
  
}

correlation_atac_motifs(correlation="Neurod2_correlation_sorted_peaks.bed", factor="neurod2", nquantiles=50, boundary_1=25, boundary_2=200)
correlation_atac_motifs(correlation="Neurod2_correlation_sorted_peaks.bed", factor="neurog2", nquantiles=50, boundary_1=25, boundary_2=200)
correlation_atac_motifs(correlation="Neurod2_correlation_sorted_peaks.bed", factor="shared", nquantiles=50, boundary_1=25, boundary_2=200)
correlation_atac_motifs(correlation="Neurod2_correlation_sorted_peaks.bed", factor="proneuralnt", nquantiles=50, boundary_1=25, boundary_2=200)

correlation_atac_motifs(correlation="Neurog2_correlation_sorted_peaks.bed", factor="neurod2", nquantiles=50, boundary_1=25, boundary_2=200)
correlation_atac_motifs(correlation="Neurog2_correlation_sorted_peaks.bed", factor="neurog2", nquantiles=50, boundary_1=25, boundary_2=200)
correlation_atac_motifs(correlation="Neurog2_correlation_sorted_peaks.bed", factor="shared", nquantiles=50, boundary_1=25, boundary_2=200)
correlation_atac_motifs(correlation="Neurog2_correlation_sorted_peaks.bed", factor="proneuralnt", nquantiles=50, boundary_1=25, boundary_2=200)

correlation_atac_motifs(correlation="Proneural_correlation_sorted_peaks.bed", factor="neurod2", nquantiles=50, boundary_1=25, boundary_2=200)
correlation_atac_motifs(correlation="Proneural_correlation_sorted_peaks.bed", factor="neurog2", nquantiles=50, boundary_1=25, boundary_2=200)
correlation_atac_motifs(correlation="Proneural_correlation_sorted_peaks.bed", factor="shared", nquantiles=50, boundary_1=25, boundary_2=200)
correlation_atac_motifs(correlation="Proneural_correlation_sorted_peaks.bed", factor="proneuralnt", nquantiles=50, boundary_1=25, boundary_2=200)


#### DINUCLEOTIDES IN CHIP-SEQ SUMMITS #
#


neurog2_extended_peaks_motifs<-readRDS("~/neurog2/neurog2_extended_peaks_motifs.Rds")
neurod2_extended_peaks_motifs<-readRDS("~/neurod2/neurod2_extended_peaks_motifs.Rds")

neurod2_peaks<-read.table("~/neurod2/chipseq/e14.5.neurod2_highconfsubpeaks_summit.bed")
neurog2_peaks<-read.table("~/neurog2/chipseq/e14.5.neurog2_peaks.subpeaks.bed")



correlation_factor_motifs<-function(correlation,factor,extended,nquantiles,boundary_1,boundary_2) {
  require(dplyr)
  require(bedr)
  require(ggpubr)
  setwd("~/openchr/noack")
  cor<-read.table(correlation) %>% distinct()
  colnames(cor)<-c("chr","start","end","orden")
  #intersecar los motivos de los factores con los picos ordenados
  interr<-bedr.join.region( bedr.sort.region(factor),bedr.sort.region(cor) ) %>% filter(chr.b !=".") %>% distinct()
  interr<-arrange(interr, as.numeric(orden) ) 
  
  #poner los ids de los quantiles
  quantil<-as.integer(length(interr$chr)/nquantiles)
  quantiles<-rep(1:nquantiles, each=quantil)
  quantiles<-append(quantiles,rep(nquantiles,length(interr$chr)-length(quantiles)))
  interr$quantil<-quantiles
  
  interr<-interr %>% select(V5,quantil) %>% distinct()
  colnames(interr)[1]<-"peak_id"
  
  #sacar por cada quantil cuantos picos hay
  npeaks<-interr %>% select(quantil,peak_id) %>% distinct() %>% group_by(quantil) %>% summarise(npeaks=n())
  
  
  #mergear la tabla de los quantiles con la de los motivos extendidos
  interr_motifs<-merge(interr, extended, by="peak_id")
  
  #filtrar los picos que no tienen ebox
  interr_motifs<-interr_motifs %>% filter(ebox_seq != ".")
  #dividir la tabla en los cercanos al summit y los lejanos
  interr_motifs$dist_summit<-abs(as.numeric(interr_motifs$summit)-as.numeric(interr_motifs$ebox_start)+5)
  interr_close<-interr_motifs %>% filter(dist_summit<boundary_1)
  interr_far<-interr_motifs %>% filter(dist_summit>boundary_2)
  
  interr_close<-interr_close %>% group_by(quantil,central_dinucleotide) %>% summarise(nmotifs=n())
  interr_far<-interr_far %>% group_by(quantil,central_dinucleotide) %>% summarise(nmotifs=n())
  
  #crear un dataframe con 0s para añadir los motivos que faltan
  indices_motivos<-as.character(rep(1:nquantiles, each=10))
  indices_motifs<-rep( c("CAT-CAT","CAT-CAG","CAG-CAG","CAG-CAC","CAT-CAC","CAC-CAC","CAA-CAA","CAA-CAG","CAA-CAT","CAA-CAC"), nquantiles )
  nmotivos<-rep(0,nquantiles)
  extraa<-data.frame("quantil"=indices_motivos,"central_dinucleotide"=indices_motifs,"nmotifs"=nmotivos)
  extraa$quantil<-as.double(extraa$quantil)
  #añadir los 0s
  interr_close<-rbind(interr_close,extraa) %>% group_by(quantil,central_dinucleotide) %>% summarise(nmotifs=sum(nmotifs))
  interr_far<-rbind(interr_far,extraa) %>% group_by(quantil,central_dinucleotide) %>% summarise(nmotifs=sum(nmotifs))
  #añadir el background y el npeaks
  interr_close$nmotifs_background<-interr_far$nmotifs
  interr_close<-arrange(interr_close,as.numeric(quantil))
  interr_close$npeaks<-rep(npeaks$npeaks, each=10)
  
  interr_close$centrality<-(interr_close$nmotifs/(boundary_1*2))/(interr_close$nmotifs_background/((400-boundary_2)*2))
  interr_close$nmotifs_peak<-interr_close$nmotifs/interr_close$npeaks
  interr_close$quantil<-as.numeric(interr_close$quantil)
  
  interr_close<-interr_close %>% select(quantil,central_dinucleotide,centrality,nmotifs_peak)
  
  #añadir una columna que sirva para escalar las demas columnas
  scalee<-data.frame("quantil"=1:nquantiles)
  scalee$centrality<-8
  scalee$nmotifs_peak<-0.4
  scalee$central_dinucleotide<-"scale"
  summ<-rbind(interr_close,scalee)
  
  summ$centrality[which(summ$centrality>8)]<-8
  summ$quantil<-as.character(summ$quantil)
  summ$quantil<-as.factor(summ$quantil)
  summ$quantil<-factor(summ$quantil, levels = 1:nquantiles)
  
  require(ggplot2)
  p<-ggballoonplot(summ, y = "central_dinucleotide", x = "quantil",
                   size = "nmotifs_peak", fill = "centrality")+ 
    gradient_fill(c("black","darkred","yellow2"))
  print(p)
}


pdf("~/neurog2peaks_cor_balloon.pdf")
correlation_factor_motifs(correlation="Neurog2_correlation_sorted_summits_extended.bed", factor=neurog2_peaks, 
                          extended=neurog2_extended_peaks_motifs,nquantiles=10, boundary_1=25, boundary_2=150 )
correlation_factor_motifs(correlation="Proneural_correlation_sorted_summits_extended.bed", factor=neurog2_peaks, 
                          extended=neurog2_extended_peaks_motifs,nquantiles=10, boundary_1=25, boundary_2=150 )
correlation_factor_motifs(correlation="Neurod2_correlation_sorted_summits_extended.bed", factor=neurog2_peaks, 
                          extended=neurog2_extended_peaks_motifs,nquantiles=10, boundary_1=25, boundary_2=150 )
dev.off()

pdf("~/neurod2peaks_cor_balloon.pdf")
correlation_factor_motifs(correlation="Neurog2_correlation_sorted_summits_extended.bed", factor=neurod2_peaks, 
                          extended=neurod2_extended_peaks_motifs,nquantiles=10, boundary_1=25, boundary_2=150 )
correlation_factor_motifs(correlation="Proneural_correlation_sorted_summits_extended.bed", factor=neurod2_peaks, 
                          extended=neurod2_extended_peaks_motifs,nquantiles=10, boundary_1=, boundary_2=150 )
correlation_factor_motifs(correlation="Neurod2_correlation_sorted_summits_extended.bed", factor=neurod2_peaks, 
                          extended=neurod2_extended_peaks_motifs,nquantiles=10, boundary_1=25, boundary_2=200 )
dev.off()



#plot del numero de motivos por pico
correlation_factor_nmotifs_peak<-function(correlation,factor,motifs,nquantiles) {
  require(dplyr)
  require(bedr)
  require(ggpubr)
  setwd("~/openchr/noack")
  cor<-read.table(correlation) %>% distinct()
  colnames(cor)<-c("chr","start","end","orden")
  #intersecar los motivos de los factores con los picos ordenados
  interr<-bedr.join.region( bedr.sort.region(factor),bedr.sort.region(cor) ) %>% filter(chr.b !=".") %>% distinct()
  interr<-arrange(interr, as.numeric(orden) ) 
  
  #poner los ids de los quantiles
  quantil<-as.integer(length(interr$chr)/nquantiles)
  quantiles<-rep(1:nquantiles, each=quantil)
  quantiles<-append(quantiles,rep(nquantiles,length(interr$chr)-length(quantiles)))
  interr$quantil<-quantiles
  
  interr<-interr %>% select(V5,quantil) %>% distinct()
  colnames(interr)[1]<-"peak_id"
  
  #sacar por cada quantil cuantos picos hay
  npeaks<-interr %>% select(quantil,peak_id) %>% distinct() %>% group_by(quantil) %>% summarise(npeaks=n())
  
  #mergear la tabla de los quantiles con la de los motivos
  interr_motifs<-merge(interr, motifs, by="peak_id")
  
  #sacar el numero de motivos por pico
  interr_motifs$motif_presence<-rep(1,length(interr_motifs$chr ))
  interr_motifs$motif_presence[which(interr_motifs$central_dinucleotide==".")]<-0
  
  
  nmotifs_peaks_quantiles<-interr_motifs %>% group_by(peak_id,quantil) %>% summarise(nmotifs=sum(motif_presence))
  nmotifs_peaks_quantiles$nmotifs<-as.numeric(nmotifs_peaks_quantiles$nmotifs)
  nmotifs_peaks_quantiles$nmotifs[which(nmotifs_peaks_quantiles$nmotifs>6)]<-7
  npeaks_with_nmotifs<-nmotifs_peaks_quantiles %>% group_by(quantil, nmotifs) %>% summarise(freq=n())
  
  #npeaks_with_nmotifs$nmotifs<-as.character(npeaks_with_nmotifs$nmotifs)
  #npeaks_with_nmotifs$nmotifs[which(npeaks_with_nmotifs$nmotifs=="10")]<-">9"
  
  #npeaks_with_nmotifs$nmotifs<-factor(npeaks_with_nmotifs$nmotifs,levels=c(">9","9","8","7","6","5","4","3","2","1","0") )
  
  npeaks_with_nmotifs$quantil<-as.character(npeaks_with_nmotifs$quantil)
  npeaks_with_nmotifs$quantil<-factor(npeaks_with_nmotifs$quantil,levels=c(1:10))
  
  pdf("~/stacked_suki.pdf", useDingbats = F)
  ggplot(npeaks_with_nmotifs, aes(fill=nmotifs, y=freq, x=quantil)) +
    geom_bar(position="fill", stat="identity", colour="black", size=0.08) + theme_classic() +
    scale_fill_continuous(low="beige", high="darkblue")
  dev.off()
}


pdf("~/neurod2peaks_cor_nmotifs_stacked.pdf")
correlation_factor_motifs ( correlation="Neurod2_correlation_sorted_summits_extended.bed", factor=neurod2_peaks, 
                            motifs=neurod2_peaks_motifs, nquantiles=10 )
dev.off()

setwd("~/neurod2/chipseq")
neurod2_peaks<-read.table( "e14.5.neurod2_highconfsubpeaks_summit.bed" )
neurod2_peaks_motifs<-readRDS("~/neurod2/neurod2_peaks.Rds")


correlation="Neurod2_correlation_sorted_summits_extended.bed"; factor=neurod2_peaks
motifs=neurod2_peaks_motifs; nquantiles=10 




correlation_factor_motifs_fc_barplot<-function(correlation,factor,extended,nquantiles,boundary_1,boundary_2) {
  require(dplyr)
  require(bedr)
  require(ggpubr)
  setwd("~/openchr/noack")
  cor<-read.table(correlation) %>% distinct()
  colnames(cor)<-c("chr","start","end","orden")
  #intersecar los motivos de los factores con los picos ordenados
  interr<-bedr.join.region( bedr.sort.region(factor),bedr.sort.region(cor) ) %>% filter(chr.b !=".") %>% distinct()
  interr<-arrange(interr, as.numeric(orden) ) 
  
  #poner los ids de los quantiles
  quantil<-as.integer(length(interr$chr)/nquantiles)
  quantiles<-rep(1:nquantiles, each=quantil)
  quantiles<-append(quantiles,rep(nquantiles,length(interr$chr)-length(quantiles)))
  interr$quantil<-quantiles
  
  interr<-interr %>% select(V5,quantil) %>% distinct()
  colnames(interr)[1]<-"peak_id"
  
  #sacar por cada quantil cuantos picos hay
  npeaks<-interr %>% select(quantil,peak_id) %>% distinct() %>% group_by(quantil) %>% summarise(npeaks=n())
  
  
  #mergear la tabla de los quantiles con la de los motivos extendidos
  interr_motifs<-merge(interr, extended, by="peak_id")
  
  #filtrar los picos que no tienen ebox
  interr_motifs<-interr_motifs %>% filter(ebox_seq != ".")
  #dividir la tabla en los cercanos al summit y los lejanos
  interr_motifs$dist_summit<-abs(as.numeric(interr_motifs$summit)-as.numeric(interr_motifs$ebox_start)+5)
  interr_close<-interr_motifs %>% filter(dist_summit<boundary_1)
  interr_far<-interr_motifs %>% filter(dist_summit>boundary_2)
  
  interr_close<-interr_close %>% group_by(quantil,central_dinucleotide) %>% summarise(nmotifs=n())
  interr_far<-interr_far %>% group_by(quantil,central_dinucleotide) %>% summarise(nmotifs=n())
  
  #crear un dataframe con 0s para añadir los motivos que faltan
  indices_motivos<-as.character(rep(1:nquantiles, each=10))
  indices_motifs<-rep( c("CAT-CAT","CAT-CAG","CAG-CAG","CAG-CAC","CAT-CAC","CAC-CAC","CAA-CAA","CAA-CAG","CAA-CAT","CAA-CAC"), nquantiles )
  nmotivos<-rep(0,nquantiles)
  extraa<-data.frame("quantil"=indices_motivos,"central_dinucleotide"=indices_motifs,"nmotifs"=nmotivos)
  extraa$quantil<-as.double(extraa$quantil)
  #añadir los 0s
  interr_close<-rbind(interr_close,extraa) %>% group_by(quantil,central_dinucleotide) %>% summarise(nmotifs=sum(nmotifs))
  interr_far<-rbind(interr_far,extraa) %>% group_by(quantil,central_dinucleotide) %>% summarise(nmotifs=sum(nmotifs))
  #añadir el background y el npeaks
  interr_close$nmotifs_background<-interr_far$nmotifs
  interr_close<-arrange(interr_close,as.numeric(quantil))
  interr_close$npeaks<-rep(npeaks$npeaks, each=10)
  
  interr_close$centrality<-(interr_close$nmotifs/(boundary_1*2))/(interr_close$nmotifs_background/((400-boundary_2)*2))
  interr_close$nmotifs_peak<-interr_close$nmotifs/interr_close$npeaks
  interr_close$quantil<-as.numeric(interr_close$quantil)
  
  interr_close<-interr_close %>% select(quantil,central_dinucleotide,centrality,nmotifs_peak)
  
  
  require(ggplot2)
  summ_start<-interr_close %>% filter(quantil==1)
  summ_end<-interr_close %>% filter(quantil==10)
  summ_fc<-summ_end
  summ_fc$fc_nmotifs<-summ_end$nmotifs_peak/summ_start$nmotifs_peak
  summ_fc$fc_centrality<-summ_end$centrality/summ_start$centrality
  
  #añadir una columna que sirva para escalar las demas columnas
  summ_fc<-summ_fc[,c(2,5,6)]
  
  scalee<-c("scalee",6,1.8)
  
  summ_fc<-rbind(summ_fc,scalee)
  summ_fc$fc_centrality<-as.numeric(summ_fc$fc_centrality)
  summ_fc$fc_nmotifs<-as.numeric(summ_fc$fc_nmotifs)
  
  summ_fc$fc_centrality[which(summ_fc$fc_centrality>1.8)]<-1.8
  summ_fc$fc_nmotifs[which(summ_fc$fc_nmotifs>6)]<-6
  
  
  p<-ggplot(summ_fc, aes(fill=fc_centrality, y=fc_nmotifs, x=central_dinucleotide)) + geom_bar(stat="identity") +
    theme_classic() + gradient_fill(c("navy","lightblue"))
  
  print(p)
}



#### FC BARPLOTSS
pdf("~/neurog2peaks_cor_fc_barplot.pdf", useDingbats = F, width = 10)
correlation_factor_motifs_fc_barplot(correlation="Neurog2_correlation_sorted_summits_extended.bed", factor=neurog2_peaks, 
                                     extended=neurog2_extended_peaks_motifs,nquantiles=10, boundary_1=25, boundary_2=150 )
correlation_factor_motifs_fc_barplot(correlation="Proneural_correlation_sorted_summits_extended.bed", factor=neurog2_peaks, 
                                     extended=neurog2_extended_peaks_motifs,nquantiles=10, boundary_1=25, boundary_2=150 )
correlation_factor_motifs_fc_barplot(correlation="Neurod2_correlation_sorted_summits_extended.bed", factor=neurog2_peaks, 
                                     extended=neurog2_extended_peaks_motifs,nquantiles=10, boundary_1=25, boundary_2=150 )
dev.off()

pdf("~/neurod2peaks_cor_fc_barplot.pdf", useDingbats = F, width = 10)
correlation_factor_motifs_fc_barplot(correlation="Neurog2_correlation_sorted_summits_extended.bed", factor=neurod2_peaks, 
                                     extended=neurod2_extended_peaks_motifs,nquantiles=10, boundary_1=25, boundary_2=150 )
correlation_factor_motifs_fc_barplot(correlation="Proneural_correlation_sorted_summits_extended.bed", factor=neurod2_peaks, 
                                     extended=neurod2_extended_peaks_motifs,nquantiles=10, boundary_1=25, boundary_2=150 )
correlation_factor_motifs_fc_barplot(correlation="Neurod2_correlation_sorted_summits_extended.bed", factor=neurod2_peaks, 
                                     extended=neurod2_extended_peaks_motifs,nquantiles=10, boundary_1=25, boundary_2=150 )
dev.off()




################## promotores / enhancers 
correlation_factor_motifs_subdivide<-function(correlation,factor,extended,nquantiles,boundary_1,boundary_2, intersect, overlap) {
  require(dplyr)
  require(bedr)
  require(ggpubr)
  setwd("~/openchr/noack")
  cor<-read.table(correlation) %>% distinct()
  colnames(cor)<-c("chr","start","end","orden")
  #intersecar los motivos de los factores con los picos ordenados
  interr<-bedr.join.region( bedr.sort.region(factor),bedr.sort.region(cor) ) %>% filter(chr.b !=".") %>% distinct()
  interr<-arrange(interr, as.numeric(orden) )
  
  #poner los ids de los quantiles
  quantil<-as.integer(length(interr$chr)/nquantiles)
  quantiles<-rep(1:nquantiles, each=quantil)
  quantiles<-append(quantiles,rep(nquantiles,length(interr$chr)-length(quantiles)))
  interr$quantil<-quantiles
  
  colnames(intersect)<-c("cromosoma","inicio","final")
  #intersecar con el set de peaks que quiera
  if (overlap=="yes") {
    interr<-bedr.join.region( bedr.sort.region(interr), bedr.sort.region(intersect) ) %>% filter(cromosoma != ".")
  }
  
  if (overlap=="no") {
    interr<-bedr.join.region( bedr.sort.region(interr), bedr.sort.region(intersect) ) %>% filter(cromosoma == ".")
  }
  
  colnames(interr)[5]<-"peak_id"
  interr<-interr %>% dplyr::select(peak_id,quantil) %>% distinct()
  
  
  #sacar por cada quantil cuantos picos hay
  npeaks<-interr %>% dplyr::select(quantil,peak_id) %>% distinct() %>% group_by(quantil) %>% summarise(npeaks=n())
  
  
  #mergear la tabla de los quantiles con la de los motivos extendidos
  interr_motifs<-merge(interr, extended, by="peak_id")
  
  #filtrar los picos que no tienen ebox
  interr_motifs<-interr_motifs %>% filter(ebox_seq != ".")
  #dividir la tabla en los cercanos al summit y los lejanos
  interr_motifs$dist_summit<-abs(as.numeric(interr_motifs$summit)-as.numeric(interr_motifs$ebox_start)+5)
  interr_close<-interr_motifs %>% filter(dist_summit<boundary_1)
  interr_far<-interr_motifs %>% filter(dist_summit>boundary_2)
  
  interr_close<-interr_close %>% group_by(quantil,central_dinucleotide) %>% summarise(nmotifs=n())
  interr_close$quantil<-as.double(interr_close$quantil)
  interr_far<-interr_far %>% group_by(quantil,central_dinucleotide) %>% summarise(nmotifs=n())
  interr_far$quantil<-as.double(interr_far$quantil)
  
  #crear un dataframe con 0s para añadir los motivos que faltan
  indices_motivos<-as.character(rep(1:nquantiles, each=10))
  indices_motifs<-rep( c("CAT-CAT","CAT-CAG","CAG-CAG","CAG-CAC","CAT-CAC","CAC-CAC","CAA-CAA","CAA-CAG","CAA-CAT","CAA-CAC"), nquantiles )
  nmotivos<-rep(0,nquantiles)
  extraa<-data.frame("quantil"=indices_motivos,"central_dinucleotide"=indices_motifs,"nmotifs"=nmotivos)
  extraa$quantil<-as.double(extraa$quantil)
  #añadir los 0s
  interr_close<-rbind(interr_close,extraa) %>% group_by(quantil,central_dinucleotide) %>% summarise(nmotifs=sum(nmotifs))
  interr_far<-rbind(interr_far,extraa) %>% group_by(quantil,central_dinucleotide) %>% summarise(nmotifs=sum(nmotifs))
  #añadir el background y el npeaks
  interr_close$nmotifs_background<-interr_far$nmotifs
  interr_close<-arrange(interr_close,as.numeric(quantil))
  interr_close$npeaks<-rep(npeaks$npeaks, each=10)
  
  interr_close$centrality<-(interr_close$nmotifs/(boundary_1*2))/(interr_close$nmotifs_background/((400-boundary_2)*2))
  interr_close$nmotifs_peak<-interr_close$nmotifs/interr_close$npeaks
  interr_close$quantil<-as.numeric(interr_close$quantil)
  
  interr_close<-interr_close %>% dplyr::select(quantil,central_dinucleotide,centrality,nmotifs_peak)
  
  #añadir una columna que sirva para escalar las demas columnas
  scalee<-data.frame("quantil"=1:nquantiles)
  scalee$centrality<-8
  scalee$nmotifs_peak<-0.5
  scalee$central_dinucleotide<-"scale"
  
  scalee2<-data.frame("quantil"=1:nquantiles)
  scalee2$centrality<-0
  scalee2$nmotifs_peak<-0.1
  scalee2$central_dinucleotide<-"scalee2"
  summ<-rbind(interr_close,scalee,scalee2)
  
  summ$centrality[which(summ$centrality>8)]<-8
  summ$quantil<-as.character(summ$quantil)
  summ$quantil<-as.factor(summ$quantil)
  summ$quantil<-factor(summ$quantil, levels = 1:nquantiles)
  summ$nmotifs_peak[which(summ$nmotifs_peak==0)]<-NA
  require(ggplot2)
  p<-ggballoonplot(summ, y = "central_dinucleotide", x = "quantil",
                   size = "nmotifs_peak", fill = "centrality")+ 
    gradient_fill(c("black","darkred","yellow2"))
  print(p)
}


enhancers<-read.table("~/openchr/noack/enhancers.bed")
promoters<-read.table("~/openchr/noack/promoters.bed")[,c(1:3)]

pdf("~/neurod2_motifs_correlation_promenh_subdivided.pdf", height = 6, width = 6)
correlation_factor_motifs_subdivide(correlation="Neurod2_correlation_sorted_summits_extended.bed", factor=neurod2_peaks, 
                                    extended=neurod2_extended_peaks_motifs,nquantiles=10, boundary_1=25, boundary_2=150, 
                                    intersect = enhancers, overlap = "yes" )
correlation_factor_motifs_subdivide(correlation="Neurod2_correlation_sorted_summits_extended.bed", factor=neurod2_peaks, 
                                    extended=neurod2_extended_peaks_motifs,nquantiles=10, boundary_1=25, boundary_2=150, 
                                    intersect = promoters, overlap = "yes" )
dev.off()







# pero primero tengo que mirar si cambian los motivos dependiendo de si hay link o no

correlation_factor_motifs_subdivide_enhancers<-function(correlation,factor,extended,nquantiles,boundary_1,boundary_2, intersect, overlap) {
  require(dplyr)
  require(bedr)
  require(ggpubr)
  setwd("~/openchr/noack")
  cor<-read.table(correlation) %>% distinct()
  colnames(cor)<-c("chr","start","end","orden")
  #intersecar los motivos de los factores con los picos ordenados
  interr<-bedr.join.region( bedr.sort.region(factor),bedr.sort.region(cor) ) %>% filter(chr.b !=".") %>% distinct()
  interr<-arrange(interr, as.numeric(orden) )
  
  #poner los ids de los quantiles
  quantil<-as.integer(length(interr$chr)/nquantiles)
  quantiles<-rep(1:nquantiles, each=quantil)
  quantiles<-append(quantiles,rep(nquantiles,length(interr$chr)-length(quantiles)))
  interr$quantil<-quantiles
  
  
  #primero intersecar con enhancers. porque todos estos de priming son enhancers
  enhancers<-read.table("~/openchr/noack/enhancers.bed")
  colnames(enhancers)<-c("chr_enhancer","start_enhancer","end_enhancer")
  
  interr<-bedr.join.region( bedr.sort.region(interr), bedr.sort.region(enhancers) ) %>% filter(chr_enhancer != ".")
  #intersecar con el set de peaks que quiera
  if (overlap=="yes") {
    interr<-bedr.join.region( bedr.sort.region(interr), bedr.sort.region(intersect) ) %>% filter(chr.bom_enh != ".")
  }
  
  if (overlap=="no") {
    interr<-bedr.join.region( bedr.sort.region(interr), bedr.sort.region(intersect) ) %>% filter(chr.bom_enh == ".")
  }
  
  colnames(interr)[5]<-"peak_id"
  interr<-interr %>% dplyr::select(peak_id,quantil) %>% distinct()
  
  
  #sacar por cada quantil cuantos picos hay
  npeaks<-interr %>% dplyr::select(quantil,peak_id) %>% distinct() %>% group_by(quantil) %>% summarise(npeaks=n())
  
  
  #mergear la tabla de los quantiles con la de los motivos extendidos
  interr_motifs<-merge(interr, extended, by="peak_id")
  
  #filtrar los picos que no tienen ebox
  interr_motifs<-interr_motifs %>% filter(ebox_seq != ".")
  #dividir la tabla en los cercanos al summit y los lejanos
  interr_motifs$dist_summit<-abs(as.numeric(interr_motifs$summit)-as.numeric(interr_motifs$ebox_start)+5)
  interr_close<-interr_motifs %>% filter(dist_summit<boundary_1)
  interr_far<-interr_motifs %>% filter(dist_summit>boundary_2)
  
  interr_close<-interr_close %>% group_by(quantil,central_dinucleotide) %>% summarise(nmotifs=n())
  interr_close$quantil<-as.double(interr_close$quantil)
  interr_far<-interr_far %>% group_by(quantil,central_dinucleotide) %>% summarise(nmotifs=n())
  interr_far$quantil<-as.double(interr_far$quantil)
  
  #crear un dataframe con 0s para añadir los motivos que faltan
  indices_motivos<-as.character(rep(1:nquantiles, each=10))
  indices_motifs<-rep( c("CAT-CAT","CAT-CAG","CAG-CAG","CAG-CAC","CAT-CAC","CAC-CAC","CAA-CAA","CAA-CAG","CAA-CAT","CAA-CAC"), nquantiles )
  nmotivos<-rep(0,nquantiles)
  extraa<-data.frame("quantil"=indices_motivos,"central_dinucleotide"=indices_motifs,"nmotifs"=nmotivos)
  extraa$quantil<-as.double(extraa$quantil)
  #añadir los 0s
  interr_close<-rbind(interr_close,extraa) %>% group_by(quantil,central_dinucleotide) %>% summarise(nmotifs=sum(nmotifs))
  interr_far<-rbind(interr_far,extraa) %>% group_by(quantil,central_dinucleotide) %>% summarise(nmotifs=sum(nmotifs))
  #añadir el background y el npeaks
  interr_close$nmotifs_background<-interr_far$nmotifs
  interr_close<-arrange(interr_close,as.numeric(quantil))
  interr_close$npeaks<-rep(npeaks$npeaks, each=10)
  
  interr_close$centrality<-(interr_close$nmotifs/(boundary_1*2))/(interr_close$nmotifs_background/((400-boundary_2)*2))
  interr_close$nmotifs_peak<-interr_close$nmotifs/interr_close$npeaks
  interr_close$quantil<-as.numeric(interr_close$quantil)
  
  interr_close<-interr_close %>% dplyr::select(quantil,central_dinucleotide,centrality,nmotifs_peak)
  
  #añadir una columna que sirva para escalar las demas columnas
  scalee<-data.frame("quantil"=1:nquantiles)
  scalee$centrality<-8
  scalee$nmotifs_peak<-0.5
  scalee$central_dinucleotide<-"scale"
  
  scalee2<-data.frame("quantil"=1:nquantiles)
  scalee2$centrality<-0
  scalee2$nmotifs_peak<-0.1
  scalee2$central_dinucleotide<-"scalee2"
  summ<-rbind(interr_close,scalee,scalee2)
  
  summ$centrality[which(summ$centrality>8)]<-8
  summ$quantil<-as.character(summ$quantil)
  summ$quantil<-as.factor(summ$quantil)
  summ$quantil<-factor(summ$quantil, levels = 1:nquantiles)
  summ$nmotifs_peak[which(summ$nmotifs_peak==0)]<-NA
  require(ggplot2)
  p<-ggballoonplot(summ, y = "central_dinucleotide", x = "quantil",
                   size = "nmotifs_peak", fill = "centrality")+ 
    gradient_fill(c("black","darkred","yellow2"))
  print(p)
}

correlation="Neurod2_correlation_sorted_summits_extended.bed"; factor=neurod2_peaks
extended=neurod2_extended_peaks_motifs; nquantiles=10; boundary_1=25; boundary_2=150; intersect = positive_links; overlap = "yes"

pdf("~/neurod2_motifs_correlation_sorted_link_subdivided.pdf", height = 6, width = 6)
correlation_factor_motifs_subdivide_enhancers(correlation="Neurod2_correlation_sorted_summits_extended.bed", factor=neurod2_peaks, 
                                              extended=neurod2_extended_peaks_motifs,nquantiles=6, boundary_1=25, boundary_2=150, 
                                              intersect = positive_links, overlap = "no" )
correlation_factor_motifs_subdivide_enhancers(correlation="Neurod2_correlation_sorted_summits_extended.bed", factor=neurod2_peaks, 
                                              extended=neurod2_extended_peaks_motifs,nquantiles=6, boundary_1=25, boundary_2=150, 
                                              intersect = positive_links, overlap = "yes" )
dev.off()

correlation="Neurod2_correlation_sorted_summits_extended.bed"; factor=neurod2_peaks;
extended=neurod2_extended_peaks_motifs; nquantiles=6;  boundary_1=25;  boundary_2=150; 
intersect = positive_links; overlap = "no" 


##### BOXPLOTS CHIPSEQ SIGNAL NEUROD2 #######
#############################################

# Intersect the NEUROD2 peaks sorted by correlation with the bedGraph files showing the genomic coverage of NEUROD2 ChIP-seq.
system("
  module load BEDTools
  cd ~/p/openchr/noack
  cat Neurod2_correlation_sorted_summits_index.bed | awk '{print $1,$2-250,$3+250,$4}' | tr " " "\t" > Neurod2_correlation_sorted_summits_index_exteded_500bp.bed
  bedtools map -a <( sort -k1,1 -k2,2n Neurod2_correlation_sorted_summits_index_exteded_500bp.bed) -b ~/p/neurod2/e14.5.neurod2.1_treat_pileup.bdg -c 4 -o mean > neurod2_pileup_1_Neurod2_correlation_sorted_summits_index_exteded_500bp.bed
  bedtools map -a <( sort -k1,1 -k2,2n Neurod2_correlation_sorted_summits_index_exteded_500bp.bed) -b ~/p/neurod2/e14.5.neurod2.1_control_lambda.bdg -c 4 -o mean > control_lambda_1_Neurod2_correlation_sorted_summits_index_exteded_500bp.bed
  bedtools map -a <( sort -k1,1 -k2,2n Neurod2_correlation_sorted_summits_index_exteded_500bp.bed) -b ~/p/neurod2/e14.5.neurod2.2_treat_pileup.bdg -c 4 -o mean > neurod2_pileup_2_Neurod2_correlation_sorted_summits_index_exteded_500bp.bed
  bedtools map -a <( sort -k1,1 -k2,2n Neurod2_correlation_sorted_summits_index_exteded_500bp.bed) -b ~/p/neurod2/e14.5.neurod2.2_control_lambda.bdg -c 4 -o mean > control_lambda_2_Neurod2_correlation_sorted_summits_index_exteded_500bp.bed
 
  cat proneural_table.bed | grep -v neurog2 | grep -v proneuralnt > proneural_table_neurod2.bed
  cat proneural_table.bed | grep proneuralnt > proneural_table_proneuralnt.bed
  bedtools intersect -wa -a neurod2_pileup_1_Neurod2_correlation_sorted_summits_index_exteded_500bp.bed -b proneural_table_neurod2.bed > neurod2_neurod2_pileup_1_Neurod2_correlation_sorted_summits_index_exteded_500bp.bed
  bedtools intersect -wa -a neurod2_pileup_1_Neurod2_correlation_sorted_summits_index_exteded_500bp.bed -b proneural_table_proneuralnt.bed > proneuralnt_neurod2_pileup_1_Neurod2_correlation_sorted_summits_index_exteded_500bp.bed
  bedtools intersect -wa -a control_lambda_1_Neurod2_correlation_sorted_summits_index_exteded_500bp.bed -b proneural_table_neurod2.bed > neurod2_control_lambda_1_Neurod2_correlation_sorted_summits_index_exteded_500bp.bed
  bedtools intersect -wa -a control_lambda_1_Neurod2_correlation_sorted_summits_index_exteded_500bp.bed -b proneural_table_proneuralnt.bed > proneuralnt_control_lambda_1_Neurod2_correlation_sorted_summits_index_exteded_500bp.bed
  
  bedtools intersect -wa -a neurod2_pileup_2_Neurod2_correlation_sorted_summits_index_exteded_500bp.bed -b proneural_table_neurod2.bed > neurod2_neurod2_pileup_2_Neurod2_correlation_sorted_summits_index_exteded_500bp.bed
  bedtools intersect -wa -a neurod2_pileup_2_Neurod2_correlation_sorted_summits_index_exteded_500bp.bed -b proneural_table_proneuralnt.bed > proneuralnt_neurod2_pileup_2_Neurod2_correlation_sorted_summits_index_exteded_500bp.bed
  bedtools intersect -wa -a control_lambda_2_Neurod2_correlation_sorted_summits_index_exteded_500bp.bed -b proneural_table_neurod2.bed > neurod2_control_lambda_2_Neurod2_correlation_sorted_summits_index_exteded_500bp.bed
  bedtools intersect -wa -a control_lambda_2_Neurod2_correlation_sorted_summits_index_exteded_500bp.bed -b proneural_table_proneuralnt.bed > proneuralnt_control_lambda_2_Neurod2_correlation_sorted_summits_index_exteded_500bp.bed
  
  cd ~/openchr/noack
")
setwd("~/openchr/noack")


library(dplyr)
neurod2_neurod2_neurod2<-read.table("neurod2_neurod2_pileup_1_Neurod2_correlation_sorted_summits_index_exteded_500bp.bed") %>% select(V4,V5)
colnames(neurod2_neurod2_neurod2)<-c("peak_id","signal_1")
neurod2_neurod2_neurod2$factor<-"neurod2 pileup"
neurod2_lambda_neurod2<-read.table("neurod2_control_lambda_1_Neurod2_correlation_sorted_summits_index_exteded_500bp.bed") %>% select(V4,V5)
colnames(neurod2_lambda_neurod2)<-c("peak_id","signal_1")
neurod2_lambda_neurod2$factor<-"control lambda"
join_1<-rbind(neurod2_neurod2_neurod2,neurod2_lambda_neurod2)

neurod2_neurod2_neurod2<-read.table("neurod2_neurod2_pileup_2_Neurod2_correlation_sorted_summits_index_exteded_500bp.bed") %>% select(V4,V5)
colnames(neurod2_neurod2_neurod2)<-c("peak_id","signal_2")
neurod2_neurod2_neurod2$factor<-"neurod2 pileup"
neurod2_lambda_neurod2<-read.table("neurod2_control_lambda_2_Neurod2_correlation_sorted_summits_index_exteded_500bp.bed") %>% select(V4,V5)
colnames(neurod2_lambda_neurod2)<-c("peak_id","signal_2")
neurod2_lambda_neurod2$factor<-"control lambda"
join_2<-rbind(neurod2_neurod2_neurod2,neurod2_lambda_neurod2)

join<-merge(join_1,join_2, by=c("peak_id","factor") )

join$decil<-NA
join$signal<-rowMeans(cbind(join$signal_1,join$signal_2) )

npeaks<-length(read.table("Neurod2_correlation_sorted_summits_index_exteded_500bp.bed")[,1])
decil<-npeaks/10
i<-1
j<-1
while (i<npeaks) {
  print(i)
  join$decil[which(join$peak_id %in% c( i:(i+decil) ) ) ]<-j
  j<-j+1
  i<-i+decil
}



library(ggplot2)
join$decil<-as.character(join$decil)
join$decil<-factor(join$decil, levels= c (1:10))
pdf("~/neurod2_peaks_signal.pdf", useDingbats = F, height = 4, width = 7)
ggplot(join, aes(x=decil, y=signal, fill=factor)) +
  geom_boxplot() + theme_classic() + ylim(0,60) + ggtitle("Neurod2 peaks")
dev.off()

