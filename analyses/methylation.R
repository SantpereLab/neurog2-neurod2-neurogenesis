# METHYLATION DELTAS BARPLOTS ########################
#####################################################


# Obtain the bedgraph files from the ATAC-seq bam files to do map the ATAC coverage to the genomic regions.

system ("
  cd ~/p/openchr/noack
  cd cluster_bams_imputation

  module load SAMtools
  # Merge proj
  sbatch --partition=long --wrap="samtools merge clusterPN.bam clusterPN1.bam clusterPN2.bam clusterPN3.bam"
  module load Miniconda3/4.9.2
  sbatch --partition=normal --mem=150000 --wrap="samtools index clusterPN.bam; wait; bamCoverage -b clusterPN.bam -o clusterPN.bam.bw --normalizeUsing RPKM"
  
  cd ~/openchr/noack
  import_shiva p/openchr/noack/cluster_bams_imputation/clusterPN.bam.bw .
  bigWigToBedGraph clusterNSC.bam.bw clusterNSC.bam.bdg
  bigWigToBedGraph clusterIPC.bam.bw clusterIPC.bam.bdg 
  bigWigToBedGraph clusterPN.bam.bw clusterPN.bam.bdg
  cat clusterNSC.bam.bdg | grep -v JH | grep -v chrUn | grep -v GL | grep -v random > clean.clusterNSC.bam.bdg
  cat clusterIPC.bam.bdg | grep -v JH | grep -v chrUn | grep -v GL | grep -v random > clean.clusterIPC.bam.bdg
  cat clusterPN.bam.bdg | grep -v JH | grep -v chrUn | grep -v GL | grep -v random > clean.clusterPN.bam.bdg
")

# Fix the bigwig files to be integer.

nsc_bdg<-read.table("~/openchr/noack/clean.clusterNSC.bam.bdg")
ipc_bdg<-read.table("~/openchr/noack/clean.clusterIPC.bam.bdg")
pn_bdg<-read.table("~/openchr/noack/clean.clusterPN.bam.bdg")

nsc_bdg$V4<-as.integer(nsc_bdg$V4*1000)
ipc_bdg$V4<-as.integer(ipc_bdg$V4*1000)
pn_bdg$V4<-as.integer(pn_bdg$V4*1000)

write.table (nsc_bdg, "~/openchr/noack/clean.clusterNSC_fixed.bam.bdg", col.names = F, row.names = F, quote = F, sep = "\t" )
write.table (ipc_bdg, "~/openchr/noack/clean.clusterIPC_fixed.bam.bdg", col.names = F, row.names = F, quote = F, sep = "\t" )
write.table (pn_bdg, "~/openchr/noack/clean.clusterPN_fixed.bam.bdg", col.names = F, row.names = F, quote = F, sep = "\t" )


# Subdivide the motifs by the intersection with the different proneural factors.

system ("
  cd ~/openchr/noack
  for file in private_neurod2.bed private_neurog2.bed shared.bed
  do
    for file2 in CA[ACTG]-CA[ACTG].bed
    do
      bedtools intersect -wa -a $file2 -b $file > $file2.$file
    done
  done
")


# Print the E-boxes that are present in ATAC-seq peaks but that do not overlap with any proneural factor.
atac_motifs<-readRDS("~/openchr/noack/atac_peaks_motifs.Rds")
atac_overlap_proneural<-read.table("~/openchr/noack/proneural_table.bed")
library(bedr)
motifs<-bedr.join.region(bedr.sort.region(atac_motifs),bedr.sort.region(atac_overlap_proneural) )
atac_motifs_proneuralnt<-motifs %>% filter(V4=="1.proneuralnt")
options(scipen = 999)
for (motif in unique(atac_motifs_proneuralnt$central_dinucleotide)) {
  motifs<-atac_motifs_proneuralnt %>% dplyr::filter(central_dinucleotide==motif)
  motifs$ebox_end<-as.numeric(motifs$ebox_start)+5
  motifs$ebox_start<-motifs$ebox_end-1
  motifs<-motifs %>% dplyr::select(chr,ebox_start,ebox_end)
  write.table(motifs,paste0("~/openchr/noack/proneuralnt",motif,".bed"),sep="\t",col.names = F,row.names = F, quote = F, )
}


# Calculate the deltas of the motifs.

system("
  # Deltas of the methylaiton PN - NSC.
  module load BEDTools
  for file in CA*.bed.private_neurod2.bed proneuralntCA*bed
  do
    echo $file
    # Extend a window around the motifs and compute the mean of the methylation around that window.
    cat $file | awk '{print $1,$2-25,$3+25}' | tr ' ' '\t' | bedtools map -a - -b ~/p/methylation/noack/2021/merged_NSC_CpGmeth.bedGraph -c 4 -o mean > NSC_methylation_${file}
    cat $file | awk '{print $1,$2-25,$3+25}' | tr ' ' '\t' | bedtools map -a - -b ~/p/methylation/noack/2021/merged_PN_CpGmeth.bedGraph -c 4 -o mean > PN_methylation_${file}
  done

  # Deltas of the ATAC.
  for file in CA*.bed.private_neurod2.bed proneuralntCA*bed
  do
    echo $file
    # Extend a window around the motifs and compute the mean of the ATAC-seq around that window.
    cat $file | awk '{print $1,$2-25,$3+25}' | tr ' ' '\t' | bedtools map -a - -b clean.clusterNSC.bam.bdg -c 4 -o mean > NSC_atac_${file}
    cat $file | awk '{print $1,$2-25,$3+25}' | tr ' ' '\t' | bedtools map -a - -b clean.clusterPN.bam.bdg -c 4 -o mean > PN_atac_${file}
  done
")

# Read bed files with mapped values and compute deltas.
delta_privii<-c()
delta_neuralnt<-c()
motivos<-c("CAT-CAT","CAT-CAG","CAG-CAG","CAG-CAC","CAT-CAC","CAC-CAC","CAA-CAA","CAA-CAG","CAA-CAT","CAA-CAC")
for (motif in motivos) {
  nsc_methylation_motif_privii<-read.table( paste0 ("~/openchr/noack/NSC_methylation_",motif,".bed.private_neurod2.bed") ) %>% filter(V4!=".")
  pn_methylation_motif_privii<-read.table( paste0 ("~/openchr/noack/PN_methylation_",motif,".bed.private_neurod2.bed") ) %>% filter(V4!=".")
  delta_motif_privii<-mean(as.numeric(pn_methylation_motif_privii$V4) )-mean(as.numeric(nsc_methylation_motif_privii$V4) )
  delta_privii<-c(delta_privii,delta_motif_privii)
  
  nsc_methylation_motif_neuralnt<-read.table( paste0 ("~/openchr/noack/NSC_methylation_proneuralnt",motif,".bed") ) %>% filter(V4!=".")
  pn_methylation_motif_neuralnt<-read.table( paste0 ("~/openchr/noack/PN_methylation_proneuralnt",motif,".bed") ) %>% filter(V4!=".")
  delta_motif_neuralnt<-mean(as.numeric(pn_methylation_motif_neuralnt$V4) )-mean(as.numeric(nsc_methylation_motif_neuralnt$V4) )
  delta_neuralnt<-c(delta_neuralnt,delta_motif_neuralnt)
}
names(delta_privii)<-motivos
names(delta_neuralnt)<-motivos


# The same with the ATAC.

delta_privii_atac<-c()
delta_neuralnt_atac<-c()
motivos<-c("CAT-CAT","CAT-CAG","CAG-CAG","CAG-CAC","CAT-CAC","CAC-CAC","CAA-CAA","CAA-CAG","CAA-CAT","CAA-CAC")
for (motif in motivos) {
  nsc_methylation_motif_privii<-read.table( paste0 ("~/openchr/noack/NSC_atac_",motif,".bed.private_neurod2.bed") ) %>% filter(V4!=".")
  pn_methylation_motif_privii<-read.table( paste0 ("~/openchr/noack/PN_atac_",motif,".bed.private_neurod2.bed") ) %>% filter(V4!=".")
  delta_motif_privii<-mean(as.numeric(pn_methylation_motif_privii$V4) )-mean(as.numeric(nsc_methylation_motif_privii$V4) )
  delta_privii_atac<-c(delta_privii_atac,delta_motif_privii)
  
  nsc_methylation_motif_neuralnt<-read.table( paste0 ("~/openchr/noack/NSC_atac_proneuralnt",motif,".bed") ) %>% filter(V4!=".")
  pn_methylation_motif_neuralnt<-read.table( paste0 ("~/openchr/noack/PN_atac_proneuralnt",motif,".bed") ) %>% filter(V4!=".")
  delta_motif_neuralnt<-mean(as.numeric(pn_methylation_motif_neuralnt$V4) )-mean(as.numeric(nsc_methylation_motif_neuralnt$V4) )
  delta_neuralnt_atac<-c(delta_neuralnt_atac,delta_motif_neuralnt)
}
names(delta_privii_atac)<-motivos
names(delta_neuralnt_atac)<-motivos

pdf("delta_barplots.pdf",height = 5, width = 14)
barplot(delta_privii, ylim = c(-40,0) )
barplot(delta_neuralnt, ylim = c(-40,0) )

barplot(delta_privii_atac, ylim = c(-5,30) )
barplot(delta_neuralnt_atac, ylim = c(-5,30) )
dev.off()

library(GenomicFeatures)
library(seqplots)


# The same but subdividing in promoters and enhancers.

system ("

  # Subdivide the motifs in promoters and enhancers. 
  module load BEDTools
  for file in CA*.bed.private_neurod2.bed proneuralntCA*bed
  do
    bedtools intersect -a $file -b ~/p/genomes/gencode.vM24.annotation.filtered.promoters.bed | sort -k1,1 -k2,2n > promoters.$file
    bedtools intersect -v -a $file -b ~/p/genomes/gencode.vM24.annotation.filtered.promoters.bed | sort -k1,1 -k2,2n > enhancers.$file
  done
  
  for file in promoters.CA*private_neurod2.bed enhancers.CA*private_neurod2.bed enhancers.proneuralntCA*bed promoters.proneuralntCA*bed
  do
    echo $file
    # Extend a window around the motifs and compute the mean of the ATAC-seq around that window.
    cat $file | awk '{print $1,$2-25,$3+25}' | tr ' ' '\t' | bedtools map -a - -b ~/p/methylation/noack/2021/merged_NSC_CpGmeth.bedGraph -c 4 -o mean > NSC_methylation_${file}
    cat $file | awk '{print $1,$2-25,$3+25}' | tr ' ' '\t' | bedtools map -a - -b ~/p/methylation/noack/2021/merged_PN_CpGmeth.bedGraph -c 4 -o mean > PN_methylation_${file}
  done
  
  
  # The same with the ATAC.
  cd ~/openchr/noack
  for file in promoters.CA*private_neurod2.bed enhancers.CA*private_neurod2.bed enhancers.proneuralntCA*bed promoters.proneuralntCA*bed
  do
    echo $file
    # Extend a window around the motifs and compute the mean of the ATAC-seq around that window.
    cat $file | awk '{print $1,$2-25,$3+25}' | tr ' ' '\t' | bedtools map -a - -b clean.clusterNSC.bam.bdg -c 4 -o mean > NSC_atac_${file}
    cat $file | awk '{print $1,$2-25,$3+25}' | tr ' ' '\t' | bedtools map -a - -b clean.clusterPN.bam.bdg -c 4 -o mean > PN_atac_${file}
  done

")


######## los plots

delta_privii<-c()
delta_neuralnt<-c()
motivos<-c("CAT-CAT","CAT-CAG","CAG-CAG","CAG-CAC","CAT-CAC","CAC-CAC","CAA-CAA","CAA-CAG","CAA-CAT","CAA-CAC")
for (motif in motivos) {
  nsc_methylation_motif_privii<-read.table( paste0 ("~/openchr/noack/NSC_methylation_promoters.",motif,".bed.private_neurod2.bed") ) %>% filter(V4!=".")
  pn_methylation_motif_privii<-read.table( paste0 ("~/openchr/noack/PN_methylation_promoters.",motif,".bed.private_neurod2.bed") ) %>% filter(V4!=".")
  delta_motif_privii<-mean(as.numeric(pn_methylation_motif_privii$V4) )-mean(as.numeric(nsc_methylation_motif_privii$V4) )
  delta_privii<-c(delta_privii,delta_motif_privii)
  
  nsc_methylation_motif_neuralnt<-read.table( paste0 ("~/openchr/noack/NSC_methylation_promoters.proneuralnt",motif,".bed") ) %>% filter(V4!=".")
  pn_methylation_motif_neuralnt<-read.table( paste0 ("~/openchr/noack/PN_methylation_promoters.proneuralnt",motif,".bed") ) %>% filter(V4!=".")
  delta_motif_neuralnt<-mean(as.numeric(pn_methylation_motif_neuralnt$V4) )-mean(as.numeric(nsc_methylation_motif_neuralnt$V4) )
  delta_neuralnt<-c(delta_neuralnt,delta_motif_neuralnt)
}
names(delta_privii)<-motivos
names(delta_neuralnt)<-motivos

#lo mesmo con el atac

delta_privii_atac<-c()
delta_neuralnt_atac<-c()
motivos<-c("CAT-CAT","CAT-CAG","CAG-CAG","CAG-CAC","CAT-CAC","CAC-CAC","CAA-CAA","CAA-CAG","CAA-CAT","CAA-CAC")
for (motif in motivos) {
  nsc_methylation_motif_privii<-read.table( paste0 ("~/openchr/noack/NSC_atac_promoters.",motif,".bed.private_neurod2.bed") ) %>% filter(V4!=".")
  pn_methylation_motif_privii<-read.table( paste0 ("~/openchr/noack/PN_atac_promoters.",motif,".bed.private_neurod2.bed") ) %>% filter(V4!=".")
  delta_motif_privii<-mean(as.numeric(pn_methylation_motif_privii$V4) )-mean(as.numeric(nsc_methylation_motif_privii$V4) )
  delta_privii_atac<-c(delta_privii_atac,delta_motif_privii)
  
  nsc_methylation_motif_neuralnt<-read.table( paste0 ("~/openchr/noack/NSC_atac_promoters.proneuralnt",motif,".bed") ) %>% filter(V4!=".")
  pn_methylation_motif_neuralnt<-read.table( paste0 ("~/openchr/noack/PN_atac_promoters.proneuralnt",motif,".bed") ) %>% filter(V4!=".")
  delta_motif_neuralnt<-mean(as.numeric(pn_methylation_motif_neuralnt$V4) )-mean(as.numeric(nsc_methylation_motif_neuralnt$V4) )
  delta_neuralnt_atac<-c(delta_neuralnt_atac,delta_motif_neuralnt)
}
names(delta_privii_atac)<-motivos
names(delta_neuralnt_atac)<-motivos

pdf("delta_barplots_promoters.pdf",height = 5, width = 14)
barplot(delta_privii, ylim = c(-40,0), main = "Neurod2 private meth")
barplot(delta_neuralnt, ylim = c(-40,0), main = "No proneural meth" )

barplot(delta_privii_atac, ylim = c(-5,40), main = "Neurod2 private atac" )
barplot(delta_neuralnt_atac, ylim = c(-5,40), main = "No proneural atac"  )
dev.off()



############ METHYLATION #####################
##############################################


# Methylation around the motifs with the Noack et al 2022 data.
tracks <- paste0( c('GSM4710397_NSC','GSM4710400_IPC','GSM4710403_PN'), '_rep1_CpGmeth.bedGraph.bw' )

plotmeth<-function(meth_matrix,main) {
  return(plotAverage(plotset=meth_matrix, labels = c('NSC','IPC','PN'), xlim = c(-750,750),
                     ylim = c(0,100), main = main, xlab = "Dist to motif (bp)", ylab = '% CpG Methylation',
                     plotScale = "linear", type = "full", error.estimates = T, yaxs='i',xaxs='i',
                     legend = TRUE, legend_ext = FALSE, legend_pos = 'bottomright', colvec = c("cadetblue1","deepskyblue3","navy"),
                     legend_ext_pos = "topleft", cex.axis = 14, cex.lab = 16,xaxt='n',cex.legend = 12, pointsize = 12) )
}





# Function to print a lineplot of the average methylation around a genomic point.
plotmeth<-function(meth_matrix,main) {
  return(plotAverage(plotset=meth_matrix, labels = c('NSC','IPC','PN'), xlim = c(-750,750),
                     ylim = c(0,100000), main = main, xlab = "Dist to motif (bp)", ylab = 'ATAC-seq RPKM',
                     plotScale = "linear", type = "full", error.estimates = T, yaxs='i',xaxs='i',
                     legend = TRUE, legend_ext = FALSE, legend_pos = 'bottomright', colvec = c("cadetblue1","deepskyblue3","navy"),
                     legend_ext_pos = "topleft", cex.axis = 14, cex.lab = 16,xaxt='n',cex.legend = 12, pointsize = 12) )
}


# NEUROD2 private motifs.
for (motif in c( "CAT-CAT", "CAT-CAG", "CAG-CAG", "CAG-CAC", "CAT-CAC" )) {
  meth<-getPlotSetArray(tracks=tracks,features=paste0("~/openchr/noack/",motif,".bed.private_neurod2.bed"),
                        refgenome='mm10',type = 'mf',add_heatmap=F,xmin=750,xmax=750,bin = 10)
  pdf( paste0(motif,"_neurod2_meth.pdf"), useDingbats = F, width = 8, height = 5)
  plotmeth(meth,main=motif)
  dev.off()
}


# Methylation in genomic regions stratified by intersection with proneural factors.

# Stratify the summits of the ATAC-seq peaks by intersection with different proneural factors. 
system("
  cd ~/openchr/noack
  awk '{print $1,$2+250,$3-250,$4}' proneural_table.bed | tr " " "\t" > proneural_table_summits.bed
  cat proneural_table_summits.bed | grep 2.private_neurog2 > proneural_table_summits_neurog2.bed
  cat proneural_table_summits.bed | grep 1.proneuralnt > proneural_table_summits_proneuralnt.bed
  cat proneural_table_summits.bed | grep 3.shared > proneural_table_summits_shared.bed
  cat proneural_table_summits.bed | grep 4.private_neurod2 > proneural_table_summits_neurod2.bed
")


plotmeth<-function(meth_matrix,main) {
  return(plotAverage(plotset=meth_matrix, labels = c('NSC','IPC','PN'), xlim = c(-2000,2000),
                     ylim = c(0,100), main = main, xlab = "Dist to summit", ylab = '% CpG Methylation',
                     plotScale = "linear", type = "full", error.estimates = T, yaxs='i',xaxs='i',
                     legend = TRUE, legend_ext = FALSE, legend_pos = 'bottomright', colvec = c("cadetblue1","deepskyblue3","navy"),
                     legend_ext_pos = "topleft", cex.axis = 14, cex.lab = 16,xaxt='n',cex.legend = 12, pointsize = 12) )
}

pdf("~/atac_centered_methylation_seqplots.pdf")
plotmeth(getPlotSetArray(tracks=tracks,features='~/openchr/noack/proneural_table_summits_neurog2.bed',refgenome='mm10',
                         type = 'mf',add_heatmap=F,xmin=2000,xmax=2000,bin = 50),main="Neurog2 private")
plotmeth(getPlotSetArray(tracks=tracks,features='~/openchr/noack/proneural_table_summits_proneuralnt.bed',refgenome='mm10',
                         type = 'mf',add_heatmap=F,xmin=2000,xmax=2000,bin = 50),main="No Neurog2/Neurod2")
plotmeth(getPlotSetArray(tracks=tracks,features='~/openchr/noack/proneural_table_summits_shared.bed',refgenome='mm10',
                         type = 'mf',add_heatmap=F,xmin=2000,xmax=750,bin = 50),main="Shared")
plotmeth(getPlotSetArray(tracks=tracks,features='~/openchr/noack/proneural_table_summits_neurod2.bed',refgenome='mm10',
                         type = 'mf',add_heatmap=F,xmin=2000,xmax=2000,bin = 50),main="Neurod2 private")
dev.off()


# METHYLATION IN THE TRAJECTORIES ##########
############################################

# Stratify the trajectories based on the intersection with NEUROD2 peaks.

system ("
  for file in enhancer_early enhancer_invariant enhancer_neuronal enhancer_transient
  do
    cat $file.neurod2_proneural_table.bed $file.sharedd_proneural_table.bed | sort -k1,1 -k2,2n > $file.all_neurod2.bed
  done
  
  for file in promoter_early promoter_invariant promoter_neuronal promoter_transient
  do
    cat $file.neurod2_proneural_table.bed $file.sharedd_proneural_table.bed | sort -k1,1 -k2,2n > $file.all_neurod2.bed
  done
  
  for file in promoter_early promoter_invariant promoter_neuronal promoter_transient
  do
    cat $file.proneuralnt_proneural_table.bed $file.neurog2_proneural_table.bed | sort -k1,1 -k2,2n > $file.neurod2nt.bed
  done
  
  for file in enhancer_early enhancer_invariant enhancer_neuronal enhancer_transient
  do
    cat $file.proneuralnt_proneural_table.bed $file.neurog2_proneural_table.bed | sort -k1,1 -k2,2n > $file.neurod2nt.bed
  done
")

# Print the lineplots of the different subsets of regions.

pdf("~/methylation_supplementary.pdf")
plotmeth(getPlotSetArray(tracks=tracks,features='~/openchr/noack/promoter_invariant.all_neurod2.bed',refgenome='mm10',
                         type = 'mf',add_heatmap=F,xmin=2000,xmax=2000,bin = 50),main="Neurog2 private" )

plotmeth(getPlotSetArray(tracks=tracks,features='~/openchr/noack/promoter_early.all_neurod2.bed',refgenome='mm10',
                         type = 'mf',add_heatmap=F,xmin=2000,xmax=2000,bin = 50),main="Neurog2 private" )

plotmeth(getPlotSetArray(tracks=tracks,features='~/openchr/noack/promoter_transient.all_neurod2.bed',refgenome='mm10',
                         type = 'mf',add_heatmap=F,xmin=2000,xmax=2000,bin = 50),main="Neurog2 private" )

plotmeth(getPlotSetArray(tracks=tracks,features='~/openchr/noack/promoter_neuronal.all_neurod2.bed',refgenome='mm10',
                         type = 'mf',add_heatmap=F,xmin=2000,xmax=2000,bin = 50),main="Neurog2 private" )

plotmeth(getPlotSetArray(tracks=tracks,features='~/openchr/noack/enhancer_invariant.all_neurod2.bed',refgenome='mm10',
                         type = 'mf',add_heatmap=F,xmin=2000,xmax=2000,bin = 50),main="Neurog2 private" )

plotmeth(getPlotSetArray(tracks=tracks,features='~/openchr/noack/enhancer_early.all_neurod2.bed',refgenome='mm10',
                         type = 'mf',add_heatmap=F,xmin=2000,xmax=2000,bin = 50),main="Neurog2 private" )

plotmeth(getPlotSetArray(tracks=tracks,features='~/openchr/noack/enhancer_transient.all_neurod2.bed',refgenome='mm10',
                         type = 'mf',add_heatmap=F,xmin=2000,xmax=2000,bin = 50),main="Neurog2 private" )

plotmeth(getPlotSetArray(tracks=tracks,features='~/openchr/noack/enhancer_neuronal.all_neurod2.bed',refgenome='mm10',
                         type = 'mf',add_heatmap=F,xmin=2000,xmax=2000,bin = 50),main="Neurog2 private" )

plotmeth(getPlotSetArray(tracks=tracks,features='~/openchr/noack/promoter_invariant.all_neurod2.bed',refgenome='mm10',
                         type = 'mf',add_heatmap=F,xmin=2000,xmax=2000,bin = 50),main="Neurog2 private" )

plotmeth(getPlotSetArray(tracks=tracks,features='~/openchr/noack/promoter_early.neurod2nt.bed',refgenome='mm10',
                         type = 'mf',add_heatmap=F,xmin=2000,xmax=2000,bin = 50),main="Neurog2 private" )

plotmeth(getPlotSetArray(tracks=tracks,features='~/openchr/noack/promoter_transient.neurod2nt.bed',refgenome='mm10',
                         type = 'mf',add_heatmap=F,xmin=2000,xmax=2000,bin = 50),main="Neurog2 private" )

plotmeth(getPlotSetArray(tracks=tracks,features='~/openchr/noack/promoter_neuronal.neurod2nt.bed',refgenome='mm10',
                         type = 'mf',add_heatmap=F,xmin=2000,xmax=2000,bin = 50),main="Neurog2 private" )

plotmeth(getPlotSetArray(tracks=tracks,features='~/openchr/noack/enhancer_invariant.neurod2nt.bed',refgenome='mm10',
                         type = 'mf',add_heatmap=F,xmin=2000,xmax=2000,bin = 50),main="Neurog2 private" )

plotmeth(getPlotSetArray(tracks=tracks,features='~/openchr/noack/enhancer_early.neurod2nt.bed',refgenome='mm10',
                         type = 'mf',add_heatmap=F,xmin=2000,xmax=2000,bin = 50),main="Neurog2 private" )

plotmeth(getPlotSetArray(tracks=tracks,features='~/openchr/noack/enhancer_transient.neurod2nt.bed',refgenome='mm10',
                         type = 'mf',add_heatmap=F,xmin=2000,xmax=2000,bin = 50),main="Neurog2 private" )

plotmeth(getPlotSetArray(tracks=tracks,features='~/openchr/noack/enhancer_neuronal.neurod2nt.bed',refgenome='mm10',
                         type = 'mf',add_heatmap=F,xmin=2000,xmax=2000,bin = 50),main="Neurog2 private" )
dev.off()




############# METHYLATION HEATMAPS #####
########################################


# Merge the replicates to do the heatmap.
system ("
  module load BEDTools
  cat GSM4710397_NSC_rep1_CpGmeth.bedGraph GSM4710399_NSC_rep3_CpGmeth.bedGraph GSM4710398_NSC_rep2_CpGmeth.bedGraph | sort -k1,1 -k2,2n |\
  bedtools merge -c 4 -o mean -i - | sort -k1,1 -k2,2n > merged_NSC_CpGmeth.bedGraph
  cat GSM4710400_IPC_rep1_CpGmeth.bedGraph GSM4710401_IPC_rep2_CpGmeth.bedGraph GSM4710402_IPC_rep3_CpGmeth.bedGraph | sort -k1,1 -k2,2n |\
  bedtools merge -c 4 -o mean -i - | sort -k1,1 -k2,2n > merged_IPC_CpGmeth.bedGraph
  cat GSM4710403_PN_rep1_CpGmeth.bedGraph GSM4710404_PN_rep2_CpGmeth.bedGraph GSM4710405_PN_rep3_CpGmeth.bedGraph | sort -k1,1 -k2,2n |\
  bedtools merge -c 4 -o mean -i - | sort -k1,1 -k2,2n > merged_PN_CpGmeth.bedGraph
")


system ("
  for file in extended*bdg
  do
    echo $file
    bedGraphToBigWig $file ~/genomes/mm10.genome $file.bw
  done
")


# Compute matrices of genomic coverage files around genomic regions, where each row is a genomic region.
matt<-function(methylation_file,tss) {
  library(GenomicRanges)
  library(dplyr)
  meth<-read.table(methylation_file) %>% select(V1,V2,V3,V4)
  colnames(meth)<-c("chr","start","end","methylation")
  meth<-makeGRangesFromDataFrame(meth,keep.extra.columns=T)
  
  # If there are a lot of rows in the matrix, downsample, because the output heatmap gives an error.
  if (length(tss$start)>50000) {
    rand_indices<-sample(1:length(tss$start), 50000, replace=F)
    tss<-tss[rand_indices,]
  }
  
  tss<-makeGRangesFromDataFrame(tss)
  library(EnrichedHeatmap)
  # Use the parameters suggested in the vignettes of the package: https://bioconductor.org/packages/devel/bioc/vignettes/EnrichedHeatmap/inst/doc/EnrichedHeatmap.html
  mat2 = normalizeToMatrix(meth, tss, value_column = "methylation", mean_mode = "absolute",
                           extend = 2000, w = 50, background = NA, smooth = T)
  return(mat2)
}

# Output a heatmap from a matrix input.
plotttt<-function(mat2, orderrr) {
  library(circlize)
  meth_col_fun = colorRamp2(c(0,50,100), c("white", "lightblue","cornflowerblue"))
  
  axisname=c("-2kb", "summit","","2kb")
  legend_params=list(col_fun = meth_col_fun, title = NULL, border = "black")
  p1<-EnrichedHeatmap(mat2, row_order = orderrr, col = meth_col_fun, top_annotation=NULL, pos_line=F, 
                      axis_name= axisname, heatmap_legend_param = legend_params, use_raster=F)
  return(p1)
}

# Output a heatmap from a matrix input, but splitting the heatmap in promoters and enhancers.
plotttt_split<-function(mat2, orderrr, spliteo) {
  library(circlize)
  meth_col_fun = colorRamp2(c(0,50,100), c("white", "lightblue","cornflowerblue"))
  
  axisname=c("-2kb", "summit","","2kb")
  legend_params=list(col_fun = meth_col_fun, title = NULL, border = "black")
  p1<-EnrichedHeatmap(mat2, row_order = orderrr, col = meth_col_fun, top_annotation=NULL, pos_line=F, 
                      axis_name= axisname, row_split=spliteo, heatmap_legend_param = legend_params, use_raster=F)
  return(p1)
}

# Function to plot together the heatmaps of the methylation in the 3 cell populations, without sorting the regions.
heatmethh_nosort<-function(summits_file) {
  tss<-read.table(summits_file) %>% select(V1,V2,V3)
  colnames(tss)<-c("chr","start","end")
  
  # Compute the matrices.
  mat1<-matt("~/p/methylation/noack/2021/GSM4710397_NSC_rep1_CpGmeth.bedGraph", tss)
  mat2<-matt("~/p/methylation/noack/2021/GSM4710400_IPC_rep1_CpGmeth.bedGraph", tss)
  mat3<-matt("~/p/methylation/noack/2021/GSM4710403_PN_rep1_CpGmeth.bedGraph", tss)
  
  # Dont change the order of the rows.
  order_default<-c(1:length(tss[,1]))
  p1<-plotttt(mat1,order_default)
  p2<-plotttt(mat2,order_default)
  p3<-plotttt(mat3,order_default)
  p<-p1+p2+p3
  pdf(paste0(summits_file,"_meth_w50.pdf"), width = 7, height = 15)
  print(p)
  dev.off()
}

# Function to plot together the heatmaps of the methylation in the 3 cell populations, subdivided by promoters and enhancers.
heatmethh_proportional<-function(promoters,enhancers, all) {
  require(dplyr)
  tss<-read.table(all) %>% select(V1,V2,V3)
  colnames(tss)<-c("chr","start","end")
  
  promoters<-read.table(promoters)
  enhancers<-read.table(enhancers)
  promenh<-c( rep("a promoter",length(promoters[,1]) ), rep("enhancer",length(enhancers[,1]) )  )
  
  # If there are more than 50000 regions, pick random 50000 regions. Because the heatmap would give an error.
  if (length(tss$start)>50000) {
    rand_indices<-sample(1:length(tss$start), 50000, replace=F)
    tss<-tss[rand_indices,]
    promenh<-promenh[rand_indices]
  }
  
  # Compute the matrices with the different methylation files.
  mat1<-matt("~/p/methylation/noack/2021/GSM4710397_NSC_rep1_CpGmeth.bedGraph", tss)
  mat2<-matt("~/p/methylation/noack/2021/GSM4710400_IPC_rep1_CpGmeth.bedGraph", tss)
  mat3<-matt("~/p/methylation/noack/2021/GSM4710403_PN_rep1_CpGmeth.bedGraph", tss)
  
  # Sort the peaks with the methylation averaged across the three replicates.
  order1<-enriched_score(mat1)
  order2<-enriched_score(mat2)
  order3<-enriched_score(mat3)
  mean_order<-order( (order1+order2+order3)/3, decreasing=F)
  
  p1<-plotttt_split(mat1, mean_order, promenh)
  p2<-plotttt_split(mat2, mean_order, promenh)
  p3<-plotttt_split(mat3, mean_order, promenh)
  p<-p1+p2+p3
  
  pdf(paste0(all,"_meth_w50.pdf"), width = 7, height = 7)
  print(p)
  dev.off()
}



# Subdivide the regions in promoters and enhancers.

system("
  module load BEDTools
  cd ~/p/openchr/noack
  bedtools intersect -a proneural_table_summits_neurod2.bed -b footprinting/promoters.bed > promoters_proneural_table_summits_neurod2.bed
  bedtools intersect -a proneural_table_summits_neurod2.bed -b footprinting/enhancers.bed > enhancers_proneural_table_summits_neurod2.bed
  bedtools intersect -a proneural_table_summits_neurog2.bed -b footprinting/promoters.bed > promoters_proneural_table_summits_neurog2.bed
  bedtools intersect -a proneural_table_summits_neurog2.bed -b footprinting/enhancers.bed > enhancers_proneural_table_summits_neurog2.bed
  bedtools intersect -a proneural_table_summits_shared.bed -b footprinting/promoters.bed > promoters_proneural_table_summits_shared.bed
  bedtools intersect -a proneural_table_summits_shared.bed -b footprinting/enhancers.bed > enhancers_proneural_table_summits_shared.bed
  bedtools intersect -a proneural_table_summits_proneuralnt.bed -b footprinting/promoters.bed > promoters_proneural_table_summits_proneuralnt.bed
  bedtools intersect -a proneural_table_summits_proneuralnt.bed -b footprinting/enhancers.bed > enhancers_proneural_table_summits_proneuralnt.bed
  
  cat promoters_proneural_table_summits_neurod2.bed enhancers_proneural_table_summits_neurod2.bed > all_proneural_table_summits_neurod2.bed
  cat promoters_proneural_table_summits_neurog2.bed enhancers_proneural_table_summits_neurog2.bed > all_proneural_table_summits_neurog2.bed
  cat promoters_proneural_table_summits_shared.bed enhancers_proneural_table_summits_shared.bed > all_proneural_table_summits_shared.bed
  cat promoters_proneural_table_summits_proneuralnt.bed enhancers_proneural_table_summits_proneuralnt.bed > all_proneural_table_summits_proneuralnt.bed
")

heatmethh_proportional(promoters = "~/p/openchr/noack/promoters_proneural_table_summits_neurod2.bed", 
                       enhancers = "~/p/openchr/noack/enhancers_proneural_table_summits_neurod2.bed" , all = "~/p/openchr/noack/all_proneural_table_summits_neurod2.bed")
heatmethh_proportional(promoters = "~/p/openchr/noack/promoters_proneural_table_summits_neurog2.bed", 
                       enhancers = "~/p/openchr/noack/enhancers_proneural_table_summits_neurog2.bed", all = "~/p/openchr/noack/all_proneural_table_summits_neurog2.bed")
heatmethh_proportional(promoters = "~/p/openchr/noack/promoters_proneural_table_summits_shared.bed", 
                       enhancers = "~/p/openchr/noack/enhancers_proneural_table_summits_shared.bed", all = "~/p/openchr/noack/all_proneural_table_summits_shared.bed")
heatmethh_proportional(promoters = "~/p/openchr/noack/promoters_proneural_table_summits_proneuralnt.bed", 
                       enhancers = "~/p/openchr/noack/enhancers_proneural_table_summits_proneuralnt.bed", all = "~/p/openchr/noack/all_proneural_table_summits_proneuralnt.bed")


# Heatmap of NEUROD2 peaks sorted according to the correlation of their chromatin accessibility with the expression of Neurod2.
system("
  bedtools intersect -wa -a Neurod2_correlation_sorted_summits.bed -b ~/p/neurod2/chipseq/e14.5.neurod2_highconfsubpeaks_summit.bed > neurod2intersect_Neurod2_correlation_sorted_summits.bed
")

heatmethh_nosort("~/p/openchr/noack/neurod2intersect_Neurod2_correlation_sorted_summits.bed")



####### POSTNATAL METHYLATION ########
######################################

# Liu 2021
#convertir esos tsv en bedgraph y luego en bw files
system("
cd ~/methylation/liu
for file in *tsv; do cat $file | awk '{print $1,$2-1,$2,$5/$6*100}' | sort -k1,1 -k2,2n | 
tr ' ' '\t' | tr ',' '\t' | cut -f1,2,3,4 > $file.bdg;
bedGraphToBigWig $file.bdg ~/genomes/mm10.genome $file.bw; done")

system("
  cat proneural_table_summits.bed | grep proneuralnt > proneuralnt_summits.bed
  cat proneural_table_summits.bed | grep -v proneuralnt | grep -v neurog2 > neurod2_atac_summits.bed
")


tracks <- c("~/methylation/liu/MajorType.IT-L23.CGN-Merge.allc.tsv.bw","~/methylation/liu/MajorType.IT-L4.CGN-Merge.allc.tsv.bw","~/methylation/liu/MajorType.IT-L5.CGN-Merge.allc.tsv.bw","~/methylation/liu/MajorType.NP-L6.CGN-Merge.allc.tsv.bw",
            "~/methylation/liu/MajorType.IT-L6.CGN-Merge.allc.tsv.bw","~/methylation/liu/MajorType.CT-L6.CGN-Merge.allc.tsv.bw","~/methylation/liu/MajorType.L6b.CGN-Merge.allc.tsv.bw","~/methylation/liu/MajorType.PT-L5.CGN-Merge.allc.tsv.bw")

plotmeth_liu<-function(meth_matrix,main) {
  print(plotAverage(plotset=meth_matrix, labels = c('IT-L23','IT-L4','IT-L5','NP-L6','IT-L6','CT-L6','L6b','PT-L5'), xlim = c(-2000,2000),
                    ylim = c(0,100), main = main, xlab = "Dist to peak summit (bp)", ylab = '% CpG Methylation',
                    plotScale = "linear", type = "full", error.estimates = T, yaxs='i',xaxs='i',
                    legend = TRUE, legend_ext = FALSE, legend_pos = 'bottomright',
                    legend_ext_pos = "topleft", cex.axis = 14, cex.lab = 16,xaxt='n',cex.legend = 12, pointsize = 12) )
}

pdf("postnatal_meth_atac_seqs_ummits_no_instersect_proneural.pdf", height = 7, width = 10)
proneuralnt <- getPlotSetArray(tracks=tracks,features="~/openchr/noack/proneuralnt_summits.bed",refgenome='mm10',type = 'mf',add_heatmap=F,xmin=2000,xmax=2000,bin = 100)
plotmeth_liu(proneuralnt,main="ATAC-seq peaks without proneural")
dev.off()

pdf("postnatal_meth_atac_seq_summits_instersect_neurod2.pdf", height = 7, width = 10)
neurod2 <- getPlotSetArray(tracks=tracks,features="~/openchr/noack/neurod2_atac_summits.bed",refgenome='mm10',type = 'mf',add_heatmap=F,xmin=2000,xmax=2000,bin = 100)
plotmeth_liu(neurod2,main="ATAC-seq peaks that intersect with Neurod2")
dev.off()

pdf("postnatal_meth_e14.5_neurod2_summits.pdf", height = 7, width = 10)
neurod2_summits <- getPlotSetArray(tracks=tracks,features="~/neurod2/chipseq/neurod2_peaks_summits.bed",refgenome='mm10',type = 'mf',add_heatmap=F,xmin=2000,xmax=2000,bin = 100)
plotmeth_liu(neurod2_summits,main="E14.5 Neurod2 summits")
dev.off()

pdf("postnatal_meth_p0_neurod2_summits.pdf", height = 7, width = 10)
neurod2_summits <- getPlotSetArray(tracks=tracks,features="~/neurod2/chipseq/p0.neurod2_summits.bed",refgenome='mm10',type = 'mf',add_heatmap=F,xmin=2000,xmax=2000,bin = 100)
plotmeth_liu(neurod2_summits,main="p0 Neurod2 summits")
dev.off()

#supdividir en promoters y enhancers
####


system("
  cd ~/openchr/noack
  bedtools intersect -wa -a proneural_table.bed -b ~/neurod2/chipseq/p0.neurod2_highconfsubpeaks_summit.bed | uniq | awk '{print $1,$2+250,$3-250}' | tr ' ' '\t' | sort -k1,1 -k2,2n > p0neurod2_proneural_table.bed
  
")





system("
  cd ~/openchr/noack
  bedtools intersect -a proneuralnt_summits.bed -b <(awk '{print $1,$2-1000,$3+1000}' promoters.bed | tr ' ' '\t') > promoters_proneuralnt_summits.bed
  bedtools intersect -a proneuralnt_summits.bed -b <(awk '{print $1,$2-1000,$3+1000}' enhancers.bed | tr ' ' '\t') > enhancers_proneuralnt_summits.bed
  bedtools intersect -a neurod2_atac_summits.bed -b <(awk '{print $1,$2-1000,$3+1000}' promoters.bed | tr ' ' '\t') > promoters_neurod2_atac_summits.bed
  bedtools intersect -a neurod2_atac_summits.bed -b <(awk '{print $1,$2-1000,$3+1000}' enhancers.bed | tr ' ' '\t') > enhancers_neurod2_atac_summits.bed
  
  
  bedtools intersect -a p0neurod2_proneural_table.bed -b <(awk '{print $1,$2-1000,$3+1000}' ~/openchr/noack/promoters.bed | tr ' ' '\t') > promoters_p0neurod2_proneural_table.bed
  bedtools intersect -a p0neurod2_proneural_table.bed -b <(awk '{print $1,$2-1000,$3+1000}' ~/openchr/noack/enhancers.bed | tr ' ' '\t') > enhancers_p0neurod2_proneural_table.bed
")

pdf("postnatal_meth.pdf", height = 7, width = 10)
proneuralnt <- getPlotSetArray(tracks=tracks,features="~/openchr/noack/promoters_proneuralnt_summits.bed",refgenome='mm10',type = 'mf',add_heatmap=F,xmin=2000,xmax=2000,bin = 50)
plotmeth_liu(proneuralnt,main="promoter ATAC-seq peaks without proneural")

neurod2 <- getPlotSetArray(tracks=tracks,features="~/openchr/noack/promoters_neurod2_atac_summits.bed",refgenome='mm10',type = 'mf',add_heatmap=F,xmin=2000,xmax=2000,bin = 50)
plotmeth_liu(neurod2,main="promoter ATAC-seq peaks that intersect with Neurod2")

p0neurod2 <- getPlotSetArray(tracks=tracks,features="~/openchr/noack/promoters_p0neurod2_proneural_table.bed",refgenome='mm10',type = 'mf',add_heatmap=F,xmin=2000,xmax=2000,bin = 50)
plotmeth_liu(neurod2,main="promoter ATAC-seq peaks that intersect with p0 Neurod2")


proneuralnt <- getPlotSetArray(tracks=tracks,features="~/openchr/noack/enhancers_proneuralnt_summits.bed",refgenome='mm10',type = 'mf',add_heatmap=F,xmin=2000,xmax=2000,bin = 50)
plotmeth_liu(proneuralnt,main="enhancer ATAC-seq peaks without proneural")

neurod2 <- getPlotSetArray(tracks=tracks,features="~/openchr/noack/enhancers_neurod2_atac_summits.bed",refgenome='mm10',type = 'mf',add_heatmap=F,xmin=2000,xmax=2000,bin = 50)
plotmeth_liu(neurod2,main="enhancer ATAC-seq peaks that intersect with Neurod2")

neurod2 <- getPlotSetArray(tracks=tracks,features="~/openchr/noack/enhancers_p0neurod2_proneural_table.bed",refgenome='mm10',type = 'mf',add_heatmap=F,xmin=2000,xmax=2000,bin = 50)
plotmeth_liu(neurod2,main="enhancer ATAC-seq peaks that intersect with p0 Neurod2")

dev.off()




