##############################
###  SUBSET BAMS #############
##############################

#subsetear cada celltype a 350 celulas (que son las celulas de la PN3) para hacer los heatmaps de la figura final
metadata<-noack_obj@meta.data
barcodes_nsc<-rownames(metadata)[metadata$predicted.id=="NSC"]
barcodes_ipc<-rownames(metadata)[metadata$predicted.id=="IPC"]
barcodes_pn1<-rownames(metadata)[metadata$predicted.id=="PN1"]
barcodes_pn2<-rownames(metadata)[metadata$predicted.id=="PN2"]
barcodes_pn3<-rownames(metadata)[metadata$predicted.id=="PN3"]
barcodes_nsc<-barcodes_nsc[sample(1:length(barcodes_nsc),350)]
barcodes_ipc<-barcodes_ipc[sample(1:length(barcodes_ipc),350)]
barcodes_pn1<-barcodes_pn1[sample(1:length(barcodes_pn1),350)]
barcodes_pn2<-barcodes_pn2[sample(1:length(barcodes_pn2),350)]
barcodes_pn3<-barcodes_pn3[sample(1:length(barcodes_pn3),350)]

metadata<-metadata[c(barcodes_nsc,barcodes_ipc,barcodes_pn1,barcodes_pn2,barcodes_pn3),]
metadata$barcoud<-rownames(metadata)
metadata$barcoud<-paste0(substr(metadata$barcoud,1,nchar(metadata$barcoud)-2),"-1")
barcodes_rep1<-metadata %>% filter(orig.ident=="E14_rep1") %>% select(barcoud,predicted.id)
write.table(barcodes_rep1,"rep1_barcodes_imputation.csv",quote = F, sep = ",", row.names = F, col.names = F)
barcodes_rep2<-metadata %>% filter(orig.ident=="E14_rep2") %>% select(barcoud,predicted.id)
write.table(barcodes_rep2,"rep2_barcodes_imputation.csv",quote = F, sep = ",", row.names = F, col.names = F)

system('
  module load Python
  python rep1_barcodes_imputation.csv SRR12385873_alignment/outs/possorted_bam.bam rep1
  python rep2_barcodes_imputation.csv SRR12385874_SRR12385875_merged_possorted_bam.bam rep2
  
  module load SAMtools
  #mergear las bam files
  samtools merge clusterNSC.bam rep1clusterNSC.bam rep2clusterNSC.bam
  samtools merge clusterIPC.bam rep1clusterIPC.bam rep2clusterIPC.bam
  samtools merge clusterPN1.bam rep1clusterPN1.bam rep2clusterPN1.bam
  samtools merge clusterPN2.bam rep1clusterPN2.bam rep2clusterPN2.bam
  samtools merge clusterPN2.bam rep1clusterPN2.bam rep2clusterPN2.bam
  samtools merge clusterPN3.bam rep1clusterPN3.bam rep2clusterPN3.bam
')


# Subset the bams for the non-excitatory cell types.

barcodes_rep1<-metadata %>% filter(orig.ident=="E14_rep1") %>% filter(celltype_atac %in% c("CR","IN","MG+Mural")) %>% select(barcoud,celltype_atac)
write.table(barcodes_rep1,"rep1_barcodes_other_cells.csv",quote = F, sep = ",", row.names = F, col.names = F)
barcodes_rep2<-metadata %>% filter(orig.ident=="E14_rep2") %>% filter(celltype_atac %in% c("CR","IN","MG+Mural")) %>% select(barcoud,celltype_atac)
write.table(barcodes_rep2,"rep2_barcodes_other_cells.csv",quote = F, sep = ",", row.names = F, col.names = F)

system ('
  python subsetbam.py rep1_barcodes_other_cells.csv SRR12385873_alignment/outs/possorted_bam.bam rep1
  python subsetbam.py rep2_barcodes_other_cells.csv SRR12385874_SRR12385875_merged_possorted_bam.bam rep2
')

system('
  samtools merge clusterIN.bam rep1clusterIN.bam rep2clusterIN.bam
  samtools merge clusterCR.bam rep1clusterCR.bam rep2clusterCR.bam
  samtools merge clusterMG+Mural.bam rep1clusterMG+Mural.bam rep2clusterMG+Mural.bam
')




#### PEAK CALLING ###############
#################################

# Merge the bam files derived from the cell-type clusters to do the peak calling and obtain the summits.
system ('
  module load SAMtools
  module load Miniconda3
  samtools merge merged_clusters.bam clusterNSC.bam clusterIPC.bam clusterPN1.bam clusterPN2.bam clusterPN3.bam
  samtools index merged_clusters.bam
  bamCoverage -b merged_clusters.bam -o merged_clusters.bdg --outFileFormat bedgraph --normalizeUsing RPKM
  
  macs2 callpeak -f BAM -t merged_clusters.bam --bdg -g mm --outdir . -n macs2_merged
  
  # Split the peaks obtained with macs2 in subpeaks
  java -jar ~/opt/PeakSplitter_Java/PeakSplitter.jar -p macs2_merged_peaks.narrowPeak -w macs2_merged*.bdg -o . -x 1 -f peaksplit_merged_macs2
')

#### COMPUTE MATRIX ###############
###################################

# Compute the count matrix with the called peaks and the fragments

library(dplyr)
library(Seurat)
library(Signac)
noack_obj<-readRDS("merged_scATAC_integrated_cicero.RDS")

fragments_rep1<-CreateFragmentObject( path = "GSM4710395_scATAC_E14_rep1_fragments.tsv.gz")
fragments_rep2<-CreateFragmentObject( path = "GSM4710396_scATAC_E14_rep2_fragments.tsv.gz")
library(genomation)
peaks<-read.table("peaksplit_merged_macs2.macs2_merged_peaks.narrowPeak.clean.subpeaks.bed")

peaks<-peaks[,c(1,2,3,5)]
colnames(peaks)<-c("chr","start","end","peak_id")

library(GenomicRanges)
peaks<-makeGRangesFromDataFrame(peaks)

rep1_barcodes<-read.table("GSM4710395_scATAC_E14_rep1_barcodes.tsv")$V1
rep1_matrix<-FeatureMatrix(
  fragments = list(fragments_rep1),
  features = peaks,
  cells = rep1_barcodes)
saveRDS(rep1_matrix,"rep1_matrix.RDS")

rep2_barcodes<-read.table("GSM4710396_scATAC_E14_rep2_barcodes.tsv")$V2
rep2_matrix<-FeatureMatrix(
  fragments = list(fragments_rep2),
  features = peaks,
  cells = rep2_barcodes)
saveRDS(rep2_matrix,"rep2_matrix.RDS")

rep1_matrix<-readRDS("rep1_matrix_peaksplit.RDS")
colnames(rep1_matrix)<-stringr::str_replace(colnames(rep1_matrix),"-1","-5")
rep2_matrix<-readRDS("rep2_matrix_peaksplit.RDS")
colnames(rep2_matrix)<-stringr::str_replace(colnames(rep2_matrix),"-1","-6")

matrix<-cbind(rep1_matrix,rep2_matrix)

saveRDS(matrix,"~/p/openchr/noack/mypeaks_featurematrix.RDS")



# TRAJECTORIES RNA-SEQ IMPUTED CLUSTERS ###########
##################################################


# Functions:

# Generate pseudobulks of a given cell types, grouping the cells based in the pseudotime ordering. And aggregate the counts within each pseudobulk.
generate_pseudobulks<-function(barcodes, celltype) {
  require(dplyr)
  
  #Order the cells according to the pseudotime
  metadata_celltype<-metadata[barcodes,]
  sort_barcodes_celltype<-rownames(metadata_celltype %>% arrange(pseudotime))
  counts_celltype<-counts[,sort_barcodes_celltype]
  
  i<-1
  j<-0
  rep1_barcodes<-c()
  rep1_pseudobulks<-c()
  rep2_barcodes<-c()
  rep2_pseudobulks<-c()
  while (i<length(sort_barcodes_celltype)) {
    pseudo_barcodes<-sort_barcodes_celltype[i:(i+174)]
    pseudo_metadata<-metadata[pseudo_barcodes,]
    pseudo_metadata_rep1<- pseudo_metadata %>% filter(orig.ident=="E14_rep1")
    pseudo_metadata_rep2<- pseudo_metadata %>% filter(orig.ident=="E14_rep2")
    pseudo_metadata_rep1_barcodes<- rownames(pseudo_metadata_rep1) %>% stringr::str_replace("-5", "-1")
    pseudo_metadata_rep2_barcodes<- rownames(pseudo_metadata_rep2) %>% stringr::str_replace("-6", "-1")
    
    rep1_barcodes<-c(rep1_barcodes,pseudo_metadata_rep1_barcodes)
    rep1_pseudobulks<-c(rep1_pseudobulks, rep(paste0("pseudotime_pseudobulk_",j,"_",celltype),length(pseudo_metadata_rep1_barcodes)) )
    rep2_barcodes<-c(rep2_barcodes,pseudo_metadata_rep2_barcodes)
    rep2_pseudobulks<-c(rep2_pseudobulks, rep(paste0("pseudotime_pseudobulk_",j,"_",celltype),length(pseudo_metadata_rep2_barcodes) ) )
    
    if (i==1) {
      pseudobulks_cell<-rowSums(as.matrix(counts_celltype[,i:(i+174)]))
    }
    else {
      pseudo<-rowSums(as.matrix(counts_celltype[,i:(i+174)]))
      pseudobulks_cell<-cbind(pseudobulks_cell,pseudo)
    }
    i<-i+175
    print(i)
    j<-j+1
  }
  colnames(pseudobulks_cell)<-c(paste0(celltype,"_",c(1:j)))
  
  #Create and export a dataframe with the barcodes and their corresponding pseudobulk id, to subset the bam files later.
  rep1<-cbind(rep1_barcodes,rep1_pseudobulks)
  rep2<-cbind(rep2_barcodes,rep2_pseudobulks)
  write.table(rep1,paste0("imputation",celltype,"_pseudotime_pseudobulks_rep1_barcodes.csv"),quote = F, sep = ",", row.names = F, col.names = F)
  write.table(rep2,paste0("imputation",celltype,"_pseudotime_pseudobulks_rep2_barcodes.csv"),quote = F, sep = ",", row.names = F, col.names = F)
  return(pseudobulks_cell)
}

# Generate pseudobulks within each cell type with aggregated ATAC counts.
pseudobulks_expression_idents<-function(noack_obj) {
  require(dplyr)
  require(reshape2)
  metadata<-noack_obj@meta.data
  #Remove cells that do not have pseudotime info
  metadata<-metadata[-which(is.na(metadata$pseudotime)), ]
  counts<-noack_obj@assays$MACS2peaks@counts
  #Generate pseudobulks of 175 cells in each cluster
  #There will be: 8 in NSC, 5 in IPC, 7 in PN1, 8 in PN2 and 2 in PN3
  #So take random cells in each cell type: 1400 in NSC, 875 in IPC, 1225 in PN1, 1400 in PN2, and 350 in PN3
  barcodes_nsc<-rownames(metadata)[metadata$predicted.id=="NSC"]
  barcodes_ipc<-rownames(metadata)[metadata$predicted.id=="IPC"]
  barcodes_pn1<-rownames(metadata)[metadata$predicted.id=="PN1"]
  barcodes_pn2<-rownames(metadata)[metadata$predicted.id=="PN2"]
  barcodes_pn3<-rownames(metadata)[metadata$predicted.id=="PN3"]
  barcodes_nsc<-barcodes_nsc[sample(1:length(barcodes_nsc),1400)]
  barcodes_ipc<-barcodes_ipc[sample(1:length(barcodes_ipc),875)]
  barcodes_pn1<-barcodes_pn1[sample(1:length(barcodes_pn1),1225)]
  barcodes_pn2<-barcodes_pn2[sample(1:length(barcodes_pn2),1400)]
  barcodes_pn3<-barcodes_pn3[sample(1:length(barcodes_pn3),350)]
  
  # Generate pseudobulks within each cell type
  pseudo_nsc<-generate_pseudobulks(barcodes_nsc,"NSC")
  pseudo_ipc<-generate_pseudobulks(barcodes_ipc,"IPC")
  pseudo_pn1<-generate_pseudobulks(barcodes_pn1,"PN1")
  pseudo_pn2<-generate_pseudobulks(barcodes_pn2,"PN2")
  pseudo_pn3<-generate_pseudobulks(barcodes_pn3,"PN3")
  
  pseudobulks<-cbind(pseudo_nsc,pseudo_ipc,pseudo_pn1,pseudo_pn2,pseudo_pn3)
  return(pseudobulks)
}

# Assign a trajectory to each peak.
obtain_trajectories<-function (cluster_signal, combinations) {
  require(HybridMTest)
  require(dplyr)
  cell_types<-c("NSC","IPC","PN1","PN2","PN3")
  cluster_total_reads<-colSums(cluster_signal)
  
  # CPM normalization.
  cluster_signal<-t(t(cluster_signal)*1000000/cluster_total_reads)
  cluster_signal<-log2(cluster_signal+1)
  
  # Convert the combinations in vectors to correlate them with the signal of the regions in the pseudobulks.
  combinations_repeated<-strsplit(combinations,"")
  combinations_repeated<-as.data.frame(sapply(combinations_repeated, as.numeric) )
  combinations_repeated<-apply(combinations_repeated, 2, function(x) rep(x, times=c(8,5,7,8,2) ) )
  combinations_repeated<-t(combinations_repeated)
  
  # Function that return which is the combination that correlates the most with the region.
  max_correlation<-function(x) {
    
    if (sum(x)==0) {
      return ("00000")
    }
    
    # Do an ANOVA to test whether the ATAC-seq signal is significantly different between cell-type clusters.
    signals<-data.frame("signal"=x,"cluster"= rep(cell_types,times=c(8,5,7,8,2)) )
    res.aov <- aov(signal ~ cluster, data = signals)
    
    # Extract F and p values of the ANOVA
    fval<-summary(res.aov)[[1]][1,4]
    pval<-summary(res.aov)[[1]][1,5]
    
    # If the signal is significantly different between cell types, compute which is the trajectory of the region using the correlation method.
    if (fval>1 && pval<0.000001 ) {
      # Correlate all the combinations with the row.
      correlations<-row.pearson(combinations_repeated,x)
      
      # Get the combination that gives the maximum correlation.
      maximum_cor_index<-which.max(correlations$stat)
      combination<-paste(combinations [maximum_cor_index] )
      
      # Return the pearson r and the p-value
      return(combination)
    }
    
    else {
      return("11111")
    }
    
  }
  
  # Apply to all rows
  
  patterns<-apply(cluster_signal, 1, max_correlation)
  
  peaks_cluster_df<-data.frame("peak"=rownames(cluster_signal),"cluster"=patterns)
  return(peaks_cluster_df)
}



library(dplyr)
library(Seurat)
library(Signac)

noack_obj<-readRDS("merged_scATAC_integrated_cicero.RDS")
pseudotime<-readRDS("imputation_pseudotime.RDS")
noack_obj$pseudotime<-as.vector(pseudotime@data)
noack_obj@assays$MACS2peaks@counts<-matrix


# Annotate mitotic cell types as non-mitotic
noack_obj$predicted.id[noack_obj$predicted.id=="NSC_M"]<-"NSC"
noack_obj$predicted.id[noack_obj$predicted.id=="IPC_M"]<-"IPC"

# Obtain the pseudobulks from the scATAC-seq object.
noack_pseudobulks_imputation<-pseudobulks_expression_idents(noack_obj)

# Divide the regions in promoters and enhancers.
library(GenomicFeatures)
mm10.gencode.txdb<-makeTxDbFromGFF("gencode.vM24.annotation.filtered.gtf", format="gtf")
library(Signac)
library(ChIPseeker)
peaks_granges<-StringToGRanges(rownames(noack_pseudobulks))

# Ignore downstream and 1st exon because they are very small.
options(ChIPseeker.ignore_1st_exon = TRUE)
options(ChIPseeker.ignore_downstream = TRUE)
anno <- annotatePeak(peaks_granges, TxDb=mm10.gencode.txdb)
annotation<-anno@anno$annotation
annotation[grep("Exon",anno@anno$annotation)]<-"enhancer"
annotation[grep("Intron",anno@anno$annotation)]<-"enhancer"
annotation[grep("UTR",anno@anno$annotation)]<-"enhancer"
annotation[grep("Intergenic",anno@anno$annotation)]<-"enhancer"
annotation[grep("Downstream",anno@anno$annotation)]<-"enhancer"
annotation[grep("Promoter",anno@anno$annotation)]<-"promoter"
annotation[which(abs(anno@anno$distanceToTSS)<=2000)]<-"promoter"

promoter_pseudobulks<-noack_pseudobulks[which(annotation=="promoter"), ]
enhancer_pseudobulks<-noack_pseudobulks[which(annotation=="enhancer"), ]


# Create all the possible combinations (trajectories).
combinations<-c("00001","00011","00111","01111","01000","00100","00010","01110","00110","01100","11110","11100","11000","10000")

# Obtain the trajectoreis in the promoter and enhancer pseudobulks.
trajectories_promoters<-obtain_trajectories(promoter_pseudobulks,combinations)
saveRDS(trajectories_promoters,"noack_trajectories_promoters_imputation_peaksplit.Rds")

trajectories_enhancers<-obtain_trajectories(enhancer_pseudobulks,combinations)
saveRDS(trajectories_enhancers,"noack_trajectories_enhancers_imputation_peaksplit.Rds")

trajectories_promoters<-readRDS("noack_trajectories_promoters_imputation_peaksplit.Rds")
trajectories_enhancers<-readRDS("noack_trajectories_enhancers_imputation_peaksplit.Rds")
table(trajectories_enhancers$cluster)/length(trajectories_enhancers$cluster)


# Enhancers
library(tidyr)
trajectories_enhancers$cluster<-paste0("enhancer_",trajectories_enhancers$cluster)
trajectories_2_bed<-separate(trajectories_enhancers,col="peak",sep="-", into = c("chr","start","end"))
trajectories<-unique(trajectories_2_bed$cluster)

options(scipen=999999999999999)
for (traj in trajectories) {
  trajectory_bed_midpoints<-trajectories_2_bed %>% dplyr::filter(cluster==traj)
  trajectory_bed_midpoints$end<-as.numeric(trajectory_bed_midpoints$start)+as.integer((as.numeric(trajectory_bed_midpoints$end)-as.numeric(trajectory_bed_midpoints$start))/2)
  trajectory_bed_midpoints$start<-trajectory_bed_midpoints$end-1
  write.table(trajectory_bed_midpoints, paste0(traj,"_imputation_traj_peak_peaksplit.bed"), sep = "\t", row.names = F, col.names = F, quote=F)
}

# Promoters
trajectories_promoters$cluster<-paste0("promoter_",trajectories_promoters$cluster) ;trajectories<-trajectories_promoters
trajectories_2_bed<-separate(trajectories_promoters,col="peak",sep="-", into = c("chr","start","end"))
trajectories<-unique(trajectories_2_bed$cluster)

for (traj in trajectories) {
  trajectory_bed_midpoints<-trajectories_2_bed %>% dplyr::filter(cluster==traj)
  trajectory_bed_midpoints$end<-as.numeric(trajectory_bed_midpoints$start)+as.integer((as.numeric(trajectory_bed_midpoints$end)-as.numeric(trajectory_bed_midpoints$start))/2)
  trajectory_bed_midpoints$start<-trajectory_bed_midpoints$end-1
  write.table(trajectory_bed_midpoints, paste0(traj,"_imputation_traj_peak_peaksplit.bed"), sep = "\t", row.names = F, col.names = F, quote=F)
}




