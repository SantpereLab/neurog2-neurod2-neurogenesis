## AXON GUIDANCE HEATMAP ################################
#########################################################

library(dplyr)
neurod2_peaks<-readRDS("neurod2_peaks.Rds") %>% select(chr,start,end, summit, ebox_start, peak_id, central_dinucleotide)
write.table(neurod2_peaks, "neurod2_peaks_motifs.bed", quote = F, sep = "\t", col.names = F, row.names = F)

# Intersect the motifs with the trajectories.
system("
   bedtools intersect -wao -a neurod2_peaks_motifs.bed -b ~/openchr/noack/trajs_imputation.bed |\
   bedtools intersect -wao -a - -b ~/openchr/noack/positive_pairs.bed | cut -f1,2,3,4,5,6,7,11,16 > neurod2_motifs_trajs_links.bed
")

master<-read.table("neurod2_motifs_trajs_links.bed") %>% filter(V8!=".",V9!=".")
colnames(master)<-c("chr","start","end","summit","ebox_start","peak_id","central_dinucleotide","traj","gene")

# Add axon guidance families.

library("readxl")
axon_families<-read_excel("AxonGuidance.xlsx") %>% filter(!is.na(...2))
colnames(axon_families)<-c("axon_guidance_family","gene")
# Merge axon guidance table with Neurod2 motifs trajectories table.
master<-merge(master, axon_families, by="gene", all.x=T)
neurod2_families<-master %>% filter(!is.na(axon_guidance_family)) %>% select(axon_guidance_family, gene, peak_id) %>% distinct %>%
  group_by(axon_guidance_family,gene) %>% summarise(npeaks_positive_correlation=n())
extraa<-axon_families
extraa$npeaks_positive_correlation<-0
familiess<-rbind(neurod2_families,extraa) %>% group_by(axon_guidance_family,gene) %>% summarise(npeaks_positive_correlation=sum(npeaks_positive_correlation) )
good_order<-familiess$gene

# Assign to each peak its closest gene.

gencode_tss<-read.table("~/genomes/gencode.vM24.annotation.filtered.genes.tss.bed")
colnames(gencode_tss)<-c("chr","start","end","gene")
library(GenomicRanges)
nearest_gene_neurod2<-distanceToNearest(makeGRangesFromDataFrame(neurod2_peaks),makeGRangesFromDataFrame(gencode_tss), select="all")
nearest_gene_neurod2<-as.data.frame(nearest_gene_neurod2)
nearest_gene_neurod2$genes<-gencode_tss$gene[nearest_gene_neurod2$subjectHits]
npeaks_neurod2_closest<-table(nearest_gene_neurod2$genes) 
npeaks_neurod2_df<-data.frame( "gene"=names(npeaks_neurod2_closest) , "neurod2peaks_closest"= as.vector(npeaks_neurod2_closest) )

familiess<-merge(familiess,npeaks_neurod2_df, by = "gene", all.x=T)
familiess$neurod2peaks_closest[is.na(familiess$neurod2peaks_closest)]<-0

# Add expression + pseudotime data.

library(Seurat)
expression<-readRDS("merged_scRNA_unfiltered_IDs_faye.RDS")
expression_matrix<-expression@assays$RNA@data
metadata<-expression@meta.data
metadata$celltype<-Idents(expression)
# Merge mitotic and non-mitotic cells.
metadata$celltype<-gsub("_M","",metadata$celltype)
metadata<-metadata %>% filter(!is.na(pseudotime))

j<-1
for (cell in c("NSC","IPC","PN1","PN2","PN3")) {
  print(cell)
  metadata_cell<-metadata %>% filter(celltype==cell)
  metadata_cell<-metadata_cell %>% arrange(pseudotime)
  i<-1
  # Do pseudobulks of 175 cells.
  while ((i+175)<=length(metadata_cell$pseudotime)) {
    print(j)
    print(i)
    pseudobulk_cells<-rownames(metadata_cell)[i:(i+174)]
    pseudobulk_expression<-rowSums(as.matrix(expression_matrix[,pseudobulk_cells]))
    if (j==1) { pseudo_expression<-pseudobulk_expression}
    if (j>1) { pseudo_expression<-cbind(pseudo_expression,pseudobulk_expression)}
    i<-i+175
    j<-j+1
  }
}
colnames(pseudo_expression)<-c(paste0("NSC_",1:11),paste0("IPC_",1:8),paste0("PN1_",1:9),paste0("PN2_",1:8),paste0("PN3_",1:2))

library("RColorBrewer")
library(pheatmap)

rownames(familiess)<-familiess$gene
familiess<-familiess[good_order,]
familiess<-familiess %>% filter (gene %in% rownames(pseudo_expression ) )
axon_genes<-familiess$gene
names(axon_genes)<-familiess$axon_guidance_family
axon_genes<-axon_genes[order(factor(names(axon_genes) )) ]
splitt<-names(axon_genes)

# Heatmap annotation

library(circlize)
col_fun = colorRamp2(c(0,150, 300), c("cornsilk", "salmon1", "darkorchid4"))
library(ComplexHeatmap)
ha<-rowAnnotation (empty = anno_empty(border = FALSE), 
                   foo = anno_block(gp = gpar(fill = 7),labels_rot=0, width=unit(20, "mm"),
                                    labels = unique(splitt), labels_gp=gpar(fontsize=10 ) ) )



pdf("~/axon_guidance_heatmap.pdf", useDingbats = F, height = 17, width = 12)

# Print annotation of the heatmap with the names of the families of the axon guidance genes.
h<-Heatmap ( pseudo_expression[axon_genes,], cluster_columns = F, cluster_rows = F, col = col_fun, 
             row_split = splitt, row_gap = unit(0, "mm"),
             right_annotation=ha )
print(h)

# Annotate the number of peaks positively linked to the gene.
h2<-Heatmap ( pseudo_expression[axon_genes,], cluster_columns = F, cluster_rows = F, col = col_fun, 
              row_split = splitt, row_gap = unit(0, "mm"),
              right_annotation= ha2)
print(h2)

# Annotate the number of peaks as a color.
h3<-Heatmap ( pseudo_expression[axon_genes,], cluster_columns = F, cluster_rows = F, col = col_fun, 
              row_split = splitt, row_gap = unit(0, "mm"),
              right_annotation= ha3)
print(h3)

# Annotate the number of peaks that are closest to the gene.
h4<-Heatmap ( pseudo_expression[axon_genes,], cluster_columns = F, cluster_rows = F, col = col_fun, 
              row_split = splitt, row_gap = unit(0, "mm"),
              right_annotation= ha4)
print(h4)

# Put the color
h5<-Heatmap ( pseudo_expression[axon_genes,], cluster_columns = F, cluster_rows = F, col = col_fun, 
              row_split = splitt, row_gap = unit(0, "mm"),
              right_annotation= ha5)
print(h5)

# Print the heatmap with the expression of the genes. 
p<-pheatmap(pseudo_expression[axon_genes,] ,cluster_cols = F, cluster_rows = F, 
            color=colorRampPalette(c( "cornsilk","salmon1","darkorchid4"))(300))
dev.off()







