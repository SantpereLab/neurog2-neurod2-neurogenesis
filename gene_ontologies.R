# GENE ONTOLOGY ###########
################################################


###
## gene ontologies cogiendo solo los picos que correlacionan positivamente con la expresion
####


#hacer primero con todos los peaks de atac en general para comprobar que funciona

positive_links<-read.table("positive_pairs.tsv", header = T)
positive_links<-positive_links[,c(2,3,4,1,5,6,7)]
positive_links$peak<-paste(positive_links$chrom_enh,positive_links$start_enh,positive_links$end_enh,sep = "-")

trajs_imputation<-rbind(readRDS("noack_trajectories_promoters_imputation_peaksplit.Rds"),(readRDS("noack_trajectories_enhancers_imputation_peaksplit.Rds")) )
library(stringr)
trajs_imputation<-cbind(str_split_fixed(trajs_imputation$peak, "-", 3),trajectories_2$cluster)
write.table(trajs_imputation,"trajs_imputation.bed", sep = "\t", quote = F, col.names = F, row.names = F)

system("bedtools intersect -wao -b trajs_imputation.bed -a neural_atac_peaks_neurod2_intersect.bed |cut -f1,2,3,8 > trajs_imputation_neurod2_neural_intersect.bed")

atac_peaks_intersect_neurod2<-read.table("/home/xabier/openchr/noack/trajs_imputation_neurod2_neural_intersect.bed")
atac_peaks_intersect_neurod2$peak<-paste(atac_peaks_intersect_neurod2$V1,atac_peaks_intersect_neurod2$V2,atac_peaks_intersect_neurod2$V3,sep = "-")
neurod2_positive_links_trajectories<-merge(atac_peaks_intersect_neurod2,positive_links,by="peak")
neurod2_positive_links_trajectories<-neurod2_positive_links_trajectories %>% filter(V4 !="." )
neurod2_positive_links_trajectories$V4<-str_replace(neurod2_positive_links_trajectories$V4, "11110|11100|11000|10000", "early")
neurod2_positive_links_trajectories$V4<-str_replace(neurod2_positive_links_trajectories$V4, "01110|01100|00110|01000|00100|00010", "transcient")
neurod2_positive_links_trajectories$V4<-str_replace(neurod2_positive_links_trajectories$V4, "01111|00111|00011|00001", "neuronal")
neurod2_positive_links_trajectories$V4<-str_replace(neurod2_positive_links_trajectories$V4, "11111", "invariant")

#no distingo entre enhancer y enhancer de momento. pq ellos usan anotacion refseq y yo gencode del chipseeker.
#deberia cambiarlo quizas. aunq de momento asi va bien pq son analisis independientes. el mio, de la señal y motivos.
# y el suyo, de correlacion con la expresion
jj<-neurod2_positive_links_trajectories %>% select(V4,gene_name) %>% distinct() %>% group_by(gene_name) %>% summarise(n=n()) %>% filter(n==1)
unique_genes<-jj$gene_name

neurod2_trajectories_genes<-split(neurod2_positive_links_trajectories$gene_name,neurod2_positive_links_trajectories$V4)

#mirar como overlapan los genes de las distintas trayectorias
#funcion para hacer el heatmap del overlap de los genes de las trayectorias
heatmap_overlap_genes<-function(trajectories_genes) {
  for(traj in names(trajectories_genes)) { 
    trajectories_genes[[traj]]<-unique(unlist(strsplit(trajectories_genes[[traj]],","))) 
  }
  unique<-c()
  for(traj in names(trajectories_genes)) { 
    uni<-length(which(trajectories_genes[[traj]] %in% unique_genes))
    unique<-c(unique,uni)
  }
  matrix<-matrix(nrow = length(names(trajectories_genes)), ncol =  length(names(trajectories_genes)))
  for (i in 1:length(names(trajectories_genes))) {
    for (j in 1:length(names(trajectories_genes))) {
      matrix[i,j]<-length(intersect(trajectories_genes[[i]], trajectories_genes[[j]]))
    }
  }
  colnames(matrix)<-names(trajectories_genes)
  rownames(matrix)<-names(trajectories_genes)
  matrix<-cbind(matrix,unique)
  total_genes_trajectory<-unlist(lapply(trajectories_genes, function(x) {length(x)} ) )
  genes_proportion_intersect<-matrix/total_genes_trajectory*100
  print(pheatmap::pheatmap(genes_proportion_intersect, cluster_rows =  F, cluster_cols = F,
                           display_numbers = round(matrix,2)))
}
heatmap_overlap_genes(neurod2_trajectories_genes)

#barplot los genes unicos de cada trayectoria
uniq_length<-c()
for (traj in names(neurod2_trajectories_genes)) {
  length<-length(which(unique_genes %in% neurod2_trajectories_genes[[traj]]))
  uniq_length<-c(uniq_length,length)
}
names(uniq_length)<-names(neurod2_trajectories_genes)













#los picos de neurod2 de p0
p0_neurod2_intersect<-read.table("neural_atac_peaks_p0_neurod2_intersect.bed")
p0_neurod2_intersect$peak<-paste0(p0_neurod2_intersect$V1,"-",p0_neurod2_intersect$V2,"-",p0_neurod2_intersect$V3)
#los picos de atac de p56
p56_intersect<-read.table("neural_atac_peaks_p56_li_intersect.bed")
p56_intersect$peak<-paste0(p56_intersect$V1,"-",p56_intersect$V2,"-",p56_intersect$V3)


trajectories_genes<-neurod2_trajectories_genes


load("~/onts.RData")

#sacar todos los genes linkeados para usarlos como background
links_genes<-unique(positive_links$gene_name)

heatmap_ontologies_neurod2_pcts<-function(trajectories_genes) {
  
  i<-1 
  print("Number of genes in trajectories:")
  for(traj in names(trajectories_genes)) {
    #quitar las trayectorias que tengan menos de 10 genes unicos. de las que quito, el maximo es 9
    # de las que quedan, el minimo es 50. asi que es buen cutoff
    print(traj)
    traj_genes<-unique(unlist(trajectories_genes[traj]))
    #traj_genes<-traj_genes[which(traj_genes %in% unique_genes)]
    print(length(traj_genes))
    traj_genes_terms<-gost(traj_genes,
                           organism = "mmusculus",user_threshold=1,sources = "GO:BP", evcodes = T)$result
    traj_genes_terms$traj<-rep(traj,length(traj_genes_terms[,1]))
    print(dim(traj_genes_terms))
    if (i==1) {
      traj_genes_terms_all<-traj_genes_terms
    }
    else { 
      traj_genes_terms_all<-rbind(traj_genes_terms_all,traj_genes_terms)    
    }
    i<-i+1
  }
  
}
#sacar los terms siginificativos en al menos una trayectoria
significant_terms<-unique(traj_genes_terms_all$term_name[which(traj_genes_terms_all$p_value<0.05)])
traj_genes_terms_all<-traj_genes_terms_all %>% filter(term_name %in% significant_terms)
saveRDS(traj_genes_terms_all,"~/neurod2/traj_genes_terms_all.Rds")

genes_ontologies<-traj_genes_terms_all %>% dplyr::select(term_name,intersection)
#hacer el enrichment mejor q el minus log2 pval. pq con el minus log2 pval las trayectorias con mas genes salen con valores mas altos
traj_genes_terms_all$odds_ratio<-(traj_genes_terms_all$intersection_size/traj_genes_terms_all$query_size)/(traj_genes_terms_all$term_size/traj_genes_terms_all$effective_domain_size)

traj_genes_terms_all_select<-traj_genes_terms_all %>% dplyr::filter(term_size>100) %>% dplyr::select(term_name, odds_ratio, traj) %>% distinct()
traj_genes_terms_significant_matrix<-reshape(traj_genes_terms_all_select, idvar = "traj", timevar = "term_name", direction = "wide")
traj_genes_terms_significant_matrix[is.na(traj_genes_terms_significant_matrix)]<-0
rownames(traj_genes_terms_significant_matrix)<-traj_genes_terms_significant_matrix$traj
colnames(traj_genes_terms_significant_matrix)<-colnames(traj_genes_terms_significant_matrix) %>% stringr::str_replace("odds_ratio.", "")
traj_genes_terms_significant_matrix<-as.matrix(traj_genes_terms_significant_matrix[,-1])
traj_genes_terms_significant_matrix_scaled<-scale(traj_genes_terms_significant_matrix)


#primero hacer con todo el pheatmap
library(pheatmap)
p<-pheatmap(traj_genes_terms_significant_matrix_scaled,
            treeheight_col = 0, cluster_rows = F, show_colnames = F, cluster_cols = T,
            treeheight_row = 0, color = hcl.colors(50, "Mint"))

pdf( "pheatmap.pdf", useDingbats = F, width = 20, height = 5)
print(p)
dev.off()
traj_genes_terms_significant_matrix<-t(traj_genes_terms_significant_matrix[c("early","transcient","neuronal"),])
#distintos numeros para cada trayectoria, para sacar 15 significativos en cada
numbers<-c(100,25,15)
i<-1
for (traj in colnames(traj_genes_terms_significant_matrix)) {
  onts_max_traj<-names(tail(sort(traj_genes_terms_significant_matrix_scaled[traj,]),numbers[i]))
  pdf( paste0(traj,"_ontologies.pdf"), useDingbats = F, width = 1, height = 4)
  pheatmap(traj_genes_terms_significant_matrix[names(sort(traj_genes_terms_significant_matrix[onts_max_traj[1:numbers[i]],traj], decreasing = T)) ,],
           breaks = seq(from=0, to=10, by=0.1), show_rownames=F, legend = F, main = traj,
           treeheight_col = 0, show_colnames = T,cluster_rows = F, cluster_cols = F, treeheight_row = 0, color = hcl.colors(100, "Mint"))
  dev.off()
  pdf( paste0(traj,"_ontologies_names.pdf"), useDingbats = F, width = 8, height = 4)
  pheatmap(traj_genes_terms_significant_matrix[names(sort(traj_genes_terms_significant_matrix[onts_max_traj[1:numbers[i]],traj], decreasing = T)) ,],
           breaks = seq(from=0, to=10, by=0.1), show_rownames=T, legend = F, main = traj,
           treeheight_col = 0, show_colnames = T,cluster_rows = F, cluster_cols = F, treeheight_row = 0, color = hcl.colors(100, "Mint"))
  dev.off()
  i<-i+1
}

onts_all_ordered<-colnames(p$gtable$grobs[[1]]$children$GRID.rect.271$gp$fill)
#asi se printeaba el heatmap en el pdf:

pdf("pheatmap1.pdf", useDingbats = F, width = 2, height = 7)
pheatmap(t(traj_genes_terms_significant_matrix[,names(sort(traj_genes_terms_significant_matrix["Descending",onts_max_descending[1:15]], decreasing = T)) ]),
         breaks = seq(from=0, to=10, by=0.1), show_rownames=F, legend = T,
         treeheight_col = 0, show_colnames = T,cluster_rows = F, cluster_cols = F, treeheight_row = 0, color = hcl.colors(100, "Mint"))
dev.off()

save.image("~/onts.RData")

load("~/onts.RData")
#las proporciones con las trayectorias reducidas
top_terms<-gprofiler2::gost(links_genes, organism = "mmusculus", user_threshold = 1, sources = "GO:BP", evcodes = T)$result
props_early<-c()
props_transcient<-c()
props_neuronal<-c()
for (term in top_terms$term_name) {
  term_genes<-unlist(strsplit(top_terms$intersection[which(top_terms$term_name==term)],","))
  term_peaks<-neurod2_positive_links_trajectories$peak[neurod2_positive_links_trajectories$gene_name %in% term_genes]
  trajectories_2_term<-neurod2_positive_links_trajectories %>% filter(peak %in% term_peaks)
  proportions<-table(trajectories_2_term$V4)/length(trajectories_2_term$V4)
  props_early<-c(props_early,proportions["early"])
  props_transcient<-c(props_transcient,proportions["transcient"])
  props_neuronal<-c(props_neuronal,proportions["neuronal"])
}
proportions_terms<-data.frame("props_early"=props_early,"props_transcient"=props_transcient,"props_neuronal"=props_neuronal)

rownames(proportions_terms)<-top_terms$term_name
proportions_terms<-as.matrix(proportions_terms)
colnames(proportions_terms)<-c("early","transcient","neuronal")
proportions_terms[is.na(proportions_terms)]<-0
proportions_terms<-t(proportions_terms)

#el heatmap general:
pdf( "proportions_pheatmap.pdf", useDingbats = F, width = 20, height = 5)
pheatmap(as.matrix(t(proportions_terms[onts_all_ordered,])), treeheight_col = 0, color = hcl.colors(100, "Lajolla"), show_rownames = T, legend=F,
         breaks = seq(from=0, to=1, by=0.01), cluster_cols=F,cluster_rows = F, show_colnames = F, treeheight_row = 0, legend = T)
dev.off()

numbers<-c(100,25,15)
i<-1
for (traj in c("early","transcient","neuronal")) {
  onts_max_traj<-names(tail(sort(traj_genes_terms_significant_matrix_scaled[traj,]),numbers[i]))
  pdf( paste0(traj,"_ontologies_proportions.pdf"), useDingbats = F, width = 8, height = 4)
  pheatmap(as.matrix(t(proportions_terms[,names(sort(traj_genes_terms_significant_matrix[onts_max_traj[1:numbers[i]],traj], decreasing = T))])), treeheight_col = 0, color = hcl.colors(100, "Lajolla"), 
           show_rownames = T,legend=F,breaks = seq(from=0, to=1, by=0.01), cluster_cols=F,cluster_rows = F, show_colnames = T, treeheight_row = 0, main = traj)
  dev.off()
  i<-i+1
}


#sacar para cada trayectoria los terms que tienen mas proporcion
#dentro de terms que sean significativos y tengan un minomo de size 100
positive_links_neurod2<-positive_links %>% filter (peak %in% neurod2_intersect$peak)
all_genes_neurod2_trajs<-unique(positive_links_neurod2$gene_name)
sig_terms<-gprofiler2::gost(all_genes_neurod2_trajs, organism = "mmusculus", user_threshold = 0.05, sources = "GO:BP", evcodes = T)$result
saveRDS(sig_terms,"sig_terms.Rds")
proportions_terms<-proportions_terms[,sig_terms$term_name]
for (traj in c("early","transcient","neuronal")) {
  onts_max_traj<-names(head(sort(proportions_terms[traj,], decreasing = T),45))
  pheatmap(as.matrix(t(proportions_terms[,onts_max_traj])), treeheight_col = 0, color = hcl.colors(100, "Lajolla"), show_rownames = T ,
           breaks = seq(from=0, to=1, by=0.01), cluster_cols=F,cluster_rows = F, show_colnames = T, treeheight_row = 0, legend = T, main = traj)
}



### sacar el porcentaje de picos de neurod2 en los picos asociados a las ontologies
#### primero sacar todos los genes con link que estan dentro de las ontologies
links_genes<-unique(positive_links$gene_name)
proportion_neurod2_peaks<-function(traj_terms, peaks) {
  neurod2_pct_terms<-c()
  significant<-c()
  i<-1
  for (term in traj_terms){
    print(i)
    i<-i+1
    term_genes<-unlist(strsplit(top_terms$intersection[which(top_terms$term_name==term)],","))
    term_peaks<-positive_links$peak[which(positive_links$gene_name %in% term_genes)]
    neurod2_pct<-length(which(term_peaks %in% peaks$peak))/length(term_peaks)
    neurod2_pct_terms<-c(neurod2_pct_terms,neurod2_pct)
    #hacer el fisher test
    neurod2_peaks_linked<-peaks$peak[which(peaks$peak %in% positive_links$peak)]
    dat <- data.frame(
      "term" = c(length(which(term_peaks %in% peaks$peak)), length(which(!term_peaks %in% peaks$peak))),
      "no_term" = c(length(which(!neurod2_peaks_linked %in% term_peaks)), length(which(!positive_links$peak %in% unique(c(peaks$peak,term_peaks))))),
      row.names = c("neurod2", "no_neurod2")
    )
    colnames(dat) <- c("term", "no_term")
    print(term)
    print(dat)
    fisher<-fisher.test(dat)
    if (fisher$estimate>1 && fisher$p.value<0.05) {sig<-1}
    else {sig<-0}
    significant<-c(significant,sig)
  }
  names(neurod2_pct_terms)<-traj_terms
  names(significant)<-traj_terms
  return(cbind(neurod2_pct_terms,significant))
}
pcts_all<-proportion_neurod2_peaks(onts_all_ordered,neurod2_intersect)
pcts_all_p0<-proportion_neurod2_peaks(onts_all_ordered,p0_neurod2_intersect)
pcts_all_p56<-proportion_neurod2_peaks(onts_all_ordered,p56_intersect)
pct_matrix<-cbind(pcts_all[,1],pcts_all_p0[,1],pcts_all_p56[,1])
sig_matrix<-cbind(pcts_all[,2],pcts_all_p0[,2],pcts_all_p56[,2])
colnames(pct_matrix)<-c("E14.5 Neurod2 peaks","P0 Neurod2 peaks","P56 ATAC-seq peaks")
library(circlize)
library(ComplexHeatmap)

pdf( "fisher_heatmap.pdf", useDingbats = F, width = 20, height = 5)
Heatmap(t(pct_matrix), col =colorRamp2(c(0,0.33,0.66,1), c("white","skyblue","darkblue","black")),
        row_names_max_width = unit(18, "cm"), rect_gp = gpar(col = "gray", lwd = 1),
        column_split = NULL, show_row_names = T, show_heatmap_legend = F,
        row_names_gp = grid::gpar(fontsize = 12),  heatmap_legend_param = list(title = "proportion"), show_column_names = F,
        cluster_rows = F, cluster_columns = F ) 
dev.off()


#he probado a coger 60 en progenitors para obtener ontologies con overlap significativo de fisher con neurod2
onts_max_progenitor_ordered<-names(tail(sort(traj_genes_terms_significant_matrix_scaled["early",]),100))
onts_max_progenitor_ordered<-names(sort(traj_genes_terms_significant_matrix[onts_max_progenitor_ordered[1:100],"early"], decreasing = T))
saveRDS(onts_max_progenitor_ordered,"~/neurod2/onts_max_progenitor_ordered.Rds")

neurod2_intersect$peak<-paste0(neurod2_intersect$V1,"-",neurod2_intersect$V2,"-",neurod2_intersect$V3)
pcts_progenitor<-proportion_neurod2_peaks(onts_max_progenitor_ordered,neurod2_intersect)
pcts_progenitor_p0<-proportion_neurod2_peaks(onts_max_progenitor_ordered,p0_neurod2_intersect)
pcts_progenitor_p56<-proportion_neurod2_peaks(onts_max_progenitor_ordered,p56_intersect)
pct_matrix<-cbind(pcts_progenitor[,1],pcts_progenitor_p0[,1],pcts_progenitor_p56[,1])
sig_matrix<-cbind(pcts_progenitor[,2],pcts_progenitor_p0[,2],pcts_progenitor_p56[,2])
colnames(pct_matrix)<-c("E14.5 Neurod2 peaks","P0 Neurod2 peaks","P56 ATAC-seq peaks")
library(circlize)
library(ComplexHeatmap)

pdf( "fisher_heatmap1.pdf", useDingbats = F, width = 8, height = 14)
Heatmap(pct_matrix, col =colorRamp2(c(0,0.33,0.66,1), c("white","skyblue","darkblue","black")),
        row_names_max_width = unit(18, "cm"), rect_gp = gpar(col = "gray", lwd = 1),
        column_names_gp = grid::gpar(fontsize = 11),
        column_split = NULL, show_row_names = T, show_heatmap_legend = F,
        row_names_gp = grid::gpar(fontsize = 12),  heatmap_legend_param = list(title = "proportion"), show_column_names = T,
        cluster_rows = F, cluster_columns = F, cell_fun = function(j, i, x, y, w, h, fill) {
          if(sig_matrix[i, j]==1) {
            grid.text("*", x, y)
          }
        })
dev.off()

onts_max_transcient_ordered<-names(tail(sort(traj_genes_terms_significant_matrix_scaled["transcient",]),25))
onts_max_transcient_ordered<-names(sort(traj_genes_terms_significant_matrix[onts_max_transcient_ordered[1:25],"transcient"], decreasing = T))
saveRDS(onts_max_transcient_ordered,"~/neurod2/onts_max_transcient_ordered.Rds")


pcts_transcient<-proportion_neurod2_peaks(onts_max_transcient_ordered,neurod2_intersect)
pcts_transcient_p0<-proportion_neurod2_peaks(onts_max_transcient_ordered,p0_neurod2_intersect)
pcts_transcient_p56<-proportion_neurod2_peaks(onts_max_transcient_ordered,p56_intersect)
pct_matrix<-cbind(pcts_transcient[,1],pcts_transcient_p0[,1],pcts_transcient_p56[,1])
sig_matrix<-cbind(pcts_transcient[,2],pcts_transcient_p0[,2],pcts_transcient_p56[,2])
colnames(pct_matrix)<-c("E14.5 Neurod2 peaks","P0 Neurod2 peaks","P56 ATAC-seq peaks")
library(circlize)
library(ComplexHeatmap)

pdf( "fisher_heatmap2.pdf", useDingbats = F, width = 8, height = 14)
Heatmap(pct_matrix[1:25,], col =colorRamp2(c(0,0.33,0.66,1), c("white","skyblue","darkblue","black")),
        row_names_max_width = unit(18, "cm"), rect_gp = gpar(col = "gray", lwd = 1),
        column_names_gp = grid::gpar(fontsize = 11),
        column_split = NULL, show_row_names = T, show_heatmap_legend = F,
        row_names_gp = grid::gpar(fontsize = 12),  heatmap_legend_param = list(title = "proportion"), show_column_names = T,
        cluster_rows = F, cluster_columns = F, cell_fun = function(j, i, x, y, w, h, fill) {
          if(sig_matrix[i, j]==1) {
            grid.text("*", x, y)
          } 
        }) 
dev.off()

onts_max_neuronal_ordered<-names(tail(sort(traj_genes_terms_significant_matrix_scaled["neuronal",]),15))
onts_max_neuronal_ordered<-names(sort(traj_genes_terms_significant_matrix[onts_max_neuronal_ordered[1:15],"neuronal"], decreasing = T))
saveRDS(onts_max_neuronal_ordered,"~/neurod2/onts_max_neuronal_ordered.Rds")

pcts_neuronal<-proportion_neurod2_peaks(onts_max_neuronal_ordered,neurod2_intersect)
pcts_neuronal_p0<-proportion_neurod2_peaks(onts_max_neuronal_ordered,p0_neurod2_intersect)
pcts_neuronal_p56<-proportion_neurod2_peaks(onts_max_neuronal_ordered,p56_intersect)
pct_matrix<-cbind(pcts_neuronal[,1],pcts_neuronal_p0[,1],pcts_neuronal_p56[,1])
sig_matrix<-cbind(pcts_neuronal[,2],pcts_neuronal_p0[,2],pcts_neuronal_p56[,2])
colnames(pct_matrix)<-c("E14.5 Neurod2 peaks","P0 Neurod2 peaks","P56 ATAC-seq peaks")
library(circlize)
library(ComplexHeatmap)

pdf( "fisher_heatmap3.pdf", useDingbats = F, width = 8, height = 7)
Heatmap(pct_matrix[1:15,], col =colorRamp2(c(0,0.33,0.66,1), c("white","skyblue","darkblue","black")),
        row_names_max_width = unit(18, "cm"), rect_gp = gpar(col = "gray", lwd = 1),
        column_names_gp = grid::gpar(fontsize = 11),
        column_split = NULL, show_row_names = T, show_heatmap_legend = F,
        row_names_gp = grid::gpar(fontsize = 12),  heatmap_legend_param = list(title = "proportion"), show_column_names = T,
        cluster_rows = F, cluster_columns = F, cell_fun = function(j, i, x, y, w, h, fill) {
          if(sig_matrix[i, j]==1) {
            grid.text("*", x, y)
          } 
        }) 
dev.off()







