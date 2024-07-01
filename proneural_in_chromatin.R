## NEUROG2 / NEUROD2 OVERLAP ########
#############################################
library(eulerr)
library(bedr)
library(dplyr)
neurog2_peaks<-read.table("~/neurog2/chipseq/e14.5.neurog2_peaks.subpeaks.bed")
neurod2_peaks<-read.table("~/neurod2/chipseq/e14.5.neurod2_highconfsubpeaks_summit.bed")
join<-bedr.join.region(neurog2_peaks,neurod2_peaks)
colnames(join)[6]<-c("chr_ng2")
interseccion<-length((join[,c(1:6)] %>% filter(chr_ng2 != "."))[,1])
private_neurog2<-length((join[,c(1:6)] %>% filter(chr_ng2 == "."))[,1])
join2<-bedr.join.region(neurod2_peaks,neurog2_peaks)
colnames(join2)[6]<-c("chr_ng2")
private_neurod2<-length((join2[,c(1:6)] %>% filter(chr_ng2 == "."))[,1])

plot(euler(c("A" = private_neurod2, "B" = private_neurog2,
             "A&B" = interseccion) ), list(type = "percent")  )

plot(euler(c("A" = 6, "B" = 10, "C" = 10,
             "A&B" = 4, "A&C" = 4, "B&C" = 4,
             "A&B&C" = 2)), list(type = "percent") )





# TRAJECTORIES ##########################
################################################


# Pie chart proportions trajectories all ATAC-seq peaks.

system("
  cd ~/openchr/noack/trajectories_beds/definitive
  cat enhancer_invariant_imputation_traj_peak_peaksplit.bed promoter_invariant_imputation_traj_peak_peaksplit.bed | sort -k1,1 -k2,2n > invariant.bed
  cat enhancer_neuronal_imputation_traj_peak_peaksplit.bed promoter_neuronal_imputation_traj_peak_peaksplit.bed | sort -k1,1 -k2,2n > neuronal.bed
  cat enhancer_transcient_imputation_traj_peak_peaksplit.bed promoter_transcient_imputation_traj_peak_peaksplit.bed | sort -k1,1 -k2,2n > transcient.bed
  cat enhancer_early_imputation_traj_peak_peaksplit.bed promoter_early_imputation_traj_peak_peaksplit.bed | sort -k1,1 -k2,2n > early.bed

  bedtools intersect -wa -a ../../peaksplit_merged_macs2.macs2_merged_peaks.narrowPeak.clean.subpeaks.bed -b invariant.bed > invariant_peaks.bed
  bedtools intersect -wa -a ../../peaksplit_merged_macs2.macs2_merged_peaks.narrowPeak.clean.subpeaks.bed -b early.bed > early_peaks.bed
  bedtools intersect -wa -a ../../peaksplit_merged_macs2.macs2_merged_peaks.narrowPeak.clean.subpeaks.bed -b transcient.bed > transcient_peaks.bed
  bedtools intersect -wa -a ../../peaksplit_merged_macs2.macs2_merged_peaks.narrowPeak.clean.subpeaks.bed -b neuronal.bed > neuronal_peaks.bed
  
  bedtools intersect -wa -a ../../proneural_table.bed -b invariant_peaks.bed | sort -k1,1 -k2,2n | uniq > invariant_proneural_table.bed
  bedtools intersect -wa -a ../../proneural_table.bed -b early_peaks.bed | sort -k1,1 -k2,2n | uniq > early_proneural_table.bed
  bedtools intersect -wa -a ../../proneural_table.bed -b transcient_peaks.bed | sort -k1,1 -k2,2n | uniq > transcient_proneural_table.bed
  bedtools intersect -wa -a ../../proneural_table.bed -b neuronal_peaks.bed | sort -k1,1 -k2,2n | uniq > neuronal_proneural_table.bed
")

invariant_proneural_table<-read.table("~/openchr/noack/trajectories_beds/definitive/invariant_proneural_table.bed")[,4]
early_proneural_table<-read.table("~/openchr/noack/trajectories_beds/definitive/early_proneural_table.bed")[,4]
transcient_proneural_table<-read.table("~/openchr/noack/trajectories_beds/definitive/transcient_proneural_table.bed")[,4]
neuronal_proneural_table<-read.table("~/openchr/noack/trajectories_beds/definitive/neuronal_proneural_table.bed")[,4]

pdf("~/venns.pdf")
pie(table(invariant_proneural_table))
pie(table(early_proneural_table))
pie(table(transcient_proneural_table))
pie(table(neuronal_proneural_table))
dev.off()

# Pie chart proportions trajectories all ATAC-seq peaks.

invariant_proneural_table_neurod2<-read.table("~/openchr/noack/trajectories_beds/definitive/invariant_proneural_table.bed")[,4]
early_proneural_table_neurod2<-read.table("~/openchr/noack/trajectories_beds/definitive/early_proneural_table.bed")[,4]
transcient_proneural_table_neurod2<-read.table("~/openchr/noack/trajectories_beds/definitive/transcient_proneural_table.bed")[,4]
neuronal_proneural_table_neurod2<-read.table("~/openchr/noack/trajectories_beds/definitive/neuronal_proneural_table.bed")[,4]

invariant_proneural_table_neurog2<-read.table("~/openchr/noack/trajectories_beds/definitive/invariant_proneural_table.bed")[,4]
early_proneural_table_neurog2<-read.table("~/openchr/noack/trajectories_beds/definitive/early_proneural_table.bed")[,4]
transcient_proneural_table_neurog2<-read.table("~/openchr/noack/trajectories_beds/definitive/transcient_proneural_table.bed")[,4]
neuronal_proneural_table_neurog2<-read.table("~/openchr/noack/trajectories_beds/definitive/neuronal_proneural_table.bed")[,4]

# Fisher tests intersections of the different factors with the trajectories.

all_atac_peaks<-length( c(invariant_proneural_table,early_proneural_table,transcient_proneural_table,neuronal_proneural_table) )
all_neurod2_peaks_in_atac<-length ( which( c(invariant_proneural_table,early_proneural_table,transcient_proneural_table,neuronal_proneural_table) =="4.private_neurod2" ) )
all_neurog2_peaks_in_atac<-length ( which( c(invariant_proneural_table,early_proneural_table,transcient_proneural_table,neuronal_proneural_table) =="2.private_neurog2" ) )
all_shared_peaks_in_atac<-length ( which( c(invariant_proneural_table,early_proneural_table,transcient_proneural_table,neuronal_proneural_table) =="3.shared" ) )

all_transcient<-length(transcient_proneural_table)
all_early<-length(early_proneural_table)
all_neuronal<-length(neuronal_proneural_table)

early_neurog2<-length (which(early_proneural_table_neurog2=="2.private_neurog2") )
transcient_neurod2<-length (which(transcient_proneural_table_neurod2=="4.private_neurod2") )
neuronal_neurod2<-length (which(neuronal_proneural_table_neurod2=="4.private_neurod2") )
transcient_shared<-length (which(transcient_proneural_table_shared=="3.shared") )


jj<-fisher.test( matrix ( c(early_neurog2 ,all_early-early_neurog2, 
                            all_neurog2_peaks_in_atac-early_neurog2, all_atac_peaks-all_early-all_neurog2_peaks_in_atac ), nrow = 2 )   )
jj$p.value
jj$estimate


jj<-fisher.test( matrix ( c(transcient_neurod2 ,all_transcient-transcient_neurod2, 
                            all_neurod2_peaks_in_atac-transcient_neurod2, all_atac_peaks-all_transcient-all_neurod2_peaks_in_atac ), nrow = 2 )   )
jj$p.value
jj$estimate

jj<-fisher.test( matrix ( c(neuronal_neurod2 ,all_neuronal-neuronal_neurod2, 
                            all_neurod2_peaks_in_atac-neuronal_neurod2, all_atac_peaks-all_neuronal-all_neurod2_peaks_in_atac ), nrow = 2 )   )
jj$p.value
jj$estimate

jj<-fisher.test( matrix ( c(transcient_shared ,all_transcient-transcient_shared, 
                            all_shared_peaks_in_atac-transcient_shared, all_atac_peaks-all_transcient-all_shared_peaks_in_atac ), nrow = 2 )   )
jj$p.value
jj$estimate


# Pie chart of all the trajectories, without the summarization.
setwd("~/openchr/noack")
trajectories_2<-rbind(readRDS("noack_trajectories_promoters_imputation_peaksplit.Rds"),(readRDS("noack_trajectories_enhancers_imputation_peaksplit.Rds")) )
proportions_trajectories<-as.vector(table(trajectories_2$cluster))
names(proportions_trajectories)<-names(table(trajectories_2$cluster))
proportions_trajectories<-proportions_trajectories[c("11111","10000","11000","11100","11110",
                                                     "01000","01100","00100","01110","00110","00010",
                                                     "01111","00111","00011","00001")]
colors<-c(
  "papayawhip",
              rep("darkgoldenrod1",4),rep("paleturquoise1",6),
              rep("lightgreen",4)
)
pie(proportions_trajectories, main = "Proportions of trajectories" ,col = colors )


# The same but only with the peaks that intersect with Neurod2.

neurdo2_peaks_trajectories<-read.table("neurod_peak_travectories.tsv")
neurdo2_peaks_trajectories$V2<-str_split_fixed(neurdo2_peaks_trajectories$V2, "_", 2)[,1]
neurdo2_peaks_trajectories<-neurdo2_peaks_trajectories %>% filter(V2!=".") %>% distinct()
prop_neurod2<-table(neurdo2_peaks_trajectories$V2)
prop_neurod2<-prop_neurod2[c("11111","10000","11000","11100","11110",
                             "01000","01100","00100","01110","00110","00010",
                             "01111","00111","00011","00001")]
pie(prop_neurod2, main = "Proportions of trajectories" ,col = colors )

# The same but only with the peaks that intersect with Neurog2.

trajectories_2<-rbind(readRDS("noack_trajectories_promoters_imputation_peaksplit.Rds"), (readRDS("noack_trajectories_enhancers_imputation_peaksplit.Rds")) )
write.table(trajectories_2, file = "~/openchr/noack/trajs.txt", sep = "\t", col.names=F, quote = F, row.names = F)
system('
       cat ~/openchr/noack/trajs.txt | tr "-" "\t" > ~/openchr/noack/trajs.bed
       bedtools intersect -wao -a ~/neurog2/chipseq/e14.5.neurog2_peaks.subpeaks.bed -b ~/openchr/noack/trajs.bed | cut -f5,9 > ~/openchr/noack/neurog2_trajs.tsv
')


# Pie charts with the subset of peaks that intersect with neurog2.

intersect<-bedr.join.region(bedr.sort.region(neurog2_peaks),bedr.sort.region( neurod2_peaks))
colnames(intersect)[5]<-"peak_id"
colnames(intersect)[6]<-"chr.b"
private<-intersect$peak_id[intersect$chr.b=="."]
shared<-intersect$peak_id[intersect$chr.b!="."]

neurg2_peaks_trajectories<-read.table("~/openchr/noack/neurdo2_trajs.tsv")
neurg2_peaks_trajectories$V2<-str_split_fixed(neurg2_peaks_trajectories$V2, "_", 2)[,1]
neurg2_peaks_trajectories<-neurg2_peaks_trajectories %>% filter(V2!=".") %>% distinct()  %>% filter(V1 %in% shared)
prop_neurog2<-table(neurg2_peaks_trajectories$V2)
orderr<-c("11111","10000","11000","11100","11110",
          "01000","01100","00100","01110","00110","00010",
          "01111","00111","00011","00001")
ordered<-orderr[which(orderr %in% names(prop_neurog2))]
prop_neurog2<-prop_neurog2[ordered]
colors<-c(
  "papayawhip",
              rep("darkgoldenrod1",4),rep("paleturquoise1",5),
              rep("lightgreen",3)
)
pie(prop_neurog2, main = "Proportions of trajectories" ,col = colors )


# Pie charts with the subset of peaks that intersect with neurod2.

intersect<-bedr.join.region(bedr.sort.region(neurod2_peaks),bedr.sort.region( neurog2_peaks))
colnames(intersect)[5]<-"peak_id"
colnames(intersect)[6]<-"chr.b"
private<-intersect$peak_id[intersect$chr.b=="."]
shared<-intersect$peak_id[intersect$chr.b!="."]

neurg2_peaks_trajectories<-read.table("~/openchr/noack/neurod_peak_travectories.tsv")
neurg2_peaks_trajectories$V2<-str_split_fixed(neurg2_peaks_trajectories$V2, "_", 2)[,1]
neurg2_peaks_trajectories<-neurg2_peaks_trajectories %>% filter(V2!=".") %>% distinct()  %>% filter(V1 %in% private)
prop_neurog2<-table(neurg2_peaks_trajectories$V2)
orderr<-c("11111","10000","11000","11100","11110",
          "01000","01100","00100","01110","00110","00010",
          "01111","00111","00011","00001")
ordered<-orderr[which(orderr %in% names(prop_neurog2))]
prop_neurog2<-prop_neurog2[ordered]
colors<-c(
  "papayawhip",
              rep("darkgoldenrod1",4),rep("paleturquoise1",6),
              rep("lightgreen",4)
)
pie(prop_neurog2, main = "Proportions of trajectories" ,col = colors )


par(family = "arial")


# Stacked barplot of the overlap of NEUROD2 and NEUROG2 with the ATAC-seq peaks grouped by trajectories.

setwd("~/openchr/noack")
trajs_enhancer<-readRDS("noack_trajectories_enhancers_imputation_peaksplit.Rds")
trajs_promoter<-readRDS("noack_trajectories_promoters_imputation_peaksplit.Rds")
trajs_enhancer$cluster<-paste0(trajs_enhancer$cluster,"_enhancer")
trajs_promoter$cluster<-paste0(trajs_promoter$cluster,"_promoter")
trajs_imputation<-rbind(trajs_enhancer,trajs_promoter)
write.table(trajs_imputation,"trajs_imputation.bed", sep = "\t", quote = F, col.names = F, row.names = F)

# Intersect the E14.5 NEUROD2 peaks with the trajectories.
system("
       bedtools intersect -wao -b <(cat trajs_imputation.bed | tr '-' '\t') -a ~/neurog2/chipseq/e14.5.neurog2_peaks.subpeaks.bed | cut -f5,9 | uniq > neurog2_trajs.tsv
")

# Get the ones that do not overlap.
system("
  bedtools intersect -v -a <(cat trajs_imputation.bed | tr '-' '\t') -b ~/neurog2/chipseq/e14.5.neurog2_peaks.subpeaks.bed > no_neurog2_trajs_imputation.bed
")


neurog2_peaks_trajectories<-read.table("~/openchr/noack/neurog2_trajs.tsv")
neurog2_peaks_trajectories<-neurog2_peaks_trajectories %>% dplyr::filter(V2 != ".")
library(stringr)
neurog2_peaks_trajectories<-str_split_fixed(neurog2_peaks_trajectories$V2, "_", 2)
neurog2_peaks_trajectories[,2]<-paste0(neurog2_peaks_trajectories[,2],"_neurog2")

no_neurog2_peaks_trajectories<-read.table("~/openchr/noack/no_neurog2_trajs_imputation.bed")
no_neurog2_peaks_trajectories<-str_split_fixed(no_neurog2_peaks_trajectories$V4, "_", 2)
no_neurog2_peaks_trajectories[,2]<-paste0(no_neurog2_peaks_trajectories[,2],"_no_neurod2")

peaks_trajectories<-as.data.frame(rbind(neurog2_peaks_trajectories,no_neurog2_peaks_trajectories))
colnames(peaks_trajectories)<-c("trajectory","type")
peaks_trajectories<-peaks_trajectories %>% group_by(trajectory,type) %>% summarise(n=n())


peaks_trajectories$trajectory <- factor(peaks_trajectories$trajectory, levels =orderr )
peaks_trajectories_promoter<-peaks_trajectories[grep("promoter",peaks_trajectories$type),]
peaks_trajectories_enhancer<-peaks_trajectories[grep("enhancer",peaks_trajectories$type),]


# Grouped
library(ggplot2)

ggplot(peaks_trajectories_promoter, aes(fill=type, y=n, x=trajectory)) + geom_bar(position="fill", stat="identity") +
  theme_classic() + scale_fill_manual(values = c("coral1","brown3","darkslategray3","deepskyblue4"))


ggplot(peaks_trajectories_enhancer, aes(fill=type, y=n, x=trajectory)) + geom_bar(position="fill", stat="identity") +
  theme_classic() + scale_fill_manual(values = c("darkslategray3","deepskyblue4"))

# The same but with the summarised trajectories.

peaks_trajectories<-as.data.frame(rbind(peaks_trajectories_promoter,peaks_trajectories_enhancer))
colnames(peaks_trajectories)<-c("trajectory","type","n")

peaks_trajectories$trajectory<-str_replace(peaks_trajectories$trajectory, "11110|11100|11000|10000", "early")
peaks_trajectories$trajectory<-str_replace(peaks_trajectories$trajectory, "01110|01100|00110|01000|00100|00010", "transcient")
peaks_trajectories$trajectory<-str_replace(peaks_trajectories$trajectory, "01111|00111|00011|00001", "neuronal")
peaks_trajectories$trajectory<-str_replace(peaks_trajectories$trajectory, "11111", "invariant")
peaks_trajectories$type<-str_replace(peaks_trajectories$type, "enhancer_", "")
peaks_trajectories$type<-str_replace(peaks_trajectories$type, "promoter_", "")
peaks_trajectories<-na.omit(peaks_trajectories) %>% group_by(trajectory,type) %>% summarise(n2=sum(n))
for (traj in unique(peaks_trajectories$trajectory)){
  jj<-peaks_trajectories %>% filter(trajectory==traj)
  print(pie(table(jj$type), main = traj))
}

peaks_trajectories$trajectory <- factor(peaks_trajectories$trajectory, levels =c("invariant","early","transcient","neuronal"))

ggplot(peaks_trajectories, aes(fill=type, y=n2, x=trajectory)) + geom_bar(position="fill", stat="identity") +
  theme_classic() + scale_fill_manual(values = c("coral1","brown3"))


# Grouped
library(ggplot2)
ggplot(peaks_trajectories_promoter, aes(fill=type, y=n, x=trajectory)) + geom_bar(position="fill", stat="identity") +
  theme_classic() + scale_fill_manual(values = c("coral1","brown3","darkslategray3","deepskyblue4"))


ggplot(peaks_trajectories_enhancer, aes(fill=type, y=n, x=trajectory)) + geom_bar(position="fill", stat="identity") +
  theme_classic() + scale_fill_manual(values = c("darkslategray3","deepskyblue4"))



### Fold change between the promoter and the enhancer.

nd2<-cbind(peaks_trajectories_promoter,peaks_trajectories_enhancer)
nd2<-nd2[-(grep("_no_",nd2$type...2)),]
nd2<-as.data.frame(nd2 %>% filter(!is.na(trajectory...1)))
rownames(nd2)<-nd2$trajectory...1
nd2<-nd2[ordered,]
nd2$logratio<-log(nd2$n...6/nd2$n...3)
barplot(nd2$logratio)



nd2_2<-rbind(peaks_trajectories_promoter,peaks_trajectories_enhancer)
nd2_2<-as.data.frame(nd2_2 %>% filter(!is.na(trajectory)))
nd2_2$trajj<-paste0(nd2_2$trajectory,"_",nd2_2$type)

nd2_2$trajj<-str_replace(nd2_2$trajj, "11110|11100|11000|10000", "early")
nd2_2$trajj<-str_replace(nd2_2$trajj, "01110|01100|00110|01000|00100|00010", "transcient")
nd2_2$trajj<-str_replace(nd2_2$trajj, "01111|00111|00011|00001", "neuronal")
nd2_2$trajj<-str_replace(nd2_2$trajj, "11111", "invariant")
nd2_2<-nd2_2 %>% group_by(trajj) %>% summarise(n=sum(n) )
vec_nd2<-nd2_2$n
names(vec_nd2)<-nd2_2$trajj
pie(vec_nd2)





## SANKEY PLOT POSTNATAL REGIONS ##################
###################################################


###### NEUROD2

setwd("~/neurod2/chipseq")
neurod2_peaks<-read.table("e14.5.neurod2_highconfsubpeaks_summit.bed")
colnames(neurod2_peaks)<-c("chr","start","end","summit","peak_id")
setwd("~/openchr/noack")
neurdo2_peaks_trajectories<-read.table("neurod_peak_travectories.tsv")
neurdo2_peaks_trajectories$V2<-str_replace(neurdo2_peaks_trajectories$V2, "11110|11100|11000|10000", "early")
neurdo2_peaks_trajectories$V2<-str_replace(neurdo2_peaks_trajectories$V2, "01110|01100|00110|01000|00100|00010", "transcient")
neurdo2_peaks_trajectories$V2<-str_replace(neurdo2_peaks_trajectories$V2, "01111|00111|00011|00001", "neuronal")
neurdo2_peaks_trajectories$V2<-str_replace(neurdo2_peaks_trajectories$V2, "11111", "invariant")
colnames(neurdo2_peaks_trajectories)<-c("peak_id","traj")
peaks_trajectories<-merge(neurod2_peaks,neurdo2_peaks_trajectories,by="peak_id")
peaks_trajectories<-peaks_trajectories %>% filter(traj!=".")
peaks_trajectories<-peaks_trajectories %>% filter(traj!="00000_enhancer")
peaks_trajectories<-peaks_trajectories %>% filter(traj!="00000_promoter")
neurod2_peaks_trajectories<-peaks_trajectories[,c(-1,-5)]

p0_neurod2_peaks<-read.table("~/neurod2/chipseq/p0.neurod2_highconfsubpeaks.bed")
#sacar cada trayectoria cuantos overlaps tiene con los picos de p0
require(bedr)
intersect<-bedr.join.region(bedr.sort.region(neurod2_peaks_trajectories),bedr.sort.region(p0_neurod2_peaks))
colnames(intersect)<-c("chr","start","end","traj","chr_p0","start_p0","end_p0")
overlap_yes<-intersect %>% dplyr::filter(chr_p0 != ".")
trajs_to_p0<-table(overlap_yes$traj)
overlap_no<-intersect %>% dplyr::filter(chr_p0 == ".")

#los picos de p56
p56_glut<-read.table("~/openchr/li/all.bed")
#sacar los q no overlapan con p0 si estan abiertos en p56 o no
overlap_no_p0_intersect_p56<-bedr.join.region(bedr.sort.region(overlap_no),bedr.sort.region(p56_glut))
#y aqui sacar la tabla de los q no overlapan y de los q si overlapan
overlap_no_p0_overlap_yes_p56<-overlap_no_p0_intersect_p56 %>% filter(V1!=".")
overlap_no_p0_overlap_no_p56<-overlap_no_p0_intersect_p56 %>% filter(V1==".")
trajs_to_open<-table(overlap_no_p0_overlap_yes_p56$traj)
trajs_to_closed<-table(overlap_no_p0_overlap_no_p56$traj)

#sacar los q van de p0 a p56
p0_intersect_p56<-bedr.join.region(bedr.sort.region(p0_neurod2_peaks),bedr.sort.region(p56_glut))
colnames(p0_intersect_p56)<-c("chr_p0","start_p0","end_p0","chr_p56","start_p56","end_p56")
p0_to_closed<-length(which(p0_intersect_p56$chr_p56=="."))
p0_to_open<-length(which(p0_intersect_p56$chr_p56!="."))


# Library
library(networkD3)
library(dplyr)

# A connection data frame is a list of flows with intensity for each flow
links <- data.frame(
  source=c(names(trajs_to_p0), names(trajs_to_open), names(trajs_to_closed), "P0 Neurod2","P0 Neurod2"), 
  target=c(rep("P0 Neurod2",length(names(trajs_to_p0))), rep("Open in P56 glutamatergic neurons",length(names(trajs_to_open))),
           rep("Closed in P56 glutamatergic neurons",length(names(trajs_to_closed))), "Open in P56 glutamatergic neurons",
           "Closed in P56 glutamatergic neurons"), 
  value=c(trajs_to_p0,trajs_to_open,trajs_to_closed,p0_to_open,p0_to_closed)
)

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

color_scale <- 
  "d3.scaleOrdinal()
     .domain(['invariant_promoter', 'invariant_enhancer', 'early_promoter', 'early_enhancer', 
              'transcient_promoter','transcient_enhancer','neuronal_promoter','neuronal_enhancer','P0 Neurod2'])
     .range(['papayawhip', 'papayawhip', 'orange', 'orange', 
     'lightblue', 'lightblue', 'lightgreen','lightgreen','black']);
  "
# Make the Network
p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=T, colourScale = color_scale)
p




####### BARPLOT POSTNATAL CUSANOVICH ###############
####################################################

intersect_neurod2_e14.5=c()
intersect_neurod2_p0=c()
total_cell=c()
celltypes<-c("10.1","1.1","11.1","11.2","11.3","11.4","11.5","12.1","12.2","12.3","12.4","12.5","13.1",
             "14.1","14.2","14.3","15.1","15.2","15.3","16.2","16.3","17.1","17.2","17.3","17.4","17.5","18.1","18.2","18.3",
             "18.4","18.5","19.1","19.2","19.3","19.4","20.1","2.1","21.1","21.2","22.1","22.2","22.3","22.4","23.1","23.2","24.1","24.2",
             "25.1","25.2","25.3","26.1","26.2","26.3","27.1","27.2","27.3","28.1","28.2","29.1","30.1","30.2","30.3","30.4","3.1","4.1","4.2",
             "4.3","4.4","5.1","5.2","5.3","5.4","5.5","5.6","6.1","7.1","7.2","8.1","8.2","9.2","9.3")

for (cell in celltypes) {
  print(cell)
  cell_ranges<-read.table(paste0("~/openchr/sub.clusters.specific.peaks.beds/mm10.",cell,".specific.regions"))
  intersect<-bedr.join.region(bedr.sort.region(cell_ranges),bedr.sort.region(neurod2_peaks))
  colnames(intersect)<-1:9
  number_intersect_neurod2<-length(which(intersect$`5`!="."))
  intersect_p0<-bedr.join.region(bedr.sort.region(cell_ranges),bedr.sort.region(p0_neurod2_peaks))
  colnames(intersect_p0)<-1:9
  number_intersect_neurod2_p0<-length(which(intersect_p0$`5`!="."))
  intersect_neurod2_e14.5<-c(intersect_neurod2_e14.5, number_intersect_neurod2)
  intersect_neurod2_p0<-c(intersect_neurod2_p0, number_intersect_neurod2_p0)
  total_cell<-c(total_cell, length(cell_ranges$V1))
}
neurod2_p56_cells_intersect<-as.data.frame(cbind(intersect_neurod2_e14.5,intersect_neurod2_p0,total_cell) )
rownames(neurod2_p56_cells_intersect)<-celltypes
neurod2_p56_cells_intersect$pct_e14.5<-neurod2_p56_cells_intersect$intersect_neurod2_e14.5/neurod2_p56_cells_intersect$total_cell
neurod2_p56_cells_intersect$pct_p0<-neurod2_p56_cells_intersect$intersect_neurod2_p0/neurod2_p56_cells_intersect$total_cell

library(reshape2)
li_intersect_df_melt<-setNames(melt(as.matrix(neurod2_p56_cells_intersect[,c(4,5)])), c('cluster', 'neurod2_dataset', 'percentage'))

ggplot(li_intersect_df_melt, aes(x=as.character(cluster),y=percentage, fill=neurod2_dataset)) + geom_bar(position="dodge", stat="identity") + coord_flip()

neurod2_p56_cells_intersect$jaccard_e14.5<-neurod2_p56_cells_intersect$intersect_neurod2_e14.5/(neurod2_p56_cells_intersect$total_cell+
                                                                                                  length(neurod2_peaks$V1)-
                                                                                                  neurod2_p56_cells_intersect$intersect_neurod2_e14.5)
neurod2_p56_cells_intersect$jaccard_p0<-neurod2_p56_cells_intersect$intersect_neurod2_p0/(neurod2_p56_cells_intersect$total_cell+
                                                                                            length(p0_neurod2_peaks$V1)-
                                                                                            neurod2_p56_cells_intersect$intersect_neurod2_p0)
#escalar los jaccards
neurod2_p56_cells_intersect$jaccard_e14.5<-neurod2_p56_cells_intersect$jaccard_e14.5/max(neurod2_p56_cells_intersect$jaccard_e14.5)
neurod2_p56_cells_intersect$jaccard_p0<-neurod2_p56_cells_intersect$jaccard_p0/max(neurod2_p56_cells_intersect$jaccard_p0)

li_intersect_df_melt<-setNames(melt(as.matrix(neurod2_p56_cells_intersect[,c(6,7)])), c('cluster', 'neurod2_dataset', 'jaccard'))

ggplot(li_intersect_df_melt, aes(x=reorder(cluster,jaccard,
                                           function(x)+max(x)),
                                 y=jaccard, fill=neurod2_dataset)) + 
  geom_bar(position="dodge", stat="identity") + coord_flip() + theme_classic() +  scale_fill_manual(values=c("goldenrod", "royalblue" ))



intersect_neurod2_e14.5=c()
intersect_neurod2_p0=c()
total_cell=c()
celltypes<-c("atac.sc.adult.Astrocytes.19.1.cusanovich2018.bed","atac.sc.adult.Astrocytes.19.2.cusanovich2018.bed","atac.sc.adult.Astrocytes.19.3.cusanovich2018.bed","atac.sc.adult.Astrocytes.19.4.cusanovich2018.bed",
             "atac.sc.adult.Cerebellar_granule_cells.8.1.cusanovich2018.bed","atac.sc.adult.Cerebellar_granule_cells.8.2.cusanovich2018.bed","atac.sc.adult.Ex._neurons_CPN.5.1.cusanovich2018.bed","atac.sc.adult.Ex._neurons_CThPN.5.3.cusanovich2018.bed","atac.sc.adult.Ex._neurons_CThPN.5.4.cusanovich2018.bed","atac.sc.adult.Ex._neurons_SCPN.29.1.cusanovich2018.bed",
             "atac.sc.adult.Ex._neurons_SCPN.5.2.cusanovich2018.bed","atac.sc.adult.Inhibitory_neurons.15.1.cusanovich2018.bed","atac.sc.adult.Inhibitory_neurons.15.2.cusanovich2018.bed","atac.sc.adult.Inhibitory_neurons.5.5.cusanovich2018.bed","atac.sc.adult.Microglia.16.3.cusanovich2018.bed","atac.sc.adult.Oligodendrocytes.21.1.cusanovich2018.bed","atac.sc.adult.Oligodendrocytes.21.2.cusanovich2018.bed",
             "atac.sc.adult.Purkinje_cells.27.1.cusanovich2018.bed","atac.sc.adult.SOM.plus._Interneurons.15.3.cusanovich2018.bed")

for (cell in celltypes) {
  print(cell)
  cell_ranges<-read.table(paste0("~/openchr/neural.clusters.open.regions.beds/",cell))
  intersect<-bedr.join.region(bedr.sort.region(cell_ranges),bedr.sort.region(neurod2_peaks))
  colnames(intersect)<-1:8
  number_intersect_neurod2<-length(which(intersect$`4`!="."))
  intersect_p0<-bedr.join.region(bedr.sort.region(cell_ranges),bedr.sort.region(p0_neurod2_peaks))
  colnames(intersect_p0)<-1:8
  number_intersect_neurod2_p0<-length(which(intersect_p0$`4`!="."))
  intersect_neurod2_e14.5<-c(intersect_neurod2_e14.5, number_intersect_neurod2)
  intersect_neurod2_p0<-c(intersect_neurod2_p0, number_intersect_neurod2_p0)
  total_cell<-c(total_cell, length(cell_ranges$V1))
}
neurod2_p56_cells_intersect<-as.data.frame(cbind(intersect_neurod2_e14.5,intersect_neurod2_p0,total_cell) )

rownames(neurod2_p56_cells_intersect)<-celltypes
neurod2_p56_cells_intersect$pct_e14.5<-neurod2_p56_cells_intersect$intersect_neurod2_e14.5/neurod2_p56_cells_intersect$total_cell
neurod2_p56_cells_intersect$pct_p0<-neurod2_p56_cells_intersect$intersect_neurod2_p0/neurod2_p56_cells_intersect$total_cell



li_intersect_df_melt<-setNames(melt(as.matrix(neurod2_p56_cells_intersect[,c(4,5)])), c('cluster', 'neurod2_dataset', 'percentage'))

ggplot(li_intersect_df_melt, aes(x=as.character(cluster),y=percentage, fill=neurod2_dataset)) + geom_bar(position="dodge", stat="identity") + coord_flip()



neurod2_p56_cells_intersect$jaccard_e14.5<-neurod2_p56_cells_intersect$intersect_neurod2_e14.5/(neurod2_p56_cells_intersect$total_cell+
                                                                                                  length(neurod2_peaks$V1)-
                                                                                                  neurod2_p56_cells_intersect$intersect_neurod2_e14.5)
neurod2_p56_cells_intersect$jaccard_p0<-neurod2_p56_cells_intersect$intersect_neurod2_p0/(neurod2_p56_cells_intersect$total_cell+
                                                                                            length(p0_neurod2_peaks$V1)-
                                                                                            neurod2_p56_cells_intersect$intersect_neurod2_p0)
li_intersect_df_melt<-setNames(melt(as.matrix(neurod2_p56_cells_intersect[,c(6,7)])), c('cluster', 'neurod2_dataset', 'jaccard'))

ggplot(li_intersect_df_melt, aes(x=as.character(cluster),y=jaccard, fill=neurod2_dataset)) + geom_bar(position="dodge", stat="identity") + coord_flip()




