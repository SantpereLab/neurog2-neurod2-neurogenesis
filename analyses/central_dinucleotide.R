# CENTRAL DINUCLEOTIDE ######################
##############################################

###### NEUROD2
# Detect motifs in NEUROD2 peaks
neurod2_peaks<-read.table("e14.5.neurod2_highconfsubpeaks_summit.bed")
colnames(neurod2_peaks)<-c("chr","start","end","summit","peak_id")
neurod2_peaks<-find_motifs(neurod2_peaks, peaks_file = "e14.5.neurod2_highconfsubpeaks_summit.bed")

# Detect motifs in NEUROD2 extended peaks
neurod2_extended_peaks<-read.table("e14.5.neurod2_highconfsubpeaks_summit.bed")
colnames(neurod2_extended_peaks)<-c("chr","start","end","summit","peak_id")
neurod2_extended_peaks$start<-as.numeric(neurod2_extended_peaks$summit)-400
neurod2_extended_peaks$end<-as.numeric(neurod2_extended_peaks$summit)+400

write.table(neurod2_extended_peaks, file = "e14.5.neurod2_highconfsubpeaks_extended.bed", sep = "\t", col.names=F, quote = F, row.names = F )
neurod2_extended_peaks_motifs<-rbind( find_motifs (peaks=neurod2_extended_peaks, peaks_file = "e14.5.neurod2_highconfsubpeaks_extended.bed" ) )


# Plot the distribution of the central dinucleotides. 
colors<-c("deeppink", "yellow","azure4", "seagreen3", "darkgoldenrod4", "chartreuse1","blue","turquoise2","red", "black")

nd2_distribution<-distribution_motifs(peaks=neurod2_extended_peaks_motifs,lim=length(neurod2_extended_peaks$peak_id)/15,colors=colors,binwidth=7,main="Neurod2")


###### NEUROG2
# Detect motifs.
setwd("~/neurog2/chipseq")
neurog2_peaks<-read.table("e14.5.neurog2_peaks.subpeaks.bed")
colnames(neurog2_peaks)<-c("chr","start","end","summit","peak_id")
neurog2_peaks_motifs<-find_motifs(neurog2_peaks,peaks_file = "e14.5.neurog2_peaks.subpeaks.bed")
saveRDS(neurog2_peaks_motifs,"neurog2_peaks_motifs.Rds")

# Detect motifs in the extended peaks.
neurog2_extended_peaks<-neurog2_peaks
neurog2_extended_peaks$start<-as.numeric(neurog2_peaks$summit)-400
neurog2_extended_peaks$end<-as.numeric(neurog2_peaks$summit)+400
write.table(neurog2_peaks, file = "e14.5.neurog2_peaks.subpeaks_extended.bed", sep = "\t", col.names=F, quote = F, row.names = F)

neurog2_extended_peaks_motifs<-find_motifs(peaks=neurog2_extended_peaks,peaks_file = "e14.5.neurog2_peaks.subpeaks_extended.bed")

#Plot the distribution around the summit.
colors<-c("deeppink", "yellow","azure4", "seagreen3", "darkgoldenrod4", "chartreuse1","blue", "turquoise2","red", "black")

ng2_distribution<-distribution_motifs(peaks=neurog2_extended_peaks_motifs,lim=length(neurog2_extended_peaks$peak_id)/15,colors=colors,binwidth=7, main = "Neurog2")

pdf("~/central_dinucleotide_distributino.pdf", useDingbats = F, height = 3, width = 3)
print(nd2_distribution)
print(ng2_distribution)
dev.off()


#### P0 NEUROD2
# Detect motifs.
p0_neurod2_peaks<-read.table("p0.neurod2_highconfsubpeaks_summit.bed")
p0_neurod2_peaks$summit<-p0_neurod2_peaks$V2+p0_neurod2_peaks$V10
p0_neurod2_peaks<- p0_neurod2_peaks %>% select(V1, V2, V3, V4, summit )
colnames(p0_neurod2_peaks)<-c("chr","start","end","peak_id","summit")
p0_neurod2_peaks_motifs<-find_motifs(p0_neurod2_peaks,"p0.neurod2_highconfsubpeaks_summit.bed")

# Detect motifs in extended peaks.
p0_neurod2_extended_peaks<-p0_neurod2_peaks
p0_neurod2_extended_peaks$start<-p0_neurod2_peaks$summit-400
p0_neurod2_extended_peaks$end<-p0_neurod2_peaks$summit+400
write.table(p0_neurod2_extended_peaks, file = "p0.neurod2_highconfsubpeaks_extended.bed", sep = "\t", col.names=F, quote = F, row.names = F)

p0_neurod2_extended_peaks_motifs<-find_motifs(peaks=p0_neurod2_extended_peaks,peaks_file = "p0.neurod2_highconfsubpeaks_extended.bed")

# PLot the distribution.
colors<-c("deeppink", "yellow","azure4", "seagreen3", "darkgoldenrod4", "chartreuse1","blue","turquoise2","red", "black")

distribution_motifs(peaks=p0_neurod2_extended_peaks_motifs,lim=length(p0_neurod2_peaks$peak_id)/15,colors=colors,binwidth=1, "P0 Neurod2 peaks")

çººººº
neurod2_trajs_motifs<-sumarizacion_trajs_generales(extended_motifs = neurod2_extended_peaks_motifs,  frontera1 = 50, frontera2 = 200, ids=neurod2_peaks_trajectories)

ggballoonplot(neurod2_trajs_motifs, x="traj",
              y = "central_dinucleotide", size = "nmotifs", fill = "centrality") +  gradient_fill( c("black","darkred","yellow2") )

neurog2_trajs_motifs<-sumarizacion_trajs_generales(extended_motifs = neurog2_extended_peaks_motifs,  frontera1 = 50, frontera2 = 200, ids=neurog2_peaks_trajectories)

ggballoonplot(neurog2_trajs_motifs, x="traj",
              y = "central_dinucleotide", size = "nmotifs", fill = "centrality") +  gradient_fill( c("black","darkred","yellow2") )
