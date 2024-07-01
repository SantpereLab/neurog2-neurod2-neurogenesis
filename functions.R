# FUNCTIONS #################
###############################################

## FUNCION PARA SACARTE LOS MOTIVOS EN LOS PEAKS
find_motifs<-function(peaks, peaks_file) {
  require(dplyr)
  require(bedr)
  require(pracma)
  options(scipen = 999)
  write.table(peaks,peaks_file,sep="\t",quote=F, col.names=F, row.names = F )
  system(paste0("bedtools getfasta -fi ~/genomes/mm10.fa -bed ", peaks_file, " | grep -v chr > ", peaks_file,".fa") )
  seqs<-read.table(paste0(peaks_file,".fa"), fill = TRUE , header = FALSE)
  seqs$V1<-toupper(seqs$V1)
  
  ebox<-data.frame()
  for (i in 1:length(seqs$V1)){
    print(i)
    ebox_chr<-peaks$chr[i]
    matches<-refindall(seqs$V1[i],'..CA..TG..')
    if (length(matches)==0){
      next
    }
    for (j in 1:length(matches)){
      ebox_start<-peaks$start[i]+matches[j]
      ebox_end<-ebox_start+10
      ebox_seq<-substr(seqs$V1[i],matches[j],matches[j]+9)
      ebox<-rbind(ebox,c(ebox_chr,ebox_start,ebox_end,ebox_seq))
    }
  }
  colnames(ebox)<-c("ebox_chr","ebox_start","ebox_end","ebox_seq")
  ebox$ebox_start<-as.numeric(ebox$ebox_start)
  ebox$ebox_end<-as.numeric(ebox$ebox_end)
  
  #sacar los dinucleotides
  ebox$central_dinucleotide<-substr(ebox$ebox_seq,5,6)
  #pairear los dinucleos
  paired_dinucleotides<-c("GA/TC","GG/CC","AA/TT","CA/TG","AC/GT","AG/CT","GC","CG","TA","AT")
  for (i in 1:length(ebox$central_dinucleotide)) {
    if (is.na(ebox$central_dinucleotide[i])) { next }
    ebox$central_dinucleotide[i]<-paired_dinucleotides[grep(ebox$central_dinucleotide[i],paired_dinucleotides)]
  }
  #cambiar la notacion de los dinucleotides por la buena
  ebox$central_dinucleotide[ebox$central_dinucleotide=="GC"]<-"CAG-CAG"
  ebox$central_dinucleotide[ebox$central_dinucleotide=="CG"]<-"CAC-CAC"
  ebox$central_dinucleotide[ebox$central_dinucleotide=="GA/TC"]<-"CAT-CAG"
  ebox$central_dinucleotide[ebox$central_dinucleotide=="CA/TG"]<-"CAT-CAC"
  ebox$central_dinucleotide[ebox$central_dinucleotide=="TA"]<-"CAT-CAT"
  ebox$central_dinucleotide[ebox$central_dinucleotide=="AA/TT"]<-"CAA-CAT"
  ebox$central_dinucleotide[ebox$central_dinucleotide=="AC/GT"]<-"CAA-CAG"
  ebox$central_dinucleotide[ebox$central_dinucleotide=="AG/CT"]<-"CAA-CAC"
  ebox$central_dinucleotide[ebox$central_dinucleotide=="AT"]<-"CAA-CAA"
  ebox$central_dinucleotide[ebox$central_dinucleotide=="GG/CC"]<-"CAG-CAC"
  print(length(ebox$chr))
  peaks_motifs<-bedr::bedr.join.region(bedr.sort.region(peaks),bedr.sort.region(ebox))
  system("rm *fa")
  peaks_motifs<- peaks_motifs %>% dplyr::select("chr","start","end","summit","peak_id","ebox_start.b","ebox_seq","central_dinucleotide")
  names(peaks_motifs)[6]<-"ebox_start"
  return(peaks_motifs)
}


find_spacing<-function(peaks, peaks_file, spacing) {
  require(dplyr)
  require(bedr)
  require(pracma)
  system(paste0("bedtools getfasta -fi ~/genomes/mm10.fa -bed ", peaks_file, " | grep -v chr > ", peaks_file,".fa") )
  seqs<-read.table(paste0(peaks_file,".fa"))
  seqs$V1<-toupper(seqs$V1)
  ebox<-data.frame()
  for (i in 1:length(seqs$V1)){
    print(i)
    ebox_chr<-peaks$chr[i]
    motif<-paste0("CA..TG",paste(rep(".",spacing),collapse = ""),"CA..TG")
    matches<-refindall(seqs$V1[i],motif)
    if (length(matches)==0){
      next
    }
    for (j in 1:length(matches)){
      spacing_middle<-peaks$start[i]+matches[j]+6+round(spacing/2)
      spacing_middle_2<-spacing_middle-1
    }
    
    if(length(ebox)==0){
      ebox<-data.frame("chr"=ebox_chr,"start"=spacing_middle_2,"end"=spacing_middle)
    }
    else {
      ebox<-rbind(ebox,c(ebox_chr,spacing_middle_2,spacing_middle))
    }
  }
  system("rm *fa")
  return(ebox)
}

# lo mismo solo que te devuelve el spacing entero
find_spacing_2<-function(peaks, peaks_file, spacing) {
  require(dplyr)
  require(bedr)
  require(pracma)
  system(paste0("bedtools getfasta -fi ~/genomes/mm10.fa -bed ", peaks_file, " | grep -v chr > ", peaks_file,".fa") )
  seqs<-read.table(paste0(peaks_file,".fa"))
  seqs$V1<-toupper(seqs$V1)
  
  ebox<-data.frame()
  for (i in 1:length(seqs$V1)){
    print(i)
    ebox_chr<-peaks$chr[i]
    motif<-paste0("CA..TG",paste(rep(".",spacing),collapse = ""),"CA..TG")
    matches<-refindall(seqs$V1[i],motif)
    if (length(matches)==0){
      next
    }
    for (j in 1:length(matches)){
      spacing_start<-matches[j]+6+round(spacing/2)
      spacing_end<-matches[j]+6+spacing+6
    }
    ebox<-rbind(ebox,c(ebox_chr,spacing_middle_2,spacing_middle))
  }
  system("rm *fa")
  return(ebox)
}

find_noncanonical<-function(peaks, peaks_file, motif) {
  require(dplyr)
  require(bedr)
  require(pracma)
  system(paste0("bedtools getfasta -fi ~/genomes/mm10.fa -bed ", peaks_file, " | grep -v chr > ", peaks_file,".fa") )
  seqs<-read.table(paste0(peaks_file,".fa"))
  seqs$V1<-toupper(seqs$V1)
  
  ebox<-data.frame()
  for (i in 1:8000){
    print(i)
    ebox_chr<-peaks$chr[i]
    matches<-refindall(seqs$V1[i], motif)
    if (length(matches)==0){
      next
    }
    for (j in 1:length(matches)){
      ebox_start<-peaks$start[i]+matches[j]
      ebox_end<-ebox_start+10
      ebox_seq<-substr(seqs$V1[i],matches[j],matches[j]+9)
      ebox<-rbind(ebox,c(ebox_chr,ebox_start,ebox_end,ebox_seq))
    }
  }
  colnames(ebox)<-c("ebox_chr","ebox_start","ebox_end","ebox_seq")
  ebox$ebox_start<-as.numeric(ebox$ebox_start)
  ebox$ebox_end<-as.numeric(ebox$ebox_end)
  return(ebox)
}


#FUNCION PARA SACARTE LA DISTRIBUCION DE LOS MOTIVOS AL REDEDOR DEL SUMMIT
distribution_motifs<-function(peaks,lim,colors,binwidth,main){
  require(ggplot2)
  peaks$motif_center<-as.numeric(peaks$ebox_start)+5
  peaks$motif_distance_summit<-peaks$motif_center-as.numeric(peaks$summit)
  peaks<-peaks %>% filter(central_dinucleotide != ".")
  ggplot(peaks, aes(x=motif_distance_summit, color = central_dinucleotide)) + 
    geom_freqpoly(binwidth=binwidth) + ylim(0,lim) +
    scale_color_manual(values=colors) +
    scale_x_continuous(breaks = round(seq(-400, 400, by = 400),1)) +
    theme_classic() + ggtitle(main) + xlab("") + ylab("") 
}

###### balloon plot de los central dinucleotides

sumarizacion<-function(extended_motifs, frontera1,frontera2) {
  npeaks<-length(unique(extended_motifs$peak_id))
  extended_motifs$dist_summit<-abs(as.numeric(extended_motifs$summit)-as.numeric(extended_motifs$ebox_start)+5)
  extended_motifs_summit<-extended_motifs %>% filter(dist_summit<frontera1) %>% filter(central_dinucleotide != ".") 
  extended_motifs_background<-extended_motifs %>% filter(dist_summit>frontera2) %>% filter(central_dinucleotide != ".")
  
  extended_motifs_summit<-extended_motifs_summit %>% group_by(central_dinucleotide) %>% summarise(nmotifs=n() )
  #añadir 0s donde podria no haber nada
  extra_motifs<-c("CAT-CAT","CAT-CAG","CAG-CAG","CAG-CAC","CAT-CAC","CAC-CAC","CAA-CAA","CAA-CAG","CAA-CAT","CAA-CAC")
  extra_nmotifs<-rep(0, 10 )
  extra<-data.frame("central_dinucleotide"=extra_motifs,"nmotifs"=extra_nmotifs)
  extended_motifs_summit<-rbind(extended_motifs_summit,extra) %>% group_by(central_dinucleotide) %>% summarise(nmotifs=sum(nmotifs) )
  extended_motifs_background<-extended_motifs_background %>% group_by(central_dinucleotide) %>% summarise(nmotifs=n() )
  
  extended_motifs_summit$centrality<-extended_motifs_summit$nmotifs/extended_motifs_background$nmotifs
  extended_motifs_summit$nmotifs<-extended_motifs_summit$nmotifs/npeaks
  return(extended_motifs_summit)
}

sumarizacion_trajs<-function(extended_motifs, frontera1,frontera2, ids) {
  colnames(ids)<-c("peak_id","traj")
  extended_motifs<-merge(extended_motifs,ids,by="peak_id")
  npeaks<-ids %>% group_by(traj) %>% summarise(npeaks=n())
  extended_motifs$dist_summit<-abs(as.numeric(extended_motifs$summit)-as.numeric(extended_motifs$ebox_start)+5)
  extended_motifs_summit<-extended_motifs %>% filter(dist_summit<frontera1) %>% filter(central_dinucleotide != ".") 
  extended_motifs_background<-extended_motifs %>% filter(dist_summit>frontera2) %>% filter(central_dinucleotide != ".")
  
  extended_motifs_summit<-extended_motifs_summit %>% group_by(central_dinucleotide,traj) %>% summarise(nmotifs=n() )
  #añadir 0s donde podria no haber nada
  extra_motifs<-rep(c("CAT-CAT","CAT-CAG","CAG-CAG","CAG-CAC","CAT-CAC","CAC-CAC","CAA-CAA","CAA-CAG","CAA-CAT","CAA-CAC"), 15 )
  extra_nmotifs<-rep(0, 150 )
  extra_trajs<-rep( unique(ids$traj), each=10)
  extra<-data.frame("central_dinucleotide"=extra_motifs,"nmotifs"=extra_nmotifs, "traj"=extra_trajs)
  extended_motifs_summit<-rbind(extended_motifs_summit,extra) %>% group_by(central_dinucleotide,traj) %>% summarise(nmotifs=sum(nmotifs) )
  extended_motifs_background<-extended_motifs_background %>% group_by(central_dinucleotide,traj) %>% summarise(nmotifs=n() )
  extended_motifs_background<-rbind(extended_motifs_background,extra)
  extended_motifs_background<-extended_motifs_background %>% group_by(central_dinucleotide,traj) %>% summarise(nmotifs=sum(nmotifs) )
  
  
  extended_motifs_summit$centrality<-extended_motifs_summit$nmotifs/extended_motifs_background$nmotifs
  extended_motifs_summit$centrality[is.nan(extended_motifs_summit$centrality)]<-0
  npeakssss<-rep(npeaks$npeaks, 10)
  extended_motifs_summit$nmotifs<-extended_motifs_summit$nmotifs/npeakssss
  return(extended_motifs_summit)
}

sumarizacion_trajs_neurog2<-function(extended_motifs, frontera1,frontera2, ids) {
  colnames(ids)<-c("peak_id","traj")
  extended_motifs<-merge(extended_motifs,ids,by="peak_id")
  npeaks<-ids %>% group_by(traj) %>% summarise(npeaks=n())
  extended_motifs$dist_summit<-abs(as.numeric(extended_motifs$summit)-as.numeric(extended_motifs$ebox_start)+5)
  extended_motifs_summit<-extended_motifs %>% filter(dist_summit<frontera1) %>% filter(central_dinucleotide != ".") 
  extended_motifs_background<-extended_motifs %>% filter(dist_summit>frontera2) %>% filter(central_dinucleotide != ".")
  
  extended_motifs_summit<-extended_motifs_summit %>% group_by(central_dinucleotide,traj) %>% summarise(nmotifs=n() )
  #añadir 0s donde podria no haber nada
  extra_motifs<-rep(c("CAT-CAT","CAT-CAG","CAG-CAG","CAG-CAC","CAT-CAC","CAC-CAC","CAA-CAA","CAA-CAG","CAA-CAT","CAA-CAC"), 13 )
  extra_nmotifs<-rep(0, 130 )
  extra_trajs<-rep( unique(ids$traj), each=10)
  extra<-data.frame("central_dinucleotide"=extra_motifs,"nmotifs"=extra_nmotifs, "traj"=extra_trajs)
  extended_motifs_summit<-rbind(extended_motifs_summit,extra) %>% group_by(central_dinucleotide,traj) %>% summarise(nmotifs=sum(nmotifs) )
  extended_motifs_background<-extended_motifs_background %>% group_by(central_dinucleotide,traj) %>% summarise(nmotifs=n() )
  extended_motifs_background<-rbind(extended_motifs_background,extra)
  extended_motifs_background<-extended_motifs_background %>% group_by(central_dinucleotide,traj) %>% summarise(nmotifs=sum(nmotifs) )
  
  
  extended_motifs_summit$centrality<-extended_motifs_summit$nmotifs/extended_motifs_background$nmotifs
  extended_motifs_summit$centrality[is.nan(extended_motifs_summit$centrality)]<-0
  npeakssss<-rep(npeaks$npeaks, 10)
  extended_motifs_summit$nmotifs<-extended_motifs_summit$nmotifs/npeakssss
  return(extended_motifs_summit)
}

sumarizacion_trajs_generales<-function(extended_motifs, frontera1,frontera2, ids) {
  colnames(ids)<-c("peak_id","traj")
  ids$traj<-str_replace(ids$traj, "11110|11100|11000|10000", "2.early")
  ids$traj<-str_replace(ids$traj, "01110|01100|00110|01000|00100|00010", "3.transcient")
  ids$traj<-str_replace(ids$traj, "01111|00111|00011|00001", "4.neuronal")
  ids$traj<-str_replace(ids$traj, "11111", "1.invariant")
  
  extended_motifs<-merge(extended_motifs,ids,by="peak_id")
  npeaks<-ids %>% group_by(traj) %>% summarise(npeaks=n())
  extended_motifs$dist_summit<-abs(as.numeric(extended_motifs$summit)-as.numeric(extended_motifs$ebox_start)+5)
  extended_motifs_summit<-extended_motifs %>% filter(dist_summit<frontera1) %>% filter(central_dinucleotide != ".") 
  extended_motifs_background<-extended_motifs %>% filter(dist_summit>frontera2) %>% filter(central_dinucleotide != ".")
  
  extended_motifs_summit<-extended_motifs_summit %>% group_by(central_dinucleotide,traj) %>% summarise(nmotifs=n() )
  #añadir 0s donde podria no haber nada
  extra_motifs<-rep(c("CAT-CAT","CAT-CAG","CAG-CAG","CAG-CAC","CAT-CAC","CAC-CAC","CAA-CAA","CAA-CAG","CAA-CAT","CAA-CAC"), 4 )
  extra_nmotifs<-rep(0,40)
  extra_trajs<-rep( unique(ids$traj), each=10)
  extra<-data.frame("central_dinucleotide"=extra_motifs,"nmotifs"=extra_nmotifs, "traj"=extra_trajs)
  extended_motifs_summit<-rbind(extended_motifs_summit,extra) %>% group_by(central_dinucleotide,traj) %>% summarise(nmotifs=sum(nmotifs) )
  extended_motifs_background<-extended_motifs_background %>% group_by(central_dinucleotide,traj) %>% summarise(nmotifs=n() )
  extended_motifs_background<-rbind(extended_motifs_background,extra)
  extended_motifs_background<-extended_motifs_background %>% group_by(central_dinucleotide,traj) %>% summarise(nmotifs=sum(nmotifs) )
  
  
  extended_motifs_summit$centrality<-extended_motifs_summit$nmotifs/extended_motifs_background$nmotifs
  extended_motifs_summit$centrality[is.nan(extended_motifs_summit$centrality)]<-0
  npeakssss<-rep(npeaks$npeaks, 10)
  extended_motifs_summit$nmotifs<-extended_motifs_summit$nmotifs/npeakssss
  return(extended_motifs_summit)
}


#FUNCION PARA SACARTE BARPLOTS DE LOS MOTIVOS
barplots_motifs<-function(peaks,colors,main,dist_summit_filter, lim){
  require(ggplot2)
  peaks$motif_center<-as.numeric(peaks$ebox_start)+5
  npeaks<-length(unique(peaks$peak_id))
  #filtrar distancia al rededor del summit
  peaks$motif_distance_summit<-peaks$motif_center-as.numeric(peaks$summit)
  peaks<-peaks %>% filter(central_dinucleotide != ".")
  peaks<-peaks %>% filter(abs(motif_distance_summit) <= dist_summit_filter)
  #contar el numero de dinucleotidos
  dinu_count<-peaks %>% group_by(central_dinucleotide) %>% summarise(n=n())
  #normalizar por el numero de picos x el numero de bp de ventana al rededor del summit
  dinu_count$n<-dinu_count$n/(npeaks)
  colnames(dinu_count)<-c("motif","counts")
  #añadir 0s si no esta el motivo dado
  añadido<-data.frame("motif"=c("CAT-CAT","CAT-CAG","CAG-CAG","CAG-CAC","CAT-CAC","CAC-CAC","CAA-CAA","CAA-CAG","CAA-CAT","CAA-CAC"),"counts"=rep(0,10))
  dinu_count<-rbind(dinu_count,añadido)
  dinu_count<-dinu_count %>% group_by(motif) %>% summarise(`motifs/peak`=sum(counts) )
  p<-ggplot(dinu_count, aes(x=motif, y=`motifs/peak`, fill=motif)) +
    geom_bar(stat="identity")+theme_minimal() + scale_fill_manual(values=colors) + ylim(0,lim) + ggtitle(main) + theme(legend.position = "none" )
  print(p)
}


#LO MISMO PERO EN VIOLIN PLOTS
violin_motifs<-function(peaks){
  peaks<-peaks %>% filter(central_dinucleotide != ".")
  peaks$motif_center<-as.numeric(peaks$ebox_start.b)+5
  peaks$motif_distance_summit<-abs(peaks$motif_center-as.numeric(peaks$summit))
  getPalette = colorRampPalette(brewer.pal(9, "Set1"))
  ggplot(peaks, aes(x=reorder(central_dinucleotide,motif_distance_summit), y=motif_distance_summit, fill=central_dinucleotide)) +
    geom_violin(trim=T) + theme_classic() + scale_fill_manual(values = getPalette(10)) 
}


#funcion para sacar los flanking dinucleotides
flanking_dinucleotides<-function(peaks){
  #hacer el revcomp para poner todos los motivos en la misma orientacion
  revcomp<-peaks$ebox_seq
  revcomp[grep("..CAGATG..",revcomp)]<-stri_reverse(chartr("ATGC","TACG",
                                                           revcomp[grep("..CAGATG..",revcomp)]) )
  revcomp[grep("..CACATG..",revcomp)]<-stri_reverse(chartr("ATGC","TACG",
                                                           revcomp[grep("..CACATG..",revcomp)]) )
  revcomp[grep("..CACCTG..",revcomp)]<-stri_reverse(chartr("ATGC","TACG",
                                                           revcomp[grep("..CACCTG..",revcomp)]) )
  #sacar los flanking
  peaks$upstream_dinucleotide<-substr(revcomp,1,2)
  peaks$downstream_dinucleotide<-substr(revcomp,9,10)
  return(peaks)
}

#balloon plot que represente el porcentaje de los dinucleotidos + el fold change sobre los porcentajes del background
flanking_motifs_balloon<-function(peaks,distance_around_summit,lim,back_dist) {
  peaks$motif_center<-as.numeric(peaks$ebox_start)+5
  peaks$motif_distance_summit<-peaks$motif_center-as.numeric(peaks$summit)
  
  #cat half-site
  cat_cat_upstream<-peaks %>% filter(central_dinucleotide=="CAT-CAT") %>% select(central_dinucleotide, motif_distance_summit, upstream_dinucleotide)
  cat_cat_downstream<-peaks %>% filter(central_dinucleotide=="CAT-CAT") %>% select(central_dinucleotide, motif_distance_summit, downstream_dinucleotide)
  #hay que hacer el revcomp en los downstream
  cat_cat_downstream$downstream_dinucleotide<-stri_reverse(chartr("ATGC","TACG", cat_cat_downstream$downstream_dinucleotide) )
  cat_cag_upstream<-peaks %>% filter(central_dinucleotide=="CAT-CAG") %>% select(central_dinucleotide, motif_distance_summit, upstream_dinucleotide)
  cat_cac_upstream<-peaks %>% filter(central_dinucleotide=="CAT-CAC") %>% select(central_dinucleotide, motif_distance_summit, upstream_dinucleotide)
  #juntarlos
  colnames(cat_cat_upstream)[3]<-"dinucleotide"
  colnames(cat_cat_downstream)[3]<-"dinucleotide"
  colnames(cat_cag_upstream)[3]<-"dinucleotide"
  colnames(cat_cac_upstream)[3]<-"dinucleotide"
  cat_flanking<-rbind(cat_cat_upstream,cat_cat_downstream,cat_cag_upstream,cat_cac_upstream)
  
  #sacar los porcentajes de los flanking dinucleotides
  pcts<-cat_flanking %>% filter(abs(motif_distance_summit)<distance_around_summit) %>% 
    group_by(central_dinucleotide,dinucleotide) %>% summarise(n=n()) %>% 
    mutate(pct=n/sum(n))
  pcts_background<-cat_flanking %>% filter(abs(motif_distance_summit)>back_dist) %>% 
    group_by(central_dinucleotide,dinucleotide) %>% summarise(n=n()) %>% 
    mutate(pct=n/sum(n))
  pcts$central_dinucleotide<-paste0("1   CAT : ",pcts$central_dinucleotide)
  pcts_background$central_dinucleotide<-paste0("1   CAT : ",pcts_background$central_dinucleotide)
  
  wonk<-pcts
  wonk_background<-pcts_background
  
  #cag half-site
  cag_cag_upstream<-peaks %>% filter(central_dinucleotide=="CAG-CAG") %>% select(central_dinucleotide, motif_distance_summit, upstream_dinucleotide)
  cag_cac_upstream<-peaks %>% filter(central_dinucleotide=="CAG-CAC") %>% select(central_dinucleotide, motif_distance_summit, upstream_dinucleotide)
  cag_cag_downstream<-peaks %>% filter(central_dinucleotide=="CAG-CAG") %>% select(central_dinucleotide, motif_distance_summit, downstream_dinucleotide)
  cat_cag_downstream<-peaks %>% filter(central_dinucleotide=="CAT-CAG") %>% select(central_dinucleotide, motif_distance_summit, downstream_dinucleotide)
  #hay que hacer el revcomp en los downstream
  cag_cag_downstream$downstream_dinucleotide<-stri_reverse(chartr("ATGC","TACG", cag_cag_downstream$downstream_dinucleotide) )
  cat_cag_downstream$downstream_dinucleotide<-stri_reverse(chartr("ATGC","TACG", cat_cag_downstream$downstream_dinucleotide) )
  #juntarlos
  colnames(cag_cag_upstream)[3]<-"dinucleotide"
  colnames(cag_cac_upstream)[3]<-"dinucleotide"
  colnames(cag_cag_downstream)[3]<-"dinucleotide"
  colnames(cat_cag_downstream)[3]<-"dinucleotide"
  cag_flanking<-rbind(cag_cag_upstream,cag_cac_upstream,cag_cag_downstream,cat_cag_downstream)
  
  #sacar los porcentajes de los flanking dinucleotides
  pcts<-cag_flanking %>% filter(abs(motif_distance_summit)<distance_around_summit) %>% 
    group_by(central_dinucleotide,dinucleotide) %>% summarise(n=n()) %>% 
    mutate(pct=n/sum(n))
  pcts_background<-cag_flanking %>% filter(abs(motif_distance_summit)>back_dist) %>% 
    group_by(central_dinucleotide,dinucleotide) %>% summarise(n=n()) %>% 
    mutate(pct=n/sum(n))
  pcts$central_dinucleotide<-paste0("2   CAG : ",pcts$central_dinucleotide)
  pcts_background$central_dinucleotide<-paste0("2   CAG : ", pcts_background$central_dinucleotide)
  
  woonk<-pcts    
  woonk_background<-pcts_background
  
  #cac half-site
  cat_cac_downstream<-peaks %>% filter(central_dinucleotide=="CAT-CAC") %>% select(central_dinucleotide, motif_distance_summit, downstream_dinucleotide)
  cag_cac_downstream<-peaks %>% filter(central_dinucleotide=="CAG-CAC") %>% select(central_dinucleotide, motif_distance_summit, downstream_dinucleotide)
  #hay que hacer el revcomp en los downstream
  cat_cac_downstream$downstream_dinucleotide<-stri_reverse(chartr("ATGC","TACG", cat_cac_downstream$downstream_dinucleotide) )
  cag_cac_downstream$downstream_dinucleotide<-stri_reverse(chartr("ATGC","TACG", cag_cac_downstream$downstream_dinucleotide) )
  colnames(cat_cac_downstream)[3]<-"dinucleotide"
  colnames(cag_cac_downstream)[3]<-"dinucleotide"
  cac_flanking<-rbind(cat_cac_downstream,cag_cac_downstream)
  
  #sacar los porcentajes de los flanking dinucleotides
  pcts<-cac_flanking %>% filter(abs(motif_distance_summit)<distance_around_summit) %>% 
    group_by(central_dinucleotide,dinucleotide) %>% summarise(n=n()) %>% 
    mutate(pct=n/sum(n))
  pcts_background<-cac_flanking %>% filter(abs(motif_distance_summit)>back_dist) %>% 
    group_by(central_dinucleotide,dinucleotide) %>% summarise(n=n()) %>% 
    mutate(pct=n/sum(n))
  pcts$central_dinucleotide<-paste0("3   CAC : ",pcts$central_dinucleotide)
  pcts_background$central_dinucleotide<-paste0("3   CAC : ",pcts_background$central_dinucleotide)
  
  wooonk<-pcts
  wooonk_background<-pcts_background
  
  
  pcts<-rbind(wonk,woonk,wooonk)
  background<-rbind(wonk_background,woonk_background,wooonk_background)
  pcts$pct[pcts$pct>0.2]<-0.2
  extra<-data.frame("pct"=rep(0,128))
  extra$n<-rep(0,128)
  extra$dinucleotide<-rep( c("AA", "AC", "AG", "AT" ,"CA" ,"CC", "CG", "CT" ,"GA", "GC", "GG", "GT", "TA", "TC" ,"TG" ,"TT"), 8)
  extra$central_dinucleotide<-rep(c("1   CAT : CAT-CAC", "1   CAT : CAT-CAG", "1   CAT : CAT-CAT", "2   CAG : CAG-CAC" ,"2   CAG : CAG-CAG", "2   CAG : CAT-CAG", "3   CAC : CAG-CAC",
                                    "3   CAC : CAT-CAC" ) , each=16 )
  pcts<-rbind(extra,pcts)
  pcts<-pcts %>% group_by(central_dinucleotide,dinucleotide) %>% summarise(pct=sum(pct) )
  background<-rbind(extra,background)
  background<-background %>% group_by(central_dinucleotide,dinucleotide) %>% summarise(pct=sum(pct))
  
  pcts$fold_change<-pcts$pct/background$pct
  
  pcts$fold_change[is.nan(pcts$fold_change)]<-0
  pcts$fold_change[which(pcts$fold_change>lim)]<-lim
  pcts$fold_change[is.na(pcts$fold_change)]<-0
  pcts$fold_change[is.nan(pcts$fold_change)]<-0
  
  ggballoonplot(pcts, x = "central_dinucleotide", y = "dinucleotide",
                size = "pct", fill = "fold_change") + 
    gradient_fill(c("black","darkred","yellow2"))
  
}

flanking_motifs_balloon_2<-function(peaks,distance_around_summit,lim,back_dist) {
  peaks$motif_center<-as.numeric(peaks$ebox_start)+5
  peaks$motif_distance_summit<-peaks$motif_center-as.numeric(peaks$summit)
  
  #cat half-site
  cat_cat_upstream<-peaks %>% filter(central_dinucleotide=="CAT-CAT") %>% select(central_dinucleotide, motif_distance_summit, upstream_dinucleotide)
  cat_cat_downstream<-peaks %>% filter(central_dinucleotide=="CAT-CAT") %>% select(central_dinucleotide, motif_distance_summit, downstream_dinucleotide)
  #hay que hacer el revcomp en los downstream
  cat_cat_downstream$downstream_dinucleotide<-stri_reverse(chartr("ATGC","TACG", cat_cat_downstream$downstream_dinucleotide) )
  cat_cag_upstream<-peaks %>% filter(central_dinucleotide=="CAT-CAG") %>% select(central_dinucleotide, motif_distance_summit, upstream_dinucleotide)
  cat_cac_upstream<-peaks %>% filter(central_dinucleotide=="CAT-CAC") %>% select(central_dinucleotide, motif_distance_summit, upstream_dinucleotide)
  #juntarlos
  colnames(cat_cat_upstream)[3]<-"dinucleotide"
  colnames(cat_cat_downstream)[3]<-"dinucleotide"
  colnames(cat_cag_upstream)[3]<-"dinucleotide"
  colnames(cat_cac_upstream)[3]<-"dinucleotide"
  cat_flanking<-rbind(cat_cat_upstream,cat_cat_downstream,cat_cag_upstream,cat_cac_upstream)
  
  #convertir a todos en CAT
  cat_flanking$central_dinucleotide<-"CAT"
  
  #sacar los porcentajes de los flanking dinucleotides
  pcts<-cat_flanking %>% filter(abs(motif_distance_summit)<distance_around_summit) %>% 
    group_by(central_dinucleotide,dinucleotide) %>% summarise(n=n()) %>% 
    mutate(pct=n/sum(n))
  pcts_background<-cat_flanking %>% filter(abs(motif_distance_summit)>back_dist) %>% 
    group_by(central_dinucleotide,dinucleotide) %>% summarise(n=n()) %>% 
    mutate(pct=n/sum(n))
  
  wonk<-pcts
  wonk_background<-pcts_background
  
  #cag half-site
  cag_cag_upstream<-peaks %>% filter(central_dinucleotide=="CAG-CAG") %>% select(central_dinucleotide, motif_distance_summit, upstream_dinucleotide)
  cag_cac_upstream<-peaks %>% filter(central_dinucleotide=="CAG-CAC") %>% select(central_dinucleotide, motif_distance_summit, upstream_dinucleotide)
  cag_cag_downstream<-peaks %>% filter(central_dinucleotide=="CAG-CAG") %>% select(central_dinucleotide, motif_distance_summit, downstream_dinucleotide)
  cat_cag_downstream<-peaks %>% filter(central_dinucleotide=="CAT-CAG") %>% select(central_dinucleotide, motif_distance_summit, downstream_dinucleotide)
  #hay que hacer el revcomp en los downstream
  cag_cag_downstream$downstream_dinucleotide<-stri_reverse(chartr("ATGC","TACG", cag_cag_downstream$downstream_dinucleotide) )
  cat_cag_downstream$downstream_dinucleotide<-stri_reverse(chartr("ATGC","TACG", cat_cag_downstream$downstream_dinucleotide) )
  #juntarlos
  colnames(cag_cag_upstream)[3]<-"dinucleotide"
  colnames(cag_cac_upstream)[3]<-"dinucleotide"
  colnames(cag_cag_downstream)[3]<-"dinucleotide"
  colnames(cat_cag_downstream)[3]<-"dinucleotide"
  cag_flanking<-rbind(cag_cag_upstream,cag_cac_upstream,cag_cag_downstream,cat_cag_downstream)
  
  #convertir a todos en CAG
  cag_flanking$central_dinucleotide<-"CAG"
  
  #sacar los porcentajes de los flanking dinucleotides
  pcts<-cag_flanking %>% filter(abs(motif_distance_summit)<distance_around_summit) %>% 
    group_by(central_dinucleotide,dinucleotide) %>% summarise(n=n()) %>% 
    mutate(pct=n/sum(n))
  pcts_background<-cag_flanking %>% filter(abs(motif_distance_summit)>back_dist) %>% 
    group_by(central_dinucleotide,dinucleotide) %>% summarise(n=n()) %>% 
    mutate(pct=n/sum(n))
  
  woonk<-pcts
  woonk_background<-pcts_background
  
  #cac half-site
  cat_cac_downstream<-peaks %>% filter(central_dinucleotide=="CAT-CAC") %>% select(central_dinucleotide, motif_distance_summit, downstream_dinucleotide)
  cag_cac_downstream<-peaks %>% filter(central_dinucleotide=="CAG-CAC") %>% select(central_dinucleotide, motif_distance_summit, downstream_dinucleotide)
  #hay que hacer el revcomp en los downstream
  cat_cac_downstream$downstream_dinucleotide<-stri_reverse(chartr("ATGC","TACG", cat_cac_downstream$downstream_dinucleotide) )
  cag_cac_downstream$downstream_dinucleotide<-stri_reverse(chartr("ATGC","TACG", cag_cac_downstream$downstream_dinucleotide) )
  colnames(cat_cac_downstream)[3]<-"dinucleotide"
  colnames(cag_cac_downstream)[3]<-"dinucleotide"
  cac_flanking<-rbind(cat_cac_downstream,cag_cac_downstream)
  
  #convertir a todos en CAC
  cac_flanking$central_dinucleotide<-"CAC"
  
  #sacar los porcentajes de los flanking dinucleotides
  pcts<-cac_flanking %>% filter(abs(motif_distance_summit)<distance_around_summit) %>% 
    group_by(central_dinucleotide,dinucleotide) %>% summarise(n=n()) %>% 
    mutate(pct=n/sum(n))
  pcts_background<-cac_flanking %>% filter(abs(motif_distance_summit)>back_dist) %>% 
    group_by(central_dinucleotide,dinucleotide) %>% summarise(n=n()) %>% 
    mutate(pct=n/sum(n))
  
  wooonk<-pcts
  wooonk_background<-pcts_background
  
  pcts<-rbind(wonk,woonk,wooonk)
  background<-rbind(wonk_background,woonk_background,wooonk_background)
  
  
  pcts$fold_change<-pcts$pct/background$pct
  
  pcts$fold_change[is.nan(pcts$fold_change)]<-0
  pcts$fold_change[which(pcts$fold_change>lim)]<-lim
  pcts$fold_change[is.na(pcts$fold_change)]<-0
  pcts$fold_change[is.nan(pcts$fold_change)]<-0
  
  scale<-data.frame("central_dinucleotide"=c("CAT","CAG","CAC"),"dinucleotide"="ZZscale", "n"=0,
                    "pct"=rep(0.2,3), "fold_change"=rep(2,3))
  pcts<-rbind(pcts,scale)
  pcts$central_dinucleotide<-as.factor(pcts$central_dinucleotide)
  pcts$central_dinucleotide<-factor(pcts$central_dinucleotide, levels=c("CAT","CAG","CAC"))
  ggballoonplot(pcts, x = "central_dinucleotide", y = "dinucleotide",
                size = "pct", fill = "fold_change") + 
    gradient_fill(c("black","darkred","yellow2"))
}


#anotacion genica en pie chart

#anotar picos
annotate_peaks<-function(regions,start,end,chr) {
  #la anotacion esta filtrada para poner solo protein moding, lncRNA y miRNA
  mm10.gencode.txdb<-makeTxDbFromGFF("~/genomes/gencode.vM24.annotation.filtered.gtf", format="gtf")
  
  granges<-makeGRangesFromDataFrame(regions,start.field = start, end.field = end, seqnames.field = chr)
  #ignoro downstream y first exon pq son muy pequeños
  options(ChIPseeker.ignore_1st_exon = TRUE)
  options(ChIPseeker.ignore_downstream = TRUE)
  anno <- annotatePeak(granges, TxDb=mm10.gencode.txdb)
  annotation<-anno@anno$annotation
  annotation[grep("Exon",anno@anno$annotation)]<-"intragenic"
  annotation[grep("Intron",anno@anno$annotation)]<-"intragenic"
  annotation[grep("UTR",anno@anno$annotation)]<-"intragenic"
  annotation[grep("Intergenic",anno@anno$annotation)]<-"intergenic"
  annotation[grep("Downstream",anno@anno$annotation)]<-"intergenic"
  annotation[which(abs(anno@anno$distanceToTSS)<=1000)]<-"promoter"
  regions$annotation<-annotation
  regions$gene<-anno@anno$geneId
  return(regions)
}

#numero de mutaciones de cada motivo comparado con la rata
conservation_motifs<-function(motifs,fasta_mafs_dir) {
  
  super<-c()
  
  for (i in 1:length(motifs[,1])){
    motif_id<-motifs[i,4]
    print(motif_id)
    info<-file.info(paste0(fasta_mafs_dir,"/",motif_id,".maf.fa"))
    # si no hay file es que no hay alineamineto y entonces poner NAs en todas las posiciones de la conservacion
    if (info$size==0){
      unaligned<-data.frame("conservation_rat"= NA ,"motif_id"=motif_id)
      super<-rbind(super,unaligned)
      next
    }
    
    #leer el fasta del alineamiento
    ali<-read.fasta(paste0(fasta_mafs_dir,"/",motif_id,".maf.fa"))
    names(ali)<-sapply(strsplit(as.character(names(ali)), split=".", fixed=T),"[", 1)
    names(ali)->n
    
    length(ali[[1]])->l
    #o sea la referencia, mm10, es el primer elemento 100pre
    ref<-1
    
    mutations<-rep(list(0), length(n))
    names(mutations)<-n
    indel<-mutations
    mutations.imp<-mutations
    indel.imp<-mutations
    
    #cojo las core positions solo pq las del centro igual pueden variar mas facilmente.
    # puedes pasar de un GA a un TA y tener una funcion parecida
    core<-c(3,4,7,8)
    
    r<-0
    for (pos in c(1:l)){
      toupper(unlist(getFrag(ali,pos,pos)))->frag
      ## Filled gaps in mm10 counts as mutations for given species
      if(frag[ref]=="-"){
        if(r %in% core){
          ind <- n[which(frag != '-')]
          indel.imp[ind]<-as.numeric(indel.imp[ind])+1
        }
        n[which(frag != '-')] ->index
        indel[index]<-as.numeric(indel[index])+1
        next
      }
      
      r<-r+1
      
      if(r %in% core){
        mut <- n[which(frag != frag[ref] & frag != '-')]
        mutations.imp[mut]<-as.numeric(mutations.imp[mut])+1
        ind <- n[which(frag != frag[ref] & frag == '-')]
        indel.imp[ind]<-as.numeric(indel.imp[ind])+1
      }
      mut <- n[which(frag != frag[ref] & frag != '-')]
      mutations[mut]<-as.numeric(mutations[mut])+1
      ind <- n[which(frag != frag[ref] & frag == '-')]
      indel[ind]<-as.numeric(indel[ind])+1
    }
    
    #solo mutaciones single. tener en cuenta el numero de mutaciones. el total
    mutations_core_rat<-mutations.imp[[2]]
    
    #si cae un indel dentro poner NA para no tener en cuenta. o sea si cae en el central dinucleotide tampoco tener en cuenta.
    #pq se esta petando la ebox
    if (indel[[2]] > 0) { mutations_core_rat<-NA }
    
    #ahora juntarlo todo
    pre.super<-data.frame("conservation_rat"=mutations_core_rat ,"motif_id"=motif_id)
    
    #juntar esta presuper con la anterior super
    super<-rbind(super,pre.super)
  }
  return(super)
}


#funcion que te clasifica los picos segun su correlacion con la expresion del gen que quieras
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

correlation_expression_peaks_neurog2<-function(gene, return_neurog2peaks)  {
  #sortear las celulas por la expresion del gen y hacer pseudobulks
  atac_counts<-atac_counts[,barcodes]
  expression_neurod2<-sct_expression[gene,]
  expression_neurod2<-sort(expression_neurod2)
  atac_counts<-atac_counts[,names(expression_neurod2)]
  gc()
  i<-1
  #loop que te hace los seudobulks
  while (i+100<length(colnames(atac_counts))) {
    print(i)
    if (i==1) { 
      pseudobulks<-rowSums(as.matrix(atac_counts[,i:(i+100)]))
      neurod2_expression_pseudos<-mean(expression_neurod2[i:(i+100)])
    }
    else {
      pseudobulk<-rowSums(as.matrix(atac_counts[,i:(i+100)]))
      pseudobulks<-cbind(pseudobulks,pseudobulk)
      neurod2_expression_pseudo<-mean(expression_neurod2[i:(i+100)])
      neurod2_expression_pseudos<-c(neurod2_expression_pseudos,neurod2_expression_pseudo)
    }
    i<-i+100 
  }
  gc()
  rownames(pseudobulks)<-rownames(atac_counts)
  #cpms
  pseudobulks<-t(t(pseudobulks)/colSums(pseudobulks)*1000000)
  require(HybridMTest)
  sperma<-row.spearman(pseudobulks, neurod2_expression_pseudos)
  gc()
  sperma$peak_set<-rep("All_peaks", length(sperma[,1]))
  #intersecar los picos con los chipseq de neurod2
  require(Signac)
  atac_peaks<-StringToGRanges(rownames(pseudobulks), sep = c("-", "-"))
  gc()
  atac_peaks<-as.data.frame(atac_peaks)
  atac_peaks$seqnames<-as.character(atac_peaks$seqnames)
  neurod2_peaks<-read.table("~/neurog2/chipseq/e14.5.neurog2_peaks.subpeaks.bed")
  neurod2_peaks$V2<-as.integer(neurod2_peaks$V2)
  neurod2_peaks$V3<-as.integer(neurod2_peaks$V3)
  require(bedr)
  intersect<-bedr.join.region(bedr.sort.region(atac_peaks),bedr.sort.region(neurod2_peaks))
  detach("package:Signac", unload=TRUE)
  library(dplyr)
  intersect<- intersect %>% filter(V2!=-1)
  intersect$atac_id<-paste0(intersect$seqnames,".",intersect$start,".",intersect$end)
  intersect_sperma<-sperma[intersect$atac_id,]
  intersect_sperma$peak_set<-rep("Neurog2_intersect", length(intersect_sperma[,1]))
  sperma_join<-rbind(sperma,intersect_sperma)
  gc()
  if (return_neurog2peaks==T){
    intersect_sperma$atac_id<-rownames(intersect_sperma)
    merged<-merge(intersect,intersect_sperma, by="atac_id")
    write.table(merged,"~/openchr/noack/correlation_peaks_neurog2.txt", sep = "\t", quote = F)
    return(merged)
  }
  require(ggplot2)
  return(sperma_join)
}

#funcion que saca barcodes ordenando por la expresion de un gen
correlation_expression_motifs<-function(gene_expression, factor)  {
  #sortear las celulas por la expresion del gen y hacer pseudobulks
  atac_counts<-atac_counts[,barcodes]
  expression_neurod2<-sct_expression[gene_expression,]
  expression_neurod2<-sort(expression_neurod2)
  #hacer pseudobulks de 50 celulas
  i<-1
  j<-1
  while ( i<length(expression_neurod2)) {
    if(i==1) {
      pseudo<-rep(paste0(factor,"_",j),50)
    }
    else if ((i+50)>length(expression_neurod2)) {
      pseudo<-c(pseudo,rep(paste0(factor,"_",j), (i+50) - length(expression_neurod2)+3 ))
    }
    else {
      pseudo<-c(pseudo,rep(paste0(factor,"_",j), 50 ))
    }
    i<-i+50
    j<-j+1
  }
  sorted_barcodes<-paste0(names(expression_neurod2),",",pseudo)
  rep1<-sorted_barcodes[grep("-5",sorted_barcodes)]
  rep2<-sorted_barcodes[grep("-6",sorted_barcodes)]
  rep1<-rep1 %>% stringr::str_replace("-5", "-1")
  rep2<-rep2 %>% stringr::str_replace("-6", "-1")
  write.table(rep1,paste0("~/openchr/noack/rep1_",factor,".csv"), col.names = F, row.names = F, quote = F )
  write.table(rep2,paste0("~/openchr/noack/rep2_",factor,".csv"), col.names = F, row.names = F, quote = F )
  gc()
}

factor_ordered_expression<-function(factor)  {
  #sortear las celulas por la expresion del gen y hacer pseudobulks
  atac_counts<-atac_counts[,barcodes]
  expression_neurod2<-sct_expression[factor,]
  expression_neurod2<-sort(expression_neurod2)
  #hacer pseudobulks de 50 celulas
  i<-1
  while ( i<length(expression_neurod2)) {
    if(i==1) {
      media<-mean(expression_neurod2[i:50])
      pseudo<-media
    }
    else if ((i+50)>length(expression_neurod2)) {
      media<-mean(expression_neurod2[i:length(expression_neurod2)])
      pseudo<-c(pseudo,media)
    }
    else {
      media<-mean(expression_neurod2[i:50])
      pseudo<-c(pseudo,media)
    }
    i<-i+50
  }
  return(pseudo)
}


motifs_correlation<-function(correlation_table, motifs, ndecils, limit, fontsize, title) {
  require(dplyr)
  correlation_table<-correlation_table %>% arrange(stat)
  npeaks<-length(correlation_table[,1])
  window<-npeaks%/%ndecils
  i<-1
  j<-1
  while ((i+window)<npeaks) {
    peaks_decil<-correlation_table$V5[i:(i+window)]
    min_rho<-round(correlation_table$stat[i],2)
    max_rho<-round(correlation_table$stat[i+window],2)
    motifs_decil<-motifs %>% filter(peak_id %in% peaks_decil)
    plotname<-paste0("plot",j)
    plot<-distribution_motifs(peaks=motifs_decil,lim=limit,colors=colors,binwidth=8,
                              main = paste0("rho: ", min_rho,"/",max_rho))
    plot<-plot + theme(legend.position = "none") + theme(text = element_text(size = fontsize))
    assign(plotname,plot)
    i<-i+window
    j<-j+1
  }
  require(ggpubr)
  return(ggarrange(plot1,plot2,plot3,plot4,plot5 + rremove("x.text"), 
                   ncol = 5, nrow = 1 ) )
}

