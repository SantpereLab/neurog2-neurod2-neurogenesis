system("
    cd ~/p/openchr/noack
    #resizear los summits todos igual, para que no pasen cosas raras. y a parte centrimo lo hace asi tambien
    #de momento hago grande la ventana. luego ya mirare a ver como filtro eso
    #tambien poner un ID
    cat Proneural_correlation_sorted_summits.bed | awk '{print $1,$2-750,$3+750}' | awk '{ print $0, NR }' | tr ' ' '\t' > Proneural_correlation_summits_extended.bed
    cat Neurog2_correlation_sorted_summits.bed | awk '{print $1,$2-750,$3+750}' | awk '{ print $0, NR }' | tr ' ' '\t' > Neurog2_correlation_summits_extended.bed
    cat Neurod2_correlation_sorted_summits.bed | awk '{print $1,$2-750,$3+750}' | awk '{ print $0, NR }' | tr ' ' '\t' > Neurod2_correlation_summits_extended.bed
    
    #splittear las bed files en quantiles
    srun --partition=normal --mpi=none --nodes=1 --pty bash -i
    module load R
    R
    setwd('~/p/openchr/noack')

    head(sorted_neurod2)
    splitt<-function(sortedd, ndeciles) {
      sorted<-read.table(sortedd)
      npeaks<-length(sorted[,1])
      veintil<-round(npeaks/ndeciles)
      #loop que te saca los indices
      j<-1
      i<-0
      while (i+veintil<=npeaks) {
        print(i)
        print(j)
        sorted_decil<-sorted[i:(i+veintil),]
        sorted_decil$decil<-j
        write.table(sorted_decil, file=paste0('quantil_',j,'_',sortedd), sep='\t', quote=F, col.names=F, row.names=F )
        j<-j+1
        i<-i+veintil
      }
    }

    # cambiar para que salgan 50
    
    splitt(sortedd='Neurod2_correlation_summits_extended.bed',ndeciles=10)
    splitt(sortedd='Neurog2_correlation_summits_extended.bed',ndeciles=10)
    splitt(sortedd='Proneural_correlation_summits_extended.bed',ndeciles=10)
    
    module load BEDTools
    cd ~/p/openchr/noack
    for file in quantil*Neurod2*bed
    do
      sbatch --mem=300000 --partition=bigmem --wrap='sort -k1,1 -k2,2n $file | bedtools intersect -wao -a - -b ~/p/genomes/mm10.archetype_motifs.v1.0.bed | cut -f1,2,3,4,7,8,9 > viestrao.$file '
    done
    
    for file in quantil*Neurog2*bed
    do
      sbatch --mem=300000 --partition=bigmem --wrap='sort -k1,1 -k2,2n $file | bedtools intersect -wao -a - -b ~/p/genomes/mm10.archetype_motifs.v1.0.bed | cut -f1,2,3,4,7,8,9 > viestrao.$file '
    done
    
    for file in quantil*Proneural*bed
    do
      sbatch --mem=300000 --partition=bigmem --wrap='sort -k1,1 -k2,2n $file | bedtools intersect -wao -a - -b ~/p/genomes/mm10.archetype_motifs.v1.0.bed | cut -f1,2,3,4,7,8,9 > viestrao.$file '
    done
")







#dividir los viestra files segun la interseccion con los distintos proneural factors
system("
  cd ~/openchr/noack
  export_shiva atac_private_neurod2.bed p/openchr/noack
  export_shiva atac_private_neurog2.bed p/openchr/noack
  export_shiva atac_shared.bed p/openchr/noack
  export_shiva proneuralnt.bed p/openchr/noack
")

system("
    cd ~/p/openchr/noack
    module load BEDTools
    for file in viestr*bed
      do
      for file2 in atac_private_neurod2 atac_private_neurog2 atac_shared proneuralnt
      do
        bedtools intersect -wa -a $file -b $file2.bed > $file2.$file
      done
    done
")

#primero hacer el intersect entre los picos de los quantiles y los de los factores para saber el npeaks
system ("
  cd ~/p/openchr/noack
  module load BEDTools
  for file in quantil*bed
    do
    for file2 in atac_private_neurod2 atac_private_neurog2 atac_shared proneuralnt
    do
      bedtools intersect -wa -a $file -b $file2.bed > $file2.$file
    done
  done
  
  #hacer tambien con todos los picos de neurod2 y de neurog2 en general
  cat atac_shared.bed atac_private_neurod2.bed > all_atac_neurod2.bed
  cat atac_shared.bed atac_private_neurog2.bed > all_atac_neurog2.bed
  
  for file in quantil*bed vi*quantil*bed
    do
    for file2 in all_atac_neurod2 all_atac_neurog2 proneuralnt
    do
      bedtools intersect -wa -a $file -b $file2.bed > $file2.$file
    done
  done
")

#subdividir en promotores y enhancers
system("
  cd ~/p/openchr/noack
  module load BEDTools
  for file in quantil*bed vi*quantil*bed
  do
    cd ~/p/openchr/noack
    bedtools intersect -wa -b footprinting/promoters.bed -a $file > promoters_${file}
    bedtools intersect -wa -b footprinting/enhancers.bed -a $file > enhancers_${file}
  done
")



# Funcion que te saca la centralidad de los motivos
centrall<-function(factorr,subdivision,nquantiles,boundary_1,boundary_2) {
  setwd("~/p/openchr/noack")
  require(dplyr)
  #hacer por cada quantil
  for (i in 1:nquantiles) {
    print(i)

    # si hay picos
    if ( file.info(paste0(subdivision,"viestrao.quantil_",i,"_",factorr,"_correlation_summits_extended.bed"))$size> 0) {

      sortt<-read.table(paste0(subdivision,"viestrao.quantil_",i,"_",factorr,"_correlation_summits_extended.bed"))
      colnames(sortt)<-c("chr","start","end","orden","motif_start","motif_end","archetype")
      sortt$summit<-as.numeric(sortt$start+750)
      sortt$motif_midpoint<-as.integer((sortt$motif_end-sortt$motif_start)/2)+sortt$motif_start
      sortt$distance_summit<-abs(sortt$summit-sortt$motif_midpoint)
      sortt_close<-sortt %>% dplyr::filter(distance_summit<boundary_1)
      sortt_far<-sortt %>% dplyr::filter(distance_summit>boundary_2)
      
      sortt_close_summ<-sortt_close %>% group_by(archetype) %>% summarise(nmotifs=n())
      sortt_far_summ<-sortt_far %>% group_by(archetype) %>% summarise(nmotifs=n())
      sortt_close_summ$cuantil<-i
      sortt_far_summ$cuantil<-i
      npeaks<-length(read.table(paste0(subdivision,"quantil_",i,"_",factorr,"_correlation_summits_extended.bed"))[,1])
      sortt_close_summ$npeaks<-npeaks
      sortt_far_summ$npeaks<-length(read.table(paste0(subdivision,"quantil_",i,"_",factorr,"_correlation_summits_extended.bed"))[,1])
      if (i==1) {
        closee<-sortt_close_summ
        faar<-sortt_far_summ
        all_motifs<-unique(c(closee$archetype,faar$archetype) )
        npeakk<-npeaks
      }
      else {
        faar<-rbind(faar,sortt_far_summ)
        closee<-rbind(closee,sortt_close_summ)
        all_motifs<-unique(unique(c(faar$archetype,closee$archetype) ) ,all_motifs )
        npeakk<-c(npeakk,npeaks)
      }

    }
    
    # si no hay picos del factor en el cuantil, pongo que hay 0 picos. luego los NaN ya los transformare en 0.
    else {
      if (i==1) {
        npeakk<-0
        all_motifs<-c()
        closee<-data.frame("nmotifs"=c(),"cuantil"=c(),"npeaks"=c())
        faar<-data.frame("nmotifs"=c(),"cuantil"=c(),"npeaks"=c())
      }
      else {
        npeakk<-c(npeakk,npeaks)
      }
    }
  }
  extraa<-data.frame("archetype"=rep(all_motifs,nquantiles))
  extraa$npeaks<-rep(npeakk,each=length(all_motifs) )
  extraa$nmotifs<-0
  extraa$cuantil<-rep(1:nquantiles, each=length(all_motifs))
  closee<-rbind(closee,extraa) %>% group_by(archetype,cuantil,npeaks) %>% summarise(nmotifs=sum(nmotifs) )
  faar<-rbind(faar,extraa) %>% group_by(archetype,cuantil,npeaks) %>% summarise(nmotifs=sum(nmotifs) )

  closee$nmotifs_corrected<-closee$nmotifs/closee$npeaks/(boundary_1*2)
  faar$nmotifs_corrected<-faar$nmotifs/faar$npeaks/((750-boundary_2)*2)
  closee$centrality<-closee$nmotifs_corrected/faar$nmotifs_corrected
  
  saveRDS(closee,paste0(subdivision,factorr,"_archetypes.RDS"))
}

#hacerlo con 10 quantiles. y los factores sin subdividirlos tanto
centrall(factorr="Neurod2",subdivision="all_atac_neurod2.",nquantiles=10,boundary_1=25,boundary_2=250)
centrall(factorr="Neurod2",subdivision="",nquantiles=10,boundary_1=25,boundary_2=250)
centrall(factorr="Neurod2",subdivision="proneuralnt.",nquantiles=10,boundary_1=25,boundary_2=250)
centrall(factorr="Neurod2",subdivision="all_atac_neurog2.",nquantiles=10,boundary_1=25,boundary_2=250)

factorr="Neurod2";subdivision="all_atac_neurog2.";nquantiles=10;boundary_1=25;boundary_2=250

#hacerlo con 10 quantiles. y los factores sin subdividirlos tanto
centrall(factorr="Neurog2",subdivision="all_atac_neurog2.",nquantiles=10,boundary_1=25,boundary_2=250)
centrall(factorr="Neurog2",subdivision="",nquantiles=10,boundary_1=25,boundary_2=250)
centrall(factorr="Neurog2",subdivision="proneuralnt.",nquantiles=10,boundary_1=25,boundary_2=250)
centrall(factorr="Neurog2",subdivision="all_atac_neurod2.",nquantiles=10,boundary_1=25,boundary_2=250)

centrall(factorr="Proneural",subdivision="atac_shared.",nquantiles=10,boundary_1=25,boundary_2=250)





## FUNCIONES
centrality_pheat<-function(summ,factor, anchura, altura, top_motifs, return_motifs_sorted) {
  #funcion de pheatmap para el factor que qiueras
  #hacer el pheatmap primero para clusterizar los motivos
  require(dplyr)
  summ$centrality[ which(is.nan(summ$centrality)) ]<-0
  summ$nmotifs_corrected[ which(is.nan(summ$nmotifs_corrected)) ]<-0
  
  #sortear los motivos segun su fold cahnge en la centralidad
  #coger los primeros y ultimos 10, para que no haya tanto random
  start_centrality<-summ %>% filter( cuantil==1 ) %>% group_by(archetype) %>% summarise(start_centrality=mean(centrality))
  end_centrality<-summ %>% filter( cuantil==10 ) %>% group_by(archetype) %>% summarise(end_centrality=mean(centrality))
  fc_centrality<-end_centrality$end_centrality/start_centrality$start_centrality
  names(fc_centrality)<-end_centrality$archetype
  fc_centrality<-sort(fc_centrality)
  motifs_sorted<-rev(names(fc_centrality) )
  motifs_sorted<-motifs_sorted[1:top_motifs]
  if (return_motifs_sorted==T) {
    return(motifs_sorted)
  }
  
  summ<- summ %>% filter(archetype %in% motifs_sorted)
  
  summ$cuantil<-as.factor(summ$cuantil)
  summ$cuantil<-factor(summ$cuantil,levels=rev(c(1:10)))
  summ$archetype<-as.factor(summ$archetype)
  summ$archetype<-factor(summ$archetype,levels=motifs_sorted)
  
  pdf( paste0( "~/", factor, "_cor_atacseq_peaks_centrality_archetypes_fc_centrality_sorted.pdf" ), useDingbats = F, width = anchura, height = altura )
    p<-ggballoonplot(summ, x = "archetype", y = "cuantil",
                  size = "nmotifs_corrected", fill = "centrality") + gradient_fill( c("black","darkred","yellow2") )
    print(p)
  dev.off()
}

#para poner todo en la misma escala, añadir una columna extra y capar
centrality_pheat_scale<-function(summ,factor, anchura, altura, top_motifs, sorted_motis) {
  #funcion de pheatmap para el factor que qiueras
  #hacer el pheatmap primero para clusterizar los motivos
  require(dplyr)
  summ$centrality[ which(is.nan(summ$centrality)) ]<-0
  summ$nmotifs_corrected[ which(is.nan(summ$nmotifs_corrected)) ]<-0
  
  sorted_motis<-sorted_motis[1:top_motifs]

  summ<- summ %>% filter(archetype %in% sorted_motis)
  
  summ$cuantil<-as.factor(summ$cuantil)
  summ$cuantil<-factor(summ$cuantil,levels=rev(c(1:10)))
  summ$archetype<-as.factor(summ$archetype)
  summ$archetype<-factor(summ$archetype,levels=sorted_motis)
  summ<-summ %>% select(archetype,cuantil,nmotifs_corrected,centrality)
  # añadir una columna que sirva para escalar las demas columnas
  scalee<-data.frame("cuantil"= unique(summ$cuantil) )
  scalee$centrality<-4
  scalee$nmotifs_corrected<-0.02
  scalee$archetype<-unique(scalee$archetype)
  summ<-rbind(summ,scalee)
  
  scalee2<-data.frame("cuantil"= unique(summ$cuantil) )
  scalee2$centrality<-0
  scalee2$nmotifs_corrected<-0.0000001
  scalee2$archetype<-unique(scalee2$archetype)
  summ<-rbind(summ,scalee2)
  
  #tambien capar la centralidad por si hay outliers en esos quantiles que no tiene casi picos
  summ$centrality[which(summ$centrality>4)]<-4
  summ$nmotifs_corrected[which(summ$nmotifs_corrected>0.02)]<-0.02
  summ$nmotifs_corrected[which(summ$nmotifs_corrected==0)]<-NA
  
  
  pdf(paste0("~/",factor,"_cor_atacseq_peaks_centrality_archetypes_fc_centrality_sorted.pdf"), useDingbats = F, width = anchura, height = altura)
  p<-ggballoonplot(summ, x = "archetype", y = "cuantil",
                   size = "nmotifs_corrected", fill = "centrality") +
    gradient_fill(c("black","darkred","yellow2"))
  print(p)
  dev.off()
}

#el lineplot que me ha dicho gabriel
centrality_lineplot<-function(summ,factor, anchura, altura, top_motifs, sorted_motis) {
  #funcion de pheatmap para el factor que qiueras
  #hacer el pheatmap primero para clusterizar los motivos
  require(dplyr)
  summ$centrality[ which(is.nan(summ$centrality)) ]<-0
  summ$nmotifs_corrected[ which(is.nan(summ$nmotifs_corrected)) ]<-0
  
  sorted_motis<-sorted_motis[1:top_motifs]
  
  summ<- summ %>% filter(archetype %in% sorted_motis)
  
  summ$cuantil<-as.factor(summ$cuantil)
  summ$cuantil<-factor(summ$cuantil,levels=c(1:10))
  summ$archetype<-as.factor(summ$archetype)
  summ$archetype<-factor(summ$archetype,levels=sorted_motis)
  summ<-summ %>% select(archetype,cuantil,nmotifs_corrected,centrality)
  summ$ebox<-as.character(summ$archetype)
  summ$ebox[-grep("Ebox",summ$ebox)] <-"no"
  
  summ$ebox<-factor(summ$ebox,levels= c("no","Ebox/CACCTG","Ebox/CACGTG/1","Ebox/CACGTG/2","Ebox/CAGATGG","Ebox/CAGCTG","Ebox/CATATG" ) )
  
  pdf(paste0("~/",factor,"_atacseq_peaks_centrality_lineplot.pdf"), useDingbats = F, width = anchura, height = altura)
    p<-ggplot(summ,  aes(x=cuantil, y=centrality, group=archetype, alpha=ebox, color=ebox ) ) +
      scale_alpha_manual(values = c(0.2, rep(1,7) ) ) +
      scale_linewidth_manual(values = c(0.01, rep(1,7) ) ) + ylim (0,5) +
      scale_color_manual(values = c("grey","chartreuse1","darkgoldenrod4","darkgoldenrod4","red","blue","black" ) ) +
      geom_line() + theme_classic()
    
    print(p)
  dev.off()
}

centrality_lineplot_fc<-function(summ,factor, anchura, altura, top_motifs, sorted_motis) {
  #funcion de pheatmap para el factor que qiueras
  #hacer el pheatmap primero para clusterizar los motivos
  require(dplyr)
  summ$centrality[ which(is.nan(summ$centrality)) ]<-0
  summ$nmotifs_corrected[ which(is.nan(summ$nmotifs_corrected)) ]<-0
  
  sorted_motis<-sorted_motis[1:top_motifs]
  
  summ<- summ %>% filter(archetype %in% sorted_motis)
  
  summ$cuantil<-as.factor(summ$cuantil)
  summ$cuantil<-factor(summ$cuantil,levels=c(1:10))
  summ$archetype<-as.factor(summ$archetype)
  summ$archetype<-factor(summ$archetype,levels=sorted_motis)
  summ<-summ %>% select(archetype,cuantil,nmotifs_corrected,centrality)

  summ_1<-summ %>% filter(cuantil==1)
  summ_10<-summ %>% filter(cuantil==10)
  summ_10$fc<-summ_10$centrality/summ_1$centrality
  summ_10<-summ_10 %>% arrange(fc)
  return(summ_10)
}



######## en 10 deciles y sin las tres intersecciones

### todos los picos

setwd("~/openchr/noack")
require(ggpubr)
neurooooo<-readRDS("Neurod2_archetypes.RDS")
neuroggg<-readRDS("Neurog2_archetypes.RDS")



#correlacion con neurod2

motifs_neurod2_correlation_sorted<-centrality_pheat(summ=neurooooo,factor="neurod2", anchura=45, altura=14, top_motifs=233,return_motifs_sorted=T)


neurod2_intersect<-readRDS( "all_atac_neurod2.Neurod2_archetypes.RDS" )

neurog2_intersect<-readRDS( "all_atac_neurog2.Neurod2_archetypes.RDS" )

proneuralnt<-readRDS( "proneuralnt.Neurod2_archetypes.RDS" )

shared_proneural<-readRDS("atac_shared.Proneural_archetypes.RDS")


centrality_pheat_scale ( summ=neurooooo, factor="all_Neurod2cor",
                         anchura=60, altura=4, top_motifs=233, sorted_motis=motifs_neurod2_correlation_sorted )

centrality_pheat_scale ( summ=neurod2_intersect, factor="neurod2_Neurod2cor",
                         anchura=60, altura=4, top_motifs=233, sorted_motis=motifs_neurod2_correlation_sorted )

centrality_pheat_scale ( summ=neurog2_intersect, factor="neurog2_Neurod2cor",
                         anchura=60, altura=4, top_motifs=233, sorted_motis=motifs_neurod2_correlation_sorted )

centrality_pheat_scale ( summ=proneuralnt, factor="proneuralnt_Neurod2cor",
                         anchura=60, altura=4, top_motifs=233, sorted_motis=motifs_neurod2_correlation_sorted )



#con los lineplots de centralidad
centrality_lineplot(summ=neurod2_intersect,factor="neurod2peaks_Neurod2cor", anchura=7, altura=4, top_motifs=233,sorted_motis=motifs_neurog2_correlation_sorted)

centrality_lineplot(summ=proneuralnt,factor="proneuralnt_Neurod2cor", anchura=7, altura=4, top_motifs=233,sorted_motis=motifs_neurog2_correlation_sorted) 

centrality_lineplot(summ=shared_proneural,factor="shared_proneuralcor", anchura=7, altura=4, top_motifs=233,sorted_motis=motifs_neurog2_correlation_sorted)

jj<-centrality_lineplot_fc(summ=neurod2_intersect,factor="neurod2peaks_Neurod2cor", anchura=7, altura=4, top_motifs=233,sorted_motis=motifs_neurog2_correlation_sorted)

tail(jj)



## lo mismo pero con la correlacion con Neurog2


motifs_neurog2_correlation_sorted<-centrality_pheat(summ=neuroggg,factor="neurog2", anchura=45, altura=14, top_motifs=233,return_motifs_sorted=T)


neurog2_intersect<-readRDS( "all_atac_neurog2.Neurog2_archetypes.RDS" )

neurod2_intersect<-readRDS( "all_atac_neurod2.Neurog2_archetypes.RDS" )

proneuralnt<-readRDS( "proneuralnt.Neurog2_archetypes.RDS" )


centrality_pheat_scale ( summ=neuroggg, factor="all_Neurog2cor",
                         anchura=60, altura=4, top_motifs=233, sorted_motis=motifs_neurog2_correlation_sorted )

centrality_pheat_scale ( summ=neurog2_intersect, factor="neurog2_Neurog2cor",
                         anchura=60, altura=4, top_motifs=233, sorted_motis=motifs_neurog2_correlation_sorted )

centrality_pheat_scale ( summ=neurod2_intersect, factor="neurod2_Neurog2cor",
                         anchura=60, altura=4, top_motifs=233, sorted_motis=motifs_neurog2_correlation_sorted )

centrality_pheat_scale ( summ=proneuralnt, factor="proneuralnt_Neurog2cor",
                         anchura=60, altura=4, top_motifs=233, sorted_motis=motifs_neurog2_correlation_sorted )



# con los lineplots de centralidad
centrality_lineplot(summ=neurog2_intersect,factor="neurog2peaks_Neurog2cor", anchura=7, altura=4, top_motifs=233,
                    sorted_motis=motifs_neurog2_correlation_sorted) 

summ=neurog2_intersect
factor="neurog2peaks_Neurog2cor"
anchura=7
altura=4
top_motifs=233
sorted_motis=motifs_neurog2_correlation_sorted


############## CON TODOS LOS TIPOS CELULARES
###################################!

#hacer el peak calling de cada tipo celular
system('
    cd ~/p/openchr/noack/cluster_bams_imputation
    for file in clusterNSC.bam clusterIPC.bam clusterPN1.bam clusterPN2.bam clusterPN3.bam
    do
      macs2 callpeak -t $file -f BAM -g mm --outdir . -n $file
    done
    # resizear al tamaño de los analisis de arriba
    for file in cluster*bed
    do
      cat $file | awk '{print $1,$2-750,$3+750}' | tr ' ' '\t' | grep chr > resized.$file
    done
    
    for file in resized*bed
    do
      sbatch --partition=bigmem --mem=400000 --wrap=" bedtools intersect -wao -a $file -b ~/p/genomes/mm10.archetype_motifs.v1.0.bed | cut -f1,2,3,5,6,7 > vierstrao.$file "
    done
    
    cd ../cluster_bams
    for file in clusterCR.bam clusterIN.bam clusterMG+Mural.bam
    do
      macs2 callpeak -t $file -f BAM -g mm --outdir . -n $file
    done
    
    for file in cluster*bed
    do
      cat $file | awk '{print $1,$2-750,$3+750}' | tr ' ' '\t' | grep chr > resized.$file
    done
    
    for file in resized*bed
    do
      sbatch --partition=bigmem --mem=200000 --wrap=" bedtools intersect -wao -a $file -b ~/p/genomes/mm10.archetype_motifs.v1.0.bed | cut -f1,2,3,5,6,7 > vierstrao.$file "
    done
    
    
')

###### MOTIVOS EN LAS TRAYECTORIAS
############################ !
system("
  cd ~/p/openchr/noack/centrimo
  
  # resizear las beds para que sean del mismo tamaño del centrimo que hacia antes
  for file in *bed
  do
    echo $file
    cat $file | awk '{print $1,$2-500,$3+500}' | tr ' ' '\t' > resized.$file
  done
  
  for file in resized.invariant_noneurod2_summit_centered.bed resized.neuronal_noneurod2_summit_centered.bed resized.transcient_neurod2_summit_centered.bed
  do
    sbatch --partition=bigmem --mem=400000 --wrap="bedtools intersect -wao -a $file -b ~/p/genomes/mm10.archetype_motifs.v1.0.bed | cut -f1,2,3,5,6,7 > vierstrao.$file"
  done
  
")

system("
    srun --partition=bigmem --mpi=none --nodes=1 --pty bash -i
    module load R
    R
    sumarization <-function(oscarsin, set) {
      require(dplyr)
      colnames(oscarsin)<-c("chr","start","end","motif_start","motif_end","archetype")
      oscarsin$peak_id<-paste0(oscarsin$chr,":",oscarsin$motif_start,"-",oscarsin$motif_end)
      npeaks<-length( unique(oscarsin$peak_id) )
      oscarsin$summit<-as.numeric(oscarsin$start+750)
      oscarsin$motif_midpoint<-as.integer((oscarsin$motif_end-oscarsin$motif_start)/2)+oscarsin$motif_start
      oscarsin$distance_summit<-abs(oscarsin$summit-oscarsin$motif_midpoint)
      oscarsin_close<-oscarsin %>% dplyr::filter(distance_summit<25)
      oscarsin_far<-oscarsin %>% dplyr::filter(distance_summit>250)
      
      oscarsin_close_summ<-oscarsin_close %>% group_by(archetype) %>% summarise(nmotifs=n())
      oscarsin_far_summ<-oscarsin_far %>% group_by(archetype) %>% summarise(nmotifs=n())
      
      # añadir motivos a oscarsin close summ si no tiene
      all_archetypes<-read.table('~/p/genomes/all_archetpyes.txt')$V1
      extra<-data.frame('archetype'=all_archetypes, nmotifs=0)
      oscarsin_close_summ<-rbind (oscarsin_close_summ, extra ) %>% group_by(archetype) %>% summarise(nmotifs=sum(nmotifs))
      
      
      oscarsin_close_summ$nmotifs_corrected<-oscarsin_close_summ$nmotifs/npeaks/(25*2)
      oscarsin_far_summ$nmotifs_corrected<-oscarsin_far_summ$nmotifs/npeaks/((750-250)*2)
      oscarsin_close_summ$centrality<-oscarsin_close_summ$nmotifs_corrected/oscarsin_far_summ$nmotifs_corrected
      oscarsin_close_summ$peakset<-set
      oscarsin_close_summ<-oscarsin_close_summ %>% select(archetype, centrality, peakset)
      return(oscarsin_close_summ)
    }  

    setwd('~/p/openchr/noack/centrimo')
    early_neurod2<-read.table('vierstrao.resized.early_neurod2_summit_centered.bed')
    early_noneurod2<-read.table('vierstrao.resized.early_noneurod2_summit_centered.bed')
    invariant_neurod2<-read.table('vierstrao.resized.invariant_neurod2_summit_centered.bed')
    invariant_noneurod2<-read.table('vierstrao.resized.invariant_noneurod2_summit_centered.bed')
    neuronal_neurod2<-read.table('vierstrao.resized.neuronal_neurod2_summit_centered.bed')
    neuronal_noneurod2<-read.table('vierstrao.resized.neuronal_noneurod2_summit_centered.bed')
    tanscient_neurod2<-read.table('vierstrao.resized.transcient_neurod2_summit_centered.bed')
    transcient_noneurod2<-read.table('vierstrao.resized.transcient_noneurod2_summit_centered.bed')
    
    sumarization_neurod2<-rbind( sumarization(early_neurod2, "early") , sumarization(invariant_neurod2, "invariant"), sumarization(neuronal_neurod2, "neuronal"), sumarization(tanscient_neurod2, "transcient")   )
    saveRDS(sumarization_neurod2, "sumarization_neurod2_trajs.RDS")
    sumarization_noneurod2<-rbind( sumarization(early_noneurod2, "early") , sumarization(invariant_noneurod2, "invariant"), sumarization(neuronal_noneurod2, "neuronal"), sumarization(transcient_noneurod2, "transcient")   )
    saveRDS(sumarization_noneurod2, "sumarization_noneurod2_trajs.RDS")
    
    setwd('~/p/openchr/noack/cluster_bams_imputation')
    clusterIPC<-read.table('vierstrao.resized.clusterIPC.bam_summits.bed')
    clusterNSC<-read.table('vierstrao.resized.clusterNSC.bam_summits.bed')
    clusterPN1<-read.table('vierstrao.resized.clusterPN1.bam_summits.bed')
    clusterPN2<-read.table('vierstrao.resized.clusterPN2.bam_summits.bed')
    clusterPN3<-read.table('vierstrao.resized.clusterPN3.bam_summits.bed')
    
    setwd('~/p/openchr/noack/cluster_bams')
    clusterCR<-read.table('vierstrao.resized.clusterCR.bam_summits.bed')
    clusterIN<-read.table('vierstrao.resized.clusterIN.bam_summits.bed')
    clusterMGMural<-read.table('vierstrao.resized.clusterMG+Mural.bam_summits.bed')
    
    sumarization_celltype_summits<-rbind( sumarization(clusterIPC, "IPC") , sumarization(clusterNSC, "NSC"), sumarization(clusterPN1, "PN1"), sumarization(clusterPN2, "PN2"), 
    sumarization(clusterPN3, "PN3"), sumarization(clusterCR, "CR"), sumarization(clusterIN, "IN"), sumarization(clusterMGMural, "MGMural")  )
    saveRDS(sumarization_celltype_summits, "sumarization_celltype_summits.RDS")
")

system("
    tunnel_shiva
    import_shiva p/openchr/noack/centrimo/\*RDS .
    import_shiva p/openchr/noack/cluster_bams/\*RDS .
")

sumarization_neurod2_trajs<-readRDS("sumarization_neurod2_trajs.RDS")
sumarization_noneurod2_trajs<-readRDS("sumarization_noneurod2_trajs.RDS")

####### sacar para cada trayectoria un pheatmap de los motivos mas enriquecidos

mapacalor<-function(summ,traj) {
  pivote<-as.matrix(tidyr::pivot_wider(summ, names_from = archetype, values_from = centrality)[,-1])
  rownames(pivote)<-unique(summ$peakset)
  
  #sacar cuales son los motivos mas enriquecidos en X trayectoria
  centrality_traj<-pivote[traj,]
  #centrality_notraj<-colMeans( pivote[which(rownames(pivote)!=traj) ,] )
  #foldchange_traj<-centrality_traj/centrality_notraj
  top_motifs_traj<-names( rev( sort(centrality_traj) ) ) [1:15]
  pivote_traj<-pivote[ , top_motifs_traj ]
  pivote_traj<-pivote_traj[ c("invariant","early","transcient","neuronal") ,  ]
  return(pivote_traj)
}

pdf("~/vierstramotifs_atacseq_peaks_summits_neurod2intersect_trajectories.pdf", height = 2, width = 15)
  j<-mapacalor(sumarization_neurod2_trajs, "invariant")
  jj<-mapacalor(sumarization_neurod2_trajs, "early")
  jjj<-mapacalor(sumarization_neurod2_trajs, "transcient")
  jjjj<-mapacalor(sumarization_neurod2_trajs, "neuronal")
  aj<-mapacalor(sumarization_noneurod2_trajs, "invariant")
  ajj<-mapacalor(sumarization_noneurod2_trajs, "early")
  ajjj<-mapacalor(sumarization_noneurod2_trajs, "transcient")
  ajjjj<-mapacalor(sumarization_noneurod2_trajs, "neuronal")
  jiji<-cbind(j,jj,jjj,jjjj,aj,ajj,ajjj,ajjjj)
  pheat<-pheatmap::pheatmap(jiji, cluster_cols=F, cluster_rows = F, fontsize = 5, treeheight_col = 0, color = rev(hcl.colors(100, "Inferno") ) )
  print(pheat)
dev.off()

mapacalor_celltypes<-function(summ,traj) {
  pivote<-as.matrix(tidyr::pivot_wider(summ, names_from = archetype, values_from = centrality)[,-1])
  rownames(pivote)<-unique(summ$peakset)
  #sacar cuales son los motivos mas enriquecidos en X trayectoria
  centrality_traj<-pivote[traj,]
  #centrality_notraj<-colMeans( pivote[which(rownames(pivote)!=traj) ,] )
  #foldchange_traj<-centrality_traj/centrality_notraj
  top_motifs_traj<-names( rev( sort(centrality_traj) ) ) [1:15]
  pivote_traj<-pivote[ , top_motifs_traj ]
  pivote_traj<-pivote_traj[ c( "NSC","IPC","PN1","PN2","PN3","CR","IN","MGMural") ,  ]

  return(pivote_traj)
}


sumarization_celltype_summits<-readRDS("sumarization_celltype_summits.RDS")

#heatmaps separados
pdf("~/vierstramotifs_atacseq_peaks_summits_celltypes.pdf", height = 2, width = 15)

  a<-mapacalor_celltypes(sumarization_celltype_summits,"NSC")
  aa<-mapacalor_celltypes(sumarization_celltype_summits,"IPC")
  aaa<-mapacalor_celltypes(sumarization_celltype_summits,"PN1")
  aaaa<-mapacalor_celltypes(sumarization_celltype_summits,"PN2")
  aaaaa<-mapacalor_celltypes(sumarization_celltype_summits,"PN3")
  aaaaaa<-mapacalor_celltypes(sumarization_celltype_summits,"CR")
  aaaaaaa<-mapacalor_celltypes(sumarization_celltype_summits,"IN")
  aaaaaaaa<-mapacalor_celltypes(sumarization_celltype_summits,"MGMural")
  jiji<-cbind(a,aa,aaa,aaaa,aaaaa,aaaaaa,aaaaaaa,aaaaaaaa)
  pheat<-pheatmap::pheatmap(jiji, cluster_cols=F, cluster_rows = F, fontsize = 5, treeheight_col = 0, color = hcl.colors(100, "Inferno")  )
  
dev.off()




