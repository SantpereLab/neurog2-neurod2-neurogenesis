# FLANKING DINUCLEOTIDE ######################
##############################################


##### NEUROG2
# Detect flanking dinucleotide.
neurog2_extended_peaks_motifs<-flanking_dinucleotides(neurog2_extended_peaks_motifs)
neurod2_extended_peaks_motifs<-flanking_dinucleotides(neurod2_extended_peaks_motifs)

# Filter to only select the functional motifs.
neurog2_extended_peaks_motifs<-neurog2_extended_peaks_motifs %>% filter(central_dinucleotide %in% c("CAT-CAT","CAT-CAG","CAT-CAC","CAG-CAG","CAG-CAC"))
neurod2_extended_peaks_motifs<-neurod2_extended_peaks_motifs %>% filter(central_dinucleotide %in% c("CAT-CAT","CAT-CAG","CAT-CAC","CAG-CAG","CAG-CAC"))

# Plot the flanking dinucleotides of each half-site.
neurog2_balloon_flank<-flanking_motifs_balloon(peaks = neurog2_extended_peaks_motifs, distance_around_summit = 50, back_dist=100, lim = 2)
neurod2_balloon_flank<-flanking_motifs_balloon(peaks = neurod2_extended_peaks_motifs, distance_around_summit = 50, back_dist=100, lim=2)




pdf("~/flanking_balloon.pdf", useDingbats = F, width = 5)
print(neurog2_balloon_flank)
print(neurod2_balloon_flank)
dev.off()


#mirando solo el half site, no el half site de al lado
neurog2_balloon_flank_summarized<-flanking_motifs_balloon_2(peaks = neurog2_extended_peaks_motifs, distance_around_summit = 50, back_dist=100 , lim=2)
neurod2_balloon_flank_summarized<-flanking_motifs_balloon_2(peaks = neurod2_extended_peaks_motifs, distance_around_summit = 50, back_dist=100, lim=2)

pdf("~/flanking_balloon_summarised.pdf", useDingbats = F, width = 3)
print(neurog2_balloon_flank_summarized)
print(neurod2_balloon_flank_summarized)
dev.off()


