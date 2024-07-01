#!/usr/bin/bash


##########################
##### Download fastqs ####
##########################

# E14.5 NEUROD2
module load SRA-Toolkit
prefetch SRR1951389
fasterq-dump SRR1951389 -o e14.5.neurod2.1.fastq
prefetch SRR1951390
fasterq-dump SRR1951390 -o e14.5.neurod2.2.fastq
prefetch SRR1951392
fasterq-dump SRR1951392 -o e14.5.gfp.1.fastq
prefetch SRR1951393
fasterq-dump SRR1951393 -o e14.5.gfp.2.fastq

# P0 NEUROD2
prefetch SRR3955796
fasterq-dump SRR3955796 -o p0.neurod2.1.fastq
prefetch SRR3955797
fasterq-dump SRR3955797 -o p0.neurod2.2.fastq
prefetch SRR3955799
fasterq-dump SRR3955799 -o p0.gfp.fastq

# E14.5 NEUROG2
prefetch SRR1662317
fasterq-dump SRR1662317 -o neurog2.fastq
prefetch SRR1662315
fasterq-dump SRR1662315 -o input.fastq


#############################################
########## ALIGNMENT ########################
#############################################

module load BWA
module load SAMtools

#align
for file in *fastq
do
  bwa mem -t 20 ~/mm10/mm10.fa $file | samtools rmdup -s - - | samtools view -h > $file.sam; \
  cat $file.sam | grep -v -e 'XA:Z:' -e 'SA:Z:' > $file.uniqmap; cat $file.uniqmap | samtools view -bS - > $file.bam
done

rm *sam*


##############################################################
################### PEAK CALLING ############################
##############################################################

#NEUROD2 E14.5
macs3 callpeak -g mm -f BAM -t e14.5.neurod2.1.fastq.bam -c e14.5.gfp.1.fastq.bam --bw 200 --outdir . --bdg -n e14.5.neurod2.1
macs3 callpeak -g mm -f BAM -t e14.5.neurod2.2.fastq.bam -c e14.5.gfp.2.fastq.bam --bw 200 --outdir . --bdg -n e14.5.neurod2.2

#NEUROD2 P0
macs3 callpeak -g mm -f BAM -t p0.neurod2.1.fastq.sam.bam -c p0.gfp.fastq.sam.bam --bw 200 --outdir . --bdg -n p0.neurod2.1
macs3 callpeak -g mm -f BAM -t p0.neurod2.2.fastq.sam.bam -c p0.gfp.fastq.sam.bam --bw 200 --outdir . --bdg -n p0.neurod2.2

#NEUROG2
macs3 callpeak -g mm -f BAM -t neurog2.bam -c input.bam --bw 200 --outdir . --bdg -n e14.5.neurog2


#Add header to the files
sed -i '1i\track type=bedGraph' e14.5.neurod2.1_treat_pileup.bdg
sed -i '1i\track type=bedGraph' e14.5.neurod2.2_treat_pileup.bdg

sed -i '1i\track type=bedGraph' p0.neurod2.1_treat_pileup.bdg
sed -i '1i\track type=bedGraph' p0.neurod2.2_treat_pileup.bdg

sed -i '1i\track type=bedGraph' e14.5.neurog2_treat_pileup.bdg


## Split peaks in subpeaks
##########################

java -jar ~/opt/PeakSplitter_Java/PeakSplitter.jar -p e14.5.neurod2.1_peaks.narrowPeak -w e14.5.neurod2.1_treat_pileup.bdg -o . -x 1 -f false
java -jar ~/opt/PeakSplitter_Java/PeakSplitter.jar -p e14.5.neurod2.2_peaks.narrowPeak -w e14.5.neurod2.2_treat_pileup.bdg -o . -x 2 -f false

java -jar ~/opt/PeakSplitter_Java/PeakSplitter.jar -p p0.neurod2.1_peaks.narrowPeak -w p0.neurod2.1_treat_pileup.bdg -o . -x 1 -f false
java -jar ~/opt/PeakSplitter_Java/PeakSplitter.jar -p p0.neurod2.2_peaks.narrowPeak -w p0.neurod2.2_treat_pileup.bdg -o . -x 2 -f false

java -jar ~/opt/PeakSplitter_Java/PeakSplitter.jar -p e14.5.neurog2_peaks.narrowPeak -w e14.5.neurog2_treat_pileup.bdg -o . -x 1 -f false


sed -i 1d 1.e14.5.neurod2.1_peaks.subpeaks.narrowPeak
sed -i 1d 2.e14.5.neurod2.2_peaks.subpeaks.narrowPeak

sed -i 1d 1.p0.neurod2.1_peaks.subpeaks.narrowPeak
sed -i 1d 2.p0.neurod2.2_peaks.subpeaks.narrowPeak

sed -i 1d 1.e14.5.neurog2_peaks.subpeaks.narrowPeak


cat 1.e14.5.neurog2_peaks.subpeaks.narrowPeak | cut -f1,2,3,5 | awk '{printf "peak%d\t%s\n", NR, $0}' | awk '{print $2,$3,$4,$5,$1}' | tr " " "\t" > e14.5.neurog2_peaks.subpeaks.bed

# Obtain high confidence peaks that are present in more than one sample
#######################################################################

intersectBed -a 1.e14.5.neurod2.1_peaks.subpeaks.narrowPeak -b 2.e14.5.neurod2.2_peaks.subpeaks.narrowPeak | cut -f1,2,3 | egrep -v "GL|JH|chrM" > e14.5.neurod2_highconfsubpeaks.bed
intersectBed -a 1.p0.neurod2.1_peaks.subpeaks.narrowPeak -b 2.p0.neurod2.2_peaks.subpeaks.narrowPeak | cut -f1,2,3 | egrep -v "GL|JH|chrM" > p0.neurod2_highconfsubpeaks.bed


# Require that the two summits fall within the high confidence peaks

# NEUROD2

awk '{print $1,$5-1,$5}' 1.e14.5.neurod2.1_peaks.subpeaks.narrowPeak | tr " " "\t" > summits_1.tmp
awk '{print $1,$5-1,$5}' 2.e14.5.neurod2.2_peaks.subpeaks.narrowPeak | tr " " "\t" > summits_2.tmp

paste e14.5.neurod2_highconfsubpeaks.bed <(intersectBed -wao -a e14.5.neurod2_highconfsubpeaks.bed -b summits_1.tmp | cut -f6) \
<(intersectBed -wao -a e14.5.neurod2_highconfsubpeaks.bed -b summits_2.tmp | cut -f6) | grep -v "-" | \
awk '{print $1,$2,$3,int(($5-$4)/2)+$4}' | tr " " "\t" | awk '{printf "peak%d\t%s\n", NR, $0}' | awk '{print $2,$3,$4,$5,$1}' | tr " " "\t" > e14.5.neurod2_highconfsubpeaks_summit.bed

# NEUROG2

awk '{print $1,$5-1,$5}' 1.p0.neurod2.1_peaks.subpeaks.narrowPeak | tr " " "\t" > summits_1.tmp
awk '{print $1,$5-1,$5}' 2.p0.neurod2.2_peaks.subpeaks.narrowPeak | tr " " "\t" > summits_2.tmp

paste p0.neurod2_highconfsubpeaks.bed <(intersectBed -wao -a p0.neurod2_highconfsubpeaks.bed -b summits_1.tmp | cut -f6) \
<(intersectBed -wao -a p0.neurod2_highconfsubpeaks.bed -b summits_2.tmp | cut -f6) | grep -v "-" | \
awk '{print $1,$2,$3,int(($5-$4)/2)+$4}' | tr " " "\t" | awk '{printf "peak%d\t%s\n", NR, $0}' | awk '{print $2,$3,$4,$5,$1}' | tr " " "\t" > p0.neurod2_highconfsubpeaks_summit.bed

rm *tmp




