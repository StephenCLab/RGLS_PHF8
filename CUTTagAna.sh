#!/bin/bash

# Usage ： bash CUTTAGAna.sh LZandM_20260313 1.rawdata _1.fq.gz _2.fq.gz mouse LZandM
cd ./CUTTag
echo "Prepare..."

projectname=$1
mkdir $projectname
cd $projectname

mkdir 1_rawdata 2_clean 3_sam 3_bam 4_rmdup 5_bw 6_callpeak

rawdatapath=$2
fileendr1=$3
fileendr2=$4

echo "copy data..."
ls $rawdatapath/*/*$fileendr1 |while read id ; do 
    sample_name=$(basename "$id" $fileendr1)
    cp "$id" "./1_rawdata/${sample_name}_rep_R1.fq.gz";
done
ls $rawdatapath/*/*$fileendr2 |while read id ; do 
    sample_name=$(basename "$id" $fileendr2)
    cp "$id" "./1_rawdata/${sample_name}_rep_R2.fq.gz";
done
echo "Finish Copy & Prepare..."

#human/mouse
datatype=$5
pre=$6

# 1. clean data : 
source activate RNA-seq
ls 1_rawdata/*_R1.fq.gz |while read id ; do trim_galore -j 50 -q 30 --phred33 --length 120 -e 0.1  --paired -o 2_clean/ $id ${id%_R1.*gz}_R2.fq.gz;done
# 2.map: 
cd 2_clean
if [[ "$datatype" == "mouse" ]]; then
    callpeakgenome='mm'
    anno_db='mm10'
    echo "Running bowtie2 for MOUSE (mm10 reference)..."
    ls *R1*.fq.gz  | while read id ; do bowtie2  -p 50  -x /DB/bowtie2_mouse/mm10  --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -1 ${id%_R1*}_R1_val_1.fq.gz -2 ${id%_R1*}_R2_val_2.fq.gz  -S ../3_sam/${id%_R1*}.sam;done
elif [ "$datatype" == "human" ]; then
    callpeakgenome='hs'
    anno_db='hg38'
    echo "Running bowtie2 for HUMAN (GRCh38 reference)..."
    ls *R1*.fq.gz  | while read id ; do bowtie2  -p 50  -x /DB/bowtie2_human/GRCh38/GRCh38  --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -1 ${id%_R1*}_R1_val_1.fq.gz -2 ${id%_R1*}_R2_val_2.fq.gz  -S ../3_sam/${id%_R1*}.sam;done
fi
# 3.samtobam:
cd ../3_sam
ls *.sam |while read id ;do samtools view -@ 50 -bF 12 -q 30 -S $id > ../3_bam/${id%.sam}.map.q30.F12.bam;done
# 3-2.remove duplication: 
cd ../3_bam
ls *.bam |while read id ;do /home/star/anaconda3/envs/sambamba/bin/sambamba markdup -r $id ../4_rmdup/${id%.bam}.rmdup.bam ;done

# 4 sort bam and make index
cd ../4_rmdup
ls *.rmdup.bam |while read id ;do samtools sort -@ 50 $id -o ${id%.rmdup.bam}.rmdup.sorted.bam;done
ls *.rmdup.sorted.bam |while read id ;do samtools index $id;done
# 5.bam to bw: 
ls *.rmdup.sorted.bam |while read id ;do  bamCoverage -b $id --normalizeUsing RPKM -p 50 -o ../5_bw/${id%.map.q30.F12.rmdup.sorted.bam}.bw ; done
# 6.call peak(adjust for transcrip factor/ Histone modification type)：
ls *.rmdup.sorted.bam |while read id ;do  macs2  callpeak -t $id --bdg -p 1e-5 -g $callpeakgenome -n ../6_callpeak/${id%.map.q30.F12.rmdup.sorted.bam}.peak; done

for p_cut in 1e-2 1e-3 1e-4 1e-5; do
    cd ../4_rmdup
    echo "Processing: $p_cut"
    pre2=$pre$p_cut
    # 6.call peak
    source activate macs2
    macs2  callpeak -t *.map.q30.F12.rmdup.sorted.bam --bdg -p $p_cut -g $callpeakgenome -n '../6_callpeak/'$pre2'.peak'

    # 8.annotation
    annotatePeaks.pl '../6_callpeak/'$pre2'.peak_summits.bed' $anno_db > '../6_callpeak/'$pre2'_peak.anno.txt'
    tail +2 '../6_callpeak/'$pre2'_peak.anno.txt' |cut -f 2,3,4 > '../6_callpeak/'$pre2'_peak.txt'
    # 7.visualization
    source activate ChIP-seq
    cd ../6_callpeak
    computeMatrix  reference-point -p 20 -R $pre2'.peak_summits.bed' -a 3000 -b 3000 -S  ../5_bw/*.bw --skipZeros  -out ./$pre2'.computeMatrix.gz'

    plotHeatmap -m $pre2'.computeMatrix.gz' -o $pre2'.computeMatrix.gz.pdf' --colorMap RdBu --zMin -3 --zMax 3
done

################## merge PHF8、H3K27me1、H3K9me2、H4K20me1，LZ、M samples ，Perform differential analysis after calculating peak values##################
# PHF8

conda activate macs2
macs2  callpeak -t ../4_rmdup/rawbam/M_PHF8*.map.q30.F12.rmdup.sorted.bam --bdg -p 1e-5 -g mm -n M_PHF8
macs2  callpeak -t ../4_rmdup/rawbam/LZ_PHF8*.map.q30.F12.rmdup.sorted.bam --bdg -p 1e-5 -g mm -n LZ_PHF8
conda activate RNA-seq

cut -f1-3 M_PHF8_peaks.narrowPeak | sort -k1,1 -k2,2n > M_PHF8.sorted.bed
cut -f1-3 LZ_PHF8_peaks.narrowPeak | sort -k1,1 -k2,2n > LZ_PHF8.sorted.bed

cat M_PHF8.sorted.bed LZ_PHF8.sorted.bed | sort -k1,1 -k2,2n | bedtools merge -i - > PHF8_union_peaks.bed

bedtools multiinter -header -names M_PHF8 LZ_PHF8 -i M_PHF8.sorted.bed LZ_PHF8.sorted.bed > peak_overlap.txt

bedtools multicov -bams ../4_rmdup/rawbam/M_PHF8.map.q30.F12.rmdup.sorted.bam ../4_rmdup/rawbam/M_PHF8_rep.map.q30.F12.rmdup.sorted.bam ../4_rmdup/rawbam/LZ_PHF8.map.q30.F12.rmdup.sorted.bam ../4_rmdup/rawbam/LZ_PHF8_rep.map.q30.F12.rmdup.sorted.bam -bed PHF8_union_peaks.bed > PHF8_peak_counts.tsv

wc -l M_PHF8_summits.bed
wc -l LZ_PHF8_summits.bed
wc -l M_PHF8_peaks.narrowPeak
wc -l LZ_PHF8_peaks.narrowPeak
cat M_PHF8_summits.bed LZ_PHF8_summits.bed > PHF8_summits.bed
cat M_PHF8_peaks.narrowPeak LZ_PHF8_peaks.narrowPeak > PHF8_peaks.narrowPeak


Rscript --vanilla << 'EOF' 
library(openxlsx)
library(DESeq2)

PHF8_cts <- read.table("rawdata/PHF8_peak_counts.tsv", header=FALSE)
rownames(PHF8_cts) <- paste(PHF8_cts$V1, PHF8_cts$V2, PHF8_cts$V3, sep="_")
PHF8_countData <- PHF8_cts[,4:ncol(PHF8_cts)]

PHF8_coldata <- data.frame(
  row.names = colnames(PHF8_countData),
  condition = c("PHF8_M","PHF8_M","PHF8_LZ","PHF8_LZ")
)
PHF8_dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(PHF8_countData)),
  colData = PHF8_coldata,
  design = ~ condition
)
PHF8_dds <- DESeq(PHF8_dds)
PHF8_res <- results(PHF8_dds, contrast=c("condition","PHF8_LZ","PHF8_M"))
PHF8_res <- as.data.frame(PHF8_res)
PHF8_res_new <- PHF8_res %>%
  tibble::rownames_to_column("region") %>%
   extract(region, 
          into = c("chr", "start", "end"),
          regex = "^(.*)_(\\d+)_(\\d+)$",
          remove = TRUE,
          convert = TRUE)

rownames(PHF8_res_new)=rownames(PHF8_res)
write.xlsx(PHF8_res_new,'PHF8_res.xlsx',row.names=TRUE)
PHF8_summits=read.table("rawdata/PHF8_summits.bed", header=FALSE)
PHF8_summits=PHF8_summits[,c(1:4)]
colnames(PHF8_summits)=c('SummitChr','SummitStart','SummitEnd','PeakNames')
head(PHF8_summits)
PHF8_narrowPeak=read.table("rawdata/PHF8_peaks.narrowPeak", header=FALSE)
PHF8_narrowPeak=PHF8_narrowPeak[,c(1:4)]
colnames(PHF8_narrowPeak)=c('narrowPeakChr','narrowPeakStart','narrowPeakEnd','PeakNames')
head(PHF8_narrowPeak)

library(dplyr)
PHF8_narrowPeak_unique <- PHF8_narrowPeak %>%
  distinct(narrowPeakChr, narrowPeakStart, narrowPeakEnd, .keep_all = TRUE)
rownames(PHF8_narrowPeak_unique) <- paste(
  PHF8_narrowPeak_unique$narrowPeakChr, 
  PHF8_narrowPeak_unique$narrowPeakStart, 
  PHF8_narrowPeak_unique$narrowPeakEnd, 
  sep = "_"
)
PHF8_narrowPeak_unique$rownames=rownames(PHF8_narrowPeak_unique)
head(PHF8_narrowPeak_unique)

PHF8_deseq.bedup <- PHF8_res_new[which(PHF8_res_new$padj < 0.05 & PHF8_res_new$log2FoldChange > 0), c("chr", "start", "end")]
PHF8_deseq.beddown <- PHF8_res_new[which(PHF8_res_new$padj < 0.05 & PHF8_res_new$log2FoldChange < 0), c("chr", "start", "end")]
PHF8_remaining <- !(PHF8_res_new$padj < 0.05 & PHF8_res_new$log2FoldChange > 0 | PHF8_res_new$padj < 0.05 & PHF8_res_new$log2FoldChange < 0)

PHF8_deseq.bedcommon <- PHF8_res_new[PHF8_remaining,  c("chr", "start", "end")]
write.table(PHF8_deseq.bedup, file="PHF8_deseq.bedup.bed", sep="\t",  quote=F, row.names=F, col.names=F)
write.table(PHF8_deseq.beddown, file="PHF8_deseq.beddown.bed", sep="\t",  quote=F, row.names=F, col.names=F)
write.table(PHF8_deseq.bedcommon, file="PHF8_deseq.bedcommon.bed", sep="\t",  quote=F, row.names=F, col.names=F)
# up
length(rownames(PHF8_deseq.bedup))
PHF8_deseq.bedup$rownames=rownames(PHF8_deseq.bedup)
head(PHF8_deseq.bedup)
PHF8_narrowPeak_up=merge(PHF8_narrowPeak_unique,PHF8_deseq.bedup,by='rownames')
length(rownames(PHF8_narrowPeak_up))
length(rownames(PHF8_deseq.bedup))
PHF8_narrowPeak_up=PHF8_narrowPeak_up[,c(2:5)]
PHF8_summits_up=merge(PHF8_narrowPeak_up,PHF8_summits,by='PeakNames')
PHF8_summits_up_final=PHF8_summits_up[,c('SummitChr','SummitStart','SummitEnd')]
write.table(PHF8_summits_up_final, file="PHF8_summits_up_LZvsM.bed", sep="\t",  quote=F, row.names=F, col.names=F)
PHF8_narrowPeak_up_final=PHF8_summits_up[,c('narrowPeakChr','narrowPeakStart','narrowPeakEnd')]
write.table(PHF8_narrowPeak_up_final, file="PHF8_narrowPeak_up_LZvsM.bed", sep="\t",  quote=F, row.names=F, col.names=F)
# down
length(rownames(PHF8_deseq.beddown))
PHF8_deseq.beddown$rownames=rownames(PHF8_deseq.beddown)
head(PHF8_deseq.beddown)
PHF8_narrowPeak_down=merge(PHF8_narrowPeak_unique,PHF8_deseq.beddown,by='rownames')
length(rownames(PHF8_narrowPeak_down))
length(rownames(PHF8_deseq.beddown))
PHF8_narrowPeak_down=PHF8_narrowPeak_down[,c(2:5)]
PHF8_summits_down=merge(PHF8_narrowPeak_down,PHF8_summits,by='PeakNames')
PHF8_summits_down_final=PHF8_summits_down[,c('SummitChr','SummitStart','SummitEnd')]
write.table(PHF8_summits_down_final, file="PHF8_summits_down_LZvsM.bed", sep="\t",  quote=F, row.names=F, col.names=F)
PHF8_narrowPeak_down_final=PHF8_summits_down[,c('narrowPeakChr','narrowPeakStart','narrowPeakEnd')]
write.table(PHF8_narrowPeak_down_final, file="PHF8_narrowPeak_down_LZvsM.bed", sep="\t",  quote=F, row.names=F, col.names=F)
# common
length(rownames(PHF8_deseq.bedcommon))
PHF8_deseq.bedcommon$rownames=rownames(PHF8_deseq.bedcommon)
head(PHF8_deseq.bedcommon)
PHF8_narrowPeak_common=merge(PHF8_narrowPeak_unique,PHF8_deseq.bedcommon,by='rownames')
length(rownames(PHF8_narrowPeak_common))
length(rownames(PHF8_deseq.bedcommon))
PHF8_narrowPeak_common=PHF8_narrowPeak_common[,c(2:5)]
PHF8_summits_common=merge(PHF8_narrowPeak_common,PHF8_summits,by='PeakNames')
PHF8_summits_common_final=PHF8_summits_common[,c('SummitChr','SummitStart','SummitEnd')]
write.table(PHF8_summits_common_final, file="PHF8_summits_common_LZvsM.bed", sep="\t",  quote=F, row.names=F, col.names=F)
PHF8_narrowPeak_common_final=PHF8_summits_common[,c('narrowPeakChr','narrowPeakStart','narrowPeakEnd')]
write.table(PHF8_narrowPeak_common_final, file="PHF8_narrowPeak_common_LZvsM.bed", sep="\t",  quote=F, row.names=F, col.names=F)
save.image('PHF8_deg.RData')

EOF


computeMatrix reference-point  -p 56 -S M_PHF8.bw LZ_PHF8.bw -R PHF8_summits_down_LZvsM.bed PHF8_summits_common_LZvsM.bed PHF8_summits_up_LZvsM.bed -b 3000 -a 3000 -o PHF8.matrix.gz
plotHeatmap -m PHF8.matrix.gz -out PHF8.peak_heatmap.pdf --regionsLabel "Model" "Common" "LZ" --samplesLabel "PHF8_Model" "PHF8_LZ"

#H3K27me1
conda activate macs2

macs2  callpeak -t ../4_rmdup/rawbam/M_H3K27*.map.q30.F12.rmdup.sorted.bam --bdg -p 1e-5 -g mm -n M_H3K27
macs2  callpeak -t ../4_rmdup/rawbam/LZ_H3K27*.map.q30.F12.rmdup.sorted.bam --bdg -p 1e-5 -g mm -n LZ_H3K27

conda activate RNA-seq

cut -f1-3 M_H3K27_peaks.narrowPeak | sort -k1,1 -k2,2n > M_H3K27.sorted.bed
cut -f1-3 LZ_H3K27_peaks.narrowPeak | sort -k1,1 -k2,2n > LZ_H3K27.sorted.bed

cat M_H3K27.sorted.bed LZ_H3K27.sorted.bed | sort -k1,1 -k2,2n | bedtools merge -i - > H3K27_union_peaks.bed

bedtools multiinter -header -names M_H3K27 LZ_H3K27 -i M_H3K27.sorted.bed LZ_H3K27.sorted.bed > peak_overlap.txt

bedtools multicov -bams ../4_rmdup/rawbam/M_H3K27.map.q30.F12.rmdup.sorted.bam ../4_rmdup/rawbam/M_H3K27_rep.map.q30.F12.rmdup.sorted.bam ../4_rmdup/rawbam/LZ_H3K27.map.q30.F12.rmdup.sorted.bam ../4_rmdup/rawbam/LZ_H3K27_rep.map.q30.F12.rmdup.sorted.bam -bed H3K27_union_peaks.bed > H3K27_peak_counts.tsv

wc -l M_H3K27_summits.bed
wc -l LZ_H3K27_summits.bed
wc -l M_H3K27_peaks.narrowPeak
wc -l LZ_H3K27_peaks.narrowPeak
cat M_H3K27_summits.bed LZ_H3K27_summits.bed > H3K27_summits.bed
cat M_H3K27_peaks.narrowPeak LZ_H3K27_peaks.narrowPeak > H3K27_peaks.narrowPeak

Rscript --vanilla << 'EOF'
library(openxlsx)
library(DESeq2)

H3K27_cts <- read.table("rawdata/H3K27_peak_counts.tsv", header=FALSE)
rownames(H3K27_cts) <- paste(H3K27_cts$V1, H3K27_cts$V2, H3K27_cts$V3, sep="_")
H3K27_countData <- H3K27_cts[,4:ncol(H3K27_cts)]

H3K27_coldata <- data.frame(
  row.names = colnames(H3K27_countData),
  condition = c("H3K27_M","H3K27_M","H3K27_LZ","H3K27_LZ")
)
H3K27_dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(H3K27_countData)),
  colData = H3K27_coldata,
  design = ~ condition
)
H3K27_dds <- DESeq(H3K27_dds)
H3K27_res <- results(H3K27_dds, contrast=c("condition","H3K27_LZ","H3K27_M"))
H3K27_res <- as.data.frame(H3K27_res)

H3K27_res_new <- H3K27_res %>%
  tibble::rownames_to_column("region") %>%
   extract(region, 
          into = c("chr", "start", "end"),
          regex = "^(.*)_(\\d+)_(\\d+)$",
          remove = TRUE,
          convert = TRUE)

rownames(H3K27_res_new)=rownames(H3K27_res)

write.xlsx(H3K27_res_new,'H3K27_res.xlsx',row.names=TRUE)


H3K27_summits=read.table("rawdata/H3K27_summits.bed", header=FALSE)
H3K27_summits=H3K27_summits[,c(1:4)]
colnames(H3K27_summits)=c('SummitChr','SummitStart','SummitEnd','PeakNames')
head(H3K27_summits)
H3K27_narrowPeak=read.table("rawdata/H3K27_peaks.narrowPeak", header=FALSE)
H3K27_narrowPeak=H3K27_narrowPeak[,c(1:4)]
colnames(H3K27_narrowPeak)=c('narrowPeakChr','narrowPeakStart','narrowPeakEnd','PeakNames')
head(H3K27_narrowPeak)

library(dplyr)

H3K27_narrowPeak_unique <- H3K27_narrowPeak %>%
  distinct(narrowPeakChr, narrowPeakStart, narrowPeakEnd, .keep_all = TRUE)
rownames(H3K27_narrowPeak_unique) <- paste(
  H3K27_narrowPeak_unique$narrowPeakChr, 
  H3K27_narrowPeak_unique$narrowPeakStart, 
  H3K27_narrowPeak_unique$narrowPeakEnd, 
  sep = "_"
)
H3K27_narrowPeak_unique$rownames=rownames(H3K27_narrowPeak_unique)
head(H3K27_narrowPeak_unique)

H3K27_deseq.bedup <- H3K27_res_new[which(H3K27_res_new$padj < 0.05 & H3K27_res_new$log2FoldChange > 0), c("chr", "start", "end")]

H3K27_deseq.beddown <- H3K27_res_new[which(H3K27_res_new$padj < 0.05 & H3K27_res_new$log2FoldChange < 0), c("chr", "start", "end")]

H3K27_deseq.bedcommon <- H3K27_res_new[which(H3K27_res_new$padj > 0.05), c("chr", "start", "end")]
H3K27_remaining <- !(H3K27_res_new$padj < 0.05 & H3K27_res_new$log2FoldChange > 0 | H3K27_res_new$padj < 0.05 & H3K27_res_new$log2FoldChange < 0)

H3K27_deseq.bedcommon <- H3K27_res_new[H3K27_remaining,  c("chr", "start", "end")]
write.table(H3K27_deseq.bedup, file="H3K27_deseq.bedup.bed", sep="\t",  quote=F, row.names=F, col.names=F)
write.table(H3K27_deseq.beddown, file="H3K27_deseq.beddown.bed", sep="\t",  quote=F, row.names=F, col.names=F)
write.table(H3K27_deseq.bedcommon, file="H3K27_deseq.bedcommon.bed", sep="\t",  quote=F, row.names=F, col.names=F)

# up
length(rownames(H3K27_deseq.bedup))
H3K27_deseq.bedup$rownames=rownames(H3K27_deseq.bedup)
head(H3K27_deseq.bedup)
H3K27_narrowPeak_up=merge(H3K27_narrowPeak_unique,H3K27_deseq.bedup,by='rownames')
length(rownames(H3K27_narrowPeak_up))
length(rownames(H3K27_deseq.bedup))
H3K27_narrowPeak_up=H3K27_narrowPeak_up[,c(2:5)]
H3K27_summits_up=merge(H3K27_narrowPeak_up,H3K27_summits,by='PeakNames')
H3K27_summits_up_final=H3K27_summits_up[,c('SummitChr','SummitStart','SummitEnd')]
write.table(H3K27_summits_up_final, file="H3K27_summits_up_LZvsM.bed", sep="\t",  quote=F, row.names=F, col.names=F)
H3K27_narrowPeak_up_final=H3K27_summits_up[,c('narrowPeakChr','narrowPeakStart','narrowPeakEnd')]
write.table(H3K27_narrowPeak_up_final, file="H3K27_narrowPeak_up_LZvsM.bed", sep="\t",  quote=F, row.names=F, col.names=F)

# down
length(rownames(H3K27_deseq.beddown))
H3K27_deseq.beddown$rownames=rownames(H3K27_deseq.beddown)
head(H3K27_deseq.beddown)
H3K27_narrowPeak_down=merge(H3K27_narrowPeak_unique,H3K27_deseq.beddown,by='rownames')
length(rownames(H3K27_narrowPeak_down))
length(rownames(H3K27_deseq.beddown))
H3K27_narrowPeak_down=H3K27_narrowPeak_down[,c(2:5)]
H3K27_summits_down=merge(H3K27_narrowPeak_down,H3K27_summits,by='PeakNames')
H3K27_summits_down_final=H3K27_summits_down[,c('SummitChr','SummitStart','SummitEnd')]
write.table(H3K27_summits_down_final, file="H3K27_summits_down_LZvsM.bed", sep="\t",  quote=F, row.names=F, col.names=F)
H3K27_narrowPeak_down_final=H3K27_summits_down[,c('narrowPeakChr','narrowPeakStart','narrowPeakEnd')]
write.table(H3K27_narrowPeak_down_final, file="H3K27_narrowPeak_down_LZvsM.bed", sep="\t",  quote=F, row.names=F, col.names=F)

# common
length(rownames(H3K27_deseq.bedcommon))
H3K27_deseq.bedcommon$rownames=rownames(H3K27_deseq.bedcommon)
head(H3K27_deseq.bedcommon)
H3K27_narrowPeak_common=merge(H3K27_narrowPeak_unique,H3K27_deseq.bedcommon,by='rownames')
length(rownames(H3K27_narrowPeak_common))
length(rownames(H3K27_deseq.bedcommon))
H3K27_narrowPeak_common=H3K27_narrowPeak_common[,c(2:5)]
H3K27_summits_common=merge(H3K27_narrowPeak_common,H3K27_summits,by='PeakNames')
H3K27_summits_common_final=H3K27_summits_common[,c('SummitChr','SummitStart','SummitEnd')]
write.table(H3K27_summits_common_final, file="H3K27_summits_common_LZvsM.bed", sep="\t",  quote=F, row.names=F, col.names=F)
H3K27_narrowPeak_common_final=H3K27_summits_common[,c('narrowPeakChr','narrowPeakStart','narrowPeakEnd')]
write.table(H3K27_narrowPeak_common_final, file="H3K27_narrowPeak_common_LZvsM.bed", sep="\t",  quote=F, row.names=F, col.names=F)

save.image('H3K27_deg.RData')
EOF

computeMatrix reference-point  -p 56 -S M_H3K27.bw LZ_H3K27.bw -R H3K27_summits_down_LZvsM.bed H3K27_summits_common_LZvsM.bed H3K27_summits_up_LZvsM.bed -b 3000 -a 3000 -o H3K27.matrix.gz
plotHeatmap -m H3K27.matrix.gz -out H3K27.peak_heatmap.pdf --regionsLabel "Model" "Common" "LZ" --samplesLabel "H3K27_Model" "H3K27_LZ"

# H3K9me2
conda activate macs2

macs2  callpeak -t ../4_rmdup/rawbam/M_H3K9*.map.q30.F12.rmdup.sorted.bam --bdg -p 1e-5 -g mm -n M_H3K9
macs2  callpeak -t ../4_rmdup/rawbam/LZ_H3K9*.map.q30.F12.rmdup.sorted.bam --bdg -p 1e-5 -g mm -n LZ_H3K9

conda activate RNA-seq

cut -f1-3 M_H3K9_peaks.narrowPeak | sort -k1,1 -k2,2n > M_H3K9.sorted.bed
cut -f1-3 LZ_H3K9_peaks.narrowPeak | sort -k1,1 -k2,2n > LZ_H3K9.sorted.bed

cat M_H3K9.sorted.bed LZ_H3K9.sorted.bed | sort -k1,1 -k2,2n | bedtools merge -i - > H3K9_union_peaks.bed

bedtools multiinter -header -names M_H3K9 LZ_H3K9 -i M_H3K9.sorted.bed LZ_H3K9.sorted.bed > peak_overlap.txt

bedtools multicov -bams ../4_rmdup/rawbam/M_H3K9.map.q30.F12.rmdup.sorted.bam ../4_rmdup/rawbam/M_H3K9_rep.map.q30.F12.rmdup.sorted.bam ../4_rmdup/rawbam/LZ_H3K9.map.q30.F12.rmdup.sorted.bam ../4_rmdup/rawbam/LZ_H3K9_rep.map.q30.F12.rmdup.sorted.bam -bed H3K9_union_peaks.bed > H3K9_peak_counts.tsv

wc -l M_H3K9_summits.bed
wc -l LZ_H3K9_summits.bed
wc -l M_H3K9_peaks.narrowPeak
wc -l LZ_H3K9_peaks.narrowPeak
cat M_H3K9_summits.bed LZ_H3K9_summits.bed > H3K9_summits.bed
cat M_H3K9_peaks.narrowPeak LZ_H3K9_peaks.narrowPeak > H3K9_peaks.narrowPeak

Rscript --vanilla << 'EOF'
library(openxlsx)
library(DESeq2)

H3K9_cts <- read.table("rawdata/H3K9_peak_counts.tsv", header=FALSE)
rownames(H3K9_cts) <- paste(H3K9_cts$V1, H3K9_cts$V2, H3K9_cts$V3, sep="_")
H3K9_countData <- H3K9_cts[,4:ncol(H3K9_cts)]

H3K9_coldata <- data.frame(
  row.names = colnames(H3K9_countData),
  condition = c("H3K9_M","H3K9_M","H3K9_LZ","H3K9_LZ")
)
H3K9_dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(H3K9_countData)),
  colData = H3K9_coldata,
  design = ~ condition
)
H3K9_dds <- DESeq(H3K9_dds)
H3K9_res <- results(H3K9_dds, contrast=c("condition","H3K9_LZ","H3K9_M"))
H3K9_res <- as.data.frame(H3K9_res)

H3K9_res_new <- H3K9_res %>%
  tibble::rownames_to_column("region") %>%
   extract(region, 
          into = c("chr", "start", "end"),
          regex = "^(.*)_(\\d+)_(\\d+)$",
          remove = TRUE,
          convert = TRUE)

rownames(H3K9_res_new)=rownames(H3K9_res)

write.xlsx(H3K9_res_new,'H3K9_res.xlsx',row.names=TRUE)

H3K9_summits=read.table("rawdata/H3K9_summits.bed", header=FALSE)
H3K9_summits=H3K9_summits[,c(1:4)]
colnames(H3K9_summits)=c('SummitChr','SummitStart','SummitEnd','PeakNames')
head(H3K9_summits)
H3K9_narrowPeak=read.table("rawdata/H3K9_peaks.narrowPeak", header=FALSE)
H3K9_narrowPeak=H3K9_narrowPeak[,c(1:4)]
colnames(H3K9_narrowPeak)=c('narrowPeakChr','narrowPeakStart','narrowPeakEnd','PeakNames')
head(H3K9_narrowPeak)

library(dplyr)

H3K9_narrowPeak_unique <- H3K9_narrowPeak %>%
  distinct(narrowPeakChr, narrowPeakStart, narrowPeakEnd, .keep_all = TRUE)
rownames(H3K9_narrowPeak_unique) <- paste(
  H3K9_narrowPeak_unique$narrowPeakChr, 
  H3K9_narrowPeak_unique$narrowPeakStart, 
  H3K9_narrowPeak_unique$narrowPeakEnd, 
  sep = "_"
)
H3K9_narrowPeak_unique$rownames=rownames(H3K9_narrowPeak_unique)
head(H3K9_narrowPeak_unique)

H3K9_deseq.bedup <- H3K9_res_new[which(H3K9_res_new$padj < 0.05 & H3K9_res_new$log2FoldChange > 0), c("chr", "start", "end")]

H3K9_deseq.beddown <- H3K9_res_new[which(H3K9_res_new$padj < 0.05 & H3K9_res_new$log2FoldChange < 0), c("chr", "start", "end")]

H3K9_deseq.bedcommon <- H3K9_res_new[which(H3K9_res_new$padj > 0.05), c("chr", "start", "end")]
H3K9_remaining <- !(H3K9_res_new$padj < 0.05 & H3K9_res_new$log2FoldChange > 0 | H3K9_res_new$padj < 0.05 & H3K9_res_new$log2FoldChange < 0)

H3K9_deseq.bedcommon <- H3K9_res_new[H3K9_remaining,  c("chr", "start", "end")]
write.table(H3K9_deseq.bedup, file="H3K9_deseq.bedup.bed", sep="\t",  quote=F, row.names=F, col.names=F)
write.table(H3K9_deseq.beddown, file="H3K9_deseq.beddown.bed", sep="\t",  quote=F, row.names=F, col.names=F)
write.table(H3K9_deseq.bedcommon, file="H3K9_deseq.bedcommon.bed", sep="\t",  quote=F, row.names=F, col.names=F)

# up
length(rownames(H3K9_deseq.bedup))
H3K9_deseq.bedup$rownames=rownames(H3K9_deseq.bedup)
head(H3K9_deseq.bedup)
H3K9_narrowPeak_up=merge(H3K9_narrowPeak_unique,H3K9_deseq.bedup,by='rownames')
length(rownames(H3K9_narrowPeak_up))
length(rownames(H3K9_deseq.bedup))
H3K9_narrowPeak_up=H3K9_narrowPeak_up[,c(2:5)]
H3K9_summits_up=merge(H3K9_narrowPeak_up,H3K9_summits,by='PeakNames')
H3K9_summits_up_final=H3K9_summits_up[,c('SummitChr','SummitStart','SummitEnd')]
write.table(H3K9_summits_up_final, file="H3K9_summits_up_LZvsM.bed", sep="\t",  quote=F, row.names=F, col.names=F)
H3K9_narrowPeak_up_final=H3K9_summits_up[,c('narrowPeakChr','narrowPeakStart','narrowPeakEnd')]
write.table(H3K9_narrowPeak_up_final, file="H3K9_narrowPeak_up_LZvsM.bed", sep="\t",  quote=F, row.names=F, col.names=F)

# down
length(rownames(H3K9_deseq.beddown))
H3K9_deseq.beddown$rownames=rownames(H3K9_deseq.beddown)
head(H3K9_deseq.beddown)
H3K9_narrowPeak_down=merge(H3K9_narrowPeak_unique,H3K9_deseq.beddown,by='rownames')
length(rownames(H3K9_narrowPeak_down))
length(rownames(H3K9_deseq.beddown))
H3K9_narrowPeak_down=H3K9_narrowPeak_down[,c(2:5)]
H3K9_summits_down=merge(H3K9_narrowPeak_down,H3K9_summits,by='PeakNames')
H3K9_summits_down_final=H3K9_summits_down[,c('SummitChr','SummitStart','SummitEnd')]
write.table(H3K9_summits_down_final, file="H3K9_summits_down_LZvsM.bed", sep="\t",  quote=F, row.names=F, col.names=F)
H3K9_narrowPeak_down_final=H3K9_summits_down[,c('narrowPeakChr','narrowPeakStart','narrowPeakEnd')]
write.table(H3K9_narrowPeak_down_final, file="H3K9_narrowPeak_down_LZvsM.bed", sep="\t",  quote=F, row.names=F, col.names=F)

# common
length(rownames(H3K9_deseq.bedcommon))
H3K9_deseq.bedcommon$rownames=rownames(H3K9_deseq.bedcommon)
head(H3K9_deseq.bedcommon)
H3K9_narrowPeak_common=merge(H3K9_narrowPeak_unique,H3K9_deseq.bedcommon,by='rownames')
length(rownames(H3K9_narrowPeak_common))
length(rownames(H3K9_deseq.bedcommon))
H3K9_narrowPeak_common=H3K9_narrowPeak_common[,c(2:5)]
H3K9_summits_common=merge(H3K9_narrowPeak_common,H3K9_summits,by='PeakNames')
H3K9_summits_common_final=H3K9_summits_common[,c('SummitChr','SummitStart','SummitEnd')]
write.table(H3K9_summits_common_final, file="H3K9_summits_common_LZvsM.bed", sep="\t",  quote=F, row.names=F, col.names=F)
H3K9_narrowPeak_common_final=H3K9_summits_common[,c('narrowPeakChr','narrowPeakStart','narrowPeakEnd')]
write.table(H3K9_narrowPeak_common_final, file="H3K9_narrowPeak_common_LZvsM.bed", sep="\t",  quote=F, row.names=F, col.names=F)

save.image('H3K9_deg.RData')
EOF

computeMatrix reference-point  -p 56 -S M_H3K9.bw LZ_H3K9.bw -R H3K9_summits_down_LZvsM.bed H3K9_summits_common_LZvsM.bed H3K9_summits_up_LZvsM.bed -b 3000 -a 3000 -o H3K9.matrix.gz
plotHeatmap -m H3K9.matrix.gz -out H3K9.peak_heatmap.pdf --regionsLabel "Model" "Common" "LZ" --samplesLabel "H3K9_Model" "H3K9_LZ"

# H4K20me1
conda activate macs2

macs2  callpeak -t ../4_rmdup/rawbam/M_H4K20*.map.q30.F12.rmdup.sorted.bam --bdg -p 1e-5 -g mm -n M_H4K20
macs2  callpeak -t ../4_rmdup/rawbam/LZ_H4K20*.map.q30.F12.rmdup.sorted.bam --bdg -p 1e-5 -g mm -n LZ_H4K20

conda activate RNA-seq

cut -f1-3 M_H4K20_peaks.narrowPeak | sort -k1,1 -k2,2n > M_H4K20.sorted.bed
cut -f1-3 LZ_H4K20_peaks.narrowPeak | sort -k1,1 -k2,2n > LZ_H4K20.sorted.bed

cat M_H4K20.sorted.bed LZ_H4K20.sorted.bed | sort -k1,1 -k2,2n | bedtools merge -i - > H4K20_union_peaks.bed

bedtools multiinter -header -names M_H4K20 LZ_H4K20 -i M_H4K20.sorted.bed LZ_H4K20.sorted.bed > peak_overlap.txt

bedtools multicov -bams ../4_rmdup/rawbam/M_H4K20.map.q30.F12.rmdup.sorted.bam ../4_rmdup/rawbam/M_H4K20_rep.map.q30.F12.rmdup.sorted.bam ../4_rmdup/rawbam/LZ_H4K20.map.q30.F12.rmdup.sorted.bam ../4_rmdup/rawbam/LZ_H4K20_rep.map.q30.F12.rmdup.sorted.bam -bed H4K20_union_peaks.bed > H4K20_peak_counts.tsv

wc -l M_H4K20_summits.bed
wc -l LZ_H4K20_summits.bed
wc -l M_H4K20_peaks.narrowPeak
wc -l LZ_H4K20_peaks.narrowPeak
cat M_H4K20_summits.bed LZ_H4K20_summits.bed > H4K20_summits.bed
cat M_H4K20_peaks.narrowPeak LZ_H4K20_peaks.narrowPeak > H4K20_peaks.narrowPeak

Rscript --vanilla << 'EOF'
library(openxlsx)
library(DESeq2)

H4K20_cts <- read.table("rawdata/H4K20_peak_counts.tsv", header=FALSE)
rownames(H4K20_cts) <- paste(H4K20_cts$V1, H4K20_cts$V2, H4K20_cts$V3, sep="_")
H4K20_countData <- H4K20_cts[,4:ncol(H4K20_cts)]

H4K20_coldata <- data.frame(
  row.names = colnames(H4K20_countData),
  condition = c("H4K20_M","H4K20_M","H4K20_LZ","H4K20_LZ")
)
H4K20_dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(H4K20_countData)),
  colData = H4K20_coldata,
  design = ~ condition
)
H4K20_dds <- DESeq(H4K20_dds)
H4K20_res <- results(H4K20_dds, contrast=c("condition","H4K20_LZ","H4K20_M"))
H4K20_res <- as.data.frame(H4K20_res)

H4K20_res_new <- H4K20_res %>%
  tibble::rownames_to_column("region") %>%
   extract(region, 
          into = c("chr", "start", "end"),
          regex = "^(.*)_(\\d+)_(\\d+)$",
          remove = TRUE,
          convert = TRUE)

rownames(H4K20_res_new)=rownames(H4K20_res)

write.xlsx(H4K20_res_new,'H4K20_res.xlsx',row.names=TRUE)


H4K20_summits=read.table("rawdata/H4K20_summits.bed", header=FALSE)
H4K20_summits=H4K20_summits[,c(1:4)]
colnames(H4K20_summits)=c('SummitChr','SummitStart','SummitEnd','PeakNames')
head(H4K20_summits)
H4K20_narrowPeak=read.table("rawdata/H4K20_peaks.narrowPeak", header=FALSE)
H4K20_narrowPeak=H4K20_narrowPeak[,c(1:4)]
colnames(H4K20_narrowPeak)=c('narrowPeakChr','narrowPeakStart','narrowPeakEnd','PeakNames')
head(H4K20_narrowPeak)

library(dplyr)

H4K20_narrowPeak_unique <- H4K20_narrowPeak %>%
  distinct(narrowPeakChr, narrowPeakStart, narrowPeakEnd, .keep_all = TRUE)
rownames(H4K20_narrowPeak_unique) <- paste(
  H4K20_narrowPeak_unique$narrowPeakChr, 
  H4K20_narrowPeak_unique$narrowPeakStart, 
  H4K20_narrowPeak_unique$narrowPeakEnd, 
  sep = "_"
)
H4K20_narrowPeak_unique$rownames=rownames(H4K20_narrowPeak_unique)
head(H4K20_narrowPeak_unique)

H4K20_deseq.bedup <- H4K20_res_new[which(H4K20_res_new$padj < 0.05 & H4K20_res_new$log2FoldChange > 0), c("chr", "start", "end")]

H4K20_deseq.beddown <- H4K20_res_new[which(H4K20_res_new$padj < 0.05 & H4K20_res_new$log2FoldChange < 0), c("chr", "start", "end")]

H4K20_deseq.bedcommon <- H4K20_res_new[which(H4K20_res_new$padj > 0.05), c("chr", "start", "end")]
H4K20_remaining <- !(H4K20_res_new$padj < 0.05 & H4K20_res_new$log2FoldChange > 0 | H4K20_res_new$padj < 0.05 & H4K20_res_new$log2FoldChange < 0)

H4K20_deseq.bedcommon <- H4K20_res_new[H4K20_remaining,  c("chr", "start", "end")]
write.table(H4K20_deseq.bedup, file="H4K20_deseq.bedup.bed", sep="\t",  quote=F, row.names=F, col.names=F)
write.table(H4K20_deseq.beddown, file="H4K20_deseq.beddown.bed", sep="\t",  quote=F, row.names=F, col.names=F)
write.table(H4K20_deseq.bedcommon, file="H4K20_deseq.bedcommon.bed", sep="\t",  quote=F, row.names=F, col.names=F)

 
# up
length(rownames(H4K20_deseq.bedup))
H4K20_deseq.bedup$rownames=rownames(H4K20_deseq.bedup)
head(H4K20_deseq.bedup)
H4K20_narrowPeak_up=merge(H4K20_narrowPeak_unique,H4K20_deseq.bedup,by='rownames')
length(rownames(H4K20_narrowPeak_up))
length(rownames(H4K20_deseq.bedup))
H4K20_narrowPeak_up=H4K20_narrowPeak_up[,c(2:5)]
H4K20_summits_up=merge(H4K20_narrowPeak_up,H4K20_summits,by='PeakNames')
H4K20_summits_up_final=H4K20_summits_up[,c('SummitChr','SummitStart','SummitEnd')]
write.table(H4K20_summits_up_final, file="H4K20_summits_up_LZvsM.bed", sep="\t",  quote=F, row.names=F, col.names=F)
H4K20_narrowPeak_up_final=H4K20_summits_up[,c('narrowPeakChr','narrowPeakStart','narrowPeakEnd')]
write.table(H4K20_narrowPeak_up_final, file="H4K20_narrowPeak_up_LZvsM.bed", sep="\t",  quote=F, row.names=F, col.names=F)

# down
length(rownames(H4K20_deseq.beddown))
H4K20_deseq.beddown$rownames=rownames(H4K20_deseq.beddown)
head(H4K20_deseq.beddown)
H4K20_narrowPeak_down=merge(H4K20_narrowPeak_unique,H4K20_deseq.beddown,by='rownames')
length(rownames(H4K20_narrowPeak_down))
length(rownames(H4K20_deseq.beddown))
H4K20_narrowPeak_down=H4K20_narrowPeak_down[,c(2:5)]
H4K20_summits_down=merge(H4K20_narrowPeak_down,H4K20_summits,by='PeakNames')
H4K20_summits_down_final=H4K20_summits_down[,c('SummitChr','SummitStart','SummitEnd')]
write.table(H4K20_summits_down_final, file="H4K20_summits_down_LZvsM.bed", sep="\t",  quote=F, row.names=F, col.names=F)
H4K20_narrowPeak_down_final=H4K20_summits_down[,c('narrowPeakChr','narrowPeakStart','narrowPeakEnd')]
write.table(H4K20_narrowPeak_down_final, file="H4K20_narrowPeak_down_LZvsM.bed", sep="\t",  quote=F, row.names=F, col.names=F)

# common
length(rownames(H4K20_deseq.bedcommon))
H4K20_deseq.bedcommon$rownames=rownames(H4K20_deseq.bedcommon)
head(H4K20_deseq.bedcommon)
H4K20_narrowPeak_common=merge(H4K20_narrowPeak_unique,H4K20_deseq.bedcommon,by='rownames')
length(rownames(H4K20_narrowPeak_common))
length(rownames(H4K20_deseq.bedcommon))
H4K20_narrowPeak_common=H4K20_narrowPeak_common[,c(2:5)]
H4K20_summits_common=merge(H4K20_narrowPeak_common,H4K20_summits,by='PeakNames')
H4K20_summits_common_final=H4K20_summits_common[,c('SummitChr','SummitStart','SummitEnd')]
write.table(H4K20_summits_common_final, file="H4K20_summits_common_LZvsM.bed", sep="\t",  quote=F, row.names=F, col.names=F)
H4K20_narrowPeak_common_final=H4K20_summits_common[,c('narrowPeakChr','narrowPeakStart','narrowPeakEnd')]
write.table(H4K20_narrowPeak_common_final, file="H4K20_narrowPeak_common_LZvsM.bed", sep="\t",  quote=F, row.names=F, col.names=F)
save.image('H4K20_deg.RData')
EOF

computeMatrix reference-point  -p 56 -S M_H4K20.bw LZ_H4K20.bw -R H4K20_summits_down_LZvsM.bed H4K20_summits_common_LZvsM.bed H4K20_summits_up_LZvsM.bed -b 3000 -a 3000 -o H4K20.matrix.gz
plotHeatmap -m H4K20.matrix.gz -out H4K20.peak_heatmap.pdf --regionsLabel "Model" "Common" "LZ" --samplesLabel "H4K20_Model" "H4K20_LZ"

##################  get NKT bed ##################
Rscript --vanilla << 'EOF'
colors <- c("#3C5488",'#DC0000','')
my_genes <- c("Rag1", "Rag2",'Tcf7','Myc','Cd8a',
              "Cd1d1",'Cd1d2','B2m','Cstd','Lamp1',
              'Egr1',"Egr2",'Fos','Jun','Sox4',
              "Cd69",'Gata3','Id2','Klf2','Maf' )
library(rtracklayer)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
genes_gr <- genes(txdb)
genes_df <- as.data.frame(genes_gr)
library(clusterProfiler)
library(org.Mm.eg.db)
genes_df$symbol <- mapIds(org.Mm.eg.db,
                          keys = genes_df$gene_id,
                          column = "SYMBOL",
                          keytype = "ENTREZID",
                          multiVals = "first")

target_genes <- genes_df[genes_df$symbol %in% my_genes, ]
target_bed <- target_genes[, c("seqnames", "start", "end")]
write.table(target_bed, "../CUTTag/nkt_genes.bed", row.names=F, col.names=F, sep="\t",quote=F)
EOF

################## NKT TSS ±3000 ##################deeptools
computeMatrix reference-point  -p 56 -S LZ_PHF8.bw M_PHF8.bw -R nkt_genes.bed -b 3000 -a 3000 -o PHF8.matrix.gz
plotProfile -m PHF8.matrix.gz -out PHF8_TSS_profile.pdf --samplesLabel "LZ_PHF8" "M_PHF8" --perGroup

computeMatrix reference-point  -p 56 -S LZ_H3K27.bw M_H3K27.bw -R nkt_genes.bed -b 3000 -a 3000 -o H3K27.matrix.gz
plotProfile -m H3K27.matrix.gz -out H3K27_TSS_profile.pdf --samplesLabel "LZ_H3K27" "M_H3K27" --perGroup

computeMatrix reference-point  -p 56 -S LZ_H3K9.bw M_H3K9.bw -R nkt_genes.bed -b 3000 -a 3000 -o H3K9.matrix.gz
plotProfile -m H3K9.matrix.gz -out H3K9_TSS_profile.pdf --samplesLabel "LZ_H3K9" "M_H3K9" --perGroup

computeMatrix reference-point  -p 56 -S LZ_H4K20.bw M_H4K20.bw -R nkt_genes.bed -b 3000 -a 3000 -o H4K20.matrix.gz
plotProfile -m H4K20.matrix.gz -out H4K20_TSS_profile.pdf --samplesLabel "LZ_H4K20" "M_H4K20" --perGroup


####CUTTAG####
Rscript --vanilla << 'EOF'
colors <- c("#3C5488",'#DC0000')
# 
my_genes <- c("Rag1", "Rag2",'Tcf7','Myc','Cd8a',
              "Cd1d1",'Cd1d2','B2m','Cstd','Lamp1',
              'Egr1',"Egr2",'Fos','Jun','Sox4',
              "Cd69",'Gata3','Id2','Klf2','Maf' )


library(rtracklayer)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene


genes_gr <- genes(txdb)
genes_df <- as.data.frame(genes_gr)
library(clusterProfiler)
library(org.Mm.eg.db)

genes_df$symbol <- mapIds(org.Mm.eg.db,
                          keys = genes_df$gene_id,
                          column = "SYMBOL",
                          keytype = "ENTREZID",
                          multiVals = "first")

target_genes <- genes_df[genes_df$symbol %in% my_genes, ]
target_bed <- target_genes[, c("seqnames", "start", "end")]
write.table(target_bed, "../CUTTag/nkt_genes.bed", row.names=F, col.names=F, sep="\t", quote=F)
EOF

for bam in *.map.q30.F12.rmdup.sorted.bam; do
    echo "Processing $bam"
    samtools stats "$bam" > "${bam%.bam}.stats"
    grep ^IS "${bam%.bam}.stats" > "${bam%.bam}.insert_sizes.tsv"
done


library(ggplot2)
library(patchwork)  
library(dplyr)

tsv_files <- list.files(pattern = "\\.tsv$")

for (file in tsv_files) {
  sample_name <- gsub("\\.tsv$", "", file)

  data <- read.table(file, header = FALSE, row.names = 2, colClasses = c("NULL", NA, NA, NA, NA, NA))
  colnames(data) <- c("Total", "Inward", "Outward", "Other")
  
  insert_counts <- data[, "Total"]

  max_len <- min(600, length(insert_counts))
  insert_counts <- insert_counts[1:max_len]
  
  sample_name_title=gsub("\\.map.q30.F12.rmdup.sorted.insert_sizes$", "", sample_name)
  # PDF
  pdf(paste0(sample_name, "_insert_size.pdf"), width = 8, height = 5)
  barplot(insert_counts, 
          col = "steelblue", 
          border = NA,
          xlab = "Insert Size (bp)", 
          ylab = "Frequency",
          main = paste0("CUT&Tag Fragment Length Distribution - ", sample_name_title))
  abline(v = 180, col = "red", lty = 2, lwd = 2)
  abline(v = 360, col = "blue", lty = 2, lwd = 2)
  legend("topright", 
         legend = c("Mono-nucleosome (~180bp)", "Di-nucleosome (~360bp)"),
         col = c("red", "blue"), lty = 2, lwd = 2)
  dev.off()

}

for bam in *.map.q30.F12.rmdup.sorted.bam; do
    echo "Processing $bam"
    samtools flagstat "$bam" > "${bam%.bam}.flagstats"
done

library(dplyr)
library(tidyr)
flagstats_files <- list.files(pattern = "\\.flagstats$")

results_list <- list()

for (i in 1:length(flagstats_files)) {
  sample_name <- gsub("\\.map.q30.F12.rmdup.sorted.flagstats$", "", flagstats_files[i])
  
  # flagstats
  flag_lines <- readLines(flagstats_files[i])
  total_mapped_line <- flag_lines[7]
  total_mapped <- as.numeric(gsub(" .*", "", total_mapped_line))
  peaks_file <- paste0(sample_name, '.reads_in_peaks.txt')
  if (file.exists(peaks_file)) {
    peak_sum <- as.numeric(readLines(peaks_file, n = 1))
    if (!is.na(peak_sum)) {
      print(paste("peak_sum:", peak_sum))
    } else {
      peak_sum <- NA
    }
  } else {
    peak_sum <- NA
  }

  results_list[[i]] <- data.frame(
    sample = sample_name,
    total_mapped = total_mapped,
    peaks_sum = peak_sum,
    FRiP = peak_sum*100 / total_mapped
  )
}

final_result <- bind_rows(results_list)
print(final_result)

write.csv(final_result, "peak_quantification_summary.csv", row.names = FALSE)

for bam in ../4_rmdup/rawbam/*.map.q30.F12.rmdup.sorted.bam; do
    echo "Processing $bam"
    sample=${bam%.map.q30.F12.rmdup.sorted.bam}
    sample=$(basename "$sample")
    bedtools intersect -a /data/CUTTag/LZandM/4_rmdup/rawbam/$sample.map.q30.F12.rmdup.sorted.bam -b /data/CUTTag/LZandM/6_callpeak/$sample.peak_peaks.narrowPeak -bed -c | awk '{sum+=$NF} END {print sum}' > $sample.reads_in_peaks.txt
done

final_result=read.csv("peak_quantification_summary.csv")
final_result$sample=factor(final_result$sample,levels=c("LZ_PHF8_1", "LZ_PHF8_2", "LZ_H3K27_1", "LZ_H3K27_2", "LZ_H3K9_1", "LZ_H3K9_2" ,"LZ_H4K20_1","LZ_H4K20_2","M_PHF8_1", "M_PHF8_2", "M_H3K27_1", "M_H3K27_2", "M_H3K9_1", "M_H3K9_2" ,"M_H4K20_1","M_H4K20_2"  ))
p1=ggplot(final_result, aes(x = sample, y = FRiP, fill = sample))+
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(aes(label = round(FRiP, 2)), vjust = -0.3, position = position_dodge(0.9), size = 3) +
  theme_classic()+theme(axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5,color='black'),axis.text.y = element_text(color='black'),legend.position = "none")+
  scale_fill_manual(values = sci_colors)+labs(x='',y='FRiP')
ggsave(filename = "p_FRiP.pdf", plot = p1, height = 4, width = 6)









