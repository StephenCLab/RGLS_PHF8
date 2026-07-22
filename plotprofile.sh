computeMatrix reference-point  -p 56 -S LZ_PHF8.bw M_PHF8.bw -R nkt_genes.bed -b 3000 -a 3000 -o PHF8.matrix.gz
plotProfile -m PHF8.matrix.gz -out PHF8_TSS_profile.pdf --samplesLabel "LZ_PHF8" "M_PHF8" --perGroup


computeMatrix reference-point  -p 56 -S LZ_H3K27.bw M_H3K27.bw -R nkt_genes.bed -b 3000 -a 3000 -o H3K27.matrix.gz
plotProfile -m H3K27.matrix.gz -out H3K27_TSS_profile.pdf --samplesLabel "LZ_H3K27" "M_H3K27" --perGroup



computeMatrix reference-point  -p 56 -S LZ_H3K9.bw M_H3K9.bw -R nkt_genes.bed -b 3000 -a 3000 -o H3K9.matrix.gz
plotProfile -m H3K9.matrix.gz -out H3K9_TSS_profile.pdf --samplesLabel "LZ_H3K9" "M_H3K9" --perGroup




computeMatrix reference-point  -p 56 -S LZ_H4K20.bw M_H4K20.bw -R nkt_genes.bed -b 3000 -a 3000 -o H4K20.matrix.gz
plotProfile -m H4K20.matrix.gz -out H4K20_TSS_profile.pdf --samplesLabel "LZ_H4K20" "M_H4K20" --perGroup