options(stringsAsFactors = FALSE);
library(DESeq2);

#type = 'new_gene_all'
#type = 'old_gene_all'
type = 'new_exon_all'

depths <- c('0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.4', '0.5', '0.6', '0.7', '0.85', '0.9', '0.95', '0.99', '1')
N_boots <- 1:25

M <- matrix(nrow = length(depths), ncol = length(N_boots))

for (i in 1:length(depths)){

	for (j in 1:length(N_boots)){

		metadata_file_name = paste(type, '_metadata_', depths[i], '_', as.character(N_boots[j]), '.txt', sep = '')
		metadata = read.delim(file = metadata_file_name, header = T)
		cds = DESeqDataSetFromHTSeqCount(metadata, "CV",design = ~ Condition)
		print(metadata_file_name)

		cds = estimateSizeFactors(cds);
		cds = estimateDispersions(cds);
		cds = nbinomWaldTest(cds);

		de = data.frame(results(cds));
		de = de[!is.na(de$padj), ];
		gene.set = rownames(de);

		significant_genes = de[(de$padj < 0.05),];
		N_de <- length(significant_genes[,1])
		M[i,j] <- N_de

	}

}

#M_nga <- M
#M_oga <- M
M_nea <- M

#boxplot(M[] ~ depths)

png(paste("01_DE_expressed_genes_", type, ".png", sep = ''), width = 800, height = 800)
par(mar = c(8, 10, 7, 2), mgp = c(5, 1, 0))


#header = 'DE genes for exon new annotation '
#boxplot(M_nea[] ~ depths, main = header, font.main = 2, cex.main = 3, xlab = expression(bold(Depth)), ylab = expression(bold(Number~of~DE~genes)), font.lab = 2, cex.axis = 2, cex.lab = 3, ylim = c(0, 200), medcol = 'red', medlwd=6)

#header = 'DE genes for gene new annotation '
#boxplot(M_nga[] ~ depths, main = header, font.main = 2, cex.main = 3, xlab = expression(bold(Depth)), ylab = expression(bold(Number~of~DE~genes)), font.lab = 2, cex.axis = 2, cex.lab = 3, ylim = c(0, 200), medcol = 'red', medlwd=6)

header = 'DE genes for gene old annotation '
boxplot(M_oga[] ~ depths, main = header, font.main = 2, cex.main = 3, xlab = expression(bold(Depth)), ylab = expression(bold(Number~of~DE~genes)), font.lab = 2, cex.axis = 2, cex.lab = 3, ylim = c(0, 200), medcol = 'red', medlwd=6)


dev.off()

##########################
#significant.genes = de[(de$padj < .05 & abs(de$log2FoldChange) >= 2.0),];
#nrow(significant.genes);
#significant.genes = significant.genes[order(significant.genes$padj, - abs(significant.genes$log2FoldChange), decreasing=F), ];
#write.table(significant.genes, file="deseq/significant.genes.txt", quote=F, row.names=F, col.names=T, sep = "\t");
###############################















