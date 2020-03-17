options(stringsAsFactors = FALSE);
library(DESeq2);
library(ggplot2)



types = c('ER+', 'TNBC')

depths <- c('0.01','0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.85', '0.9', '0.95', '0.99', '1')
N_boots <- 1:25

ER_list <- list()
TNBC_list <- list()





for (type in types){

	depth_list <- list()

	for (i in 1:length(depths)){

		boot_list <- list()

		for (j in 1:length(N_boots)){

			metadata_file_name = paste('Metadata/', type, '_metadata_', depths[i], '_', as.character(N_boots[j]), '.txt', sep = '')
			metadata = read.delim(file = metadata_file_name, header = T)
			cds = DESeqDataSetFromHTSeqCount(metadata, "CV",design = ~ Condition)
			print(metadata_file_name)

			cds = estimateSizeFactors(cds);
			cds = estimateDispersions(cds);
			cds = nbinomWaldTest(cds);

			de = data.frame(results(cds));
			de = de[!is.na(de$padj), ];

			boot_list[[j]] <- de

		}

		depth_list[[i]] <- boot_list

	}

	if (type == 'ER+'){
		ER_list <- depth_list
	}
	if (type == 'TNBC'){
		TNBC_list <- depth_list
	}

}

save.image('subDEseq2.RData')

#############################

load('subDEseq2.RData')


pval = 0.05
fold = 'none'

de_ER_deseq <- list()
de_TNBC_deseq <- list()

for (i in 1:length(depths)){

	ER_boot_list <- list()
	TNBC_boot_list <- list()

	for (j in 1:length(N_boots)){

		de_ER <- ER_list[[i]][[j]]
		de_TNBC <- TNBC_list[[i]][[j]]

		if (fold == 'none'){
			ER_boot_list[[j]] <- rownames(de_ER[(de_ER$padj < pval),])
			TNBC_boot_list[[j]] <- rownames(de_TNBC[(de_TNBC$padj < pval),])
		}
		else {
			ER_boot_list[[j]] <- rownames(de_ER[(de_ER$padj < pval) & (de_ER$log2FoldChange > fold),])
			TNBC_boot_list[[j]] <- rownames(de_TNBC[(de_TNBC$padj < pval) & (de_TNBC$log2FoldChange > fold),])
		}

	}

	de_ER_deseq[[i]] <- ER_boot_list
	de_TNBC_deseq[[i]] <- TNBC_boot_list

}

save('de_ER_deseq', file = paste('sub_de_ER_deseq_',as.character(pval), '_', as.character(fold), sep = ''))
save('de_TNBC_deseq', file = paste('sub_de_TNBC_deseq_',as.character(pval), '_', as.character(fold), sep = ''))


########################################################

library('org.Hs.eg.db')
library('goseq');
library('GO.db');

load('DEseq2.RData')

x <- rownames(ER_list[[1]][[1]][(ER_list[[1]][[1]]$padj < 0.05)&(ER_list[[1]][[1]]$log2FoldChange > 2),])
xx <- rownames(ER_list[[1]][[1]])
xxx <- as.list(org.Hs.egALIAS2EG)

get_entrez <- function(gene_v, entr_list){

	res_v <- vector()
	counter = 1

	for (i in 1: length(gene_v)){
		if ((gene_v[i] %in% names(entr_list))== FALSE){
			#print(gene_v[i])
			next
		}

		if (length(entr_list[[gene_v[i]]]) > 1){
			next
		}

		res_v[counter] <- entr_list[[gene_v[i]]]
		counter = counter + 1
	}
	return(res_v)
}


q <- get_entrez(x, xxx)
q <- unique(q)

qq <- get_entrez(xx, xxx)
qq <- unique(qq)


qqq <- matrix(ncol = 2, nrow = length(qq))
qqq <- data.frame(qqq)
qqq[,1] <- qq

for(i in 1:length(qq)){

	if ((qq[i] %in% q) == TRUE){
		qqq[i,2] <- 1
	}
	else{
		qqq[i,2] <- 0
	}

}


names(qqq) <- c('entrezgene','diff.exprs')
gene.set <- qqq
goseq.set = as.vector(gene.set$diff.exprs);
names(goseq.set) = gene.set$entrezgene;
pwf = nullp(goseq.set, "hg19","refGene");
GO.wall = goseq(pwf,"hg19","refGene");

enriched.GO =GO.wall$category[ p.adjust(GO.wall$over_represented_pvalue, method="bonferroni")<0.01] ;
print(length(enriched.GO));



for(go in enriched.GO) { 
	print(GOTERM[[go]]); cat("------------------------------------\n");
	}






















