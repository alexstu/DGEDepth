#options(stringsAsFactors = FALSE);
library(DESeq2);
library(ggplot2)
library(doParallel)

cv_dir = 'CV/'

types = c('ER+', 'TNBC')
depths <- c('0.01','0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.85', '0.9', '0.95', '0.99', '1')
N_boots <- 1:25
ER_samples <- c('All','SRR1313092.sam','SRR1313100.sam','SRR1313105.sam','SRR1313107.sam','SRR1313124.bam','SRR1313174.bam','SRR1313175.bam','SRR1313178.bam','SRR1313197.bam','SRR1313203.bam')
TNBC_samples <- c('All', 'SRR1313141.bam','SRR1313142.bam','SRR1313147.bam','SRR1313152.bam','SRR1313164.bam','SRR1313217.bam','SRR1313220.bam','SRR1313226.bam','SRR1313227.bam', 'SRR1313229.bam')

#types = c('ER+', 'TNBC')
#depths <- c('0.01','0.05')
#N_boots <- 1:10
#ER_samples <- c('','SRR1313092.sam','SRR1313100.sam')
#TNBC_samples <- c('', 'SRR1313141.bam','SRR1313142.bam')


ER_list <- list()
TNBC_list <- list()


for (type in types){

	if (type == 'ER+'){
		samples <- ER_samples
	}
	if (type == 'TNBC'){
		samples <- TNBC_samples
	}


	for (s in 1:length(samples)){

		depth_list <- list()

		for (i in 1:length(depths)){

			boot_list <- list()

			cl <- makeCluster(13)
			registerDoParallel(cl)

			print(paste(type, samples[s], depths[i]), sep = '_')
			#for (j in 1:length(N_boots)){
			boot_list <- foreach(j = 1:length(N_boots)) %dopar% {

				library(DESeq2);
				filename <- paste(cv_dir, type, '_', as.character(depths[i]), '_', as.character(j), '.txt', sep = '')
				x <- read.table(filename, header = TRUE, row.names = 'GeneID')


				if (s == 1){
					x <- x
					if (type == 'ER+'){
						colDat = data.frame(condition = factor(c(rep(1,42), rep(2,30))))
					}
					if (type == 'TNBC'){
						colDat = data.frame(condition = factor(c(rep(1,42), rep(2,21))))
					}
				}
				if ((s > 1) & (s <= 6)){
					x <- x[, colnames(x) != samples[s]]
					if (type == 'ER+'){
						colDat = data.frame(condition = factor(c(rep(1,41), rep(2,30))))
					}
					if (type == 'TNBC'){
						colDat = data.frame(condition = factor(c(rep(1,41), rep(2,21))))
					}
				}
				if (s > 6){
					x <- x[, colnames(x) != samples[s]]
					if (type == 'ER+'){
						colDat = data.frame(condition = factor(c(rep(1,42), rep(2,29))))
					}
					if (type == 'TNBC'){
						colDat = data.frame(condition = factor(c(rep(1,42), rep(2,20))))
					}
				}


				dds<-DESeqDataSetFromMatrix(countData=x,colDat,formula(~condition))
				dds <- estimateSizeFactors(dds)
				dds <- estimateDispersions(dds)
				cds <- nbinomWaldTest(dds)

				de = data.frame(results(cds));
				de = de[!is.na(de$padj), ];

				de
			#boot_list[[j]] <- de

			}

			stopCluster(cl)

			depth_list[[i]] <- boot_list

		}

		if (type == 'ER+'){
			ER_list[[s]] <- depth_list
		}
		if (type == 'TNBC'){
			TNBC_list[[s]] <- depth_list
		}

	}
}


save.image('JKDEseq2.RData')

#############################

load('JKDEseq2.RData')


pval = 0.05
fold = 2
#fold = 'none'

#de_ER_deseq <- list(rep(list(rep(list(),length(depths))),length(samples)))
#de_TNBC_deseq <- list(rep(list(rep(list(),length(depths))),length(samples)))

de_ER_deseq <- list()
de_TNBC_deseq <- list()


for (s in 1:length(samples)){

	ER_sample_list <- list()
	TNBC_sample_list <- list()

	for (i in 1:length(depths)){

		ER_boot_list <- list()
		TNBC_boot_list <- list()

		for (j in 1:length(N_boots)){

			de_ER <- ER_list[[s]][[i]][[j]]
			de_TNBC <- TNBC_list[[s]][[i]][[j]]

			if (fold == 'none'){
				ER_boot_list[[j]] <- rownames(de_ER[(de_ER$padj < pval),])
				TNBC_boot_list[[j]] <- rownames(de_TNBC[(de_TNBC$padj < pval),])
			}
			else {
				ER_boot_list[[j]] <- rownames(de_ER[(de_ER$padj < pval) & (de_ER$log2FoldChange > fold),])
				TNBC_boot_list[[j]] <- rownames(de_TNBC[(de_TNBC$padj < pval) & (de_TNBC$log2FoldChange > fold),])
			}

		}

		ER_sample_list[[i]] <- ER_boot_list
		TNBC_sample_list[[i]] <- TNBC_boot_list

	}
	de_ER_deseq[[s]] <- ER_sample_list
	de_TNBC_deseq[[s]] <- TNBC_sample_list

}

save('de_ER_deseq', file = paste('JK_de_ER_deseq_',as.character(pval), '_', as.character(fold), sep = ''))
save('de_TNBC_deseq', file = paste('JK_de_TNBC_deseq_',as.character(pval), '_', as.character(fold), sep = ''))


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






















