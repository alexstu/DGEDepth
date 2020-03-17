options(stringsAsFactors = FALSE);
library(DESeq2);
library(ggplot2)
library(VennDiagram)
library(gridExtra)



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

save.image('DEseq2.RData')

#############################

load('DEseq2.RData')


pval = 0.05
fold = 2

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

save('de_ER_deseq', file = paste('de_ER_deseq_',as.character(pval), '_', as.character(fold), sep = ''))
save('de_TNBC_deseq', file = paste('de_TNBC_deseq_',as.character(pval), '_', as.character(fold), sep = ''))


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



###########################################################

library(samExploreR)

load('DEseq2.RData')


#genes = rownames(ER_list[[17]][[1]][with(ER_list[[17]][[1]], order(log2FoldChange,-padj)), ])[1:500]

#ER_list[[17]][[1]][with(ER_list[[17]][[1]], order(log2FoldChange,-padj)), ]

for (j in 1:length(genes[200:250])){
	#j=68
	df <- matrix(ncol = 3, nrow = 25*17*5)
	colnames(df) = c('Label', 'Variable', 'Value')
	df <- as.data.frame(df)

	for (i in 1:length(depths)){

		df[j*((i-1)*length(N_boots)+(j-1))+1,1] <- paste(as.character(j), genes[j], sep = ', ')
		df[j*((i-1)*length(N_boots)+(j-1))+1,2] <- as.numeric(depths[i])
		df[j*((i-1)*length(N_boots)+(j-1))+1,3] <- TNBC_list[[i]][[1]][genes[j], 'padj']
		}

		plotsamExplorer(df, save = TRUE, p.depth = 0.9, anova = FALSE)
		system('sleep 4')

}

#for (k in 1:(length(depths)-1)){
	k=1
#	for (j in 1:25){
	j=1
		genes = vector()
		x <- rownames(ER_list[[k]][[j]][(ER_list[[k]][[j]]$padj < 0.05)&(abs(ER_list[[k]][[j]]$log2FoldChange) > 0),])
		y <- rownames(ER_list[[17]][[1]][(ER_list[[17]][[1]]$padj < 0.05)&(abs(ER_list[[17]][[1]]$log2FoldChange) > 0),])

		z <- intersect(x,y)

		for (g in x){

			if ((g %in% z) == TRUE){
				next
			}

			if ((g %in% y) == TRUE){
				next
			}

			genes <- c(genes, g)

		}

#	}
#}
genes

xxxx <- function(g,j){

	pv <- vector()
	df <- matrix(ncol = 2, nrow = length(depths))
	df <- as.data.frame(df)
	colnames(df) <- c('f', 'p_value')

	for (i in 1:length(depths)){

		#df[j*((i-1)*length(N_boots)+(j-1))+1,1] <- paste(as.character(j), genes[j], sep = ', ')
		#df[j*((i-1)*length(N_boots)+(j-1))+1,2] <- as.numeric(depths[i])
		#pv <- c(pv,TNBC_list[[i]][[1]][genes[j], 'padj'])
		df[i,1] <- as.numeric(depths[i])
		df[i,2] <- TNBC_list[[i]][[j]][g, 'padj']
		}
	#ggplot(df, aes(x=d, y=pv)) + geom_point(aes(colour='blue',size=1, guide = FALSE)) + geom_hline(yintercept=0.05, linetype="dotted", size = 0.7)+ theme(legend.position="none") + geom_line(aes( colour = 'red' ),se = F, size = 0.7)
	ggplot(df, aes(x=f, y=p_value)) + geom_point(aes(colour='blue',size=2, guide = FALSE)) + geom_hline(yintercept=0.05, linetype="dotted", size = 0.7)+ theme(legend.position="none") + geom_line(aes( colour = 'red'),se = F, size = 1.3)+ scale_y_continuous(limits = c(0, 1)) + theme(text = element_text(size=25, face="bold")) + ggtitle(paste('DESeq2 for ', g, sep = ''))
}

#xxxx('PIK3CD',1)

j=1
pdf(paste('PIK3CD-deseq',as.character(j),'.pdf', sep = ''))
xxxx('PIK3CD',j)
dev.off()

j=1
pdf(paste('RMDN2-deseq',as.character(j),'.pdf', sep = ''))
xxxx('RMDN2',j)
dev.off()

################################################



load('DEseq2.RData')



std <- function(x) sd(x)/sqrt(length(x))

xxxxx <- function(g){

	depths <- c('0.01','0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.85', '0.9', '0.95', '0.99', '1')
	N_boots <- 1:25
	res_table = matrix(ncol = 3, nrow = length(depths))

	for (i in 1:length(depths)){

		dd = as.numeric(depths[i])
		cur_vector_pval = vector()
		for (j in N_boots){
			cur_vector_pval <- c(cur_vector_pval, TNBC_list[[i]][[j]][g, 'padj'])
		}

		res_table[i,1] <- dd
		res_table[i,2] <- (mean(cur_vector_pval))
		res_table[i,3] <- (std(cur_vector_pval))

	}

print(res_table)
epsilon = 0.02

header = paste('DESeq2 for ', g, sep = '')
#header = paste(group, ':', type, '\n', 'Expression > 1',  sep = '')
up = res_table[,2] + res_table[,3]
low = res_table[,2] - res_table[,3]
#png(paste("04_Mean_expressed_genes_", type, '_',  group, ".png", sep = ""))

par(mar=c(6, 5, 4, 2))
plot(res_table[,1], res_table[,2], main = header, font.main = 2, xlab = "f", ylab = "Mean", font.lab = 2, cex.main = 2, cex.axis = 1.5, cex.lab = 2, col = 'red', pch = 16, cex = 2, ylim = c(min(low), max(up)))
for(i in 1:length(res_table[,2])) {
	segments(res_table[i,1],low[i] , res_table[i,1], up[i], col = 'blue', lwd = 2)
	segments(res_table[i,1]-epsilon, up[i] , res_table[i,1]+epsilon, up[i], col = 'blue', lwd = 2)
	segments(res_table[i,1]-epsilon, low[i] , res_table[i,1]+epsilon, low[i], col = 'blue', lwd = 2)
}
abline(h = 0.05, lty=2, lwd=2)
#dev.off()




	#ggplot(df, aes(x=f, y=p_value)) + geom_point(aes(colour='blue',size=2, guide = FALSE)) + geom_hline(yintercept=0.05, linetype="dotted", size = 0.7)+ theme(legend.position="none") + geom_line(aes( colour = 'red'),se = F, size = 1.3)+ scale_y_continuous(limits = c(0, 1)) + theme(text = element_text(size=25, face="bold")) + ggtitle(paste('DESeq2 for ', g, ' gene', sep = ''))
}


pdf(paste('PIK3CD-deseq','all','.pdf', sep = ''))
xxxxx('PIK3CD')
dev.off()


pdf(paste('RMDN2-deseq','all','.pdf', sep = ''))
xxxxx('RMDN2')
dev.off()










