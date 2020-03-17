library('edgeR')
library(doParallel)


cv_dir = 'CV/'

types = c('ER+', 'TNBC')
depths <- c('0.01','0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.85', '0.9', '0.95', '0.99', '1')
N_boots <- 1:25
ER_samples <- c('All','SRR1313092.sam','SRR1313100.sam','SRR1313105.sam','SRR1313107.sam','SRR1313124.bam','SRR1313174.bam','SRR1313175.bam','SRR1313178.bam','SRR1313197.bam','SRR1313203.bam')
TNBC_samples <- c('All', 'SRR1313141.bam','SRR1313142.bam','SRR1313147.bam','SRR1313152.bam','SRR1313164.bam','SRR1313217.bam','SRR1313220.bam','SRR1313226.bam','SRR1313227.bam', 'SRR1313229.bam')


##types = c('ER+', 'TNBC')
#depths <- c('0.01','0.05')
#N_boots <- 1:2
#ER_samples <- c('All','SRR1313092.sam','SRR1313100.sam')
#TNBC_samples <- c('All', 'SRR1313141.bam','SRR1313142.bam')

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
		sub_names <- samples[(samples != 'All') & (samples != samples[s])]

		for (i in 1:length(depths)){

			boot_list <- list()

			cl <- makeCluster(13)
			registerDoParallel(cl)

			print(paste(type, samples[s], depths[i]), sep = '_')
			#for (j in 1:length(N_boots)){
			boot_list <- foreach(j = 1:length(N_boots)) %dopar% {

				library(edgeR);
				filename <- paste(cv_dir, type, '_', as.character(depths[i]), '_', as.character(j), '.txt', sep = '')
				x <- read.table(filename, header = TRUE, row.names = 'GeneID')
				x <- x[,sub_names]

				if (s == 1){
					if (type == 'ER+'){
						group = factor(c(rep(1,5), rep(2,5)))
					}
					if (type == 'TNBC'){
						group = factor(c(rep(1,5), rep(2,5)))
					}
				}
				if ((s > 1) & (s <= 6)){
					if (type == 'ER+'){
						group = factor(c(rep(1,4), rep(2,5)))
					}
					if (type == 'TNBC'){
						group = factor(c(rep(1,4), rep(2,5)))
					}
				}
				if (s > 6){
					if (type == 'ER+'){
						group = factor(c(rep(1,5), rep(2,4)))
					}
					if (type == 'TNBC'){
						group = factor(c(rep(1,5), rep(2,4)))
					}
				}


			y <- DGEList(counts=x,group=group)
			y <- calcNormFactors(y)
			y <- estimateCommonDisp(y)
			y <- estimateTagwiseDisp(y)
			et <- exactTest(y)
			cor <- topTags(et, n = dim(et)[1], sort.by='none')
			cor
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


save.image('sub_JKedgeR.RData')
quit(save='no')


#############################

load('sub_JKedgeR.RData')


pval = 0.05
fold = 'none'

de_ER_edger <- list()
de_TNBC_edger <- list()

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
				ER_boot_list[[j]] <- rownames(de_ER[(de_ER$table$FDR < pval),])
				TNBC_boot_list[[j]] <- rownames(de_TNBC[(de_TNBC$table$FDR < pval),])
			}
			else {
				ER_boot_list[[j]] <- rownames(de_ER[(de_ER$table$FDR < pval) & (de_ER$table$logFC > fold),])
				TNBC_boot_list[[j]] <- rownames(de_TNBC[(de_TNBC$table$FDR < pval) & (de_TNBC$table$logFC > fold),])
			}

		}

		ER_sample_list[[i]] <- ER_boot_list
		TNBC_sample_list[[i]] <- TNBC_boot_list

	}
	de_ER_edger[[s]] <- ER_sample_list
	de_TNBC_edger[[s]] <- TNBC_sample_list

}
save('de_ER_edger', file = paste('sub_JK_de_ER_edger_',as.character(pval), '_', as.character(fold), sep = ''))
save('de_TNBC_edger', file = paste('sub_JK_de_TNBC_edger_',as.character(pval), '_', as.character(fold), sep = ''))


########################################################
###########################################################

library(samExploreR)

load('edgeR.RData')


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
		df[i,2] <- as.numeric(unlist(TNBC_list[[i]][[j]][g, 'FDR'])[1])
		}
	#ggplot(df, aes(x=d, y=pv)) + geom_point(aes(colour='blue',size=1, guide = FALSE)) + geom_hline(yintercept=0.05, linetype="dotted", size = 0.7)+ theme(legend.position="none") + geom_line(aes( colour = 'red' ),se = F, size = 0.7)
	ggplot(df, aes(x=f, y=p_value)) + geom_point(aes(colour='blue',size=2, guide = FALSE)) + geom_hline(yintercept=0.05, linetype="dotted", size = 0.7)+ theme(legend.position="none") + geom_line(aes( colour = 'red'),se = F, size = 1.3)+ scale_y_continuous(limits = c(0, 1)) + theme(text = element_text(size=25, face="bold")) + ggtitle(paste('Edger for ', g, ' gene', sep = ''))
}



#pdf("PV3.pdf")
#xxxx('PIK3CD',1)
#xxxx('RMDN2',1)
#xxxx('GRR75-ASB3',1)
#dev.off()



j=1
pdf(paste('RMDN2-edger',as.character(j),'.pdf', sep = ''))
xxxx('RMDN2',j)
dev.off()

j=1
pdf(paste('PIK3CD-edger',as.character(j),'.pdf', sep = ''))
xxxx('PIK3CD',j)
dev.off()

####################################

load('edgeR.RData')



std <- function(x) sd(x)/sqrt(length(x))

xxxxx <- function(g){

	depths <- c('0.01','0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.85', '0.9', '0.95', '0.99', '1')
	N_boots <- 1:25
	res_table = matrix(ncol = 3, nrow = length(depths))

	for (i in 1:length(depths)){

		dd = as.numeric(depths[i])
		cur_vector_pval = vector()
		for (j in N_boots){
			cur_vector_pval <- c(cur_vector_pval, as.numeric(unlist(TNBC_list[[i]][[j]][g, 'FDR'])[1]))
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


pdf(paste('PIK3CD-edger','all','.pdf', sep = ''))
xxxxx('PIK3CD')
dev.off()


pdf(paste('RMDN2-edger','all','.pdf', sep = ''))
xxxxx('RMDN2')
dev.off()





