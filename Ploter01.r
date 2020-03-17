'All N(DEG)~d'
library(ggplot2)
library(samExploreR)



depths <- c('0.01','0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.85', '0.9', '0.95', '0.99', '1')
N_boots <- 1:25
filter = '2'


load(paste('JK_de_ER_deseq_0.05_',filter, sep = ''))
load(paste('JK_de_ER_ebseq_0.05_',filter, sep = ''))
load(paste('JK_de_ER_edger_0.05_',filter, sep = ''))
load(paste('JK_de_ER_noiseq_0.05_',filter, sep = ''))
load(paste('JK_de_ER_voom_0.05_',filter, sep = ''))
load(paste('JK_de_TNBC_deseq_0.05_',filter, sep = ''))
load(paste('JK_de_TNBC_ebseq_0.05_',filter, sep = ''))
load(paste('JK_de_TNBC_edger_0.05_',filter, sep = ''))
load(paste('JK_de_TNBC_noiseq_0.05_',filter, sep = ''))
load(paste('JK_de_TNBC_voom_0.05_',filter, sep = ''))

de_ER_deseq <-de_ER_deseq[[1]]
de_ER_voom <-de_ER_voom[[1]]
de_ER_edger <-de_ER_edger[[1]]
de_ER_ebseq <-de_ER_ebseq[[1]]
de_ER_noiseq <-de_ER_noiseq[[1]]
de_TNBC_deseq <-de_TNBC_deseq[[1]]
de_TNBC_voom <-de_TNBC_voom[[1]]
de_TNBC_edger <-de_TNBC_edger[[1]]
de_TNBC_ebseq <-de_TNBC_ebseq[[1]]
de_TNBC_noiseq <-de_TNBC_noiseq[[1]]


if (filter == 'none'){
	xxxx = 'no '
}
if (filter == '2'){
	xxxx = '2 '
}


#############################

ER_df <- matrix(ncol = 3, nrow = 25*17*5)
TNBC_df <- matrix(ncol = 3, nrow = 25*17*5)
colnames(TNBC_df) = c('Label', 'Variable', 'Value')
colnames(ER_df) = c('Label', 'Variable', 'Value')
ER_df <- as.data.frame(ER_df)
TNBC_df <- as.data.frame(TNBC_df)

#de_ER_deseq <- list()
#de_TNBC_deseq <- list()

for (i in 1:length(depths)){

	for (j in 1:length(N_boots)){

		ER_df[5*((i-1)*length(N_boots)+(j-1))+1,1] <- 'DESeq2'
		ER_df[5*((i-1)*length(N_boots)+(j-1))+1,2] <- as.numeric(depths[i])
		ER_df[5*((i-1)*length(N_boots)+(j-1))+1,3] <- length(de_ER_deseq[[i]][[j]])

		ER_df[5*((i-1)*length(N_boots)+(j-1))+2,1] <- 'voom'
		ER_df[5*((i-1)*length(N_boots)+(j-1))+2,2] <- as.numeric(depths[i])
		ER_df[5*((i-1)*length(N_boots)+(j-1))+2,3] <-length(de_ER_voom[[i]][[j]])

		ER_df[5*((i-1)*length(N_boots)+(j-1))+3,1] <- 'edgeR'
		ER_df[5*((i-1)*length(N_boots)+(j-1))+3,2] <- as.numeric(depths[i])
		ER_df[5*((i-1)*length(N_boots)+(j-1))+3,3] <- length(de_ER_edger[[i]][[j]])

		ER_df[5*((i-1)*length(N_boots)+(j-1))+4,1] <- 'EBSeq'
		ER_df[5*((i-1)*length(N_boots)+(j-1))+4,2] <- as.numeric(depths[i])
		ER_df[5*((i-1)*length(N_boots)+(j-1))+4,3] <- length(de_ER_ebseq[[i]][[j]])

		ER_df[5*((i-1)*length(N_boots)+(j-1))+5,1] <- 'NOISeq'
		ER_df[5*((i-1)*length(N_boots)+(j-1))+5,2] <- as.numeric(depths[i])
		ER_df[5*((i-1)*length(N_boots)+(j-1))+5,3] <- length(de_ER_noiseq[[i]][[j]])

		TNBC_df[5*((i-1)*length(N_boots)+(j-1))+1,1] <- 'DESeq2'
		TNBC_df[5*((i-1)*length(N_boots)+(j-1))+1,2] <- as.numeric(depths[i])
		TNBC_df[5*((i-1)*length(N_boots)+(j-1))+1,3] <- length(de_TNBC_deseq[[i]][[j]])

		TNBC_df[5*((i-1)*length(N_boots)+(j-1))+2,1] <- 'voom'
		TNBC_df[5*((i-1)*length(N_boots)+(j-1))+2,2] <- as.numeric(depths[i])
		TNBC_df[5*((i-1)*length(N_boots)+(j-1))+2,3] <-length(de_TNBC_voom[[i]][[j]])

		TNBC_df[5*((i-1)*length(N_boots)+(j-1))+3,1] <- 'edgeR'
		TNBC_df[5*((i-1)*length(N_boots)+(j-1))+3,2] <- as.numeric(depths[i])
		TNBC_df[5*((i-1)*length(N_boots)+(j-1))+3,3] <- length(de_TNBC_edger[[i]][[j]])

		TNBC_df[5*((i-1)*length(N_boots)+(j-1))+4,1] <- 'EBSeq'
		TNBC_df[5*((i-1)*length(N_boots)+(j-1))+4,2] <- as.numeric(depths[i])
		TNBC_df[5*((i-1)*length(N_boots)+(j-1))+4,3] <- length(de_TNBC_ebseq[[i]][[j]])

		TNBC_df[5*((i-1)*length(N_boots)+(j-1))+5,1] <- 'NOISeq'
		TNBC_df[5*((i-1)*length(N_boots)+(j-1))+5,2] <- as.numeric(depths[i])
		TNBC_df[5*((i-1)*length(N_boots)+(j-1))+5,3] <- length(de_TNBC_noiseq[[i]][[j]])

		}


}


plotsamExplorer(ER_df, save = TRUE, p.depth = 0.9, filename = paste('ER_fold_',filter, sep = ''),y.lab="Number of DEGs", leg.name = 'Methods', title = paste('DE genes for ER+ data, ',xxxx, 'fold filtering',sep=''), x.lab = 'Data used in the analysis (f)', leg.pos = 'bottom', anova=F)
plotsamExplorer(TNBC_df, save = TRUE, p.depth = 0.9, filename = paste('TNBC_fold_',filter, sep = ''),y.lab="Number of DEGs", leg.name = 'Methods', title = paste('DE genes for TNBC data, ',xxxx, 'fold filtering',sep=''), x.lab = 'Data used in the analysis (f)', leg.pos = 'bottom', anova=F)








