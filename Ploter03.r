'All N(GO)~d'
library(ggplot2)
library(samExploreR)



depths <- c('0.01','0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.85', '0.9', '0.95', '0.99', '1')
N_boots <- 1:25
filter = '2'


load(paste('go_ER_deseq_0.05_',filter, sep = ''))
load(paste('go_ER_ebseq_0.05_',filter, sep = ''))
load(paste('go_ER_edger_0.05_',filter, sep = ''))
load(paste('go_ER_noiseq_0.05_',filter, sep = ''))
load(paste('go_ER_voom_0.05_',filter, sep = ''))
load(paste('go_TNBC_deseq_0.05_',filter, sep = ''))
load(paste('go_TNBC_ebseq_0.05_',filter, sep = ''))
load(paste('go_TNBC_edger_0.05_',filter, sep = ''))
load(paste('go_TNBC_noiseq_0.05_',filter, sep = ''))
load(paste('go_TNBC_voom_0.05_',filter, sep = ''))

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
		ER_df[5*((i-1)*length(N_boots)+(j-1))+1,3] <- length(go_ER_deseq[[i]][[j]])

		ER_df[5*((i-1)*length(N_boots)+(j-1))+2,1] <- 'voom'
		ER_df[5*((i-1)*length(N_boots)+(j-1))+2,2] <- as.numeric(depths[i])
		ER_df[5*((i-1)*length(N_boots)+(j-1))+2,3] <-length(go_ER_voom[[i]][[j]])

		ER_df[5*((i-1)*length(N_boots)+(j-1))+3,1] <- 'edgeR'
		ER_df[5*((i-1)*length(N_boots)+(j-1))+3,2] <- as.numeric(depths[i])
		ER_df[5*((i-1)*length(N_boots)+(j-1))+3,3] <- length(go_ER_edger[[i]][[j]])

		ER_df[5*((i-1)*length(N_boots)+(j-1))+4,1] <- 'EBSeq'
		ER_df[5*((i-1)*length(N_boots)+(j-1))+4,2] <- as.numeric(depths[i])
		ER_df[5*((i-1)*length(N_boots)+(j-1))+4,3] <- length(go_ER_ebseq[[i]][[j]])

		ER_df[5*((i-1)*length(N_boots)+(j-1))+5,1] <- 'NOISeq'
		ER_df[5*((i-1)*length(N_boots)+(j-1))+5,2] <- as.numeric(depths[i])
		ER_df[5*((i-1)*length(N_boots)+(j-1))+5,3] <- length(go_ER_noiseq[[i]][[j]])

		TNBC_df[5*((i-1)*length(N_boots)+(j-1))+1,1] <- 'DESeq2'
		TNBC_df[5*((i-1)*length(N_boots)+(j-1))+1,2] <- as.numeric(depths[i])
		TNBC_df[5*((i-1)*length(N_boots)+(j-1))+1,3] <- length(go_TNBC_deseq[[i]][[j]])

		TNBC_df[5*((i-1)*length(N_boots)+(j-1))+2,1] <- 'voom'
		TNBC_df[5*((i-1)*length(N_boots)+(j-1))+2,2] <- as.numeric(depths[i])
		TNBC_df[5*((i-1)*length(N_boots)+(j-1))+2,3] <-length(go_TNBC_voom[[i]][[j]])

		TNBC_df[5*((i-1)*length(N_boots)+(j-1))+3,1] <- 'edgeR'
		TNBC_df[5*((i-1)*length(N_boots)+(j-1))+3,2] <- as.numeric(depths[i])
		TNBC_df[5*((i-1)*length(N_boots)+(j-1))+3,3] <- length(go_TNBC_edger[[i]][[j]])

		TNBC_df[5*((i-1)*length(N_boots)+(j-1))+4,1] <- 'EBSeq'
		TNBC_df[5*((i-1)*length(N_boots)+(j-1))+4,2] <- as.numeric(depths[i])
		TNBC_df[5*((i-1)*length(N_boots)+(j-1))+4,3] <- length(go_TNBC_ebseq[[i]][[j]])

		TNBC_df[5*((i-1)*length(N_boots)+(j-1))+5,1] <- 'NOISeq'
		TNBC_df[5*((i-1)*length(N_boots)+(j-1))+5,2] <- as.numeric(depths[i])
		TNBC_df[5*((i-1)*length(N_boots)+(j-1))+5,3] <- length(go_TNBC_noiseq[[i]][[j]])

		}
}

plotsamExplorer(ER_df, save = TRUE, p.depth = 0.9, filename = paste('GO_ER_fold_',filter, sep = ''),y.lab="Number of GO terms", leg.name = 'Methods', title = paste('GO terms for ER+ data, ',xxxx, 'fold filtering',sep=''), x.lab = 'Data used in the analysis (f)', leg.pos = 'bottom', anova=F)
plotsamExplorer(TNBC_df, save = TRUE, p.depth = 0.9, filename = paste('GO_TNBC_fold_',filter, sep = ''),y.lab="Number of GO terms", leg.name = 'Methods', title = paste('GO terms for TNBC data, ',xxxx, 'fold filtering',sep=''), x.lab = 'Data used in the analysis (f)', leg.pos = 'bottom', anova=F)





