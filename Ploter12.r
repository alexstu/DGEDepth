'-1 FDR(DEG)~d'
library(ggplot2)
library(samExploreR)



depths <- c('0.01','0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.85', '0.9', '0.95', '0.99', '1')
N_boots <- 1:25
filter = 'none'
ER_samples <- c('All','SRR1313092.sam','SRR1313100.sam','SRR1313105.sam','SRR1313107.sam','SRR1313124.bam','SRR1313174.bam','SRR1313175.bam','SRR1313178.bam','SRR1313197.bam','SRR1313203.bam')
TNBC_samples <- c('All', 'SRR1313141.bam','SRR1313142.bam','SRR1313147.bam','SRR1313152.bam','SRR1313164.bam','SRR1313217.bam','SRR1313220.bam','SRR1313226.bam','SRR1313227.bam', 'SRR1313229.bam')
methods <- c('DESeq2', 'voom', 'edgeR', 'EBSeq', 'NOISeq')

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

if (filter == 'none'){
	xxxx = 'no '
}
if (filter == '2'){
	xxxx = '2-'
}

#############################

ER_df <- matrix(ncol = 3, nrow = 25*17*11)
TNBC_df <- matrix(ncol = 3, nrow = 25*17*11)
colnames(TNBC_df) = c('Label', 'Variable', 'Value')
colnames(ER_df) = c('Label', 'Variable', 'Value')
ER_df <- as.data.frame(ER_df)
TNBC_df <- as.data.frame(TNBC_df)

ER_df_list <- list()
TNBC_df_list <- list()

ER_results <- list(de_ER_deseq, de_ER_voom, de_ER_edger, de_ER_ebseq, de_ER_noiseq)
TNBC_results <- list(de_TNBC_deseq, de_TNBC_voom, de_TNBC_edger, de_TNBC_ebseq, de_TNBC_noiseq)

#de_ER_deseq <- list()
#de_TNBC_deseq <- list()

for (m in 1:length(methods)){

	for (i in 1:length(depths)){
		for (j in 1:length(N_boots)){
			for (k in 1:11){


				ER_df[11*((i-1)*length(N_boots)+(j-1))+k,1] <- ER_samples[k]
				ER_df[11*((i-1)*length(N_boots)+(j-1))+k,2] <- as.numeric(depths[i])
				TP = length(intersect(ER_results[[m]][[k]][[i]][[j]], ER_results[[m]][[k]][[17]][[j]]))
				FN = length(ER_results[[m]][[k]][[17]][[j]]) - TP
				FP = length(ER_results[[m]][[k]][[i]][[j]]) - TP
				ER_df[11*((i-1)*length(N_boots)+(j-1))+k,3] <- (FP/(TP+FP))


				TNBC_df[11*((i-1)*length(N_boots)+(j-1))+k,1] <- TNBC_samples[k]
				TNBC_df[11*((i-1)*length(N_boots)+(j-1))+k,2] <- as.numeric(depths[i])
				TP = length(intersect(TNBC_results[[m]][[k]][[i]][[j]], TNBC_results[[m]][[k]][[17]][[j]]))
				FN = length(TNBC_results[[m]][[k]][[17]][[j]]) - TP
				FP = length(TNBC_results[[m]][[k]][[i]][[j]]) - TP
				TNBC_df[11*((i-1)*length(N_boots)+(j-1))+k,3] <- (FP/(TP+FP))

			}
		}
	}
	ER_df_list[[m]] <- ER_df
	TNBC_df_list[[m]] <- TNBC_df
}


for (m in 1:length(methods)){
	plotsamExplorer(ER_df_list[[m]], save = TRUE, p.depth = 0.9, filename = paste('-1_', methods[m], '_', 'FDR_ER_fold_',filter, sep = ''),anova = FALSE,y.lab="FDR", leg.name = 'Methods', notch = FALSE, title = paste('FDR for ER+ data, ',xxxx, 'fold filtering',sep=''))
	plotsamExplorer(TNBC_df_list[[m]], save = TRUE, p.depth = 0.9, filename = paste('-1_', methods[m], '_', 'FDR_TNBC_fold_',filter, sep = ''),anova = FALSE,y.lab="FDR", leg.name = 'Methods', notch = FALSE,title = paste('FDR for TNBC data, ',xxxx, 'fold filtering',sep=''))
}

#plotsamExplorer(ER_df_list[[m]], save = TRUE, p.depth = 0.9, filename = paste('-1_', methods[m], '_', 'FDR_ER_fold_',filter, sep = ''),anova = FALSE)#,y.lab="FDR", leg.name = 'Methods', title = paste('FDR for ER+ data, ',xxxx, 'fold filtering',sep=''))


