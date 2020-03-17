'-1 FDR~d slopes and Friedman test'
library(ggplot2)
library(samExploreR)
library(PMCMR)


depths <- c('0.01','0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.85', '0.9', '0.95', '0.99', '1')
N_boots <- 1:25
filter = '2'
#f_test <- c('0.7','0.8','0.85', '0.9', '0.95', '0.99')
f_test <- c('0.8','0.85', '0.9', '0.95', '0.99')
#f_test <- c('0.85', '0.9', '0.95', '0.99')
#f_test <- c('0.9', '0.95', '0.99')

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
f_vector <- depths %in% f_test
x_vector <- vector()
y_ER_vector <- vector()
y_TNBC_vector <- vector()


ER_results <- list(de_ER_deseq, de_ER_voom, de_ER_edger, de_ER_ebseq, de_ER_noiseq)
TNBC_results <- list(de_TNBC_deseq, de_TNBC_voom, de_TNBC_edger, de_TNBC_ebseq, de_TNBC_noiseq)

ER_slopes_list <- list()
TNBC_slopes_list <- list()

ER_FDR_list <- list()
TNBC_FDR_list <- list()

ER_df_list <- list()
TNBC_df_list <- list()

ER_df <- matrix(ncol = 4, nrow = length(methods)*length(ER_samples)*length(N_boots))
TNBC_df <- matrix(ncol = 4, nrow = length(methods)*length(TNBC_samples)*length(N_boots))
colnames(TNBC_df) = c('Software', 'Variable', 'Slope', 'Boot')
colnames(ER_df) = c('Software', 'Variable', 'Slope', 'Boot')
ER_df <- as.data.frame(ER_df)
TNBC_df <- as.data.frame(TNBC_df)



#de_ER_deseq <- list()
#de_TNBC_deseq <- list()

for (m in 1:length(methods)){

	methods_ER_list <- list()
	methods_TNBC_list <- list()

	for (k in 1:11){

		sample_ER_list <- list()
		sample_TNBC_list <- list()

		for (i in 1:length(depths)){

			depth_ER_list <- list()
			depth_TNBC_list <- list()

			for (j in 1:length(N_boots)){

				boot_ER_list <- list()
				boot_TNBC_list <- list()

				TP = length(intersect(ER_results[[m]][[k]][[i]][[j]], ER_results[[m]][[k]][[17]][[j]]))
				FN = length(ER_results[[m]][[k]][[17]][[j]]) - TP
				FP = length(ER_results[[m]][[k]][[i]][[j]]) - TP
				FDR <- (FP/(TP+FP))
				depth_ER_list[[j]] <- FDR


				TP = length(intersect(TNBC_results[[m]][[k]][[i]][[j]], TNBC_results[[m]][[k]][[17]][[j]]))
				FN = length(TNBC_results[[m]][[k]][[17]][[j]]) - TP
				FP = length(TNBC_results[[m]][[k]][[i]][[j]]) - TP
				FDR <- (FP/(TP+FP))
				depth_TNBC_list[[j]] <- FDR
			}
			sample_ER_list[[i]] <- depth_ER_list
			sample_TNBC_list[[i]] <- depth_TNBC_list
		}
	methods_ER_list[[k]] <- sample_ER_list
	methods_TNBC_list[[k]] <- sample_TNBC_list

	}
	ER_FDR_list[[m]] <- methods_ER_list
	TNBC_FDR_list[[m]] <- methods_TNBC_list

}


x_vector <- f_test

for (m in 1:length(methods)){

	methods_slope_ER <- list()
	methods_slope_TNBC <- list()

	for (k in 1:11){

		sample_slope_ER <- list()
		sample_slope_TNBC <- list()

		f_ER_FDR <- ER_FDR_list[[m]][[k]][f_vector]
		f_TNBC_FDR <- TNBC_FDR_list[[m]][[k]][f_vector]

		for (j in 1:length(N_boots)){

			y_ER_vector <- vector()
			y_TNBC_vector <- vector()

			for (i in 1:length(f_ER_FDR)){

				y_ER_vector <- c(y_ER_vector, f_ER_FDR[[i]][[j]])
				y_TNBC_vector <- c(y_TNBC_vector, f_TNBC_FDR[[i]][[j]])

			}

			fit_ER <- lm(y_ER_vector ~ x_vector)
			fit_TNBC <- lm(y_TNBC_vector ~ x_vector)

			sample_slope_ER[[j]] <- fit_ER$coefficients[[2]]
			sample_slope_TNBC[[j]] <- fit_TNBC$coefficients[[2]]
		}
		methods_slope_ER[[k]] <- sample_slope_ER
		methods_slope_TNBC[[k]] <- sample_slope_TNBC
	}
	ER_slopes_list[[m]] <- methods_slope_ER
	TNBC_slopes_list[[m]] <- methods_slope_TNBC
}



for (m in 1:length(methods)){
	for (k in 1:11){
		for (j in 1:length(N_boots)){

			ER_df[(m-1)*11*length(N_boots) + (k-1)*length(N_boots)+j,1] <- methods[m]
			ER_df[(m-1)*11*length(N_boots) + (k-1)*length(N_boots)+j,2] <- ER_samples[k]
			ER_df[(m-1)*11*length(N_boots) + (k-1)*length(N_boots)+j,3] <- ER_slopes_list[[m]][[k]][[j]]
			ER_df[(m-1)*11*length(N_boots) + (k-1)*length(N_boots)+j,4] <- j


			TNBC_df[(m-1)*11*length(N_boots) + (k-1)*length(N_boots)+j,1] <- methods[m]
			TNBC_df[(m-1)*11*length(N_boots) + (k-1)*length(N_boots)+j,2] <- TNBC_samples[k]
			TNBC_df[(m-1)*11*length(N_boots) + (k-1)*length(N_boots)+j,3] <- TNBC_slopes_list[[m]][[k]][[j]]
			TNBC_df[(m-1)*11*length(N_boots) + (k-1)*length(N_boots)+j,4] <- j
		}
	}
}


pdf(paste('TNBC_slopes_', f_test[1], '_fold_',filter, '.pdf', sep = ''))
ggplot(data=TNBC_df, aes(x=Software, y=Slope)) + geom_boxplot(fill = c('red', 'darkorange', 'green','blue','purple')) + ggtitle(paste("TNBC data -1 slopes\nfor f values: ",paste(f_test,collapse=", "), sep='')) + theme(plot.title = element_text(hjust=0.5))
dev.off()

pdf(paste('ER_slopes_', f_test[1], '_fold_',filter, '.pdf', sep = ''))
ggplot(data=ER_df, aes(x=Software, y=Slope)) + geom_boxplot(fill = c('red', 'darkorange', 'green','blue','purple')) + ggtitle(paste("ER data -1 slopes\nfor f values: ",paste(f_test,collapse=", "), sep='')) + theme(plot.title = element_text(hjust=0.5))
dev.off()

ER_df_pval <- matrix(ncol = 2, nrow = 5)
TNBC_df_pval <- matrix(ncol = 2, nrow = 5)
ER_df_pval <- as.data.frame(ER_df_pval)
TNBC_df_pval <- as.data.frame(TNBC_df_pval)
colnames(TNBC_df_pval) = c('Software', 'Pval')
colnames(ER_df_pval) = c('Software', 'Pval')


for (m in 1:length(methods)){

	TNBC_df_method <- TNBC_df[TNBC_df[,1] == methods[m],]
	ER_df_method <- ER_df[ER_df[,1] == methods[m],]

	ER_df_pval[m,1] <- methods[m]
	TNBC_df_pval[m,1] <- methods[m]
	TNBC_df_pval[m,2] <- (friedman.test(Slope ~ Variable | Boot, data=TNBC_df_method)$p.value)
	ER_df_pval[m,2] <- (friedman.test(Slope ~ Variable | Boot, data=ER_df_method)$p.value)
}
pdf(paste('TNBC_pvalues_', f_test[1], '_fold_',filter, '.pdf', sep = ''))
ggplot(data=TNBC_df_pval, aes(x=Software, y=Pval)) + geom_boxplot(fill = c('red', 'darkorange', 'green','blue','purple')) + ggtitle(paste("TNBC data -1 p-values\nfor f values: ",paste(f_test,collapse=", "), sep='')) + theme(plot.title = element_text(hjust=0.5))
dev.off()
pdf(paste('ER_pvalues_', f_test[1], '_fold_',filter, '.pdf', sep = ''))
ggplot(data=ER_df_pval, aes(x=Software, y=Pval)) + geom_boxplot(fill = c('red', 'darkorange', 'green','blue','purple')) + ggtitle(paste("ER data -1 p-values\nfor f values: ",paste(f_test,collapse=", "), sep='')) + theme(plot.title = element_text(hjust=0.5))
dev.off()
#posthoc.friedman.nemenyi.test(ER_matrix)
#posthoc.friedman.nemenyi.test(TNBC_matrix)


#posthoc.friedman.conover.test(ER_matrix)
#posthoc.friedman.conover.test(TNBC_matrix)


table(ER_df[,1])

