'-1 FDR~d slopes'
library(ggplot2)
library(samExploreR)
library(PMCMR)


depths <- c('0.01','0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.85', '0.9', '0.95', '0.99', '1')
N_boots <- 1:25
filter = '2'
f_test_070 <- c('0.7','0.8','0.85', '0.9', '0.95', '0.99')
f_test_080 <- c('0.8','0.85', '0.9', '0.95', '0.99')
f_test_085 <- c('0.85', '0.9', '0.95', '0.99')
f_test_090 <- c('0.9', '0.95', '0.99')



f_vector_list <- list(f_test_070, f_test_080, f_test_085, f_test_090)

ER_samples <- c('All','SRR1313092.sam','SRR1313100.sam','SRR1313105.sam','SRR1313107.sam','SRR1313124.bam','SRR1313174.bam','SRR1313175.bam','SRR1313178.bam','SRR1313197.bam','SRR1313203.bam')
TNBC_samples <- c('All', 'SRR1313141.bam','SRR1313142.bam','SRR1313147.bam','SRR1313152.bam','SRR1313164.bam','SRR1313217.bam','SRR1313220.bam','SRR1313226.bam','SRR1313227.bam', 'SRR1313229.bam')
#methods <- c('DESeq2', 'voom', 'edgeR', 'EBSeq', 'NOISeq')
methods <- c('DESeq2', 'EBSeq', 'voom', 'edgeR', 'NOISeq')

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

#load(paste('sub_JK_de_ER_deseq_0.05_',filter, sep = ''))
#load(paste('sub_JK_de_ER_ebseq_0.05_',filter, sep = ''))
#load(paste('sub_JK_de_ER_edger_0.05_',filter, sep = ''))
#load(paste('sub_JK_de_ER_noiseq_0.05_',filter, sep = ''))
#load(paste('sub_JK_de_ER_voom_0.05_',filter, sep = ''))
#load(paste('sub_JK_de_TNBC_deseq_0.05_',filter, sep = ''))
#load(paste('sub_JK_de_TNBC_ebseq_0.05_',filter, sep = ''))
#load(paste('sub_JK_de_TNBC_edger_0.05_',filter, sep = ''))
#load(paste('sub_JK_de_TNBC_noiseq_0.05_',filter, sep = ''))
#load(paste('sub_JK_de_TNBC_voom_0.05_',filter, sep = ''))


if (filter == 'none'){
	xxxx = 'no '
}
if (filter == '2'){
	xxxx = '2-'
}
ER_df_list <- list()
TNBC_df_list <- list()

for (f_i in 1:length(f_vector_list)){

	f_test <- f_vector_list[[f_i]] 
	f_vector <- depths %in% f_test
	x_vector <- vector()
	y_ER_vector <- vector()
	y_TNBC_vector <- vector()


#ER_results <- list(de_ER_deseq, de_ER_voom, de_ER_edger, de_ER_ebseq, de_ER_noiseq)
#TNBC_results <- list(de_TNBC_deseq, de_TNBC_voom, de_TNBC_edger, de_TNBC_ebseq, de_TNBC_noiseq)

	ER_results <- list(de_ER_deseq, de_ER_ebseq, de_ER_voom, de_ER_edger, de_ER_noiseq)
	TNBC_results <- list(de_TNBC_deseq, de_TNBC_ebseq, de_TNBC_voom, de_TNBC_edger, de_TNBC_noiseq)


	ER_slopes_list <- list()
	TNBC_slopes_list <- list()

	ER_FDR_list <- list()
	TNBC_FDR_list <- list()

	#ER_df_list <- list()
	#TNBC_df_list <- list()

	ER_df <- matrix(ncol = 3, nrow = length(methods)*length(ER_samples))
	TNBC_df <- matrix(ncol = 3, nrow = length(methods)*length(TNBC_samples))
	colnames(TNBC_df) = c('Label', 'Variable', 'Slope')
	colnames(ER_df) = c('Label', 'Variable', 'Slope')
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


	for (f in f_test){
		for (j in 1:length(N_boots)){
			x_vector <- c(x_vector, f)
		}
	}

	for (m in 1:length(methods)){

		methods_slope_ER <- list()
		methods_slope_TNBC <- list()

		sample_slope_ER <- list()
		sample_slope_TNBC <- list()

		for (k in 1:11){

			y_ER_vector <- list()
			y_TNBC_vector <- list()

			f_ER_FDR <- ER_FDR_list[[m]][[k]][f_vector]
			for (i in 1:length(f_ER_FDR)){
				for (j in 1:length(f_ER_FDR[[i]])){
					y_ER_vector[[(i-1)*length(f_ER_FDR[[i]])+j]] <- f_ER_FDR[[i]][[j]]
				}
			}
			f_TNBC_FDR <- TNBC_FDR_list[[m]][[k]][f_vector]
			for (i in 1:length(f_TNBC_FDR)){
				for (j in 1:length(f_TNBC_FDR[[i]])){
					y_TNBC_vector[[(i-1)*length(f_TNBC_FDR[[i]])+j]] <- f_TNBC_FDR[[i]][[j]]
				}
			}
			fit_ER <- lm(unlist(y_ER_vector) ~ x_vector)
			fit_TNBC <- lm(unlist(y_TNBC_vector) ~ x_vector)

			methods_slope_ER[[k]] <- fit_ER$coefficients[[2]]
			methods_slope_TNBC[[k]] <- fit_TNBC$coefficients[[2]]
		}
		ER_slopes_list[[m]] <- methods_slope_ER
		TNBC_slopes_list[[m]] <- methods_slope_TNBC
	}

	for (m in 1:length(methods)){
		for (k in 1:11){


			ER_df[11*(m-1)+k,1] <- methods[m]
			ER_df[11*(m-1)+k,2] <- ER_samples[k]
			ER_df[11*(m-1)+k,3] <- ER_slopes_list[[m]][[k]]


			TNBC_df[11*(m-1)+k,1] <- methods[m]
			TNBC_df[11*(m-1)+k,2] <- TNBC_samples[k]
			TNBC_df[11*(m-1)+k,3] <- TNBC_slopes_list[[m]][[k]]

		}
	}
	print(f_i)
	ER_df_list[[f_i]] <- ER_df
	TNBC_df_list[[f_i]] <- TNBC_df

#	ER_matrix <- matrix(ncol = length(methods), nrow = length(ER_samples))
#	colnames(ER_matrix) <- methods
#	rownames(ER_matrix) <- ER_samples
#	for (i in 1:length(methods)){
#		for (j in 1:length(ER_samples)){
#			ER_matrix[j,i] <- ER_df[((ER_df[,1]==methods[i])&(ER_df[,2]==ER_samples[j])),3]
#		}
#	}

#	TNBC_matrix <- matrix(ncol = length(methods), nrow = length(TNBC_samples))
#	colnames(TNBC_matrix) <- methods
#	rownames(TNBC_matrix) <- TNBC_samples
#	for (i in 1:length(methods)){
#		for (j in 1:length(TNBC_samples)){
#			TNBC_matrix[j,i] <- TNBC_df[((TNBC_df[,1]==methods[i])&(TNBC_df[,2]==TNBC_samples[j])),3]
#		}
#	}

}



#pdf(paste('TNBC_slopes_', f_test[1], '.pdf', sep = ''))
pdf(paste('TNBC_slopes_', '.pdf', sep = ''))
#TNBC_df$Label <- factor(TNBC_df$Label, as.character(TNBC_df$Label))
#ggplot(data=TNBC_df, aes(x=Label, y=Slope)) + geom_boxplot(fill = c('red', 'darkorange', 'purple','green','blue')) + ggtitle(paste("TNBC data -1 slopes\nfor f values: ",paste(f_test,collapse=", "), sep=''))
par(mfrow=c(2,2))
ggplot(data=TNBC_df_list[[1]], aes(x=Label, y=Slope)) + geom_boxplot(fill = c('red', 'darkorange', 'green','blue','purple')) + ggtitle(paste("TNBC data -1 slopes\nfor f values: ",paste(f_vector_list[[1]],collapse=", "), sep=''))
ggplot(data=TNBC_df_list[[2]], aes(x=Label, y=Slope)) + geom_boxplot(fill = c('red', 'darkorange', 'green','blue','purple')) + ggtitle(paste("TNBC data -1 slopes\nfor f values: ",paste(f_vector_list[[2]],collapse=", "), sep=''))
ggplot(data=TNBC_df_list[[3]], aes(x=Label, y=Slope)) + geom_boxplot(fill = c('red', 'darkorange', 'green','blue','purple')) + ggtitle(paste("TNBC data -1 slopes\nfor f values: ",paste(f_vector_list[[3]],collapse=", "), sep=''))
ggplot(data=TNBC_df_list[[4]], aes(x=Label, y=Slope)) + geom_boxplot(fill = c('red', 'darkorange', 'green','blue','purple')) + ggtitle(paste("TNBC data -1 slopes\nfor f values: ",paste(f_vector_list[[4]],collapse=", "), sep=''))
dev.off()

pdf(paste('ER_slopes_', f_test[1], '.pdf', sep = ''))
#ggplot(data=ER_df, aes(x=Label, y=Slope)) + geom_boxplot(fill = c('red', 'darkorange', 'green','blue','purple')) + ggtitle(paste("ER data -1 slopes\nfor f values: ",paste(f_test,collapse=", "), sep=''))
ER_df$Label <- factor(ER_df$Label, as.character(ER_df$Label))
ggplot(data=ER_df, aes(x=Label, y=Slope)) + geom_boxplot(fill = c('red', 'darkorange', 'purple','green','blue')) + ggtitle(paste("ER data -1 slopes\nfor f values: ",paste(f_test,collapse=", "), sep=''))
dev.off()



pdf('ER_cluster')
#d <- posthoc.friedman.nemenyi.test(ER_matrix)
d <- posthoc.friedman.conover.test(ER_matrix)

d <- d$statistic
dd <- matrix(ncol=5, nrow=5)

rownames(dd) <- methods
colnames(dd) <- methods

dd[lower.tri(dd, diag=F)] <- d[lower.tri(d, diag=T)]

for (i in 1:length(methods)){
	for (j in 1:length(methods)){
		if (i==j){
			dd[i,j] <- 0
		}
		if (is.na(dd[i,j])==T){
			dd[i,j] <- dd[j,i]
		}
	}
}
plot(hclust(dist(dd)))
dev.off()




pdf('TNBC_cluster')
#d <- posthoc.friedman.nemenyi.test(TNBC_matrix)
d <- posthoc.friedman.conover.test(TNBC_matrix)


d <- d$statistic
dd <- matrix(ncol=5, nrow=5)

rownames(dd) <- methods
colnames(dd) <- methods

dd[lower.tri(dd, diag=F)] <- d[lower.tri(d, diag=T)]

for (i in 1:length(methods)){
	for (j in 1:length(methods)){
		if (i==j){
			dd[i,j] <- 0
		}
		if (is.na(dd[i,j])==T){
			dd[i,j] <- dd[j,i]
		}
	}
}
plot(hclust(dist(dd)))
dev.off()










