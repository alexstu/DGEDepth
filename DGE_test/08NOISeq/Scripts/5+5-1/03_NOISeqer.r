library("NOISeq")
library(doParallel)

cv_dir = 'CV/'

types = c('ER+', 'TNBC')
depths <- c('0.01','0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.85', '0.9', '0.95', '0.99', '1')
N_boots <- 1:25
ER_samples <- c('All','SRR1313092.sam','SRR1313100.sam','SRR1313105.sam','SRR1313107.sam','SRR1313124.bam','SRR1313174.bam','SRR1313175.bam','SRR1313178.bam','SRR1313197.bam','SRR1313203.bam')
TNBC_samples <- c('All', 'SRR1313141.bam','SRR1313142.bam','SRR1313147.bam','SRR1313152.bam','SRR1313164.bam','SRR1313217.bam','SRR1313220.bam','SRR1313226.bam','SRR1313227.bam', 'SRR1313229.bam')


#types = c('ER+', 'TNBC')
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

				library("NOISeq");
				filename <- paste(cv_dir, type, '_', as.character(depths[i]), '_', as.character(j), '.txt', sep = '')
				x <- read.table(filename, header = TRUE, row.names = 'GeneID')
				x <- x[,sub_names]

				if (s == 1){
					x <- x
					if (type == 'ER+'){
						y = data.frame(Group = factor(c(rep(1,5), rep(0,5))))
					}
					if (type == 'TNBC'){
						y = data.frame(Group = factor(c(rep(1,5), rep(0,5))))
					}
				}
				if ((s > 1) & (s <= 6)){
					x <- x[, colnames(x) != samples[s]]
					if (type == 'ER+'){
						y = data.frame(Group = factor(c(rep(1,4), rep(0,5))))
					}
					if (type == 'TNBC'){
						y = data.frame(Group = factor(c(rep(1,4), rep(0,5))))
					}
				}
				if (s > 6){
					x <- x[, colnames(x) != samples[s]]
					if (type == 'ER+'){
						y = data.frame(Group = factor(c(rep(1,5), rep(0,4))))
					}
					if (type == 'TNBC'){
						y = data.frame(Group = factor(c(rep(1,5), rep(0,4))))
					}
				}


				data <- readData(data = x, factors = y)
				 #uqua(assayData(data)$exprs, long = 1000, lc = 0)
				res <- noiseqbio(data, norm = 'uqua', factor = 'Group', conditions = c(1,0), lc = 0)
				#boot_list[[j]] <-degenes(res)
				res

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


save.image('sub_JKNOISeq.RData')
#########################
#############################

load('sub_JKNOISeq.RData')


de_ER_noiseq_no_fold <- list()
de_TNBC_noiseq_no_fold <- list()
de_ER_noiseq_fold <- list()
de_TNBC_noiseq_fold <- list()

for (s in 1:length(samples)){

	s_de_ER_noiseq_no_fold <- list()
	s_de_TNBC_noiseq_no_fold <- list()
	s_de_ER_noiseq_fold <- list()
	s_de_TNBC_noiseq_fold <- list()

	for (i in 1:length(depths)){

		ER_boot_list_no_fold <- list()
		TNBC_boot_list_no_fold <- list()
		ER_boot_list_fold <- list()
		TNBC_boot_list_fold <- list()

		for (j in 1:length(N_boots)){

			ER_res <- degenes(ER_list[[s]][[i]][[j]])
			TNBC_res <- degenes(TNBC_list[[s]][[i]][[j]])

			ER_boot_list_no_fold[[j]] <- rownames(ER_res)
			TNBC_boot_list_no_fold[[j]] <- rownames(TNBC_res)

			ER_boot_list_fold[[j]] <- rownames(ER_res[ER_res$log2FC > 2,])
			TNBC_boot_list_fold[[j]] <-rownames(TNBC_res[TNBC_res$log2FC > 2,])

		}

		s_de_ER_noiseq_no_fold[[i]] <- ER_boot_list_no_fold
		s_de_TNBC_noiseq_no_fold[[i]] <- TNBC_boot_list_no_fold
		s_de_ER_noiseq_fold[[i]] <- ER_boot_list_fold
		s_de_TNBC_noiseq_fold[[i]] <- TNBC_boot_list_fold

	}

	de_ER_noiseq_no_fold[[s]] <- s_de_ER_noiseq_no_fold
	de_TNBC_noiseq_no_fold[[s]] <- s_de_TNBC_noiseq_no_fold
	de_ER_noiseq_fold[[s]] <- s_de_ER_noiseq_fold
	de_TNBC_noiseq_fold[[s]] <- s_de_TNBC_noiseq_fold

}

de_ER_noiseq <- list()
de_TNBC_noiseq <- list()


de_ER_noiseq <- de_ER_noiseq_no_fold
save('de_ER_noiseq', file = 'sub_JK_de_ER_noiseq_0.05_none')

de_TNBC_noiseq <- de_TNBC_noiseq_no_fold
save('de_TNBC_noiseq', file = 'sub_JK_de_TNBC_noiseq_0.05_none')

de_ER_noiseq <- de_ER_noiseq_fold
save('de_ER_noiseq', file = 'sub_JK_de_ER_noiseq_0.05_2')

de_TNBC_noiseq <- de_TNBC_noiseq_fold
save('de_TNBC_noiseq', file = 'sub_JK_de_TNBC_noiseq_0.05_2')














