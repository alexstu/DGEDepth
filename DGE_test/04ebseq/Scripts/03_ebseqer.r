library("EBSeq")
types = c('ER+', 'TNBC')

depths <- c('0.01','0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.85', '0.9', '0.95', '0.99', '1')
N_boots <- 1:25

#depths = c(0.01, 0.05, 0.15)
#N_boots <- 1:2

ER_list <- list()
TNBC_list <- list()





for (type in types){

	depth_list <- list()

	if (type == 'ER+'){
		type_factor <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
	}
	if (type == 'TNBC'){
		type_factor <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
	}



	for (i in 1:length(depths)){

		boot_list <- list()

		for (j in 1:length(N_boots)){

			file = paste('CV/', type, '_', depths[i], '_', as.character(N_boots[j]), '.txt', sep = '')
			print(file)
			x <- read.table(file,header=TRUE, row.names = 'GeneID')
			x <- as.matrix(x)
			group <- factor(type_factor)

			Sizes=MedianNorm(x)
			EBOut=EBTest(Data=x, Conditions=group,sizeFactors=Sizes, maxround=5)
			#EBDERes=GetDEResults(EBOut, FDR=0.05)

			#boot_list[[j]] <- EBDERes
			boot_list[[j]] <- EBOut

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


save.image('ebseq.RData')

#############################

load('edgeR.RData')


pval = 0.001
fold = 3

de_ER_edger <- list()
de_TNBC_edger <- list()

for (i in 1:length(depths)){

	ER_boot_list <- list()
	TNBC_boot_list <- list()

	for (j in 1:length(N_boots)){

		de_ER <- ER_list[[i]][[j]]
		de_TNBC <- TNBC_list[[i]][[j]]

		if (fold == 'none'){
			ER_boot_list[[j]] <- rownames(de_ER[(de_ER$table$FDR < pval),])
			TNBC_boot_list[[j]] <- rownames(de_TNBC[(de_TNBC$table$FDR < pval),])
		}
		else {
			ER_boot_list[[j]] <- rownames(de_ER[(de_ER$table$FDR < pval) & (de_ER$table$logFC > fold),])
			TNBC_boot_list[[j]] <- rownames(de_TNBC[(de_TNBC$table$FDR < pval) & (de_TNBC$table$logFC > fold),])
		}

	}

	de_ER_edger[[i]] <- ER_boot_list
	de_TNBC_edger[[i]] <- TNBC_boot_list

}

save('de_ER_edger', file = paste('de_ER_edger_',as.character(pval), '_', as.character(fold), sep = ''))
save('de_TNBC_edger', file = paste('de_TNBC_edger_',as.character(pval), '_', as.character(fold), sep = ''))


########################################################


de_ER_ebseq <- list()
de_TNBC_ebseq <- list()



for (i in 1:length(depths)){

	ER_boot_list <- list()
	TNBC_boot_list <- list()

	for (j in 1:length(N_boots)){

		de_ER <- ER_list[[i]][[j]]
		de_TNBC <- TNBC_list[[i]][[j]]

		ER_boot_list[[j]] <- de_ER$DEfound
		TNBC_boot_list[[j]] <- de_TNBC$DEfound

		}

	de_ER_ebseq[[i]] <- ER_boot_list
	de_TNBC_ebseq[[i]] <- TNBC_boot_list

}

save('de_ER_ebseq', file = 'de_ER_ebseq')
save('de_TNBC_ebseq', file = 'de_TNBC_ebseq')













