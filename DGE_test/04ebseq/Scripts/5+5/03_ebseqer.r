library("EBSeq")
library(foreach)
library(doParallel)

types = c('ER+', 'TNBC')
depths <- c('0.01','0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.85', '0.9', '0.95', '0.99', '1')
N_boots <- 1:25

#types = c('ER+')
#depths = c(0.01,0.99)
#N_boots <- 1:2

ER_list_de_fold <- list()
TNBC_list_de_fold <- list()

ER_list_de_all <- list()
TNBC_list_de_all <- list()

#ER_list_go <- list()
#TNBC_list_go <- list()

EBDERes_all <- list()
EBDERes_fold <- list()


for (type in types){

	depth_list_all <- list()
	depth_list_fold <- list()

	if (type == 'ER+'){
		type_factor <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
	}
	if (type == 'TNBC'){
		type_factor <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
	}



	for (i in 1:length(depths)){

		boot_list <- list()
		print(paste(type, '_', depths[i], sep = ''))
		#for (j in 1:length(N_boots)){

		cl<-makeCluster(7)
		registerDoParallel(cl)

		boot_list <- foreach (j = 1:length(N_boots)) %dopar% {
			library("EBSeq")

			file = paste('CV/', type, '_', depths[i], '_', as.character(N_boots[j]), '.txt', sep = '')
			##print(file)
			x <- read.table(file,header=TRUE, row.names = 'GeneID')
			x <- as.matrix(x)
			group <- factor(type_factor)

			Sizes=MedianNorm(x)
			EBOut=EBTest(Data=x, Conditions=group,sizeFactors=Sizes, maxround=5)

			#boot_list[[j]] <- EBDERes
			##boot_list[[j]] <- EBOut
			EBOut

		}

		stopCluster(cl)
		for (j in 1:length(N_boots)){

			FC <- PostFC(EBOut)
			FC2 <- names(FC[[1]][log2(FC[[1]])>2])

			EBDERes_fold[[j]]= intersect(GetDEResults(boot_list[[j]], FDR=0.05, Threshold_FCRatio=0)$DEfound,FC2)
			EBDERes_all[[j]]=GetDEResults(boot_list[[j]], FDR=0.05, Threshold_FCRatio=0)$DEfound

		}

		depth_list_all[[i]] <- EBDERes_all
		depth_list_fold[[i]] <- EBDERes_fold

	}

	if (type == 'ER+'){
		ER_list_de_all <- depth_list_all
		ER_list_de_fold <- depth_list_fold
	}
	if (type == 'TNBC'){
		TNBC_list_de_all <- depth_list_all
		TNBC_list_de_fold <- depth_list_fold
	}

}

save.image('ebseq.RData')

de_ER_ebseq <- list()
de_TNBC_ebseq <- list()

de_ER_ebseq <- ER_list_de_all
save('de_ER_ebseq', file = 'de_ER_ebseq_0.05_none')

de_TNBC_ebseq <- TNBC_list_de_all
save('de_TNBC_ebseq', file = 'de_TNBC_ebseq_0.05_none')

de_ER_ebseq <- ER_list_de_fold
save('de_ER_ebseq', file = 'de_ER_ebseq_0.05_2')

de_TNBC_ebseq <- TNBC_list_de_fold
save('de_TNBC_ebseq', file = 'de_TNBC_ebseq_0.05_2')





#############################

load('ebseq.RData')

de_ER_ebseq_all <- list()
de_TNBC_ebseq_all <- list()
de_ER_ebseq_fold <- list()
de_TNBC_ebseq_fold <- list()


for (i in 1:length(depths)){

	ER_boot_list_all <- list()
	TNBC_boot_list_all <- list()
	ER_boot_list_fold <- list()
	TNBC_boot_list_fold <- list()

	for (j in 1:length(N_boots)){

		de_ER_all <- ER_list_de_all[[i]][[j]]$DEfound
		de_TNBC_all <- TNBC_list_de_all[[i]][[j]]$DEfound
		de_ER_no_fold <- ER_list_de_fold[[i]][[j]]$DEfound
		de_TNBC_fold <- TNBC_list_de_fold[[i]][[j]]$DEfound

		ER_boot_list_all[[j]] <- de_ER_all
		TNBC_boot_list_all[[j]] <- de_TNBC_all
		ER_boot_list_fold[[j]] <- de_ER_no_fold
		TNBC_boot_list_fold[[j]] <- de_TNBC_fold

	}

	de_ER_ebseq_all[[i]] <- ER_boot_list_all
	de_TNBC_ebseq_all[[i]] <- TNBC_boot_list_all
	de_ER_ebseq_fold[[i]] <- ER_boot_list_fold
	de_TNBC_ebseq_fold[[i]] <- TNBC_boot_list_fold

}



de_ER_ebseq <- list()

de_ER_ebseq <- de_ER_ebseq_all
save('de_ER_ebseq', file = 'de_ER_ebseq_0.05_none')

de_TNBC_ebseq <- de_TNBC_ebseq_all
save('de_TNBC_ebseq', file = 'de_TNBC_ebseq_0.05_none')

de_ER_ebseq <- de_ER_ebseq_fold
save('de_ER_ebseq', file = 'de_ER_ebseq_0.05_2')

de_TNBC_ebseq <- de_TNBC_ebseq_fold
save('de_TNBC_ebseq', file = 'de_TNBC_ebseq_0.05_2')


########################################################


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


library('org.Hs.eg.db')
library('goseq');
library('GO.db');

load('de_ER_ebseq_0.05_2')
load('de_TNBC_ebseq_0.05_2')

depths <- c('0.01','0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.85', '0.9', '0.95', '0.99', '1')
N_boots <- 1:25


g <- read.table('ER+_0.1_1.txt', header = TRUE)
gg <- g[,1]

zzz <- as.list(org.Hs.egALIAS2EG)
xx <- gg
yy <- gg
pp <- get_entrez(xx, zzz)
pp <- unique(pp)
qq <- get_entrez(yy, zzz)
qq <- unique(qq)


go_ER_ebseq <- list()
go_TNBC_ebseq <- list()

for (i in 1:length(depths)){

	ER_boot_list <- list()
	TNBC_boot_list <- list()

	for (j in 1:length(N_boots)){

		print(paste(as.character(i), '_', as.character(j), sep = ''))

		eg_ER <- de_ER_ebseq[[i]][[j]]
		eg_TNBC <- de_TNBC_ebseq[[i]][[j]]


		x <- eg_ER
		y <- eg_TNBC

		p <- get_entrez(x, zzz)
		p <- unique(p)

		q <- get_entrez(y, zzz)
		q <- unique(q)

		ppp <- matrix(ncol = 2, nrow = length(pp))
		ppp <- data.frame(ppp)
		ppp[,1] <- pp
		for(k in 1:length(pp)){
			if ((pp[k] %in% p) == TRUE){
				ppp[k,2] <- 1
			}
			else{
				ppp[k,2] <- 0
			}
		}
		names(ppp) <- c('entrezgene','diff.exprs')
		gene.set <- ppp
		goseq.set = as.vector(gene.set$diff.exprs);
		names(goseq.set) = gene.set$entrezgene;
		pwf = nullp(goseq.set, "hg19","knownGene");
		GO.wall = goseq(pwf,"hg19","knownGene");
		go_ER =GO.wall$category[ p.adjust(GO.wall$over_represented_pvalue, method="bonferroni")<0.05] ;

		qqq <- matrix(ncol = 2, nrow = length(qq))
		qqq <- data.frame(qqq)
		qqq[,1] <- qq
		for(k in 1:length(qq)){
			if ((qq[k] %in% q) == TRUE){
				qqq[k,2] <- 1
			}
			else{
				qqq[k,2] <- 0
			}
		}
		names(qqq) <- c('entrezgene','diff.exprs')
		gene.set <- qqq
		goseq.set = as.vector(gene.set$diff.exprs);
		names(goseq.set) = gene.set$entrezgene;
		pwf = nullp(goseq.set, "hg19","knownGene");
		GO.wall = goseq(pwf,"hg19","knownGene");
		go_TNBC =GO.wall$category[ p.adjust(GO.wall$over_represented_pvalue, method="bonferroni")<0.05] ;


		ER_boot_list[[j]] <- go_ER
		TNBC_boot_list[[j]] <- go_TNBC
	}
	go_ER_ebseq[[i]] <- ER_boot_list
	go_TNBC_ebseq[[i]] <- TNBC_boot_list
}

save('go_ER_ebseq', file = 'go_ER_ebseq_0.05_2')
save('go_TNBC_ebseq', file = 'go_TNBC_ebseq_0.05_2')








E <- GetDEResults(EBOut, FDR=0.05, Threshold_FCRatio=0.9)
x <- E$DEfound

		p <- get_entrez(x, zzz)
		p <- unique(p)

		q <- get_entrez(y, zzz)
		q <- unique(q)

		ppp <- matrix(ncol = 2, nrow = length(pp))
		ppp <- data.frame(ppp)
		ppp[,1] <- pp
		for(k in 1:length(pp)){
			if ((pp[k] %in% p) == TRUE){
				ppp[k,2] <- 1
			}
			else{
				ppp[k,2] <- 0
			}
		}
		names(ppp) <- c('entrezgene','diff.exprs')
		gene.set <- ppp
		goseq.set = as.vector(gene.set$diff.exprs);
		names(goseq.set) = gene.set$entrezgene;
		pwf = nullp(goseq.set, "hg19","knownGene");
		GO.wall = goseq(pwf,"hg19","knownGene");
		go_ER =GO.wall$category[ p.adjust(GO.wall$over_represented_pvalue, method="bonferroni")<0.05] ;














