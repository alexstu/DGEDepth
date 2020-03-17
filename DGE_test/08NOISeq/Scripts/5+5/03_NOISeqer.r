library("NOISeq")


types = c('ER+', 'TNBC')
depths <- c('0.01','0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.85', '0.9', '0.95', '0.99', '1')
N_boots <- 1:25

#types = c('ER+')
#depths = c(0.01, 0.05, 0.15, 1)
#N_boots <- 1:2

ER_list <- list()
TNBC_list <- list()





for (type in types){

	depth_list <- list()

	if (type == 'ER+'){
		type_factor <- c(1,1,1,1,1,2,2,2,2,2)
		sub_names <- c('SRR1313092.sam','SRR1313100.sam','SRR1313105.sam','SRR1313107.sam','SRR1313124.bam','SRR1313174.bam','SRR1313175.bam','SRR1313178.bam','SRR1313197.bam','SRR1313203.bam')
	}
	if (type == 'TNBC'){
		type_factor <- c(1,1,1,1,1,2,2,2,2,2)
		sub_names <- c('SRR1313141.bam','SRR1313142.bam','SRR1313147.bam','SRR1313152.bam','SRR1313164.bam','SRR1313217.bam','SRR1313220.bam','SRR1313226.bam','SRR1313227.bam','SRR1313229.bam')
	}



	for (i in 1:length(depths)){

		boot_list <- list()

		for (j in 1:length(N_boots)){

			file = paste('../../../02voomlimma/Scripts/5+5/CV/', type, '_', depths[i], '_', as.character(N_boots[j]), '.txt', sep = '')
			
			print(file)
			x <- read.table(file,header=TRUE, row.names = 'GeneID')

			x <- x[,sub_names]

			x <- as.matrix(x)
			y = data.frame(Group = type_factor)

			data <- readData(data = x, factors = y)

			 #uqua(assayData(data)$exprs, long = 1000, lc = 0)
			res <- noiseqbio(data, norm = 'uqua', factor = 'Group', conditions = c(1,0), lc = 0)

			#boot_list[[j]] <-degenes(res)
			boot_list[[j]] <- res

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


#save.image('NOISeq.RData')

#############################

load('NOISeq.RData')
library('NOISeq')

pval = 0.05
fold = 2

de_ER_noiseq_all <- list()
de_TNBC_noiseq_all <- list()
de_ER_noiseq_fold <- list()
de_TNBC_noiseq_fold <- list()


for (i in 1:length(depths)){

	ER_boot_list_all <- list()
	TNBC_boot_list_all <- list()
	ER_boot_list_fold <- list()
	TNBC_boot_list_fold <- list()

	for (j in 1:length(N_boots)){

		ER_res <- degenes(ER_list[[i]][[j]], q = 0.95)
		TNBC_res <- degenes(TNBC_list[[i]][[j]], q = 0.95)

		ER_all <- rownames(ER_res)
		TNBC_all <- rownames(TNBC_res)

		ER_FC <- rownames(ER_res[abs(ER_res$log2FC) > 2 & !is.na(ER_res$log2FC),])
		TNBC_FC <- rownames(TNBC_res[abs(TNBC_res$log2FC) > 2 & !is.na(TNBC_res$log2FC),])

		ER_boot_list_all[[j]] <- ER_all
		TNBC_boot_list_all[[j]] <- TNBC_all
		ER_boot_list_fold[[j]] <- ER_FC
		TNBC_boot_list_fold[[j]] <- TNBC_FC

	}

		de_ER_noiseq_all[[i]] <- ER_boot_list_all
		de_TNBC_noiseq_all[[i]] <- TNBC_boot_list_all
		de_ER_noiseq_fold[[i]] <- ER_boot_list_fold
		de_TNBC_noiseq_fold[[i]] <- TNBC_boot_list_fold

}


de_ER_noiseq <- list()
de_TNBC_noiseq <- list()

de_ER_noiseq <- de_ER_noiseq_all
save('de_ER_noiseq', file = 'sub_de_ER_noiseq_0.05_none')

de_TNBC_noiseq <- de_TNBC_noiseq_all
save('de_TNBC_noiseq', file = 'sub_de_TNBC_noiseq_0.05_none')

de_ER_noiseq <- de_ER_noiseq_fold
save('de_ER_noiseq', file = 'sub_de_ER_noiseq_0.05_2')

de_TNBC_noiseq <- de_TNBC_noiseq_fold
save('de_TNBC_noiseq', file = 'sub_de_TNBC_noiseq_0.05_2')


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

load('de_ER_noiseq_0.05_2')
load('de_TNBC_noiseq_0.05_2')

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


go_ER_noiseq <- list()
go_TNBC_noiseq <- list()

for (i in 1:length(depths)){

	ER_boot_list <- list()
	TNBC_boot_list <- list()

	for (j in 1:length(N_boots)){

		print(paste(as.character(i), '_', as.character(j), sep = ''))

		eg_ER <- de_ER_noiseq[[i]][[j]]
		eg_TNBC <- de_TNBC_noiseq[[i]][[j]]


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
	go_ER_noiseq[[i]] <- ER_boot_list
	go_TNBC_noiseq[[i]] <- TNBC_boot_list
}

save('go_ER_noiseq', file = 'go_ER_noiseq_0.05_2')
save('go_TNBC_noiseq', file = 'go_TNBC_noiseq_0.05_2')















