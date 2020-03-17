library("limma")


types = c('ER+', 'TNBC')
depths <- c('0.01','0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.85', '0.9', '0.95', '0.99', '1')
N_boots <- 1:25

#types = c('ER+')
#depths = c(0.01, 0.05, 0.15, 1)
#N_boots <- 1:2

ER_list <- list()
TNBC_list <- list()





for (type in types){

	depth_list<- list()
	#depth_list_fold <- list()

	if (type == 'ER+'){
		type_factor <- c(0,0,0,0,0,1,1,1,1,1)
		sub_names <- c('SRR1313092.sam','SRR1313100.sam','SRR1313105.sam','SRR1313107.sam','SRR1313124.bam','SRR1313174.bam','SRR1313175.bam','SRR1313178.bam','SRR1313197.bam','SRR1313203.bam')
	}
	if (type == 'TNBC'){
		type_factor <- c(0,0,0,0,0,1,1,1,1,1)
		sub_names <- c('SRR1313141.bam','SRR1313142.bam','SRR1313147.bam','SRR1313152.bam','SRR1313164.bam','SRR1313217.bam','SRR1313220.bam','SRR1313226.bam','SRR1313227.bam','SRR1313229.bam')
	}



	for (i in 1:length(depths)){

		boot_list <- list()

		for (j in 1:length(N_boots)){

			design = matrix(ncol = 2, nrow = length(sub_names))
			design <- as.data.frame(design)
			row.names(design) <- sub_names
			design[,1] <- rep(1,length(type_factor))
			design[,2] <- type_factor

			file = paste('../../../02voomlimma/Scripts/5+5/CV/', type, '_', depths[i], '_', as.character(N_boots[j]), '.txt', sep = '')
			
			print(file)
			x <- read.table(file,header=TRUE, row.names = 'GeneID')

			x <- x[,sub_names]

			x <- as.matrix(x)

			v <- voom(x,design,plot=TRUE,normalize="quantile")
			fit <- lmFit(v,design)
			fit <- eBayes(fit)
			

			boot_list[[j]] <- topTable(fit, n = dim(fit)[1],coef = 'V2', sort.by = 'none', adjust="BH")

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


save.image('voom.RData')

#############################

load('voom.RData')


pval = 0.05
#fold = 'none'
fold = 2

de_ER_voom <- list()
de_TNBC_voom <- list()

for (i in 1:length(depths)){

	ER_boot_list <- list()
	TNBC_boot_list <- list()

	for (j in 1:length(N_boots)){

		de_ER <- ER_list[[i]][[j]]
		de_TNBC <- TNBC_list[[i]][[j]]

		if (fold == 'none'){
			ER_boot_list[[j]] <- rownames(de_ER[(de_ER$adj.P.Val < pval),])
			TNBC_boot_list[[j]] <- rownames(de_TNBC[(de_TNBC$adj.P.Val < pval),])
		}
		else {
			ER_boot_list[[j]] <- rownames(de_ER[(de_ER$adj.P.Val < pval) & (abs(de_ER$logFC) > fold),])
			TNBC_boot_list[[j]] <- rownames(de_TNBC[(de_TNBC$adj.P.Val < pval) & (abs(de_TNBC$logFC) > fold),])
		}

	}

	de_ER_voom[[i]] <- ER_boot_list
	de_TNBC_voom[[i]] <- TNBC_boot_list

}

save('de_ER_voom', file = paste('sub_de_ER_voom_',as.character(pval), '_', as.character(fold), sep = ''))
save('de_TNBC_voom', file = paste('sub_de_TNBC_voom_',as.character(pval), '_', as.character(fold), sep = ''))


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

load('voom.RData')

pval = 0.05
fold = 2


zzz <- as.list(org.Hs.egALIAS2EG)
xx <- rownames(ER_list[[1]][[1]])
yy <- rownames(TNBC_list[[1]][[1]])
pp <- get_entrez(xx, zzz)
pp <- unique(pp)
qq <- get_entrez(yy, zzz)
qq <- unique(qq)


go_ER_voom <- list()
go_TNBC_voom <- list()

for (i in 1:length(depths)){

	ER_boot_list <- list()
	TNBC_boot_list <- list()

	for (j in 1:length(N_boots)){

		print(paste(as.character(i), '_', as.character(j), sep = ''))

		eg_ER <- ER_list[[i]][[j]]
		eg_TNBC <- TNBC_list[[i]][[j]]


		x <- rownames(eg_ER[(eg_ER$adj.P.Val < pval) & (abs(eg_ER$logFC) > fold),])
		y <- rownames(eg_TNBC[(eg_TNBC$adj.P.Val < pval) & (abs(eg_TNBC$logFC) > fold),])

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
		pwf = nullp(goseq.set, "hg19","knownGene",plot.fit=FALSE);
		pwff <- pwf[!is.na(pwf[,3]),]
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
	go_ER_voom[[i]] <- ER_boot_list
	go_TNBC_voom[[i]] <- TNBC_boot_list
}

save('go_ER_voom', file = paste('go_ER_voom_',as.character(pval), '_', as.character(fold), sep = ''))
save('go_TNBC_voom', file = paste('go_TNBC_voom_',as.character(pval), '_', as.character(fold), sep = ''))

###############################



