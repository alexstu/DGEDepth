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

	depth_list <- list()

	if (type == 'ER+'){
		type_factor <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
		pattern_files_list <- c('SRR1313090.sam','SRR1313091.sam','SRR1313092.sam','SRR1313093.sam','SRR1313094.sam','SRR1313095.sam','SRR1313096.sam','SRR1313097.sam','SRR1313098.sam','SRR1313099.sam','SRR1313100.sam','SRR1313101.sam','SRR1313102.sam','SRR1313103.sam','SRR1313104.sam','SRR1313105.sam','SRR1313106.sam','SRR1313107.sam','SRR1313108.sam','SRR1313109.sam','SRR1313110.sam','SRR1313111.sam','SRR1313112.sam','SRR1313113.sam','SRR1313114.sam','SRR1313115.sam','SRR1313116.sam','SRR1313117.sam','SRR1313118.sam','SRR1313119.sam','SRR1313120.bam','SRR1313121.bam','SRR1313122.bam','SRR1313123.bam','SRR1313124.bam','SRR1313125.bam','SRR1313126.bam','SRR1313127.bam','SRR1313128.bam','SRR1313129.bam','SRR1313130.bam','SRR1313131.bam','SRR1313174.bam','SRR1313175.bam','SRR1313176.bam','SRR1313177.bam','SRR1313178.bam','SRR1313179.bam','SRR1313180.bam','SRR1313181.bam','SRR1313182.bam','SRR1313183.bam','SRR1313184.bam','SRR1313185.bam','SRR1313186.bam','SRR1313187.bam','SRR1313188.bam','SRR1313189.bam','SRR1313190.bam','SRR1313191.bam','SRR1313192.bam','SRR1313193.bam','SRR1313194.bam','SRR1313195.bam','SRR1313196.bam','SRR1313197.bam','SRR1313198.bam','SRR1313199.bam','SRR1313200.bam','SRR1313201.bam','SRR1313202.bam','SRR1313203.bam')
	}
	if (type == 'TNBC'){
		type_factor <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
		pattern_files_list <- c('SRR1313132.bam','SRR1313133.bam','SRR1313134.bam','SRR1313135.bam','SRR1313136.bam','SRR1313137.bam','SRR1313138.bam','SRR1313139.bam','SRR1313140.bam','SRR1313141.bam','SRR1313142.bam','SRR1313143.bam','SRR1313144.bam','SRR1313145.bam','SRR1313146.bam','SRR1313147.bam','SRR1313148.bam','SRR1313149.bam','SRR1313150.bam','SRR1313151.bam','SRR1313152.bam','SRR1313153.bam','SRR1313154.bam','SRR1313155.bam','SRR1313156.bam','SRR1313157.bam','SRR1313158.bam','SRR1313159.bam','SRR1313160.bam','SRR1313161.bam','SRR1313162.bam','SRR1313163.bam','SRR1313164.bam','SRR1313165.bam','SRR1313166.bam','SRR1313167.bam','SRR1313168.bam','SRR1313169.bam','SRR1313170.bam','SRR1313171.bam','SRR1313172.bam','SRR1313173.bam','SRR1313209.bam','SRR1313210.bam','SRR1313211.bam','SRR1313212.bam','SRR1313213.bam','SRR1313214.bam','SRR1313215.bam','SRR1313216.bam','SRR1313217.bam','SRR1313218.bam','SRR1313219.bam','SRR1313220.bam','SRR1313221.bam','SRR1313222.bam','SRR1313223.bam','SRR1313224.bam','SRR1313225.bam','SRR1313226.bam','SRR1313227.bam','SRR1313228.bam','SRR1313229.bam')
	}



	for (i in 1:length(depths)){

		boot_list <- list()

		for (j in 1:length(N_boots)){

			design = matrix(ncol = 2, nrow = length(pattern_files_list))
			design <- as.data.frame(design)
			row.names(design) <- pattern_files_list
			design[,1] <- rep(1,length(type_factor))
			design[,2] <- type_factor

			file = paste('CV/', type, '_', depths[i], '_', as.character(N_boots[j]), '.txt', sep = '')
			
			print(file)
			x <- read.table(file,header=TRUE, row.names = 'GeneID')
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
fold = 3

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
			ER_boot_list[[j]] <- rownames(de_ER[(de_ER$adj.P.Val < pval) & (de_ER$logFC > fold),])
			TNBC_boot_list[[j]] <- rownames(de_TNBC[(de_TNBC$adj.P.Val < pval) & (de_TNBC$logFC > fold),])
		}

	}

	de_ER_voom[[i]] <- ER_boot_list
	de_TNBC_voom[[i]] <- TNBC_boot_list

}

save('de_ER_voom', file = paste('de_ER_voom_',as.character(pval), '_', as.character(fold), sep = ''))
save('de_TNBC_voom', file = paste('de_TNBC_voom_',as.character(pval), '_', as.character(fold), sep = ''))


########################################################
###########################################################

library(samExploreR)

load('voom.RData')


#genes = rownames(ER_list[[17]][[1]][with(ER_list[[17]][[1]], order(log2FoldChange,-padj)), ])[1:500]

#ER_list[[17]][[1]][with(ER_list[[17]][[1]], order(log2FoldChange,-padj)), ]

for (j in 1:length(genes[200:250])){
	#j=68
	df <- matrix(ncol = 3, nrow = 25*17*5)
	colnames(df) = c('Label', 'Variable', 'Value')
	df <- as.data.frame(df)

	for (i in 1:length(depths)){

		df[j*((i-1)*length(N_boots)+(j-1))+1,1] <- paste(as.character(j), genes[j], sep = ', ')
		df[j*((i-1)*length(N_boots)+(j-1))+1,2] <- as.numeric(depths[i])
		df[j*((i-1)*length(N_boots)+(j-1))+1,3] <- TNBC_list[[i]][[1]][genes[j], 'padj']
		}

		plotsamExplorer(df, save = TRUE, p.depth = 0.9, anova = FALSE)
		system('sleep 4')

}

#for (k in 1:(length(depths)-1)){
	k=1
#	for (j in 1:25){
	j=1
		genes = vector()
		x <- rownames(ER_list[[k]][[j]][(ER_list[[k]][[j]]$padj < 0.05)&(abs(ER_list[[k]][[j]]$log2FoldChange) > 0),])
		y <- rownames(ER_list[[17]][[1]][(ER_list[[17]][[1]]$padj < 0.05)&(abs(ER_list[[17]][[1]]$log2FoldChange) > 0),])

		z <- intersect(x,y)

		for (g in x){

			if ((g %in% z) == TRUE){
				next
			}

			if ((g %in% y) == TRUE){
				next
			}

			genes <- c(genes, g)

		}

#	}
#}
genes


xxxx <- function(g,j){

	pv <- vector()
	df <- matrix(ncol = 2, nrow = length(depths))
	df <- as.data.frame(df)
	colnames(df) <- c('f', 'p_value')

	for (i in 1:length(depths)){

		#df[j*((i-1)*length(N_boots)+(j-1))+1,1] <- paste(as.character(j), genes[j], sep = ', ')
		#df[j*((i-1)*length(N_boots)+(j-1))+1,2] <- as.numeric(depths[i])
		#pv <- c(pv,TNBC_list[[i]][[1]][genes[j], 'padj'])
		df[i,1] <- as.numeric(depths[i])
		df[i,2] <- TNBC_list[[i]][[j]][g, 'adj.P.Val']
		}
	#ggplot(df, aes(x=d, y=pv)) + geom_point(aes(colour='blue',size=1, guide = FALSE)) + geom_hline(yintercept=0.05, linetype="dotted", size = 0.7)+ theme(legend.position="none") + geom_line(aes( colour = 'red' ),se = F, size = 0.7)
	ggplot(df, aes(x=f, y=p_value)) + geom_point(aes(colour='blue',size=2, guide = FALSE)) + geom_hline(yintercept=0.05, linetype="dotted", size = 0.7)+ theme(legend.position="none") + geom_line(aes( colour = 'red'),se = F, size = 1.3)+ scale_y_continuous(limits = c(0, 1)) + theme(text = element_text(size=25, face="bold")) + ggtitle(paste('Voom for ', g, ' gene', sep = ''))
}








#xxxx('RMDN2',1)
#xxxx('GRR75-ASB3',1)




j=1
pdf(paste('RMDN2-voom',as.character(j),'.pdf', sep = ''))
xxxx('RMDN2',j)
dev.off()

j=15
pdf(paste('PIK3CD-voom',as.character(j),'.pdf', sep = ''))
xxxx('PIK3CD',j)
dev.off()

################################################



load('voom.RData')



std <- function(x) sd(x)/sqrt(length(x))

xxxxx <- function(g){

	depths <- c('0.01','0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.85', '0.9', '0.95', '0.99', '1')
	N_boots <- 1:25
	res_table = matrix(ncol = 3, nrow = length(depths))

	for (i in 1:length(depths)){

		dd = as.numeric(depths[i])
		cur_vector_pval = vector()
		for (j in N_boots){
			cur_vector_pval <- c(cur_vector_pval, TNBC_list[[i]][[j]][g, 'adj.P.Val'])
		}

		res_table[i,1] <- dd
		res_table[i,2] <- (mean(cur_vector_pval))
		res_table[i,3] <- (std(cur_vector_pval))

	}

print(res_table)
epsilon = 0.02

header = paste('Voom for ', g, sep = '')
#header = paste(group, ':', type, '\n', 'Expression > 1',  sep = '')
up = res_table[,2] + res_table[,3]
low = res_table[,2] - res_table[,3]
#png(paste("04_Mean_expressed_genes_", type, '_',  group, ".png", sep = ""))

par(mar=c(6, 5, 4, 2))
plot(res_table[,1], res_table[,2], main = header, font.main = 2, xlab = "f", ylab = "Mean", font.lab = 2, cex.main = 2, cex.axis = 1.5, cex.lab = 2, col = 'red', pch = 16, cex = 2, ylim = c(min(low), max(up)))
for(i in 1:length(res_table[,2])) {
	segments(res_table[i,1],low[i] , res_table[i,1], up[i], col = 'blue', lwd = 2)
	segments(res_table[i,1]-epsilon, up[i] , res_table[i,1]+epsilon, up[i], col = 'blue', lwd = 2)
	segments(res_table[i,1]-epsilon, low[i] , res_table[i,1]+epsilon, low[i], col = 'blue', lwd = 2)
}
abline(h = 0.05, lty=2, lwd=2)
#dev.off()




	#ggplot(df, aes(x=f, y=p_value)) + geom_point(aes(colour='blue',size=2, guide = FALSE)) + geom_hline(yintercept=0.05, linetype="dotted", size = 0.7)+ theme(legend.position="none") + geom_line(aes( colour = 'red'),se = F, size = 1.3)+ scale_y_continuous(limits = c(0, 1)) + theme(text = element_text(size=25, face="bold")) + ggtitle(paste('DESeq2 for ', g, ' gene', sep = ''))
}


pdf(paste('PIK3CD-voom','all','.pdf', sep = ''))
xxxxx('PIK3CD')
dev.off()


pdf(paste('RMDN2-voom','all','.pdf', sep = ''))
xxxxx('RMDN2')
dev.off()





