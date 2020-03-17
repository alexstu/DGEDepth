'All N(VENN)~d'
library(ggplot2)
library(VennDiagram)



depths <- c('0.01','0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.85', '0.9', '0.95', '0.99', '1')
N_boots <- 1:25
filter = '2'
j = 1

if(filter == 'none'){
	m = 'No fold filtering, f = '
}
if(filter == '2'){
	m = '2-fold filtering, f = '
}


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

de_ER_deseq <-de_ER_deseq[[1]]
de_ER_voom <-de_ER_voom[[1]]
de_ER_edger <-de_ER_edger[[1]]
de_ER_ebseq <-de_ER_ebseq[[1]]
de_ER_noiseq <-de_ER_noiseq[[1]]
de_TNBC_deseq <-de_TNBC_deseq[[1]]
de_TNBC_voom <-de_TNBC_voom[[1]]
de_TNBC_edger <-de_TNBC_edger[[1]]
de_TNBC_ebseq <-de_TNBC_ebseq[[1]]
de_TNBC_noiseq <-de_TNBC_noiseq[[1]]

#############################



venn_depth = c('0.01','0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.85', '0.9', '0.95', '0.99', '1')


for (d in venn_depth){

	i = match(d,depths)
	#print(j)
	venn_name = paste('VENN_ER_fold_',filter, '_',d, '_',as.character(j),'.tiff', sep = '')
	col <- c('red', 'purple', 'green','darkorange','blue')
venn.diagram(list(DESeq=de_ER_deseq[[i]][[j]],VOOM=de_ER_voom[[i]][[j]],EdgeR=de_ER_edger[[i]][[j]],EBSeq=de_ER_ebseq[[i]][[j]],NOISeq=de_ER_noiseq[[i]][[j]]),fill = col, cat.cex = rep(1.5,5),main.cex = 3, main = paste(m,d,'\n '), filename = venn_name,rotation.degree = 029,cat.col = col,fontface = rep('bold', 31),cat.fontface = rep('bold',5),cat.just=list(c(1.9,3) , c(0.2,6.9) , c(1.25,-1.6) , c(1,1) , c(0.25,0)))

}


for (d in venn_depth){

	i = match(d,depths)
	#print(j)
	venn_name = paste('VENN_TNBC_fold_',filter, '_',d, '_',as.character(j),'.tiff', sep = '')
	col <- c('red', 'purple', 'green','darkorange','blue')
venn.diagram(list(DESeq=de_TNBC_deseq[[i]][[j]],VOOM=de_TNBC_voom[[i]][[j]],EdgeR=de_TNBC_edger[[i]][[j]],EBSeq=de_TNBC_ebseq[[i]][[j]],NOISeq=de_TNBC_noiseq[[i]][[j]]),fill = col, cat.cex = rep(1.5,5),main.cex = 3, main = paste(m,d,'\n '), filename = venn_name,rotation.degree = 029,cat.col = col,fontface = rep('bold', 31),cat.fontface = rep('bold',5),cat.just=list(c(1.9,3) , c(0.2,6.9) , c(1.25,-1.6) , c(1,1) , c(0.25,0)))

}





