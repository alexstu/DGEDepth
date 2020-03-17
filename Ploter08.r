'GO N(VENN)~d'
library(ggplot2)
library(VennDiagram)



depths <- c('0.01','0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.85', '0.9', '0.95', '0.99', '1')
N_boots <- 1:25

filter = '2'
j = 1


load(paste('go_ER_deseq_0.05_',filter, sep = ''))
load(paste('go_ER_ebseq_0.05_',filter, sep = ''))
load(paste('go_ER_edger_0.05_',filter, sep = ''))
load(paste('go_ER_noiseq_0.05_',filter, sep = ''))
load(paste('go_ER_voom_0.05_',filter, sep = ''))
load(paste('go_TNBC_deseq_0.05_',filter, sep = ''))
load(paste('go_TNBC_ebseq_0.05_',filter, sep = ''))
load(paste('go_TNBC_edger_0.05_',filter, sep = ''))
load(paste('go_TNBC_noiseq_0.05_',filter, sep = ''))
load(paste('go_TNBC_voom_0.05_',filter, sep = ''))

#############################



venn_depth = c('0.01','0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.85', '0.9', '0.95', '0.99', '1')

for (d in venn_depth){

	i = match(d,depths)
	#print(j)
	venn_name = paste('GO_VENN_ER_fold_',filter, '_',d, '_',as.character(j),'.tiff', sep = '')
	col <- c('red', 'purple', 'green','darkorange','blue')
	m = 'ER+ dataset, f = '
	venn.diagram(list(DESeq=go_ER_deseq[[i]][[j]],VOOM=go_ER_voom[[i]][[j]],EdgeR=go_ER_edger[[i]][[j]],EBSeq=go_ER_ebseq[[i]][[j]],NOISeq=go_ER_noiseq[[i]][[j]]),fill = col, cat.cex = rep(1.5,5),main.cex = 3, main = paste(m,d,'\n '), filename = venn_name,rotation.degree = 029,cat.col = col,fontface = rep('bold', 31),cat.fontface = rep('bold',5),cat.just=list(c(1.9,3) , c(0.2,6.9) , c(1.25,-1.6) , c(1,1) , c(0.25,0)))

}


for (d in venn_depth){

	i = match(d,depths)
	venn_name = paste('GO_VENN_TNBC_fold_',filter, '_',d, '_',as.character(j),'.tiff', sep = '')
	col <- c('red', 'purple', 'green','darkorange','blue')
	m = 'TNBC dataset, f = '
	venn.diagram(list(DESeq=go_TNBC_deseq[[i]][[j]],VOOM=go_TNBC_voom[[i]][[j]],EdgeR=go_TNBC_edger[[i]][[j]],EBSeq=go_TNBC_ebseq[[i]][[j]],NOISeq=go_TNBC_noiseq[[i]][[j]]),fill = col, cat.cex = rep(1.5,5),main.cex = 3, main = paste(m,d,'\n '), filename = venn_name,rotation.degree = 029,cat.col = col,fontface = rep('bold', 31),cat.fontface = rep('bold',5),cat.just=list(c(1.9,3) , c(0.2,6.9) , c(1.25,-1.6) , c(1,1) , c(0.25,0)))

}






