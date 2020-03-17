#N of total expressed genes ~ depth (average for)

library('e1071')
library('xtable')

group <- 'TNBC-free'
#group <- 'TNBC'

type <- 'Normalized'
#type <- 'Raw'


templates_1 <- c('SRR1313137_short.sam', 'SRR1313135_short.sam', 'SRR1313134_short.sam', 'SRR1313133_short.sam')
templates_2 <- c('SRR1313211_short.sam', 'SRR1313214_short.sam', 'SRR1313219_short.sam', 'SRR1313220_short.sam')







if (group == 'TNBC'){
templates <- templates_1
}
if (group == 'TNBC-free'){
templates <- templates_2
}



if (type == 'Normalized'){

	y_lim_M <- c(40,45)
	y_lim_K <- c(0,2500)
	y_lim_V <- c(15000,35000)
	y_lim_S <- c(10,35)

}

if (type == 'Raw'){

	y_lim_M <- c(0,2000)
	y_lim_K <- c(0,2500)
	y_lim_V <- c(0,6e7)
	y_lim_S <- c(10,35)

}


resample_D <- 10
res_table = matrix(ncol = 9, nrow = resample_D)
colnames(res_table) <- c('Depth', 'Mean', 'Mean_SD', 'Variance', 'Variance_SD', 'Skewness', 'Skewness_SD', 'Kurtosis', 'Kurtosis_SD')

setwd('CM')
for (dd in 1:(resample_D-1)){

	cur_vector_mean = vector()
	cur_vector_var = vector()
	cur_vector_sk = vector()
	cur_vector_kurt = vector()
	d <- dd/10

	for (i in 1:length(templates)){
		for (N in 1:24){

			cur_table_name <- paste(templates[i], '_', as.character(d), '_', as.character(N), '.txt', sep = '')
			cur_table <- read.table(cur_table_name, header = TRUE)

			if (type == 'Normalized'){
				cur_table[,2] <- cur_table[,2]/sum(cur_table[,2])*1000000
			}

			cur_vector_mean <- c(cur_vector_mean, mean(cur_table[,2]))
			cur_vector_var <- c(cur_vector_var, var(cur_table[,2]))
			cur_vector_sk <- c(cur_vector_sk, skewness(cur_table[,2]))
			cur_vector_kurt <- c(cur_vector_kurt, kurtosis(cur_table[,2]))

			#setEPS()
			#postscript(paste("02_All_expressed_genes_", cur_table_name, "_.eps", sep = ""))
		}
	}

	res_table[dd,1] <- as.numeric(d)
	res_table[dd,2] <- (mean(cur_vector_mean))
	res_table[dd,3] <- (sd(cur_vector_mean))
	res_table[dd,4] <- (mean(cur_vector_var))
	res_table[dd,5] <- (sd(cur_vector_var))
	res_table[dd,6] <- (mean(cur_vector_sk))
	res_table[dd,7] <- (sd(cur_vector_sk))
	res_table[dd,8] <- (mean(cur_vector_kurt))
	res_table[dd,9] <- (sd(cur_vector_kurt))

}

setwd("../")
setwd("Result")
epsilon = 0.02


#header = paste("Mean_for_",type, '_', group,  sep = '')
header = paste(group, ':', type, '\n', 'Expression > 1',  sep = '')
up = res_table[,2] + res_table[,3]
low = res_table[,2] - res_table[,3]
png(paste("04_Mean_expressed_genes_", type, '_',  group, ".png", sep = ""))
par(mar=c(6, 5, 4, 2))
plot(res_table[,1], res_table[,2], main = header, font.main = 2, xlab = "Depth", ylab = "Mean", ylim = y_lim_M, font.lab = 2, cex.main = 2, cex.axis = 1.5, cex.lab = 2, col = 'red', pch = 16, cex = 2)
for(i in 1:length(res_table[,2])) {
	segments(res_table[i,1],low[i] , res_table[i,1], up[i], col = 'blue', lwd = 2)
	segments(res_table[i,1]-epsilon, up[i] , res_table[i,1]+epsilon, up[i], col = 'blue', lwd = 2)
	segments(res_table[i,1]-epsilon, low[i] , res_table[i,1]+epsilon, low[i], col = 'blue', lwd = 2)
}
dev.off()

#header = paste("Variance_for_",type, '_', group,  sep = '')
header = paste(group, ':', type, '\n', 'Expression > 1',  sep = '')
up = res_table[,4] + res_table[,5]
low = res_table[,4] - res_table[,5]
png(paste("04_Variance_expressed_genes_", type, '_', group, ".png", sep = ""))
par(mar=c(6, 5, 4, 2))
plot(res_table[,1], res_table[,4], main = header, font.main = 2, xlab = "Depth", ylab = "Variance", ylim = y_lim_V, font.lab = 2, cex.main = 2, cex.axis = 1.5, cex.lab = 2, col = 'red', pch = 16, cex = 2)
for(i in 1:length(res_table[,2])) {
	segments(res_table[i,1],low[i] , res_table[i,1], up[i], col = 'blue', lwd = 2)
	segments(res_table[i,1]-epsilon, up[i] , res_table[i,1]+epsilon, up[i], col = 'blue', lwd = 2)
	segments(res_table[i,1]-epsilon, low[i] , res_table[i,1]+epsilon, low[i], col = 'blue', lwd = 2)
}
dev.off()



#header = paste("Skewness_for_",type, '_', group,  sep = '')
header = paste(group, ':', type, '\n', 'Expression > 1',  sep = '')
up = res_table[,6] + res_table[,7]
low = res_table[,6] - res_table[,7]
png(paste("04_Skewness_expressed_genes_", type, '_', group, ".png", sep = ""))
par(mar=c(6, 5, 4, 2))
plot(res_table[,1], res_table[,6], main = header, font.main = 2, xlab = "Depth", ylab = "Skewness", ylim = y_lim_S, font.lab = 2, cex.main = 2, cex.axis = 1.5, cex.lab = 2, col = 'red', pch = 16, cex = 2)
for(i in 1:length(res_table[,2])) {
	segments(res_table[i,1],low[i] , res_table[i,1], up[i], col = 'blue', lwd = 2)
	segments(res_table[i,1]-epsilon, up[i] , res_table[i,1]+epsilon, up[i], col = 'blue', lwd = 2)
	segments(res_table[i,1]-epsilon, low[i] , res_table[i,1]+epsilon, low[i], col = 'blue', lwd = 2)
}
dev.off()



#header = paste("Kurtosis_for_",type, '_', group,  sep = '')
header = paste(group, ':', type, '\n', 'Expression > 1',  sep = '')
up = res_table[,8] + res_table[,9]
low = res_table[,8] - res_table[,9]
png(paste("04_Kurtosis_expressed_genes_", type, '_',  group, ".png", sep = ""))
par(mar=c(6, 5, 4, 2))
plot(res_table[,1], res_table[,8], main = header, font.main = 2, xlab = "Depth", ylab = "Kurtosis", ylim = y_lim_K, font.lab = 2, cex.main = 2, cex.axis = 1.5, cex.lab = 2, col = 'red', pch = 16, cex = 2)
for(i in 1:length(res_table[,2])) {
	segments(res_table[i,1],low[i] , res_table[i,1], up[i], col = 'blue', lwd = 2)
	segments(res_table[i,1]-epsilon, up[i] , res_table[i,1]+epsilon, up[i], col = 'blue', lwd = 2)
	segments(res_table[i,1]-epsilon, low[i] , res_table[i,1]+epsilon, low[i], col = 'blue', lwd = 2)
}
dev.off()




setwd("../")







