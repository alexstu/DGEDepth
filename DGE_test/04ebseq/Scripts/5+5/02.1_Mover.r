
samples_er <- c('SRR1313090.sam','SRR1313091.sam','SRR1313092.sam','SRR1313093.sam','SRR1313094.sam','SRR1313095.sam','SRR1313096.sam','SRR1313097.sam','SRR1313098.sam','SRR1313099.sam','SRR1313100.sam','SRR1313101.sam','SRR1313102.sam','SRR1313103.sam','SRR1313104.sam','SRR1313105.sam','SRR1313106.sam','SRR1313107.sam','SRR1313108.sam','SRR1313109.sam','SRR1313110.sam','SRR1313111.sam','SRR1313112.sam','SRR1313113.sam','SRR1313114.sam','SRR1313115.sam','SRR1313116.sam','SRR1313117.sam','SRR1313118.sam','SRR1313119.sam','SRR1313120.bam','SRR1313121.bam','SRR1313122.bam','SRR1313123.bam','SRR1313124.bam','SRR1313125.bam','SRR1313126.bam','SRR1313127.bam','SRR1313128.bam','SRR1313129.bam','SRR1313130.bam','SRR1313131.bam')

samples_er_free <- c('SRR1313174.bam','SRR1313175.bam','SRR1313176.bam','SRR1313177.bam','SRR1313178.bam','SRR1313179.bam','SRR1313180.bam','SRR1313181.bam','SRR1313182.bam','SRR1313183.bam','SRR1313184.bam','SRR1313185.bam','SRR1313186.bam','SRR1313187.bam','SRR1313188.bam','SRR1313189.bam','SRR1313190.bam','SRR1313191.bam','SRR1313192.bam','SRR1313193.bam','SRR1313194.bam','SRR1313195.bam','SRR1313196.bam','SRR1313197.bam','SRR1313198.bam','SRR1313199.bam','SRR1313200.bam','SRR1313201.bam','SRR1313202.bam','SRR1313203.bam')

samples_tnbc = c('SRR1313132.bam','SRR1313133.bam','SRR1313134.bam','SRR1313135.bam','SRR1313136.bam','SRR1313137.bam','SRR1313138.bam','SRR1313139.bam','SRR1313140.bam','SRR1313141.bam','SRR1313142.bam','SRR1313143.bam','SRR1313144.bam','SRR1313145.bam','SRR1313146.bam','SRR1313147.bam','SRR1313148.bam','SRR1313149.bam','SRR1313150.bam','SRR1313151.bam','SRR1313152.bam','SRR1313153.bam','SRR1313154.bam','SRR1313155.bam','SRR1313156.bam','SRR1313157.bam','SRR1313158.bam','SRR1313159.bam','SRR1313160.bam','SRR1313161.bam','SRR1313162.bam','SRR1313163.bam','SRR1313164.bam','SRR1313165.bam','SRR1313166.bam','SRR1313167.bam','SRR1313168.bam','SRR1313169.bam','SRR1313170.bam','SRR1313171.bam','SRR1313172.bam','SRR1313173.bam')


samples_tnbc_free = c('SRR1313209.bam','SRR1313210.bam','SRR1313211.bam','SRR1313212.bam','SRR1313213.bam','SRR1313214.bam','SRR1313215.bam','SRR1313216.bam','SRR1313217.bam','SRR1313218.bam','SRR1313219.bam','SRR1313220.bam','SRR1313221.bam','SRR1313222.bam','SRR1313223.bam','SRR1313224.bam','SRR1313225.bam','SRR1313226.bam','SRR1313227.bam','SRR1313228.bam','SRR1313229.bam')







depths <- c('0.01','0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8','0.85', '0.9', '0.95', '0.99', '1')
N_boots <- 1:25
N_genes <- 1:57341

all_pattern_er <- c(samples_er, samples_er_free)
all_pattern_tnbc <- c(samples_tnbc, samples_tnbc_free)




for (i in 1:length(depths)){
	for (j in 1:length(N_boots)){
		dta <- matrix(ncol = (length(all_pattern_er)+1), nrow = length(N_genes))
		dta <- as.data.frame(dta)
		dim(dta) <- c()
		cnames <- vector()
		for (k in 1:length(all_pattern_er)){
			file <- paste(all_pattern_er[k], '_', as.character(depths[i]), '_', as.character(N_boots[j]), '.txt', sep = '')
			cv <- read.table(file, header = TRUE)
			#file.remove(file)
			dta[,1] <- cv[,1]
			dta[,k+1] <- cv[,2]
			cnames[1] <- 'GeneID'
			cnames[k+1] <- all_pattern_er[k]
		}
		colnames(dta) <- cnames
		write.table(dta, paste('ER+', '_', as.character(depths[i]), '_', as.character(N_boots[j]), '.txt', sep = ''), quote = FALSE, col.names = TRUE, row.names = FALSE)
	}
}



for (i in 1:length(depths)){
	for (j in 1:length(N_boots)){
		dta <- matrix(ncol = (length(all_pattern_tnbc)+1), nrow = length(N_genes))
		dta <- as.data.frame(dta)
		dim(dta) <- c()
		cnames <- vector()
		for (k in 1:length(all_pattern_tnbc)){
			file <- paste(all_pattern_tnbc[k], '_', as.character(depths[i]), '_', as.character(N_boots[j]), '.txt', sep = '')
			cv <- read.table(file, header = TRUE)
			#file.remove(file)
			dta[,1] <- cv[,1]
			dta[,k+1] <- cv[,2]
			cnames[1] <- 'GeneID'
			cnames[k+1] <- all_pattern_tnbc[k]
		}
		colnames(dta) <- cnames
		write.table(dta, paste('TNBC', '_', as.character(depths[i]), '_', as.character(N_boots[j]), '.txt', sep = ''), quote = FALSE, col.names = TRUE, row.names = FALSE)
	}
}





