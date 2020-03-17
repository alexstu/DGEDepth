typ = 'ER+'


header = 'Sample Name\tFile\tCondition'
samples = [
['GSM1401676','SRR1313090.sam','ER+'],
['GSM1401677','SRR1313091.sam','ER+'],
['GSM1401678','SRR1313092.sam','ER+'],
['GSM1401679','SRR1313093.sam','ER+'],
['GSM1401680','SRR1313094.sam','ER+'],
['GSM1401681','SRR1313095.sam','ER+'],
['GSM1401682','SRR1313096.sam','ER+'],
['GSM1401683','SRR1313097.sam','ER+'],
['GSM1401684','SRR1313098.sam','ER+'],
['GSM1401685','SRR1313099.sam','ER+'],
['GSM1401686','SRR1313100.sam','ER+'],
['GSM1401687','SRR1313101.sam','ER+'],
['GSM1401688','SRR1313102.sam','ER+'],
['GSM1401689','SRR1313103.sam','ER+'],
['GSM1401690','SRR1313104.sam','ER+'],
['GSM1401691','SRR1313105.sam','ER+'],
['GSM1401692','SRR1313106.sam','ER+'],
['GSM1401693','SRR1313107.sam','ER+'],
['GSM1401694','SRR1313108.sam','ER+'],
['GSM1401695','SRR1313109.sam','ER+'],
['GSM1401696','SRR1313110.sam','ER+'],
['GSM1401697','SRR1313111.sam','ER+'],
['GSM1401698','SRR1313112.sam','ER+'],
['GSM1401699','SRR1313113.sam','ER+'],
['GSM1401700','SRR1313114.sam','ER+'],
['GSM1401701','SRR1313115.sam','ER+'],
['GSM1401702','SRR1313116.sam','ER+'],
['GSM1401703','SRR1313117.sam','ER+'],
['GSM1401704','SRR1313118.sam','ER+'],
['GSM1401705','SRR1313119.sam','ER+'],
['GSM1401706','SRR1313120.bam','ER+'],
['GSM1401707','SRR1313121.bam','ER+'],
['GSM1401708','SRR1313122.bam','ER+'],
['GSM1401709','SRR1313123.bam','ER+'],
['GSM1401710','SRR1313124.bam','ER+'],
['GSM1401711','SRR1313125.bam','ER+'],
['GSM1401712','SRR1313126.bam','ER+'],
['GSM1401713','SRR1313127.bam','ER+'],
['GSM1401714','SRR1313128.bam','ER+'],
['GSM1401715','SRR1313129.bam','ER+'],
['GSM1401716','SRR1313130.bam','ER+'],
['GSM1401717','SRR1313131.bam','ER+'],
['GSM1401760','SRR1313174.bam','ER+-free'],
['GSM1401761','SRR1313175.bam','ER+-free'],
['GSM1401762','SRR1313176.bam','ER+-free'],
['GSM1401763','SRR1313177.bam','ER+-free'],
['GSM1401764','SRR1313178.bam','ER+-free'],
['GSM1401765','SRR1313179.bam','ER+-free'],
['GSM1401766','SRR1313180.bam','ER+-free'],
['GSM1401767','SRR1313181.bam','ER+-free'],
['GSM1401768','SRR1313182.bam','ER+-free'],
['GSM1401769','SRR1313183.bam','ER+-free'],
['GSM1401770','SRR1313184.bam','ER+-free'],
['GSM1401771','SRR1313185.bam','ER+-free'],
['GSM1401772','SRR1313186.bam','ER+-free'],
['GSM1401773','SRR1313187.bam','ER+-free'],
['GSM1401774','SRR1313188.bam','ER+-free'],
['GSM1401775','SRR1313189.bam','ER+-free'],
['GSM1401776','SRR1313190.bam','ER+-free'],
['GSM1401777','SRR1313191.bam','ER+-free'],
['GSM1401778','SRR1313192.bam','ER+-free'],
['GSM1401779','SRR1313193.bam','ER+-free'],
['GSM1401780','SRR1313194.bam','ER+-free'],
['GSM1401781','SRR1313195.bam','ER+-free'],
['GSM1401782','SRR1313196.bam','ER+-free'],
['GSM1401783','SRR1313197.bam','ER+-free'],
['GSM1401784','SRR1313198.bam','ER+-free'],
['GSM1401785','SRR1313199.bam','ER+-free'],
['GSM1401786','SRR1313200.bam','ER+-free'],
['GSM1401787','SRR1313201.bam','ER+-free'],
['GSM1401788','SRR1313202.bam','ER+-free'],
['GSM1401789','SRR1313203.bam','ER+-free'],

]

depths = ['0.01','0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8','0.85', '0.9', '0.95', '0.99', '1']

for n in range(25):
	N = n+1
	for d in depths:
		file_name = typ + '_metadata_' + d + '_' + str(N) + '.txt'
		fd = open(file_name, 'wb')
		fd.write('Sample Name\tFile\tCondition\n')
		for i in range(len(samples)):
			line = samples[i][0] + '\t' + samples[i][1] + '_' + d + '_' + str(N) + '.txt\t' + samples[i][2] + '\n'
			fd.write(line)
		fd.close()


