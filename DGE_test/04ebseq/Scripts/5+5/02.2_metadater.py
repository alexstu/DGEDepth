typ = 'TNBC'


header = 'Sample Name\tFile\tCondition'
samples = [
['GSM1401718','SRR1313132.bam','TNBC'],
['GSM1401719','SRR1313133.bam','TNBC'],
['GSM1401720','SRR1313134.bam','TNBC'],
['GSM1401721','SRR1313135.bam','TNBC'],
['GSM1401722','SRR1313136.bam','TNBC'],
['GSM1401723','SRR1313137.bam','TNBC'],
['GSM1401724','SRR1313138.bam','TNBC'],
['GSM1401725','SRR1313139.bam','TNBC'],
['GSM1401726','SRR1313140.bam','TNBC'],
['GSM1401727','SRR1313141.bam','TNBC'],
['GSM1401728','SRR1313142.bam','TNBC'],
['GSM1401729','SRR1313143.bam','TNBC'],
['GSM1401730','SRR1313144.bam','TNBC'],
['GSM1401731','SRR1313145.bam','TNBC'],
['GSM1401732','SRR1313146.bam','TNBC'],
['GSM1401733','SRR1313147.bam','TNBC'],
['GSM1401734','SRR1313148.bam','TNBC'],
['GSM1401735','SRR1313149.bam','TNBC'],
['GSM1401736','SRR1313150.bam','TNBC'],
['GSM1401737','SRR1313151.bam','TNBC'],
['GSM1401738','SRR1313152.bam','TNBC'],
['GSM1401739','SRR1313153.bam','TNBC'],
['GSM1401740','SRR1313154.bam','TNBC'],
['GSM1401741','SRR1313155.bam','TNBC'],
['GSM1401742','SRR1313156.bam','TNBC'],
['GSM1401743','SRR1313157.bam','TNBC'],
['GSM1401744','SRR1313158.bam','TNBC'],
['GSM1401745','SRR1313159.bam','TNBC'],
['GSM1401746','SRR1313160.bam','TNBC'],
['GSM1401747','SRR1313161.bam','TNBC'],
['GSM1401748','SRR1313162.bam','TNBC'],
['GSM1401749','SRR1313163.bam','TNBC'],
['GSM1401750','SRR1313164.bam','TNBC'],
['GSM1401751','SRR1313165.bam','TNBC'],
['GSM1401752','SRR1313166.bam','TNBC'],
['GSM1401753','SRR1313167.bam','TNBC'],
['GSM1401754','SRR1313168.bam','TNBC'],
['GSM1401755','SRR1313169.bam','TNBC'],
['GSM1401756','SRR1313170.bam','TNBC'],
['GSM1401757','SRR1313171.bam','TNBC'],
['GSM1401758','SRR1313172.bam','TNBC'],
['GSM1401759','SRR1313173.bam','TNBC'],
['GSM1401795','SRR1313209.bam','TNBC-free'],
['GSM1401796','SRR1313210.bam','TNBC-free'],
['GSM1401797','SRR1313211.bam','TNBC-free'],
['GSM1401798','SRR1313212.bam','TNBC-free'],
['GSM1401799','SRR1313213.bam','TNBC-free'],
['GSM1401800','SRR1313214.bam','TNBC-free'],
['GSM1401801','SRR1313215.bam','TNBC-free'],
['GSM1401802','SRR1313216.bam','TNBC-free'],
['GSM1401803','SRR1313217.bam','TNBC-free'],
['GSM1401804','SRR1313218.bam','TNBC-free'],
['GSM1401805','SRR1313219.bam','TNBC-free'],
['GSM1401806','SRR1313220.bam','TNBC-free'],
['GSM1401807','SRR1313221.bam','TNBC-free'],
['GSM1401808','SRR1313222.bam','TNBC-free'],
['GSM1401809','SRR1313223.bam','TNBC-free'],
['GSM1401810','SRR1313224.bam','TNBC-free'],
['GSM1401811','SRR1313225.bam','TNBC-free'],
['GSM1401812','SRR1313226.bam','TNBC-free'],
['GSM1401813','SRR1313227.bam','TNBC-free'],
['GSM1401814','SRR1313228.bam','TNBC-free'],
['GSM1401815','SRR1313229.bam','TNBC-free'],

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


