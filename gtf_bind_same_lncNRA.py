gene=[]
zd={}
zdhs={}
with open('D:\\CRC_lncRNA\\filter\\lncRNA\\only_min_max_position_lncRNA.final.v2.gtf','w') as fwgtf:
	with open('D:\\CRC_lncRNA\\filter\\lncRNA\\lncRNA.final.v2.gtf') as gtf:
		for line2 in gtf.readlines():
			gene_id=line2.split("\"")[1]
			ch,start,end=line2.split("\t")[0],line2.split("\t")[3],line2.split("\t")[4]
			zdhs[gene_id]=ch

			if gene_id not in gene:
				newlis=[]
				gene.append(gene_id)
				newlis.append(start)				
				newlis.append(end)
				zd[gene_id]=newlis
			else:
				newlis.append(start)				
				newlis.append(end)
	for key in zd:
		start=min(zd[key])
		end=max(zd[key])
		fwgtf.write(zdhs[key]+'\t'+start+'\t'+end+'\t'+key+'\n')

