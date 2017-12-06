'''
from collections import defaultdict
d = defaultdict(list)

#known
with open('D:\\CRC_lncRNA\\filter\\lncRNA\\protein_coding.final.gtf') as novel:
	for line in novel.readlines()[0:]:
		#if line.split('\t')[2]=='exon':
		if line.split('\t')[2]=='transcript':
			tranid=line.split('\t')[-1].split('"')[3]
			length=int(line.split('\t')[4])-int(line.split('\t')[3])
			# if length=='NA':
			# 	length=0
			d[tranid].append(length)

#novel
with open('pc_length.txt','w') as fprint:
	for k in d.keys():
		allsum=0
		le=len(d[k])
		val=d[k]
		for i in range(le):
			allsum =allsum+int(val[i])
		fprint.write(k+'\t'+str(allsum)+'\n')
'''
gene=[]
with open('D:\\CRC_lncRNA\\filter\\lncRNA\\lncRNA.final.v2.known.gtf') as novel:
	for line in novel.readlines():
		gene_id=line.split('\t')[-1].split('"')[1]
		#length=int(line.split('\t')[4])-int(line.split('\t')[3])
		# if length=='NA':
		# 	length=0
		gene.append(gene_id)

		
unigeneid=set(gene)
with open('known_geneid2.txt','w') as fa:
	for g in unigeneid:
		fa.write(g+'\n')

		