from collections import defaultdict
dic= defaultdict(list)
filelist=['refgene_lncRNA.known.lncRNA.gtf.refmap','lncRNA_lncipedia.known.lncRNA.gtf.refmap','gencode_lncRNA.known.lncRNA.gtf.refmap']
for file in filelist:
	with open(file) as ref:
		for line in ref.readlines()[1:]:
			gene=line.split('\t')[0].strip(' ')
			cufid=line.split('\t')[3].split('|')[0]
			dic[cufid].append(gene)

# with open('cuffid_geneid.txt','w') as fw:
# 	for k in dic:
# 		fw.write(k+'\t')
# 		vaset=set(dic[k])
# 		for va in vaset:
# 			fw.write(va+',')
# 		fw.write('\n')
# 		
with open('geneid_lncRNA.txt','w') as fw:
	with open('lncRNA.final.v2.map') as mapf:
		for l in mapf.readlines():
			cuffid=l.split('\t')[1].strip('\n')
			if cuffid in dic:
				fw.write(cuffid+'\t')
				vaset=set(dic[cuffid])
				for va in vaset:
					fw.write(va+',')
				fw.write('\n')