import sys
#bedtools中存在的
hadgene=[]
with open("rep_all_novel.Amp.geneid_precent.gistic") as geneid:
	for line in geneid.readlines():
		hadgene.append(line.split(" ")[5].strip("\n"))
set_hadgene=list(set(hadgene))
print (len(set_hadgene))


#所有的lncRNA中的信息,方便取出信息
dic={}
with open("D:\\CRC_lncRNA\\cnv\\lncRNA.final.v2.novel.sorted.bed") as bed:
	for line3 in bed.readlines():
		ge=line3.split("\t")[3]
		if ge not in dic:
			dic[ge]=line3

overlap_gene=[]

#所有lncRNA，找出没有比对的
with open("novel_Amp_no_bedtools_lncRNA.bed",'w') as fw:
	with open("lncRNA.final.v2.novel.geneid.txt") as scores:
		for line2 in scores.readlines():
			allgeid=line2.strip("\n") 
			if allgeid not in set_hadgene:
				overlap_gene.append(allgeid)
				l=dic[allgeid].split("\t")
				fw.write(l[0]+"\t"+l[1]+"\t"+l[3]+"\t0\n")


set_gene=set(overlap_gene)
print (len(set_gene))


