import sys
#had cnv
hadgene=[]
with open(sys.argv[1]) as geneid:
	for line in geneid.readlines():
		hadgene.append(line.split("\t")[2])
set_hadgene=list(set(hadgene))


# dic={}
# with open("/disk/zhw/CRClncRNA/filter/lncRNA/lncRNA.final.v2.novel_known.sorted.bed") as bed:
# 	for line3 in bed.readlines():
# 		ge=line3.split("\t")[3]
# 		if ge not in dic:
# 			dic[ge]=line3

overlap_gene=[]

with open(sys.argv[2],'w') as fw:
	with open("/disk/zhw/CRClncRNA/cnv/differentgene_updown_heatmap/intersect_normal_cnv_up_down_gtf.bed") as scores:
		for line2 in scores.readlines():
			allgeid=line2.split("\t")[-1].strip("\r\n") 
			if allgeid not in set_hadgene:
				overlap_gene.append(allgeid)
				#l=dic[allgeid].split("\t")
				fw.write("chr0\t0\t"+allgeid+"\t0"+"\n")


set_gene=set(overlap_gene)
print (len(set_gene))
print (len(set_hadgene))
print (len(set_gene)+len(set_hadgene))

