
def getchpositionlog(up_down,rec_normal):
	hadgene=[]
	gene=[]
	gene_log={}
	with open(rec_normal+'_DESeq2_edgeR_res_intersect_'+up_down+'.txt') as dfgene:
		for line in dfgene.readlines()[1:]:
			genename=line.split('\t')[0]
			if (up_down=="down"):
				log=-float(line.split('\t')[2])
			else:
				log=float(line.split('\t')[2])			
			gene.append(genename)
			gene_log[genename]=log
	print ('gene_log is ojbk')

	with open("D:\\CRC_lncRNA\\cnv\\circos\\"+rec_normal+"_DESeq2_edgeR_"+up_down+"_gene_for_circos.txt",'w') as fw:

		with open('D:\\CRC_lncRNA\\filter\\lncRNA\\only_min_max_position_lncRNA.final.v2.gtf') as gtf:
			for line2 in gtf.readlines():
				gene_id=line2.split("\t")[3].strip('\n')
				if gene_id in gene:
					if gene_id not in hadgene:
						hadgene.append(gene_id)
						ch,start,end=line2.split("\t")[0],line2.split("\t")[1],line2.split("\t")[2]
						#position=int((float(end)+float(start))/2)
						fw.write("hs"+ch.strip("chr")+' '+start+' '+end+' '+str(gene_log[gene_id])+'\n')
						#fw.write("hs"+ch.strip("chr")+' '+str(position)+' '+str(position+1)+' '+str(gene_log[gene_id])+'\n')

getchpositionlog('up','rec')
getchpositionlog('down','rec')
getchpositionlog('up','normal')
getchpositionlog('down','normal')