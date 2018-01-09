from collections import defaultdict

def protein_coding():
#最开始从gtf文件中的蛋白质编码基因中筛选
	with open('gencode.v27.onlygrep_protein_coding.gtf') as f:
		with open ('gencode.v27.protein_coding2.gtf','w') as fw:
			for line in f.readlines():
				if line.split('\t')[2]=='gene':
					fw.write(line)
				else:
					if (line.split('\t')[8].split('"')[9])=="protein_coding":
						fw.write(line)


#lncRNA 转录起始位点前后2kb
def getlncRNAposition(inputfile,outputfile):
	lncRNAname=[]
	with open (inputfile) as intersetlncRNA:
		for line in intersetlncRNA.readlines():
			lncRNAname.append(line.strip())

	with open(outputfile,'w') as fw:
		with open ('D:\\CRC_lncRNA\\filter\\lncRNA\\only_min_max_position_lncRNA2.final.v2.gtf') as alllncRNA:
			for line2 in alllncRNA.readlines():
				if line2.split('\t')[3] in lncRNAname:
					zf=line2.split('\t')[-1].strip('\n').strip(' ')
					if zf=='+':
						start=int(line2.split('\t')[1])-2000
						end=start+4000
						fw.write(line2.split('\t')[0]+'\t'+str(start)+'\t'+str(end)+'\t'+line2.split('\t')[3]+'\r')
					if zf=='-':
						start=int(line2.split('\t')[2])-2000
						end=start+4000
						fw.write(line2.split('\t')[0]+'\t'+str(start)+'\t'+str(end)+'\t'+line2.split('\t')[3]+'\r')

#getlncRNAposition("D:\\CRC_lncRNA\\diffexp\\lncRNA_TF.txt","D:\\CRC_lncRNA\\diffexp\\5151_lncRNA_TF.bed")


#lncRNA 转录起始位点前后2kb,每5bp 一行
def getlncRNApositionper5bp(inputfile,outputfile):
	lncRNAname=[]
	#with open ('D:\\CRC_lncRNA\\diffexp\\exprmorethan2_lncRNA_genename.txt') as intersetlncRNA:
	with open (inputfile) as intersetlncRNA:
		for line in intersetlncRNA.readlines()[0:]:
			lncRNAname.append(line.strip())

	#with open("D:\\CRC_lncRNA\\chipseq\\exp_more_than2_lncRNA.bed",'w') as fw:
	with open(outputfile,'w') as fw:
		with open ('D:\\CRC_lncRNA\\filter\\lncRNA\\only_min_max_position_lncRNA2.final.v2.gtf') as alllncRNA:
			for line2 in alllncRNA.readlines():
				if line2.split('\t')[3] in lncRNAname:
					zf=line2.split('\t')[-1].strip('\n').strip(' ')
					if zf=='+':
						start=int(line2.split('\t')[1])-2000
					if zf=='-':
						start=int(line2.split('\t')[2])-2000
					#方向应该相反

					end=start+4000
					se=start #1 #5
					while se <end:
						fw.write(line2.split('\t')[0]+'\t'+str(se)+'\t'+str(se+4)+'\t'+line2.split('\t')[3]+'\n')
						se=se+5 #6						

#getlncRNApositionper5bp("D:\\CRC_lncRNA\\diffexp\\lncRNA_TF.txt","D:\\CRC_lncRNA\\diffexp\\5151_lncRNA_TF_per5bp.bed")

#lncRNA 转录起始位点前后2kb,每5bp 一行
def getlncRNAdirection(inputfile,outputfile):
	lncRNAname=[]
	#with open ('D:\\CRC_lncRNA\\diffexp\\exprmorethan2_lncRNA_genename.txt') as intersetlncRNA:
	with open (inputfile) as intersetlncRNA:
		for line in intersetlncRNA.readlines()[0:]:
			lncRNAname.append(line.strip())

	#with open("D:\\CRC_lncRNA\\chipseq\\exp_more_than2_lncRNA.bed",'w') as fw:
	with open(outputfile,'w') as fw:
		with open ('D:\\CRC_lncRNA\\filter\\lncRNA\\only_min_max_position_lncRNA2.final.v2.gtf') as alllncRNA:
			for line2 in alllncRNA.readlines():
				if line2.split('\t')[3] in lncRNAname:
					zf=line2.split('\t')[-1].strip('\n').strip(' ')		
					fw.write(line2.split('\t')[3]+'\t'+line2.split('\t')[4])
			

#getlncRNAdirection("D:\\CRC_lncRNA\\diffexp\\lncRNA_TF.txt","D:\\CRC_lncRNA\\diffexp\\5151_lncRNA_TF_direc.bed")



#先从gtf中筛选出有转录活性的lncRNA
#在样本间表达量大于等于1个数大于等于2的筛选后
def gtf_lncRNA_TF(known_novel):
	#所有有转录活性的lncRNA
	tflncRNA=[]
	with open('D:\\CRC_lncRNA\\diffexp\\lncRNA.rsem.FPKM_sort_morethan1_num_2_genenames.txt') as TF_lncRNA:
		for line in TF_lncRNA.readlines():
			tflncRNA.append(line.strip())

	with open ('D:\\CRC_lncRNA\\filter\\lncRNA\\num2_lncRNA.final.v2.'+known_novel+'.gtf','w') as fwknown:
		with open('D:\\CRC_lncRNA\\filter\\lncRNA\\lncRNA.final.v2.'+known_novel+'.gtf') as knowngtf:
			for line2 in knowngtf.readlines():
				if line2.split("\"")[1] in tflncRNA:
					fwknown.write(line2)
#gtf_lncRNA_TF('known')
#gtf_lncRNA_TF('novel')


def exon_length(filename,outputfile):
	d = defaultdict(list)
	#known
	with open('D:\\CRC_lncRNA\\filter\\lncRNA\\'+filename) as novel:
		for line in novel.readlines()[0:]:
			#lncRNA
			if line.split('\t')[2]=='exon':
			#if line.split('\t')[2]=='transcript':
				tranid=line.split('\t')[-1].split('"')[3]
				length=int(line.split('\t')[4])-int(line.split('\t')[3])
				# if length=='NA':
				# 	length=0
				d[tranid].append(length)

	with open('D:\\CRC_lncRNA\\filter\\lncRNA\\'+outputfile,'w') as fprint:
		for k in d.keys():
			allsum=0
			le=len(d[k])
			val=d[k]
			for i in range(le):
				allsum =allsum+int(val[i])
			fprint.write(k+'\t'+str(allsum)+'\n') 

#exon_length('num2_lncRNA.final.v2.known.gtf','num2_lncRNA.known.length.txt')
#exon_length('num2_lncRNA.final.v2.novel.gtf','num2_lncRNA.novel.length.txt')

#gtf文件相同的基因合并，取最小最大位点
def gtf_bind_same_lncRNA():
	gene=[]
	zd={}
	zdhs={}
	zfl={}
	with open('D:\\CRC_lncRNA\\filter\\lncRNA\\only_min_max_position_lncRNA2.final.v2.gtf','w') as fwgtf:
		with open('D:\\CRC_lncRNA\\filter\\lncRNA\\lncRNA.final.v2.gtf') as gtf:
			for line2 in gtf.readlines():
				gene_id=line2.split("\"")[1]
				ch,start,end,zf=line2.split("\t")[0],line2.split("\t")[3],line2.split("\t")[4],line2.split("\t")[6]
				zdhs[gene_id]=ch
				zfl[gene_id]=zf

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
			fwgtf.write(zdhs[key]+'\t'+start+'\t'+end+'\t'+key+'\t'+zfl[key]+'\n')
#gtf_bind_same_lncRNA()



def getchpositionlog(up_down,rec_normal):
	hadgene=[]
	gene=[]
	gene_log={}
	with open('D:/CRC_lncRNA/diffexp/'+rec_normal+'_DESeq2_edgeR_res_intersect_'+up_down+'.txt') as dfgene:
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

# getchpositionlog('up','rec')
# getchpositionlog('down','rec')
# getchpositionlog('up','normal')
# getchpositionlog('down','normal')


def dfgene_location():
	alldfgene=[]
	with open('D:\\CRC_lncRNA\\diffexp\\tumor_vs_normal_DESeq2_edgeR_intersect_up_down.txt') as dfgene:
		for line in dfgene.readlines():
			genename=line.split('\n')[0]
			alldfgene.append(genename)
			#print (genename)

	up_alldfgene=alldfgene[0:685]
	down_alldfgene=alldfgene[685:len(alldfgene)]

	with open("D:\\CRC_lncRNA\\cnv\\differentgene_updown_heatmap\\intersect_normal_cnv_up_down_gtf.bed","w") as  fw:
		with open('D:\\CRC_lncRNA\\filter\\lncRNA\\only_min_max_position_lncRNA.final.v2.gtf') as gtf:
				for line2 in gtf.readlines():
					gene_id=line2.split("\t")[3].strip("\n")
					if gene_id in up_alldfgene:
						fw.write(line2)

	with open("D:\\CRC_lncRNA\\cnv\\differentgene_updown_heatmap\\intersect_normal_cnv_up_down_gtf.bed","a") as  fw:
		with open('D:\\CRC_lncRNA\\filter\\lncRNA\\only_min_max_position_lncRNA.final.v2.gtf') as gtf:
				for line2 in gtf.readlines():
					gene_id=line2.split("\t")[3].strip("\n")
					if gene_id in down_alldfgene:
						fw.write(line2)

#dfgene_location()


def getlncRNA_position():
	lncRNA=[]
	with open("D:\\CRC_lncRNA\\cnv\\percentCNV\\num2_normal_rec_0.25CNV_lncRNA.txt")as  f:
		for l in f.readlines():
			lncRNA.append(l.strip("\n"))
	with open("D:\CRC_lncRNA\TCGA_survive\surive_cnv\\num2_normal_rec_0.25CNV_lncRNA.bed","w") as  fw:
		with open('D:\\CRC_lncRNA\\filter\\lncRNA\\only_min_max_position_lncRNA.final.v2.gtf') as gtf:
				for line2 in gtf.readlines():
					gene_id=line2.split("\t")[3].strip("\n")
					if gene_id in lncRNA:
						fw.write(line2)

getlncRNA_position()