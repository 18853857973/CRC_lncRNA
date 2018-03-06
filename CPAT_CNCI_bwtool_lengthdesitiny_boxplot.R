setwd('D:\\CRC_lncRNA\\filter\\RSEM_expression')
#样本按照肿瘤和正常排序
# data=read.table('lncRNA.rsem.FPKM.txt',header = T,stringsAsFactors = F,sep = '\t',check.names = F)
# samplename=strsplit(colnames(data),split="_",fixed=T)
# rown=c()
# for (n in (1:length(samplename))){
#   g=samplename[[n]][1]
#   rown=c(rown,g)
# }
# colnames(data)=rown
# data=data[,sort(rown)]
# write.table(data,'lncRNA.rsem.FPKM_sort.txt',quote = F,sep='\t')

#筛选在各个样本中表达大于等于1个数大于等于2的
data2=read.table('D:\\CRC_lncRNA\\filter\\RSEM_expression\\lncRNA.rsem.FPKM_sort.txt',header = T,stringsAsFactors = F,sep = '\t',check.names = F)
tof=apply(data2,1,function(x) length(which(x>=1))>=2)
data4=data2[tof,]
write.table(data4,'D:\\CRC_lncRNA\\filter\\RSEM_expression\\lncRNA.rsem.FPKM_sort_morethan1_num_2.txt',quote=F,sep='\t')

data2=data4

write.table(row.names(data4),"D:\\CRC_lncRNA\\diffexp\\lncRNA.rsem.FPKM_sort_morethan1_num_2_genenames.txt",quote=F,row.names=F,col.names=F)



#data2=read.table('D:\\CRC_lncRNA\\filter\\RSEM_expression\\lncRNA.rsem.FPKM_sort_morethan1_num_2.txt',header = T,stringsAsFactors = F,sep = '\t',check.names = F)
data=matrix(as.numeric(unlist(data2)),ncol=40)
rowmean_data=data.frame(rowMeans(data[,c(1:10)]),rowMeans(data[,c(11:20)]),rowMeans(data[,c(21:30)]),rowMeans(data[,c(31:40)]))
colnames(rowmean_data)=c("A","B","C","D")
row.names(rowmean_data)=row.names(data2)

rowmean_data=log2(rowmean_data+1)

gtf=read.table('D:\\CRC_lncRNA\\filter\\lncRNA.final.v2.cp.2.txt',sep = '\t',stringsAsFactors = F)
colnames(gtf)=c("chr","known","start","end","geneid")

data2=data.frame(rownames(rowmean_data),rowmean_data)
#筛选出转录活性的lncRNA
#data2=data2[data2[,1]%in%lncRNA_TF,]

colnames(data2)=c("geneid",colnames(rowmean_data))
gtf_data=merge(gtf,data2,by="geneid",sort = F)

hs=sub('chr','hs',gtf_data[,2]) #更换chr 为hs
gtf_data=data.frame(gtf_data[,1],hs,gtf_data[,c(3:9)])

know_gtf_data=gtf_data[which(gtf_data[,3]=="known"),]
novel_gtf_data=gtf_data[which(gtf_data[,3]=="novel"),]

normal_know_gtf_data=data.frame(know_gtf_data[,c(2,4,5)],rowMeans(know_gtf_data[,c(6,7)]))
recurrence_know_gtf_data=data.frame(know_gtf_data[,c(2,4,5)],rowMeans(know_gtf_data[,c(6,8)])) 
unrecurrence_know_gtf_data=data.frame(know_gtf_data[,c(2,4,5)],rowMeans(know_gtf_data[,c(7,9)]))

normal_novel_gtf_data=data.frame(novel_gtf_data[,c(2,4,5)],rowMeans(novel_gtf_data[,c(6,7)]))
recurrence_novel_gtf_data=data.frame(novel_gtf_data[,c(2,4,5)],rowMeans(novel_gtf_data[,c(6,8)]))
unrecurrence_novel_gtf_data=data.frame(novel_gtf_data[,c(2,4,5)],rowMeans(novel_gtf_data[,c(7,9)]))

setwd('D:\\CRC_lncRNA\\filter\\RSEM_expression')
write.table(normal_know_gtf_data,'hs_normal_know_gtf_data.txt',quote = F,col.names = F,row.names = F)
write.table(recurrence_know_gtf_data,'hs_recurrence_know_gtf_data.txt',quote = F,col.names = F,row.names = F)
write.table(unrecurrence_know_gtf_data,'hs_unrecurrence_know_gtf_data.txt',quote = F,col.names = F,row.names = F)

write.table(normal_novel_gtf_data,'hs_normal_novel_gtf_data.txt',quote = F,col.names = F,row.names = F)
write.table(recurrence_novel_gtf_data,'hs_recurrence_novel_gtf_data.txt',quote = F,col.names = F,row.names = F)
write.table(unrecurrence_novel_gtf_data,'hs_unrecurrence_novel_gtf_data.txt',quote = F,col.names = F,row.names = F)


###################################################################
setwd('D:\\CRC_lncRNA\\coding_potential')

# novel_trans_id=read.table('D:\\CRC_lncRNA\\coding_potential\\novel_transcript_id_unique.txt',stringsAsFactors = F)
# novel_id=novel_trans_id[,2]
# known_trans_id=read.table('D:\\CRC_lncRNA\\coding_potential\\known_transcript_id_unique.txt',stringsAsFactors = F)
# known_id=known_trans_id[,2]
# CPAT=read.table('D:\\CRC_lncRNA\\coding_potential\\CPAT\\lncRNA',sep='\t',stringsAsFactors = F)
# CPAT_known=CPAT[which(rownames(CPAT)%in%known_id),]
# CPAT_novel=CPAT[which(rownames(CPAT)%in%novel_id),]

###CPAT###
CPAT_known=read.table('D:\\CRC_lncRNA\\coding_potential\\CPAT\\CPAT_lncRNA.final.v2.known',sep = '\t',stringsAsFactors = F)
CPAT_novel=read.table('D:\\CRC_lncRNA\\coding_potential\\CPAT\\CPAT_lncRNA.final.v2.novel',sep = '\t',stringsAsFactors = F)
CPAT_protein=read.table('D:\\CRC_lncRNA\\coding_potential\\CPAT\\CPAT_protein_coding.final',sep = '\t',stringsAsFactors = F)

library(ggplot2)
x=CPAT_known[,5]
y=ecdf(x)
df1<-data.frame(x=x,y=y(x))

x=CPAT_novel[,5]
y=ecdf(x)
df2<-data.frame(x=x,y=y(x))

x=CPAT_protein[,5]
y=ecdf(x)
df3<-data.frame(x=x,y=y(x))
pdf(file="D:\\CRC_lncRNA\\coding_potential\\CPAT\\CPAT.pdf")
sp2=ggplot(df1, aes(x,y)) +geom_line(size = 1.5,aes(color="known"))+geom_line(size = 1.5,data=df2,aes(color="novel"))+geom_line(size = 1.5,data=df3,aes(color="protein"))+labs(color="Legend text")+xlab("") + ylab("") 
sp2+theme_bw() + theme(legend.title=element_blank(),legend.position=c(0.8,0.15),panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),axis.title.x = element_text(size = 15, face = "bold"),axis.title.y= element_text(size = 15, face = "bold"))
dev.off()

###CNCI###
CNCI_known=read.table('D:\\CRC_lncRNA\\coding_potential\\CNCI\\CNCI_lncRNA.final.v2.known\\CNCI.index',stringsAsFactors = F,sep='\t')
CNCI_novel=read.table('D:\\CRC_lncRNA\\coding_potential\\CNCI\\CNCI_lncRNA.final.v2.novel\\CNCI.index',stringsAsFactors = F,sep='\t')
CNCI_pro=read.table('D:\\CRC_lncRNA\\coding_potential\\CNCI\\CNCI_protein_coding.final\\CNCI.index',stringsAsFactors = F,sep='\t')

x=CNCI_known[-1,3]
x=as.numeric(x)
y=ecdf(x)
df1<-data.frame(x=x,y=y(x))

x=CNCI_novel[-1,3]
x=as.numeric(x)
y=ecdf(x)
df2<-data.frame(x=x,y=y(x))

x=CNCI_pro[-1,3]
x=as.numeric(x)
y=ecdf(x)
df3<-data.frame(x=x,y=y(x))
pdf(file="D:\\CRC_lncRNA\\coding_potential\\CNCI\\CNCI.pdf")

sp3=ggplot(df1, aes(x,y)) +geom_line(size = 1.5,aes(color="known"))+geom_line(size = 1.5,data=df2,aes(color="novel"))+geom_line(size = 1.5,data=df3,aes(color="protein"))+labs(color="Legend text")+xlab("Protein coding protential") + ylab("Cumulative distribution") 
sp3+theme_bw() + theme(legend.title=element_blank(),legend.position=c(0.8,0.8),panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),axis.title.x = element_text(size = 15, face = "bold"),axis.title.y= element_text(size = 15, face = "bold"))
dev.off()

###bwtool###
setwd('D:\\CRC_lncRNA\\coding_potential\\bwtools')
bwtool_known=read.table('bwtool_lncRNA.final.v2.known_mean.txt',sep='\t')
bwtool_known=na.omit(bwtool_known)
x=bwtool_known[,5]
y=ecdf(x)
df1<-data.frame(x=x,y=y(x))

bwtool_novel=read.table('bwtool_lncRNA.final.v2.novel_mean.txt',sep='\t')
bwtool_novel=na.omit(bwtool_novel)
x=bwtool_novel[,5]
y=ecdf(x)
df2<-data.frame(x=x,y=y(x))

bwtool_protein=read.table('bwtool_protein_coding.final_mean.txt',sep='\t')
bwtool_protein=na.omit(bwtool_protein)
x=bwtool_protein[,5]
y=ecdf(x)
df3<-data.frame(x=x,y=y(x))
pdf(file="D:\\CRC_lncRNA\\coding_potential\\bwtools\\bwtools.pdf")
sp4=ggplot(df1, aes(x,y)) +geom_line(size = 1.5,aes(color="known"))+geom_line(size = 1.5,data=df2,aes(color="novel"))+geom_line(size = 1.5,data=df3,aes(color="protein"))+labs(color="Legend text")+xlab("") + ylab("") 
sp4+theme_bw() + theme(legend.title=element_blank(),legend.position=c(0.8,0.2),panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),axis.title.x = element_text(size = 15, face = "bold"),axis.title.y= element_text(size = 15, face = "bold"))
dev.off()


#######################################################################################
#length_desitiny
#exon_length()
setwd('D:\\CRC_lncRNA\\filter\\lncRNA')
novel_gtf=read.table('D:\\CRC_lncRNA\\filter\\lncRNA\\num2_lncRNA.novel.length.txt',sep='\t',stringsAsFactors = F)
novel_length=novel_gtf[,2]
x=novel_length

known_gtf=read.table('D:\\CRC_lncRNA\\filter\\lncRNA\\num2_lncRNA.known.length.txt',sep='\t',stringsAsFactors = F)
known_length=known_gtf[,2]
x2=known_length

coding_gtf=read.table('D:\\CRC_lncRNA\\filter\\lncRNA\\pc_length.txt',sep='\t',stringsAsFactors = F)
coding_length=coding_gtf[,2]
x3=coding_length

df1=data.frame(length=x2,type='known lncRNA')
df2=data.frame(length=x,type='novel lncRNA')
df3=data.frame(length=x3,type='protein-coding gene')

df=rbind(df1,df2,df3)

pdf(file="D:\\CRC_lncRNA\\filter\\lncRNA\\num2_length_desitiny_plot.pdf")

sp=ggplot(df,aes(x = length, fill= type,colour = type)) +geom_density()+xlim(limits = c(-1000, 90000))+labs(x="",y = " ")
sp+theme_bw() + theme(legend.title=element_blank(),legend.position=c(0.8,0.8),panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),axis.title.x = element_text(size = 15, face = "bold"),axis.title.y= element_text(size = 15, face = "bold"))
dev.off()

############################################################################################################################
#箱线图
setwd('D:\\CRC_lncRNA\\filter\\RSEM_expression')
lnc=read.table('lncRNA.rsem.FPKM_sort.txt',sep='\t',stringsAsFactors = F,header=T)
pc=read.table('pcRNA.rsem.FPKM_sort.txt',sep='\t',stringsAsFactors = F,header=T)

pc_lnc=rbind(lnc,pc)
pc_lnc2=matrix(as.numeric(unlist(pc_lnc)),ncol=40)
rownames(pc_lnc2)=rownames(pc_lnc)
colnames(pc_lnc2)=colnames(pc_lnc)
pc_lnc3<-na.omit(pc_lnc2)
#pc_lnc4=(pc_lnc3-min(pc_lnc3))/(max(pc_lnc3)-pc_lnc3)
pc_lnc5=log2(pc_lnc3+0.01)

lnc=pc_lnc5[1:35169,]
pc=pc_lnc5[-c(1:35169),]

known_gtf=read.table('D:\\CRC_lncRNA\\filter\\lncRNA\\num2_lncRNA.final.v2.known.gtf',stringsAsFactors = F,sep='"')
known_gene_id=known_gtf[,2]
known_gene_id=unique(known_gene_id)
know_data=lnc[rownames(lnc)%in%known_gene_id,]

novel_gtf=read.table('D:\\CRC_lncRNA\\filter\\lncRNA\\num2_lncRNA.final.v2.novel.gtf',stringsAsFactors = F,sep='"')
novel_gene_id=novel_gtf[,2]
novel_gene_id=unique(novel_gene_id)
novel_data=lnc[rownames(lnc)%in%novel_gene_id,]

know_normal=rowMeans(know_data[,c(1:20)])
know_rec=  rowMeans(know_data[,c(11:20)])
know_norec=rowMeans(know_data[,c(21:30)]) 

novel_normal= rowMeans(novel_data[,c(1:20)])
novel_rec=rowMeans(novel_data[,c(11:20)])
novel_norec=rowMeans(novel_data[,c(21:30)])

pc_normal=rowMeans(pc[,c(1:20)])
pc_rec=  rowMeans(pc[,c(11:20)])
pc_norec=rowMeans(pc[,c(21:30)]) 

df_know_normal=data.frame(num=know_normal,type='know_normal',RNA='know lncRNA')
df_know_rec=data.frame(num=know_rec,type='know_rec',RNA='know lncRNA')
df_know_norec=data.frame(num=know_norec,type='know_norec',RNA='know lncRNA')

df_novel_normal=data.frame(num=novel_normal,type='novel_normal',RNA='novel lncRNA')
df_novel_rec=data.frame(num=novel_rec,type='novel_rec',RNA='novel lncRNA')
df_novel_norec=data.frame(num=novel_norec,type='novel_norec',RNA='novel lncRNA')

df_pc_normal=data.frame(num=pc_normal,type='pc_normal',RNA='protein-coding gene')
df_pc_rec=data.frame(num=pc_rec,type='pc_rec',RNA='protein-coding gene')
df_pc_norec=data.frame(num=pc_norec,type='pc_norec',RNA='protein-coding gene')

df1=rbind(df_know_normal,df_novel_normal,df_pc_normal)
df1=data.frame(df1,flag='normal')

df2=rbind(df_know_rec,df_novel_rec,df_pc_rec)
df2=data.frame(df2,flag='rec')

df3=rbind(df_know_norec,df_novel_norec,df_pc_norec)
df3=data.frame(df3,flag='norec')

df=rbind(df1,df2,df3)

pdf(file="D:\\CRC_lncRNA\\filter\\RSEM_expression\\num2_three_RNA_expression_boxplot.pdf")
sp=ggplot(df)+geom_boxplot(aes(x=flag,y=num,fill=RNA))+labs(x="",y = "")+guides(fill=guide_legend(title="Legend_Title"))+ ylim(limits = c(-8, 25))
sp+theme_bw() + theme(legend.title=element_blank(),legend.position=c(0.16,0.9),panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),axis.title.x = element_text(size = 15, face = "bold"),axis.title.y= element_text(size = 15, face = "bold"))
dev.off()


