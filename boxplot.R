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


known_gtf=read.table('D:\\CRC_lncRNA\\filter\\lncRNA\\lncRNA.final.v2.known.gtf',stringsAsFactors = F,sep='"')
known_gene_id=known_gtf[,2]
known_gene_id=unique(known_gene_id)
know_data=lnc[rownames(lnc)%in%known_gene_id,]

novel_gtf=read.table('D:\\CRC_lncRNA\\filter\\lncRNA\\lncRNA.final.v2.novel.gtf',stringsAsFactors = F,sep='"')
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

library(ggplot2)

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

#y=Log2(FPKM normalized expression)
pdf(file="D:\\CRC_lncRNA\\filter\\RSEM_expression\\three_RNA_expression_boxplot.pdf")

sp=ggplot(df)+geom_boxplot(aes(x=flag,y=num,fill=RNA))+labs(x="",y = "")+guides(fill=guide_legend(title="Legend_Title"))+ ylim(limits = c(-8, 25))
sp+theme_bw() + theme(legend.title=element_blank(),legend.position=c(0.16,0.9),panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),axis.title.x = element_text(size = 15, face = "bold"),axis.title.y= element_text(size = 15, face = "bold"))
dev.off()
