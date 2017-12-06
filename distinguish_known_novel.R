setwd('D:\\CRC_lncRNA\\filter\\RSEM_expression')

data=read.table('lncRNA.rsem.FPKM.txt',header = T,stringsAsFactors = F,sep = '\t')

samplename=strsplit(colnames(data),split="_",fixed=T)
rown=c()
for (n in (1:length(samplename))){
  g=samplename[[n]][1]
  rown=c(rown,g)
}
colnames(data)=rown

data=data[,sort(rown)]
write.table(data,'lncRNA.rsem.FPKM_sort.txt',quote = F,sep='\t')


data2=read.table('lncRNA.rsem.FPKM_sort.txt',header = T,stringsAsFactors = F,sep = '\t')
data=matrix(as.numeric(unlist(data2)),ncol=40)
rowmean_data=data.frame(rowMeans(data[,c(1:10)]),rowMeans(data[,c(11:20)]),rowMeans(data[,c(21:30)]),rowMeans(data[,c(31:40)]))
colnames(rowmean_data)=c("A","B","C","D")
row.names(rowmean_data)=row.names(data2)

rowmean_data=log2(rowmean_data+1)

gtf=read.table('D:\\CRC_lncRNA\\filter\\lncRNA.final.v2.cp.2.txt',sep = '\t',stringsAsFactors = F)
colnames(gtf)=c("chr","known","start","end","geneid")

data2=data.frame(rownames(rowmean_data),rowmean_data)

colnames(data2)=c("geneid",colnames(rowmean_data))
gtf_data=merge(gtf,data2,by="geneid",sort = F)

hs=sub('chr','hs',gtf_data[,2]) #¸ü»»chr Îªhs
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

