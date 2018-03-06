setwd('D:\\CRC_lncRNA\\filter\\RSEM_expression')
data=read.table("pcRNA.rsem.FPKM_sort.txt",check.names ="F",header = T)
data2=data[,c(11:30)]

data3=data2[,-c(2,5,7,9,13,16,17,18)]

data3=data2
gene=c("METTL3","METTL14","WTAP","ALKBH5","FTO","YTHDF1","YTHDF2","YTHDF3")

data4=data3[gene,]

sampleInfo=data.frame(sample=colnames(data4),Subset=rep(c("rec","norec"),c(10,10)))

colnum=2
library('pheatmap')
library("RColorBrewer")
library("RColorBrewer")

mergeSample2=data.frame(sampleInfo$Subset)
rownames(mergeSample2)=sampleInfo[,1]
colnames(mergeSample2)=c("Subset")
#col= brewer.pal(colnum, "Paired")
col=c("purple" ,"orange")
Subset=col
names(Subset) = levels(mergeSample2$Subset)
colll=list(Subset=Subset)
tMatrix=t(data4)

tMatrix_num=matrix(as.numeric(unlist(tMatrix)),ncol=ncol(tMatrix))
rownames(tMatrix_num)=rownames(tMatrix)
colnames(tMatrix_num)=colnames(tMatrix)
mydata = scale(tMatrix_num, center = TRUE, scale = TRUE)
mydata=t(mydata)
mydata[which(mydata>2,arr.ind=T)]=2
mydata[which(mydata<(-2),arr.ind=T)]=(-2)
pheatmap(mydata,annotation_col = mergeSample2,cluster_cols = F,cluster_rows =F ,
         colorRampPalette(c("green", "black", "red"))(50),show_rownames=T,show_colnames=F,annotation_colors =colll)

mydata2=mydata[,-c(2,5,7,9,11,17,19,20)]
pheatmap(mydata2,annotation_col = mergeSample2,cluster_cols = F,cluster_rows =F ,
         colorRampPalette(c("green", "black", "red"))(50),show_rownames=T,show_colnames=F,annotation_colors =colll)
