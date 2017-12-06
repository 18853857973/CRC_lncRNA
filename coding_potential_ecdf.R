setwd('D:\\CRC_lncRNA\\coding_potential')

# novel_trans_id=read.table('D:\\CRC_lncRNA\\coding_potential\\novel_transcript_id_unique.txt',stringsAsFactors = F)
# novel_id=novel_trans_id[,2]
# known_trans_id=read.table('D:\\CRC_lncRNA\\coding_potential\\known_transcript_id_unique.txt',stringsAsFactors = F)
# known_id=known_trans_id[,2]
# CPAT=read.table('D:\\CRC_lncRNA\\coding_potential\\CPAT\\lncRNA',sep='\t',stringsAsFactors = F)
# CPAT_known=CPAT[which(rownames(CPAT)%in%known_id),]
# CPAT_novel=CPAT[which(rownames(CPAT)%in%novel_id),]

#CPAT
CPAT_known=read.table('D:\\CRC_lncRNA\\coding_potential\\CPAT\\CPAT_lncRNA.final.v2.known',sep = '\t',stringsAsFactors = F)
CPAT_novel=read.table('D:\\CRC_lncRNA\\coding_potential\\CPAT\\CPAT_lncRNA.final.v2.novel',sep = '\t',stringsAsFactors = F)
CPAT_protein=read.table('D:\\CRC_lncRNA\\coding_potential\\CPAT\\CPAT_protein_coding.final',sep = '\t',stringsAsFactors = F)

#ECDF
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
sp2+theme_bw() + theme(legend.title=element_blank(),legend.position=c(0.8,0.3),panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),axis.title.x = element_text(size = 15, face = "bold"),axis.title.y= element_text(size = 15, face = "bold"))
dev.off()


#CNCI
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


######################
#bwtool
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





