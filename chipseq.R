setwd('D:\\CRC_lncRNA\\chipseq')

getgene_chipseqnum=function(H3K){
  TF_lncRNA_per5bp=read.table(paste("D:\\CRC_lncRNA\\chipseq\\",H3K,"_5151_lncRNA_TF_per5bp.bed",sep=''),check.names = F,stringsAsFactors = F,sep='\t')
  TF_lncNRA=unique(TF_lncRNA_per5bp[,4])
  which(TF_lncRNA_per5bp[,4]==TF_lncNRA[1])
  gene_num=TF_lncRNA_per5bp[1:800,5]#length(TF_lncNRA)
  for (k in c(2:length(TF_lncNRA))){
    b=TF_lncRNA_per5bp[((800*(k-1))+1):(800*k),5]
    gene_num=rbind(gene_num,b)
    print(k)
  }
  rownames(gene_num)=TF_lncNRA
  write.table(gene_num,paste('D:\\CRC_lncRNA\\chipseq\\',H3K,'_5151_lncRNA_TF_per5bp.txt',sep=''),quote = F,col.names = F,sep = '\t')
}
getgene_chipseqnum('H3K27ac')
getgene_chipseqnum('H3K4me3')
getgene_chipseqnum('H3K4me1')

H3K4me3=read.table("D:\\CRC_lncRNA\\chipseq\\H3K4me3_5151_lncRNA_TF_per5bp.txt",sep='\t',check.names = F,stringsAsFactors = F,row.names = 1)
H3K4me1=read.table("D:\\CRC_lncRNA\\chipseq\\H3K4me1_5151_lncRNA_TF_per5bp.txt",sep='\t',check.names = F,stringsAsFactors = F,row.names = 1)
H3K27ac=read.table("D:\\CRC_lncRNA\\chipseq\\H3K27ac_5151_lncRNA_TF_per5bp.txt",sep='\t',check.names = F,stringsAsFactors = F,row.names = 1)


#m1总read:62073675
#me3:55858490

H3K4me3_normalize=H3K4me3/55858490*1000000
H3K4me1_normalize=H3K4me1/62073675*1000000
H3K27ac_noramlize=H3K27ac/42343594*1000000

#####################负链反转
lncNRA_TF_direc=read.table('D:\\CRC_lncRNA\\diffexp\\5151_lncRNA_TF_direc.bed',sep='\t',stringsAsFactors = F,check.names = F)
lncNRA_TF_direc_gene=lncNRA_TF_direc[lncNRA_TF_direc[,2]=='-',1]
lncNRA_TF_direc_gen_pos=lncNRA_TF_direc[lncNRA_TF_direc[,2]=='+',1]

H3K4me3_normalize_direc=H3K4me3_normalize[rownames(H3K4me3_normalize)%in%lncNRA_TF_direc_gene,]
H3K4me3_normalize_direc_rev=rev(H3K4me3_normalize_direc)
H3K4me3_normalize_direc_pos=H3K4me3_normalize[rownames(H3K4me3_normalize)%in%lncNRA_TF_direc_gen_pos,]
H3K4me3_normalize_direc_pos_rev=rbind(H3K4me3_normalize_direc_pos,H3K4me3_normalize_direc_rev)

H3K4me1_normalize_direc=H3K4me1_normalize[rownames(H3K4me1_normalize)%in%lncNRA_TF_direc_gene,]
H3K4me1_normalize_direc_rev=rev(H3K4me1_normalize_direc)
H3K4me1_normalize_direc_pos=H3K4me1_normalize[rownames(H3K4me1_normalize)%in%lncNRA_TF_direc_gen_pos,]
H3K4me1_normalize_direc_pos_rev=rbind(H3K4me1_normalize_direc_pos,H3K4me1_normalize_direc_rev)

H3K27ac_normalize_direc=H3K27ac_noramlize[rownames(H3K27ac_noramlize)%in%lncNRA_TF_direc_gene,]
H3K27ac_normalize_direc_rev=rev(H3K27ac_normalize_direc)
H3K27ac_normalize_direc_pos=H3K27ac_noramlize[rownames(H3K27ac_noramlize)%in%lncNRA_TF_direc_gen_pos,]
H3K27ac_normalize_direc_pos_rev=rbind(H3K27ac_normalize_direc_pos,H3K27ac_normalize_direc_rev)
######################

H3K4me3_rowsum= rowSums(H3K4me3_normalize_direc_pos_rev)
H3K4me1_rowsum= rowSums(H3K4me1_normalize_direc_pos_rev)

H3K4me3_H3K4me1_rowsum=rbind(H3K4me3_rowsum,H3K4me1_rowsum)
#区分plncRNA、elncRNA
plncRNA=colnames(H3K4me3_H3K4me1_rowsum[,H3K4me3_H3K4me1_rowsum[1,]>H3K4me3_H3K4me1_rowsum[2,]])
elncRNA=colnames(H3K4me3_H3K4me1_rowsum[,which(H3K4me3_H3K4me1_rowsum[1,]<H3K4me3_H3K4me1_rowsum[2,])])
#合并三个组蛋白
allH3K_normalize=cbind(H3K4me3_normalize_direc_pos_rev,H3K4me1_normalize_direc_pos_rev)
allH3K_normalize=cbind(allH3K_normalize,H3K27ac_normalize_direc_pos_rev)

H3K4me3_normalize_plncRNA=allH3K_normalize[rownames(allH3K_normalize)%in%plncRNA,]
H3K4me3_normalize_elncRNA=allH3K_normalize[rownames(allH3K_normalize)%in%elncRNA,]
#elncNRA在上，plncRNA在下
H3K4me3_normalize_p_elncRNA=rbind(H3K4me3_normalize_elncRNA,H3K4me3_normalize_plncRNA)

write.table(H3K4me3_normalize_p_elncRNA,'D:\\CRC_lncRNA\\chipseq\\H3K4me3_normalize_p_elncRNA.txt',quote=F,col.names = F)
library(pheatmap)
pheatmap(H3K4me3_normalize_p_elncRNA,gaps_row=1182,gaps_col=c(801,1601),cluster_cols = F,cluster_rows =F,show_rownames = F,show_colnames = F,
         color = colorRampPalette(c("white", "red"))(100))

pheatmap(H3K4me3_normalize_p_elncRNA[1:1000,1:1600],cluster_cols = F,cluster_rows =F,show_rownames = F,show_colnames = F,
         color = colorRampPalette(c("white", "red"))(1000))

test=matrix(as.numeric(unlist(H3K4me3_normalize_p_elncRNA[1:100,1:800])),ncol=100)

