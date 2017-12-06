setwd('D:\\CRC_lncRNA\\cnv\\percentCNV')
#全部lncRNA
all_novel=read.table('lncRNA.final.v2.novel.transcriptid.txt')
all_novel=all_novel[,1]
all_novel_num=length(all_novel)

all_known=read.table('lncRNA.final.v2.known.transcriptid.txt')
all_known=all_known[,1]
all_known_num=length(all_known)

#大于25的
per_novel=read.table('percentages25novel.transcriptid.txt')
per_novel=per_novel[,1]
per_novel_num=length(per_novel)

per_known=read.table('percentages25known.transcriptid.txt')
per_known=per_known[,1]
per_known_num=length(per_known)

#剩下的
novel_less=all_novel_num-per_novel_num
known_less=all_known_num-per_known_num

type=c(rep("novel",all_novel_num),rep("known",all_known_num))
num=c(rep("CNV",per_novel_num),rep("NON_CNV",novel_less),rep("CNV",per_known_num),rep("NON_CNV",known_less))
per_df=data.frame(type=type,num=num)

png(file='lncRNA_CNV_percent.png',bg="transparent")

sp=ggplot(per_df,aes(type,fill=factor(num))) + geom_bar(position='fill',width=0.5)+labs(x="",y="percent")+ggtitle("CNV_percent_in_novel_known")
sp+theme_bw() + theme(title=element_text(size=15,color="black"
),plot.title = element_text(hjust = 0.5),legend.title=element_blank(),panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),axis.title.x = element_text(size = 20, face = "bold"),axis.title.y= element_text(size = 30, face = "bold"),axis.text.x=element_text(size=25,color="black"))

dev.off()


############################
#全部lncRNA
all_novel=read.table('lncRNA.final.v2.novel.geneid.txt')
all_novel=all_novel[,1]
all_novel_num=length(all_novel)

all_known=read.table('lncRNA.final.v2.known.geneid.txt')
all_known=all_known[,1]
all_known_num=length(all_known)

#大于25的
per_novel=read.table('percentages25novel.geneid.txt')
per_novel=unique(per_novel[,1])
per_novel_num=length(per_novel)

per_known=read.table('percentages25known.geneid.txt')
per_known=unique(per_known[,1])
per_known_num=length(per_known)

#剩下的
novel_less=all_novel_num-per_novel_num
known_less=all_known_num-per_known_num

type=c(rep("novel",all_novel_num),rep("known",all_known_num))
num=c(rep("CNV",per_novel_num),rep("NON_CNV",novel_less),rep("CNV",per_known_num),rep("NON_CNV",known_less))
per_df=data.frame(type=type,num=num)

pdf(file='D:\\CRC_lncRNA\\cnv\\percentCNV\\lncRNA_CNV_percent_bar.pdf')

sp=ggplot(per_df,aes(type,fill=factor(num))) + geom_bar(position='fill',width=0.5)+labs(x="",y="percent")+ggtitle("CNV_percent_in_novel_known")
sp+theme_bw() + theme(title=element_text(size=15,color="black"
),plot.title = element_text(hjust = 0.5),legend.title=element_blank(),panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),axis.title.x = element_text(size = 20, face = "bold"),axis.title.y= element_text(size = 30, face = "bold"),axis.text.x=element_text(size=25,color="black"))

dev.off()
#####################################################












#venplot
#大于25的
per_novel_gene=read.table('percentages25novel.geneid.txt')
per_novel_gene=unique(per_novel_gene[,1])

per_known_gene=read.table('percentages25known.geneid.txt')
per_known_gene=unique(per_known_gene[,1])

normal_DESeq_edgR_intersect=read.table('D:\\CRC_lncRNA\\diffexp\\normal_DESeq_edgR_intersect_gene.txt',sep='\t')
normal_DESeq_edgR_intersect=normal_DESeq_edgR_intersect[,1]
rec_DESeq_edgR_intersect=read.table('D:\\CRC_lncRNA\\diffexp\\rec_DESeq_edgR_intersect_gene.txt',sep='\t')
rec_DESeq_edgR_intersect=rec_DESeq_edgR_intersect[,1]

union_per_novel_known_gene=union(per_novel_gene,per_known_gene)
venn.diagram(list(normal_tumor_differentlncRNA=normal_DESeq_edgR_intersect,CNVlncRNA=union_per_novel_known_gene),cat.cex=c(1,1),lwd=c(1,1),cex=2,fill=c("red","blue"),"D:\\CRC_lncRNA\\cnv\\normal_tumor_CNV_intersectgene.pdf")
venn.diagram(list(recornot_differentlncRNA=rec_DESeq_edgR_intersect,CNVlncRNA=union_per_novel_known_gene),cat.cex=c(1,1),lwd=c(1,1),cex=2,fill=c("red","blue"),"D:\\CRC_lncRNA\\cnv\\rec_ornot_CNV_intersectgene.pdf")

intersect_normal_cnv=(intersect(normal_DESeq_edgR_intersect,union_per_novel_known_gene))
#肿瘤vs正常的且有cnv
intersect_normal_cnv_up=intersect(intersect_normal_cnv,normal_intersect_up)
intersect_normal_cnv_down=intersect(intersect_normal_cnv,normal_intersect_down)

intersect_normal_cnv_up_down=c(intersect_normal_cnv_up,intersect_normal_cnv_down)
#找出差异基因对应的gtf文件为bed文件
intersect_normal_cnv_up_down=data.frame(lncRNA=intersect_normal_cnv_up_down)
lncRNA_gtf=read.table('D:\\CRC_lncRNA\\filter\\lncRNA\\only_min_max_position_lncRNA.final.v2.gtf',stringsAsFactors = F)
colnames(lncRNA_gtf)=c("chr","start","end","lncRNA")
intersect_normal_cnv_up_down_gtf=merge(intersect_normal_cnv_up_down,lncRNA_gtf,by='lncRNA',sort=F)
intersect_normal_cnv_up_down_gtf=intersect_normal_cnv_up_down_gtf[,c(2,3,4,1)]
write.table(intersect_normal_cnv_up_down_gtf,'D:\\CRC_lncRNA\\diffexp\\intersect_normal_cnv_up_down_gtf.bed',quote = F,col.names = F,row.names = F,sep = '\t')


intersect_normal_cnv_up_logFC=res_normal[rownames(res_normal)%in%intersect_normal_cnv_up,]
intersect_normal_cnv_up_logFC=intersect_normal_cnv_up_logFC[order(intersect_normal_cnv_up_logFC[,2],decreasing=T),]
intersect_normal_cnv_down_logFC=res_normal[rownames(res_normal)%in%intersect_normal_cnv_down,]
intersect_normal_cnv_down_logFC=intersect_normal_cnv_down_logFC[order(intersect_normal_cnv_down_logFC[,2],decreasing=F),]
write.table(intersect_normal_cnv_up_logFC,'intersect_normal_cnv_up_logFC.txt',quote = F)
write.table(intersect_normal_cnv_down_logFC,'intersect_normal_cnv_down_logFC.txt',quote = F)

intersect_normal_cnv_up_down_logFC=rbind(intersect_normal_cnv_up_logFC,intersect_normal_cnv_down_logFC)

#rec_to_not_rec
# intersect_rec_cnv=(intersect(rec_DESeq_edgR_intersect,union_per_novel_known_gene))
# intersect_rec_cnv_up=intersect(intersect_rec_cnv,rec_intersect_up)
# intersect_rec_cnv_down=intersect(intersect_rec_cnv,rec_intersect_down)
# intersect_normal_cnv_up_down=c(intersect_rec_cnv_up,intersect_rec_cnv_down)

intersect_rec_cnv_up_logFC=res_rec[rownames(res_rec)%in%intersect_normal_cnv_up,]
intersect_rec_cnv_up_logFC=intersect_rec_cnv_up_logFC[order(intersect_rec_cnv_up_logFC[,2],decreasing=T),]
intersect_rec_cnv_down_logFC=res_rec[rownames(res_rec)%in%intersect_normal_cnv_down,]
intersect_rec_cnv_down_logFC=intersect_rec_cnv_down_logFC[order(intersect_rec_cnv_down_logFC[,2],decreasing=T),]
intersect_rec_cnv_up_down_logFC=rbind(intersect_rec_cnv_up_logFC,intersect_rec_cnv_down_logFC)

intersect_normal_rec_cnv_up_down_logFC=cbind(-intersect_normal_cnv_up_down_logFC[,2],intersect_rec_cnv_up_down_logFC[,2])
colnames(intersect_normal_rec_cnv_up_down_logFC)=c("tumor-normal","rec-nonrec")

#heatmap 上为上调，下为下调



library(pheatmap)
pheatmap(intersect_normal_rec_cnv_up_down_logFC,cluster_cols = F,cluster_rows =F ,
         colorRampPalette(c("green", "black", "red"))(50),show_rownames=F,show_colnames=F)





#差异lncRNA_CNV的三种交集
venn.diagram(list(normal_differentlncRNA=normal_DESeq_edgR_intersect,rec_differentlncRNA=rec_DESeq_edgR_intersect,CNVlncRNA=union_per_novel_known_gene),cat.cex=c(1,1,1),lwd=c(1,1,1),cex=2,fill=c("red","blue","yellow"),"D:\\CRC_lncRNA\\cnv\\normal_rec_ornot_CNV_intersectgene.pdf")
lncRNA_CNV2=intersect(normal_DESeq_edgR_intersect,rec_DESeq_edgR_intersect)
lncRNA_CNV=intersect(lncRNA_CNV2,union_per_novel_known_gene)
write.table(lncRNA_CNV,'D:\\CRC_lncRNA\\cnv\\percentCNV\\normal_rec_0.25CNV_lncRNA.txt',quote = F,col.names = F,row.names = F)


#肿瘤相对于正常
# tumor_vs_normal_rec_0.25CNV_lncRNA_down=intersect(lncRNA_CNV,normal_intersect_up)
# tumor_vs_normal_rec_0.25CNV_lncRNA_up=intersect(lncRNA_CNV,normal_intersect_down)
# 
# rec_tumor_vs_normal_0.25CNV_lncRNA_up=intersect(lncRNA_CNV,rec_intersect_up)
# rec_tumor_vs_normal_0.25CNV_lncRNA_down=intersect(lncRNA_CNV,rec_intersect_down)



#取差异lncRNA和CNV 变异的交集
intersect_normal_dflncRNA_CNV_up=intersect(normal_intersect_up,union_per_novel_known_gene)
intersect_normal_dflncRNA_CNV_down=intersect(normal_intersect_down,union_per_novel_known_gene)
intersect_normal_dflncRNA_CNV_up_data=countData_normal[rownames(countData_normal)%in%intersect_normal_dflncRNA_CNV_up,]
intersect_normal_dflncRNA_CNV_down_data=countData_normal[rownames(countData_normal)%in%intersect_normal_dflncRNA_CNV_down,]
intersect_normal_dflncRNA_CNV_up_down_data=rbind(intersect_normal_dflncRNA_CNV_up_data,intersect_normal_dflncRNA_CNV_down_data)
write.table(intersect_normal_dflncRNA_CNV_up_down_data,'intersect_normal_dflncRNA_CNV_up_down_data.txt',quote = F)


intersect_rec_dflncRNA_CNV_up=intersect(rec_intersect_up,union_per_novel_known_gene)
intersect_rec_dflncRNA_CNV_down=intersect(rec_intersect_down,union_per_novel_known_gene)
intersect_rec_dflncRNA_CNV_up_data=countData_rec[rownames(countData_rec)%in%intersect_rec_dflncRNA_CNV_up,]
intersect_rec_dflncRNA_CNV_down_data=countData_rec[rownames(countData_rec)%in%intersect_rec_dflncRNA_CNV_down,]
intersect_rec_dflncRNA_CNV_up_down_data=rbind(intersect_rec_dflncRNA_CNV_up_data,intersect_rec_dflncRNA_CNV_down_data)
write.table(intersect_rec_dflncRNA_CNV_up_down_data,'intersect_rec_dflncRNA_CNV_up_down_data.txt',quote = F)



#heatmap
upregulateMatrix=intersect_normal_dflncRNA_CNV_up_down_data
sampleInfo=data.frame(colnames(normal_intersect_data),Subset=group_list_normal)
colnum=2
source('D:\\R\\heatmap.R')


upregulateMatrix=intersect_rec_dflncRNA_CNV_up_down_data
sampleInfo=data.frame(colnames(normal_intersect_data),Subset=group_list_rec)
colnum=2
source('D:\\R\\heatmap.R')



#############
#上下调基因和cnv 的Amp和Del，卡方检验
setwd('D:\\CRC_lncRNA\\cnv\\percentCNV')
res_up_normal=read.table("D:\\CRC_lncRNA\\diffexp\\tumor_vs_normalType_tumor_vs_normal_lfc_1_pval_0.05.deseq.up_regulate.xls",sep='\t')
res_down_normal=read.table("D:\\CRC_lncRNA\\diffexp\\tumor_vs_normalType_tumor_vs_normal_lfc_1_pval_0.05.deseq.down_regulate.xls",sep='\t')

diff_gene_edgeR_up_normal=read.csv( "D:\\CRC_lncRNA\\diffexp\\up_PValue0.05_diff_gene_edgeR_normal_vs_tumor_edgeR.csv",header=T,row.names = 1)
diff_gene_edgeR_down_normal=read.csv( "D:\\CRC_lncRNA\\diffexp\\down_PValue0.05_diff_gene_edgeR_normal_vs_tumor_edgeR.csv",header=T,row.names = 1)

normal_intersect_up=intersect(rownames(res_up_normal),rownames(diff_gene_edgeR_down_normal)) 
normal_intersect_down=intersect(rownames(res_down_normal),rownames(diff_gene_edgeR_up_normal))

normal_intersect_up_length=length(normal_intersect_up)
normal_intersect_down_length=length(normal_intersect_down)

cnv_known_novel_Amp_geneid=read.table('D:\\CRC_lncRNA\\cnv\\percentCNV\\percentages25.Amp.geneid.txt')
cnv_known_novel_Amp_geneid=cnv_known_novel_Amp_geneid[,1]

cnv_known_novel_Del_geneid=read.table('D:\\CRC_lncRNA\\cnv\\percentCNV\\percentages25.Del.geneid.txt')
cnv_known_novel_Del_geneid=cnv_known_novel_Del_geneid[,1]

cnv_known_novel_Amp_geneid_length=length(cnv_known_novel_Amp_geneid)
cnv_known_novel_Del_geneid_length=length(cnv_known_novel_Del_geneid)

x5 = matrix(c(cnv_known_novel_Amp_geneid_length,cnv_known_novel_Del_geneid_length,normal_intersect_up_length,normal_intersect_down_length),nc = 2 , byrow = T)
chisq.test(x5)

#p-value = 3.072e-13




diff_gene_edgeR_up_rec=read.csv("D:\\CRC_lncRNA\\diffexp\\up_PValue0.05_diff_gene_edgeR_rec_vs_norec_edgeR.csv")
diff_gene_edgeR_down_rec=read.csv("D:\\CRC_lncRNA\\diffexp\\down_PValue0.05_diff_gene_edgeR_rec_vs_norec_edgeR.csv")

