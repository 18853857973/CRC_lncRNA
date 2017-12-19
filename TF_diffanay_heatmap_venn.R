setwd('D:\\CRC_lncRNA\\diffexp')
library(DESeq2)
library(limma)
library(pasilla)
library(VennDiagram)
library(edgeR)
library(ggplot2)
library(ggpubr)
library(ggsignif)


#肿瘤的vs正常的：上调和下调

countData_normal=read.table('data_normal.txt',sep='\t',header = T,stringsAsFactors = F)
colData_nomal=data.frame(sample=colnames(countData_normal),Type=c(rep('normal',20),rep('tumor',20)))
colData=colData_nomal
countData_normal[is.na(countData_normal)] <- 0
keep <- rowSums(countData_normal>0) >= 0 #a Count>0 in at least 3 samples
countData_normal <- countData_normal[keep,]
type_level <- levels(colData$Type)
comb <- combn(type_level,2)
dds <- DESeqDataSetFromMatrix(countData = countData_normal,
                              colData = colData,
                              design = ~ Type)

dds2 <- DESeq(dds) 
resultsNames(dds2)
res_normal <- results(dds2)
res_normal <- as.data.frame(res_normal)
res_normal$log2FoldChange <- -res_normal$log2FoldChange 
summary(res_normal)
lfc=1
pval=0.05
args='tumor_vs_normal'

res_down_normal <- subset(res_normal,log2FoldChange >= lfc& pvalue <= pval)
res_up_normal <- subset(res_normal, log2FoldChange < -lfc&pvalue < pval)
res_normal<- rbind(res_up_normal,res_down_normal)


write.table(res_normal, file = paste('D:\\CRC_lncRNA\\diffexp\\',args,"_lfc_", lfc, "_pval_",pval, ".deseq.xls",sep = ""),sep = "\t", quote = FALSE)
write.table(res_up_normal, file = paste('D:\\CRC_lncRNA\\diffexp\\',args,"_lfc_", lfc, "_pval_",pval, ".deseq.up_regulate.xls",sep = ""),sep = "\t", quote = FALSE)
write.table(res_down_normal, file = paste('D:\\CRC_lncRNA\\diffexp\\',args,"_lfc_", lfc, "_pval_", pval, ".deseq.down_regulate.xls",sep = ""),sep = "\t", quote = FALSE)

#edgR
group_list_normal<- factor(c(rep('normal',20),rep('tumor',20)))
exprSet_normal<- countData_normal
exprSet_normal <- exprSet_normal[rowSums(cpm(exprSet_normal) > 1) >= 2,]
exprSet_normal <- DGEList(counts = exprSet_normal, group = group_list_normal)
exprSet_normal <- calcNormFactors(exprSet_normal) #对因子矫正#
exprSet_normal <- estimateCommonDisp(exprSet_normal) #估计变异系数，即估计方差；估计内部差异程度，看组间差异是否比内部差异大，如果大，可选为差异基因#
exprSet_normal <- estimateTagwiseDisp(exprSet_normal)
et <- exactTest(exprSet_normal)
tTag <- topTags(et, n=nrow(exprSet_normal))
diff_gene_edgeR_normal <- subset(tTag$table,  PValue < 0.05 & (logFC > 1 | logFC < -1))
diff_gene_edgeR_normal <- row.names(diff_gene_edgeR_normal)

#肿瘤相对正常的
diff_gene_edgeR_up_normal <- subset(tTag$table,  PValue < 0.05 & logFC > 1 )
diff_gene_edgeR_down_normal <- subset(tTag$table,  PValue < 0.05 & logFC < -1 )
diff_edgeR_normal=rbind(diff_gene_edgeR_up_normal,diff_gene_edgeR_down_normal)


write.csv(diff_gene_edgeR_up_normal,file = "up_PValue0.05_diff_gene_edgeR_tumor_vs_normal_edgeR.csv")
write.csv(diff_gene_edgeR_down_normal,file = "down_PValue0.05_diff_gene_edgeR_tumor_vs_normal_edgeR.csv")
write.csv(diff_edgeR_normal,file = "up_down_PValue0.05_diff_gene_edgeR_tumor_vs_normal_edgeR.csv")

normal_intersect_up=intersect(rownames(res_up_normal),rownames(diff_gene_edgeR_up_normal)) 
normal_intersect_down=intersect(rownames(res_down_normal),rownames(diff_gene_edgeR_down_normal))
normal_intersect_up_down=union(normal_intersect_up,normal_intersect_down)

write.table(normal_intersect_down,'D:\\CRC_lncRNA\\diffexp\\tumor_vs_normal_DESeq2_edgeR_intersect_down.txt',sep='\t',col.names = F,row.names = F,quote = F)
write.table(normal_intersect_up,'D:\\CRC_lncRNA\\diffexp\\tumor_vs_normal_DESeq2_edgeR_intersect_up.txt',sep='\t',col.names = F,row.names = F,quote = F)

normal_intersect_DESeq2_edgeR_up_data=res_normal[rownames(res_normal)%in%normal_intersect_up,]
normal_intersect_DESeq2_edgeR_down_data=res_normal[rownames(res_normal)%in%normal_intersect_down,]

write.table(normal_intersect_DESeq2_edgeR_up_data,'D:\\CRC_lncRNA\\diffexp\\df_tumor_vs_normal_DESeq2_edgeR_res_intersect_up.txt',sep='\t',quote = F)
write.table(normal_intersect_DESeq2_edgeR_down_data,'D:\\CRC_lncRNA\\diffexp\\df_tumor_vs_normal_DESeq2_edgeR_res_intersect_down.txt',sep='\t',quote = F)


normal_intersect_up_data=countData_normal[rownames(countData_normal)%in%normal_intersect_up,]
normal_intersect_down_data=countData_normal[rownames(countData_normal)%in%normal_intersect_down,]
normal_intersect_data=rbind(normal_intersect_up_data,normal_intersect_down_data)
dim(normal_intersect_up_data)
dim(normal_intersect_down_data)
dim(normal_intersect_data)

#ven图
#DESeq&edgR
#deseq_edgr_normal_tumor

venn.diagram(list(DESeq2=rownames(res_normal),edgeR=rownames(diff_edgeR_normal)),cat.cex=c(1.6,1.6),lwd=c(1,1),cex=2,fill=c("red","blue"),"D:\\CRC_lncRNA\\diffexp\\normal_tumor_DESeq_edgR.pdf")

#DESeq_edgR的基因做交集
normal_DESeq_edgR_intersect=intersect(rownames(res_normal),rownames(diff_edgeR_normal))
write(normal_DESeq_edgR_intersect,'normal_DESeq_edgR_intersect_gene.txt',sep='\t')

#heatmap
upregulateMatrix=normal_intersect_data
sampleInfo=data.frame(colnames(normal_intersect_data),Subset=group_list_normal)
colnum=2
pdf(file="D:\\CRC_lncRNA\\diffexp\\normal_tumor_heatmap.pdf")
source('D:\\R\\heatmap.R')
dev.off()

## 柱状图
novel_lncRNA=read.table('D:\\CRC_lncRNA\\filter\\lncRNA\\novel_geneid.txt',sep='\t',stringsAsFactors = F)
novel_lncRNA_genesymbol=novel_lncRNA[,1]

known_lncRNA=read.table('D:\\CRC_lncRNA\\filter\\lncRNA\\known_geneid.txt',sep='\t',stringsAsFactors = F)
known_lncRNA_genesymbol=known_lncRNA[,1]

normal_df_gene=rownames(normal_intersect_data) #3785
rec_df_gene=rownames(rec_intersect_data) #503
normal_rec_gene_interset=intersect(normal_df_gene,rec_df_gene) #199

write.table(normal_df_gene,'normal_df_gene.txt',quote = F,row.names = F,col.names = F)
write.table(rec_df_gene,'rec_df_gene.txt',quote = F,row.names = F,col.names = F)
write.table(normal_rec_gene_interset,'normal_rec_gene_interset.txt',quote = F,row.names = F,col.names = F)

#柱状图

#normal in  novel/known
normal_novel=length(intersect(normal_df_gene,novel_lncRNA_genesymbol))
normal_known=length(intersect(normal_df_gene,known_lncRNA_genesymbol))
# normal_bar=c(rep('novel',normal_novel),rep('known',normal_known))
# normal_bar_df=data.frame(normal=normal_bar)
# per_normal_novel=paste(normal_novel,length(novel_lncRNA_genesymbol),sep='/')
# per_normal_known=paste(normal_known,length(known_lncRNA_genesymbol),sep='/')
# per_normal_novel=normal_novel/length(novel_lncRNA_genesymbol)
# per_normal_known=normal_known/length(known_lncRNA_genesymbol)
no_df_novel=length(novel_lncRNA_genesymbol)-normal_novel
no_df_known=length(known_lncRNA_genesymbol)-normal_known

x1 = matrix(c(normal_novel,normal_known,no_df_novel,no_df_known),nc = 2 , byrow = T)
chisq.test(x1)

noraml_bar_df=data.frame(normal=c(rep("novel",length(novel_lncRNA_genesymbol)),rep("known",length(known_lncRNA_genesymbol))),normal_val=c(rep('DF',normal_novel),rep('NON_DF',no_df_novel),rep('DF',normal_known),rep('NON_DF',no_df_known)))
#,per=c(per_normal_novel,per_normal_known)

pdf(file='D:\\CRC_lncRNA\\diffexp\\tumorigenesis_in_novel_known_lncRNA_bar_percent.pdf')
sp=ggplot(noraml_bar_df,aes(normal,fill=factor(normal_val))) + geom_bar(position='fill',width=0.5)+labs(x="",y="percent")+ggtitle("tumorigenesis_in_novel_known")
sp+theme_bw() + theme(title=element_text(size=15,color="black"
                                         ),plot.title = element_text(hjust = 0.5),legend.title=element_blank(),panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),axis.title.x = element_text(size = 20, face = "bold"),axis.title.y= element_text(size = 30, face = "bold"),axis.text.x=element_text(size=25,color="black"))

dat <- data.frame(Group = c("known", "known", "novel", "novel"),
                  Sub   = c("DF", "NON_DF", "DF", "NON_DF"),
                  Value = c(normal_known,no_df_known,normal_novel,no_df_novel))  



# ggplot(dat, aes(Group, Value)) +
#   geom_bar(aes(fill = Sub), stat="identity", position="fill", width=.5) +geom_signif(y_position=1,comparisons = list(c("known", "novel")), 
#                                                                                      map_signif_level=TRUE)



dev.off()



#rec in  novel/known
rec_novel=length(intersect(rec_df_gene,novel_lncRNA_genesymbol))
rec_known=length(intersect(rec_df_gene,known_lncRNA_genesymbol))
# per_rec_novel=paste(rec_novel,length(novel_lncRNA_genesymbol),sep='/')
# per_rec_known=paste(rec_known,length(known_lncRNA_genesymbol),sep='/')
no_rec_df_novel=length(novel_lncRNA_genesymbol)-rec_novel
no_rec_df_known=length(known_lncRNA_genesymbol)-rec_known

x2 = matrix(c(rec_novel,rec_known,no_rec_df_novel,no_rec_df_known),nc = 2 , byrow = T)
chisq.test(x2)



rec_bar_df=data.frame(normal=c(rep("novel",length(novel_lncRNA_genesymbol)),rep("known",length(known_lncRNA_genesymbol))),normal_val=c(rep('DF',rec_novel),rep('NON_DF',no_rec_df_novel),rep('DF',rec_known),rep('NON_DF',no_rec_df_known)))
pdf(file='D:\\CRC_lncRNA\\diffexp\\metastasis_in_novel_known_lncRNA_bar_precent.pdf')
# rec_bar=c(rep('novel',rec_novel),rep('known',rec_known))
# rec_bar_df=data.frame(rec=rec_bar)
sp=ggplot(rec_bar_df,aes(normal,fill=factor(normal_val))) + geom_bar(position='fill',width=0.5)+labs(x="",y="percent")+ggtitle("metastasis_in_novel_known")#+geom_text(aes(x = 1.5, y = 1.1), label = "***", size = 5)
sp2=sp+theme_bw() + theme(title=element_text(size=15,color="black"
),plot.title = element_text(hjust = 0.5),legend.title=element_blank(),panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),axis.title.x = element_text(size = 20, face = "bold"),axis.title.y= element_text(size = 30, face = "bold"),axis.text.x=element_text(size=25,color="black"))
sp2

geom_signif(comparisons = list(c("known", "novel")))

dev.off()



#normal_rec in novel /known
png(file='tumorigenesis_metastasis_in_novel_known_lncRNA_bar_percent.png',bg="transparent")
normal_rec_novel=length(intersect(normal_rec_gene_interset,novel_lncRNA_genesymbol))
normal_rec_known=length(intersect(normal_rec_gene_interset,known_lncRNA_genesymbol))
# per_normal_rec_novel=paste(normal_rec_novel,length(novel_lncRNA_genesymbol),sep='/')
# per_normal_rec_known=paste(normal_rec_known,length(known_lncRNA_genesymbol),sep='/')
# normal_rec_bar=c(rep('novel',normal_rec_novel),rep('known',normal_rec_known))
# normal_rec_bar_df=data.frame(normal_rec=normal_rec_bar)
no_normal_rec_df_novel=length(novel_lncRNA_genesymbol)-normal_rec_novel
no_normal_rec_df_known=length(known_lncRNA_genesymbol)-normal_rec_known


x3= matrix(c(normal_rec_novel,normal_rec_known,no_normal_rec_df_novel,no_normal_rec_df_known),nc = 2 , byrow = T)
chisq.test(x3)


normal_rec_bar_df=data.frame(normal=c(rep("novel",length(novel_lncRNA_genesymbol)),rep("known",length(known_lncRNA_genesymbol))),normal_val=c(rep('DF',normal_rec_novel),rep('NON_DF',no_normal_rec_df_novel),rep('DF',normal_rec_known),rep('NON_DF',no_normal_rec_df_known)))
pdf(file='D:\\CRC_lncRNA\\diffexp\\metastasis&tumorigenesis_in_novel_known_lncRNA_bar_precent.pdf')
# rec_bar=c(rep('novel',rec_novel),rep('known',rec_known))
# rec_bar_df=data.frame(rec=rec_bar)
sp=ggplot(normal_rec_bar_df,aes(normal,fill=factor(normal_val))) + geom_bar(position='fill',width=0.5)+labs(x="",y="percent")+ggtitle("metastasis&tumorigenesis_in_novel_known")
sp+theme_bw() + theme(title=element_text(size=15,color="black"
),plot.title = element_text(hjust = 0.5),legend.title=element_blank(),panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),axis.title.x = element_text(size = 20, face = "bold"),axis.title.y= element_text(size = 30, face = "bold"),axis.text.x=element_text(size=25,color="black"))
dev.off()




#################################

dflncRNA=read.table('D:\\CRC_lncRNA\\filter\\lncRNA\\diferent_lncRNA_three_type.txt',sep='\t')
dflncRNA=dflncRNA[-1,]

normal_tumor_lncRNA=dflncRNA[,1]
gene_cluster=normal_tumor_lncRNA
#source("D:\\R\\clusterProfiler_genelist.R")




