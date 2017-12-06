#########################################################################
#复发相对于未复发:上下调：B_to_C

countData_rec=read.table('lncRNA.rsem.count_sort_rec_not.txt',sep='\t',header = T,stringsAsFactors = F)
colData=data.frame(sample=colnames(countData_rec),Type=c(rep('recu',10),rep('unrecu',10)))
countData=countData_rec
countData[is.na(countData)] <- 0
keep <- rowSums(countData>0) >= 0 #a Count>0 in at least 3 samples
countData <- countData[keep,]

type_level <- levels(colData$Type)
comb <- combn(type_level,2)

dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ Type)

dds2 <- DESeq(dds) 
resultsNames(dds2)
res_rec <- results(dds2)
res_rec <- as.data.frame(res_rec)
res_rec$log2FoldChange <- -res_rec$log2FoldChange 
summary(res_rec)
lfc=1
pval=0.05
args='rec_vs_norec'
res_up_rec <- subset(res_rec,log2FoldChange > lfc& pvalue < pval)
res_down_rec <- subset(res_rec, log2FoldChange < -lfc&pvalue < pval)
res_rec<- rbind(res_up_rec,res_down_rec)


write.table(res_rec, file = paste(args,resultsNames(dds2)[2],"_lfc_", lfc, "_pval_",pval, ".deseq.xls",sep = ""),sep = "\t", quote = FALSE)
write.table(res_up_rec, file = paste(args,resultsNames(dds2)[2],"_lfc_", lfc, "_pval_",pval, ".deseq.up_regulate.xls",sep = ""),sep = "\t", quote = FALSE)
write.table(res_down_rec, file = paste(args,resultsNames(dds2)[2],"_lfc_", lfc, "_pval_", pval, ".deseq.down_regulate.xls",sep = ""),sep = "\t", quote = FALSE)

#edgR
group_list_rec<- factor(c(rep('rec',10),rep('norec',10)))
exprSet<- countData
exprSet <- exprSet[rowSums(cpm(exprSet) > 1) >= 2,]
##设置分组信息，并做TMM标准化
exprSet <- DGEList(counts = exprSet, group = group_list_rec)
exprSet <- calcNormFactors(exprSet) #对因子矫正#
##使用qCML（quantile-adjusted conditional maximum likelihood）估计离散度（只针对单因素实验设计）
exprSet <- estimateCommonDisp(exprSet) #估计变异系数，即估计方差；估计内部差异程度，看组间差异是否比内部差异大，如果大，可选为差异基因#
exprSet <- estimateTagwiseDisp(exprSet)
##寻找差异gene(这里的exactTest函数还是基于qCML并且只针对单因素实验设计)，然后按照阈值进行筛选即可
et <- exactTest(exprSet)
tTag <- topTags(et, n=nrow(exprSet))
diff_gene_edgeR_rec <- subset(tTag$table,  PValue < 0.05 & (logFC > 1 | logFC < -1))
#diff_gene_edgeR <- subset(tTag$table,  FDR < 0.01 & (logFC > 1 | logFC < -1))#PValue < 0.05 )
diff_gene_edgeR_rec <- row.names(diff_gene_edgeR_rec)

diff_gene_edgeR_up_rec <- subset(tTag$table,  PValue < 0.05 & logFC > 1 )
diff_gene_edgeR_down_rec <- subset(tTag$table,  PValue < 0.05 & logFC < -1 )
diff_edgeR_rec=rbind(diff_gene_edgeR_up_rec,diff_gene_edgeR_down_rec)

write.csv(diff_gene_edgeR_up_rec,file = "up_PValue0.05_diff_gene_edgeR_rec_vs_norec_edgeR.csv")
write.csv(diff_gene_edgeR_down_rec,file = "down_PValue0.05_diff_gene_edgeR_rec_vs_norec_edgeR.csv")
write.csv(diff_edgeR_rec,file = "up_down_PValue0.05_diff_gene_edgeR_rec_vs_norec_edgeR.csv")


rec_intersect_up=intersect(rownames(res_up_rec),rownames(diff_gene_edgeR_up_rec))
rec_intersect_down=intersect(rownames(res_down_rec),rownames(diff_gene_edgeR_down_rec))
rec_intersect_up_down=union(rec_intersect_up,rec_intersect_down)

#复发未复发两种差异表达方法差异基因交集（上下调）
write.table(rec_intersect_up,'D:\\CRC_lncRNA\\diffexp\\rec_DESeq2_edgeR_intersect_up.txt',sep='\t',col.names = F,row.names = F,quote = F)
write.table(rec_intersect_down,'D:\\CRC_lncRNA\\diffexp\\rec_DESeq2_edgeR_intersect_down.txt',sep='\t',col.names = F,row.names = F,quote = F)

rec_intersect_DESeq2_edgeR_up_data=res_rec[rownames(res_rec)%in%rec_intersect_up,]
rec_intersect_DESeq2_edgeR_down_data=res_rec[rownames(res_rec)%in%rec_intersect_down,]

write.table(rec_intersect_DESeq2_edgeR_up_data,'D:\\CRC_lncRNA\\diffexp\\rec_DESeq2_edgeR_res_intersect_up.txt',sep='\t',quote = F)
write.table(rec_intersect_DESeq2_edgeR_down_data,'D:\\CRC_lncRNA\\diffexp\\rec_DESeq2_edgeR_res_intersect_down.txt',sep='\t',quote = F)




rec_intersect_up_data=countData_rec[rownames(countData_rec)%in%rec_intersect_up,]
rec_intersect_down_data=countData_rec[rownames(countData_rec)%in%rec_intersect_down,]
rec_intersect_data=rbind(rec_intersect_up_data,rec_intersect_down_data)


#deseq_edgr_recornot
venn.diagram(list(DESeq2=rownames(res_rec),edgeR=rownames(diff_edgeR_rec)),cat.cex=c(1.6,1.6),sub.cex = 10,lwd=c(1,1),cex=2,fill=c("red","blue"),"D:\\CRC_lncRNA\\diffexp\\rec_ornot_DESeq_edgR.pdf")


#DESeq_edgR的基因做交集
rec_DESeq_edgR_intersect=intersect(rownames(res_rec),rownames(diff_edgeR_rec))
write(rec_DESeq_edgR_intersect,'rec_DESeq_edgR_intersect_gene.txt',sep='\t')


#正常和肿瘤以及复发未复发交集
rec_normal_intersect_up=intersect(rec_intersect_up,normal_intersect_down)
rec_normal_intersect_down=intersect(rec_intersect_down,normal_intersect_up)



upregulateMatrix=rec_intersect_data
sampleInfo=data.frame(colnames(rec_intersect_data),Subset=group_list_rec)
colnum=2
pdf(file="D:\\CRC_lncRNA\\diffexp\\rec_ornot_heatmap.pdf")

source('D:\\R\\heatmap.R')

dev.off()



#正常和复发交集ven
venn.diagram(list(metastasis =rownames(rec_intersect_data),tumorigenesis =rownames(normal_intersect_data)),cat.cex=c(1,1),lwd=c(1,1),cex=2,fill=c("red","blue"),"D:\\CRC_lncRNA\\diffexp\\rec_normal.pdf")

#正常和复发上下调分别做交集
normal_tumor_interset_up=intersect(rownames(rec_intersect_up_data),rownames(normal_intersect_up_data))
normal_tumor_interset_down=intersect(rownames(rec_intersect_down_data),rownames(normal_intersect_down_data))

normal_tumor_interset_up_data=countData_rec[rownames(countData_rec)%in%normal_tumor_interset_up,]
normal_tumor_interset_down_data=countData_rec[rownames(countData_rec)%in%normal_tumor_interset_down,]
normal_tumor_interset_data=rbind(normal_tumor_interset_up_data,normal_tumor_interset_down_data)

#heatmap
upregulateMatrix=normal_tumor_interset_data
pdf(file="D:\\CRC_lncRNA\\diffexp\\rec_ornot_heatmap.pdf")

sampleInfo=data.frame(colnames(rec_intersect_data),Subset=group_list_rec)
colnum=2
source('D:\\R\\heatmap.R')
dev.off()


