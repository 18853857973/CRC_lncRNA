library(edgeR)
##跟DESeq2一样，导入数据，预处理（用了cpm函数）

# exprSet2=countData[,c(1:10,31:40,11:30)]
# group_list <- factor(c(rep('normal',20),rep('tumor',20)))

# exprSet2=countData[,c(11:30,1:10,31:40)]
# group_list <- factor(c(rep('tumor',20),rep('normal',20)))

group_list<- factor(c(rep('rec',10),rep('norec',10)))

exprSet<- countData
exprSet <- exprSet[rowSums(cpm(exprSet) > 1) >= 2,]

##设置分组信息，并做TMM标准化
exprSet <- DGEList(counts = exprSet, group = group_list)
exprSet <- calcNormFactors(exprSet) #对因子矫正#

##使用qCML（quantile-adjusted conditional maximum likelihood）估计离散度（只针对单因素实验设计）
exprSet <- estimateCommonDisp(exprSet) #估计变异系数，即估计方差；估计内部差异程度，看组间差异是否比内部差异大，如果大，可选为差异基因#
exprSet <- estimateTagwiseDisp(exprSet)

##寻找差异gene(这里的exactTest函数还是基于qCML并且只针对单因素实验设计)，然后按照阈值进行筛选即可
et <- exactTest(exprSet)
tTag <- topTags(et, n=nrow(exprSet))
diff_gene_edgeR <- subset(tTag$table,  PValue < 0.05 & (logFC > 1 | logFC < -1))
#diff_gene_edgeR <- subset(tTag$table,  FDR < 0.01 & (logFC > 1 | logFC < -1))#PValue < 0.05 )
diff_gene_edgeR <- row.names(diff_gene_edgeR)

diff_gene_edgeR_up <- subset(tTag$table,  PValue < 0.05 & logFC > 1 )
diff_gene_edgeR_down <- subset(tTag$table,  PValue < 0.05 & logFC < -1 )

write.csv(diff_gene_edgeR_up,file = "up_PValue0.05_diff_gene_edgeR_rec_vs_norec_edgeR.csv")
write.csv(diff_gene_edgeR_down,file = "down_PValue0.05_diff_gene_edgeR_rec_vs_norec_edgeR.csv")





