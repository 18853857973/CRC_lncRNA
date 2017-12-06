setwd('D:\\CRC_lncRNA\\diffexp')
# data=read.table('D:\\CRC_lncRNA\\filter\\RSEM_expression\\lncRNA.rsem.count.txt',sep='\t',stringsAsFactors=F,header=T)
# samplename=strsplit(colnames(data),split="_",fixed=T)
# rown=c()
# for (n in (1:length(samplename))){
#   g=samplename[[n]][1]
#   rown=c(rown,g)
# }
# colnames(data)=rown
# countData=data[,sort(rown)]
# write.table(countData,'lncRNA.rsem.count_sort.txt',quote = F,sep='\t')


data=read.table('lncRNA.rsem.count_sort_interger.txt',sep='\t',stringsAsFactors = F)
dim(data)
data=data[rowSums(data)>0,]
dim(data)

#colData=data.frame(sample=colnames(data),Type=c(rep('recu_normal',10),rep('recu_tumor',10),rep('norecu_tumor',10),rep('norecu_normal',10)))
#复发和未复发

#正常和肿瘤
colData=data.frame(sample=colnames(data),Type=c(rep('normal',10),rep('tumor',10),rep('tumor',10),rep('normal',10)))
colData2=colData[c(11:30,1:10,31:40),]
rownames(colData2)=c(1:40)
colData=colData2

countData_recu=data
countData2=matrix(as.integer(unlist(countData_recu)),ncol=40)
countData2=data.frame(countData2)
row.names(countData2)=row.names(countData_recu)
colnames(countData2)=colnames(countData_recu)

countData2[is.na(countData2)] <- 0
countData=countData2
save(countData,file='countData.RData')
write.table(countData,'lncRNA.rsem.count_sort_interger.txt',sep='\t',quote = F)

library("DESeq2")
library("ggplot2")

lfc = 1
pval = 0.01
args='normal_tumor'
# List all possible pared comparision
type_level <- levels(colData$Type)
comb <- combn(type_level,2)

#i = 1
for (i in 1:length(comb[1,])){
  # Extract specific info for comparision
  colData1 <- subset(colData, Type == comb[1,i] | Type == comb[2,i])
  countData1 <- countData[,colData1[,1]]
  keep <- rowSums(countData1>0) >= 3 #a Count>0 in at least 3 samples
  countData1 <- countData1[keep,]
  # DESeq
  dds <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~ Type)
  dds <- DESeq(dds)
  resultsNames(dds)
  res <- results(dds)
  res <- as.data.frame(res)
  res$log2FoldChange <- -res$log2FoldChange  # -logFC
  res_up <- subset(res, log2FoldChange > lfc & padj < pval)
  res_down <- subset(res, log2FoldChange < -lfc & padj < pval)
  write.table(res, file = paste(args,"_lfc_", lfc, "_pval_",pval, "_", resultsNames(dds)[2], ".deseq.xls",sep = ""),sep = "\t", quote = FALSE)
  write.table(res_up, file = paste(args,"_lfc_", lfc, "_pval_",pval, "_", resultsNames(dds)[2], ".deseq.up_regulate.xls",sep = ""),sep = "\t", quote = FALSE)
  write.table(res_down, file = paste(args,"_lfc_", lfc, "_pval_",pval, "_", resultsNames(dds)[2], ".deseq.down_regulate.xls",sep = ""),sep = "\t", quote = FALSE)
  res <- na.omit(res)
  # Vocano plot
  pdf(file = paste(args,"_", comb[1,i], "_vs_", comb[2,i], ".deseq.VocanoPlot.pdf", sep=""))
  par(mar = c(5, 6, 5, 5))
  tab = data.frame(logFC = res$log2FoldChange, negLogPval = -log10(res$padj)) 
  res$baseMean[res$baseMean>5000]=5000
  res$baseMean[res$baseMean<10]=10
  nosigGene = (abs(tab$logFC) < lfc | tab$negLogPval < -log10(pval))
  signGenes_up = (tab$logFC > lfc & tab$negLogPval > -log10(pval))
  signGenes_down = (tab$logFC < -lfc & tab$negLogPval > -log10(pval))
  gap = max(tab$logFC)/50
  
  
  plot(tab, pch = 21, xlab = expression(log[2]~fold~change), ylab = expression(-log[10]~pvalue), cex.lab = 1.5, col = alpha("black", 0))
  points(tab[nosigGene, ], pch = 21, xlab = expression(log[2]~fold~change), ylab = expression(-log[10]~pvalue), col = "black", bg = "grey")
  if (length(unique(signGenes_up)) > 1){
    points(tab[signGenes_up, ], pch = 21, col = "black", bg = "red") 
  }
  if (length(unique(signGenes_down)) > 1){
    points(tab[signGenes_down, ], pch = 21, col = "black", bg = "cornflowerblue") 
  }
  abline(h = -log10(pval), col = "green3", lty = 2) 
  abline(v = c(-lfc, lfc), col = "orange", lty = 2) 
  if (length(unique(signGenes_up)) > 1){
    text(tab[signGenes_up, ]$logFC, tab[signGenes_up, ]$negLogPval+gap, row.names(res[signGenes_up,]), cex = 0.5, col = "red")
  }
  if (length(unique(signGenes_down)) > 1){
    text(tab[signGenes_down, ]$logFC, tab[signGenes_down, ]$negLogPval+gap, row.names(res[signGenes_down,]), cex = 0.5, col = "blue")
  }
  mtext(paste("pval =", pval), side = 4, at = -log10(pval), cex = 0.8, line = 0.5, las = 1) 
  mtext(c(paste("-", 2, "fold"), paste("+", 2, "fold")), side = 3, at = c(-lfc, lfc), cex = 0.8, line = 0.5)
  mtext(c(comb[1,i], comb[2,i]), side = 3, at = c(3*lfc, -3*lfc), cex = 1, line=2)
  legend("top",legend = c("Upregulate","Downregulate"),pch = c(16, 16), col = c("red", "cornflowerblue"))
  
  plot(tab, pch = 16, cex = 0.6, xlab = expression(log[2]~fold~change), ylab = expression(-log[10]~pvalue), cex.lab = 1.5, col= alpha("black", 0))
  points(tab[nosigGene, ], pch = 16, cex = res$baseMean/1000, xlab = expression(log[2]~fold~change), ylab = expression(-log[10]~pvalue), col= alpha("black", 0.1))
  if (length(unique(signGenes_up)) > 1){
    points(tab[signGenes_up, ], pch = 16, cex = res$baseMean/1000, col = alpha("red", 0.4)) 
  }
  if (length(unique(signGenes_down)) > 1){
    points(tab[signGenes_down, ], pch = 16, cex = res$baseMean/1000, col = alpha("blue", 0.4)) 
  }
  abline(h = -log10(pval), col = "green3", lty = 2) 
  abline(v = c(-lfc, lfc), col = "orange", lty = 2) 
  if (length(unique(signGenes_up)) > 1){
    text(tab[signGenes_up, ]$logFC, tab[signGenes_up, ]$negLogPval+gap, row.names(res[signGenes_up,]), cex = 0.5, col = "red")
  }
  if (length(unique(signGenes_down)) > 1){
    text(tab[signGenes_down, ]$logFC, tab[signGenes_down, ]$negLogPval+gap, row.names(res[signGenes_down,]), cex = 0.5, col = "blue")
  }
  mtext(paste("pval =", pval), side = 4, at = -log10(pval), cex = 0.8, line = 0.5, las = 1) 
  mtext(c(paste("-", 2, "fold"), paste("+", 2, "fold")), side = 3, at = c(-lfc, lfc), cex = 0.8, line = 0.5)
  legend("top",legend = c("Upregulate","Downregulate"),pch = c(16, 16), col = c(alpha("red", 0.6),alpha("blue", 0.6)))
  
  plot(tab, pch = 21, xlab = expression(log[2]~fold~change), ylab = expression(-log[10]~pvalue), cex.lab = 1.5, col = alpha("black", 0))
  points(tab[nosigGene, ], pch = 21, cex = res$baseMean/1000, xlab = expression(log[2]~fold~change), ylab = expression(-log[10]~pvalue), col = "black", bg = alpha("black", 0.1))
  if (length(unique(signGenes_up)) > 1){ 
    points(tab[signGenes_up, ], pch = 21, cex = res$baseMean/1000, col = "black", bg = alpha("red", 0.3)) 
  }
  if (length(unique(signGenes_down)) > 1){
    points(tab[signGenes_down, ], pch = 21, cex = res$baseMean/1000, col = "black", bg = alpha("cornflowerblue", 0.3)) 
  }
  abline(h = -log10(pval), col = "green3", lty = 2) 
  abline(v = c(-lfc, lfc), col = "orange", lty = 2) 
  if (length(unique(signGenes_up)) > 1){
    text(tab[signGenes_up, ]$logFC, tab[signGenes_up, ]$negLogPval+gap, row.names(res[signGenes_up,]), cex = 0.5, col = "red")
  }
  if (length(unique(signGenes_down)) > 1){
    text(tab[signGenes_down, ]$logFC, tab[signGenes_down, ]$negLogPval+gap, row.names(res[signGenes_down,]), cex = 0.5, col = "blue")
  }
  mtext(paste("pval =", pval), side = 4, at = -log10(pval), cex = 0.8, line = 0.5, las = 1) 
  mtext(c(paste("-", 2, "fold"), paste("+", 2, "fold")), side = 3, at = c(-lfc, lfc), cex = 0.8, line = 0.5)
  legend("top",legend = c("Upregulate","Downregulate"),pch = c(16, 16), col = c("red", "cornflowerblue"))
  dev.off()
}


# 
# #heatmap
# upregulateMatrix=recrudesce
# sampleInfo=sampleInfo2
# colnum=2
# library('pheatmap')
# library("RColorBrewer")
# 
# mergeSample2=data.frame(sampleInfo$Subset)
# rownames(mergeSample2)=sampleInfo[,1]
# colnames(mergeSample2)=c("Subset")
# #col= brewer.pal(colnum, "Paired")
# col=c("red" ,"blue")
# Subset=col
# names(Subset) = levels(mergeSample2$Subset)
# colll=list(Subset=Subset)
# tupregulateMatrix=t(upregulateMatrix)
# mydata = scale(tupregulateMatrix, center = TRUE, scale = TRUE)
# mydata=t(mydata)
# mydata[which(mydata>2,arr.ind=T)]=2
# mydata[which(mydata<(-2),arr.ind=T)]=(-2)
# 
# pheatmap(mydata,annotation_col = mergeSample2,cluster_cols = F,cluster_rows =F ,
#          colorRampPalette(c("green", "black", "red"))(50),show_rownames=F,show_colnames=F,annotation_colors =colll)
# 
# 
# 
# 
#venn.diagram
recrudesce=read.table('fufa_not.txt',sep='\t',header=T,stringsAsFactors = F)
#number cluster 1 : 1735
#number cluster 2 : 832
normal_tumor=read.table('normal_tumor2.txt',sep='\t',header=T,stringsAsFactors = F)

# number cluster 1 : 3041
# number cluster 2 : 26

recrudesce_gene=row.names(recrudesce)
normal_tumor_gene=row.names(normal_tumor)
venn.diagram(list(X=normal_tumor_gene,Y=recrudesce_gene),fill=c("red","blue"),"venn2.tiff")

mergene=intersect(normal_tumor_gene,recrudesce_gene)
recrudesce2=data.frame(rownames(recrudesce),recrudesce)
colnames(recrudesce2)=c("gene",colnames(recrudesce2)[2:41])

normal_tumor2=data.frame(rownames(normal_tumor),normal_tumor)
colnames(normal_tumor2)=c("gene",colnames(normal_tumor2)[2:41])

#取交集
recrudesce_normal_tumor=merge(recrudesce2,normal_tumor2,by='gene')
recrudesce_normal_tumor=recrudesce_normal_tumor[,1:41]
colnames(recrudesce_normal_tumor)=c("gene",colnames(normal_tumor)[2:41])
write.table(recrudesce_normal_tumor,'recrudesce_normal_tumor.txt',sep='\t',quote = F,row.names = F)

# 






