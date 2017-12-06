known_novel="novel_known"
setwd(paste("D:\\CRC_lncRNA\\cnv\\differentgene_updown_heatmap"))

cancername=read.table("D:\\CRC_lncRNA\\cnv\\percentCNV\\cancer13.txt",stringsAsFactors = F)
cancername=cancername[,1]

COAD_amp_df=read.table(paste('COADREAD.res_',known_novel,'.Amp.geneid_precent_sorted.gistic',sep=''),sep='\t',stringsAsFactors = F)
COAD_amp=COAD_amp_df[,c(3,4)]
colnames(COAD_amp)=c("gene","COADREAD_Amp")
COAD_del_df=read.table(paste('COADREAD.res_',known_novel,'.Del.geneid_precent_sorted2.gistic',sep=''),sep='\t',stringsAsFactors = F)
COAD_del=COAD_del_df[,c(3,4)]
colnames(COAD_del)=c("gene","COADREAD_Del")
all_cancer_amp_del=merge(COAD_amp,COAD_del,by='gene',sort = F)
dim(all_cancer_amp_del)

#all_cancer_amp_del_sorted=all_cancer_amp_del[order(intersect_normal_cnv_up_down),]

for (can in cancername){
  print (can)
  COAD_amp_df=read.table(paste(can,'.res_',known_novel,'.Amp.geneid_precent_sorted.gistic',sep=''),sep='\t',stringsAsFactors = F)
  COAD_amp=COAD_amp_df[,c(3,4)]
  colnames(COAD_amp)=c("gene",paste(can,"_Amp",sep=''))
  COAD_del_df=read.table(paste(can,'.res_',known_novel,'.Del.geneid_precent_sorted2.gistic',sep=''),sep='\t',stringsAsFactors = F)
  COAD_del=COAD_del_df[,c(3,4)]
  colnames(COAD_del)=c("gene",paste(can,"_Del",sep=''))
  COAD_amp_del=merge(COAD_amp,COAD_del,by='gene',sort = F)
  dim(COAD_amp_del)
  all_cancer_amp_del=merge(all_cancer_amp_del,COAD_amp_del,by='gene',sort = F)
}

all_cancer_amp_del2=all_cancer_amp_del

all_cancer_amp_del_up=all_cancer_amp_del2[1:1098,]
all_cancer_amp_del_down=all_cancer_amp_del2[1099:nrow(all_cancer_amp_del2),]
all_cancer_amp_del_up_sorted=all_cancer_amp_del_up[order(all_cancer_amp_del_up[,2],decreasing = T),]
all_cancer_amp_del_down_sorted=all_cancer_amp_del_down[order(all_cancer_amp_del_down[,3],decreasing = F),]
all_cancer_amp_del_up_down_sorted=rbind(all_cancer_amp_del_up_sorted,all_cancer_amp_del_down_sorted)
#all_cancer_amp_del3=merge(intersect_normal_cnv_up_down,all_cancer_amp_del2,by.x='lncRNA',by.y='gene',sort=F)
dim(all_cancer_amp_del_up_down_sorted)

pdf(paste("D:\\CRC_lncRNA\\cnv\\percentCNV\\",known_novel,"_heatmapinallcancer13.pdf",sep=''),width = 2000, height = 1500)
png(paste("D:\\CRC_lncRNA\\cnv\\percentCNV\\",known_novel,"_heatmapinallcancer13.png",sep=''),width = 2000, height = 1500)
all_cancer_data=all_cancer_amp_del_up_down_sorted[,-1]
pheatmap(all_cancer_data,gaps_row=1099,cluster_cols = F,cluster_rows =F,show_rownames = F,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

dev.off()



# ACC_amp_df=read.table('ACC.res_novel.Amp.geneid_precent_sorted.gistic',sep='\t',stringsAsFactors = F)
# ACC_amp=ACC_amp_df[,c(3,4)]
# colnames(ACC_amp)=c("gene","per")
# 
# ACC_del_df=read.table('ACC.res_novel.Del.geneid_precent_sorted2.gistic',sep='\t',stringsAsFactors = F)
# ACC_del=ACC_del_df[,c(3,4)]
# colnames(ACC_del)=c("gene","per")
# ACC_amp_del=merge(ACC_amp,ACC_del,by='gene')
# dim(ACC_amp_del)
# 
# all_cancer_amp_del=merge(COAD_amp_del,ACC_amp_del,by='gene')
# 
# all_cancer_data=all_cancer_amp_del[,-1]
# pheatmap(all_cancer_data,cluster_cols = F,cluster_rows =F,show_rownames = F)
