setwd('D:\\CRC_lncRNA\\filter\\lncRNA')
novel_gtf=read.table('D:\\CRC_lncRNA\\filter\\lncRNA\\novel_length.txt',sep='\t',stringsAsFactors = F)
novel_length=novel_gtf[,2]
x=novel_length
#y=density(x,na.rm=T) 

known_gtf=read.table('D:\\CRC_lncRNA\\filter\\lncRNA\\known_length.txt',sep='\t',stringsAsFactors = F)
known_length=known_gtf[,2]
x2=known_length
#y2=density(x2,na.rm=T) 

coding_gtf=read.table('D:\\CRC_lncRNA\\filter\\lncRNA\\pc_length.txt',sep='\t',stringsAsFactors = F)
coding_length=coding_gtf[,2]
x3=coding_length
#y3=density(x3,na.rm=T) 

# plot(y, main = "density",lwd=2,col='red', ylim = c(0,0.0008))
# lines(y2,col='blue',lwd=2)
# lines(y3,col='black',lwd=2)

df1=data.frame(length=x2,type='known lncRNA')
df2=data.frame(length=x,type='novel lncRNA')
df3=data.frame(length=x3,type='protein-coding gene')

df=rbind(df1,df2,df3)

library(ggplot2)
pdf(file="D:\\CRC_lncRNA\\filter\\lncRNA\\length_desitiny_plot.pdf")

sp=ggplot(df,aes(x = length, fill= type,colour = type)) +geom_density()+xlim(limits = c(-1000, 90000))+labs(x="",y = " ")
sp+theme_bw() + theme(legend.title=element_blank(),legend.position=c(0.8,0.8),panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),axis.title.x = element_text(size = 15, face = "bold"),axis.title.y= element_text(size = 15, face = "bold"))
dev.off()

#删除legend.tittle\修改位置、去除背景和格子、修改坐标轴字体大小加粗



