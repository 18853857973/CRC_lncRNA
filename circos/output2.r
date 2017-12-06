setwd('D:\\CRC_lncRNA\\CPAT')
load("Human_logitModel.RData")
test <- read.table(file="output2.dat",sep="\t",col.names=c("ID","mRNA","ORF","Fickett","Hexamer"))
test$prob <- predict(mylogit,newdata=test,type="response")
attach(test)
output <- cbind("mRNA_size"=mRNA,"ORF_size"=ORF,"Fickett_score"=Fickett,"Hexamer_score"=Hexamer,"coding_prob"=test$prob)
write.table(output,file="output2",quote=F,sep="\t",row.names=ID)


library(rgl)
plot3d(test[,3],test[,4],test[,5],col="red", size=5)
