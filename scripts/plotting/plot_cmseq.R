library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

singlerun=read.table(args[1],header=FALSE,stringsAsFactor=FALSE)
coassembly=read.table(args[2],header=FALSE,fill=TRUE,stringsAsFactor=FALSE)
colnames(singlerun)=c("path","strain")
colnames(coassembly)=c("path","strain")

coassembly$status="Coassembly"
singlerun$status="Single run"
coassembly$strain=as.numeric(as.vector(coassembly$strain))
singlerun$strain=as.numeric(as.vector(singlerun$strain))
coassembly=coassembly[complete.cases(coassembly$strain),]
singlerun=singlerun[complete.cases(singlerun$strain),]

all=rbind(coassembly,singlerun)
all$`Assembly approach`=all$status
ggplot(all, aes(x=status, y=strain, fill=`Assembly approach`)) + 
  geom_boxplot(outlier.shape = NA)+theme_classic()+ylab("Strain heterogeneity") + xlab("Assembly approach")+scale_fill_manual(breaks=c("Coassembly","Single run"), values=c("#a8ddb5","#c994c7"))+ylim(0,12)
ggsave("data/figures/cmseq_plot.png",width=5,height=5)
