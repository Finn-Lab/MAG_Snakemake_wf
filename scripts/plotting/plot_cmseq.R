library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

singlerun=read.table(args[1],header=FALSE,stringsAsFactor=FALSE)
coassembly=read.table(args[2],header=FALSE,fill=TRUE,stringsAsFactor=FALSE)
colnames(singlerun)=c("path","strain")
colnames(coassembly)=c("path","strain")

coassembly$status="Co-assembly"
singlerun$status="Single run"
coassembly$strain=as.numeric(as.vector(coassembly$strain))
singlerun$strain=as.numeric(as.vector(singlerun$strain))
coassembly=coassembly[complete.cases(coassembly$strain),]
singlerun=singlerun[complete.cases(singlerun$strain),]

all=rbind(coassembly,singlerun)
all$status=factor(all$status, levels=c("Single run", "Co-assembly"))
ggplot(all, aes(x=status, y=strain, fill=status)) +geom_boxplot(outlier.shape=NA)+theme_classic()+ylab("Strain Heterogeneity") + xlab("Approach")+scale_fill_manual(breaks=c("Single run","Co-assembly"), values=c("#BBBBBB","#4477AA"))+guides(color=guide_legend(title="Approach"))+ylim(0,12)+theme(legend.position="none")
ggsave("data/figures/cmseq_plot.png",width=5,height=5)
