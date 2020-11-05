library(ggplot2)


args <- commandArgs(trailingOnly = TRUE)

checkm_sr=read.csv(args[1],header=TRUE,stringsAsFactor=FALSE)
checkm_coas=read.csv(args[2],header=TRUE,stringsAsFactor=FALSE)

checkm_sr$Status2="Single run"
checkm_coas$Status2="Co-assembly"

checkm=rbind.data.frame(checkm_sr,checkm_coas)
checkm$Status=factor(checkm$Status, levels=c("Single run", "Co-assembly"))

setwd("data/figures/")
ggplot(checkm, aes(x=Status, y=completeness, fill=Status)) +
  geom_boxplot(outlier.shape=NA)+theme_classic()+ylab("Completeness") + xlab("Assembly approach")+scale_fill_manual(breaks=c("Single run","Co-assembly"), values=c("#BBBBBB","#4477AA"))+guides(color=guide_legend(title="Approach"))+theme(legend.position="none")
ggsave("checkm_completeness.png",plot = last_plot())

ggplot(checkm, aes(x=Status, y=contamination, fill=Status))+geom_boxplot(outlier.shape=NA)+theme_classic()+ylab("Contamination") + xlab("Assembly approach")+scale_fill_manual(breaks=c("Single run","Co-assembly"), values=c("#BBBBBB","#4477AA"))+theme(legend.position="none")+ylim(0,10)
ggsave("checkm_contam.png",width=5,height=5)
