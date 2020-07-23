library(ggplot2)


args <- commandArgs(trailingOnly = TRUE)

checkm_sr=read.csv(args[1],header=TRUE,stringsAsFactor=FALSE)
checkm_coas=read.csv(args[2],header=TRUE,stringsAsFactor=FALSE)

checkm_sr$Status="Single run"
checkm_coas$Status="Coassembly"

checkm=rbind.data.frame(checkm_sr,checkm_coas)
setwd("data/figures/")
ggplot(checkm, aes(x=Status, y=completeness, fill=Status)) + 
  geom_boxplot()+theme_classic()+ylab("Completeness") + xlab("Assembly approach")+scale_fill_manual(breaks=c("Coassembly","Single run"), values=c("#a8ddb5","#c994c7"))
ggsave("checkm_completeness.png",plot = last_plot())

ggplot(checkm, aes(x=Status, y=contamination, fill=Status))+geom_boxplot()+theme_classic()+ylab("Contamination") + xlab("Assembly approach")+scale_fill_manual(breaks=c("Coassembly","Single run"), values=c("#a8ddb5","#c994c7"))
ggsave("checkm_contam.png",width=5,height=5)
