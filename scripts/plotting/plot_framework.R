library(ggplot2)
library(gridExtra)

args <- commandArgs(trailingOnly = TRUE)

readcounts=read.table(args[1],header=TRUE,stringsAsFactor=FALSE)
df_flag=read.table(args[2],header=FALSE,stringsAsFactor=FALSE)
bwa.counts_sr=read.table(args[3],header=FALSE,stringsAsFactor=FALSE)
bwa.counts_coas=read.table(args[4],header=FALSE,stringsAsFactor=FALSE)

readcounts$Run=gsub("data/00_preprocessing/processed/singlerun/","",readcounts$Run)
readcounts=readcounts[!grepl("_2.fastq", readcounts$Run),]
readcounts$Run=gsub("_1.fastq","",readcounts$Run)
readcounts$Readcount=2*readcounts$Readcount

colnames(df_flag)=c("Run","Assembly")
df_comb=merge(df_flag,readcounts,by="Run")
df_comb$percassemb=df_comb$Assembly/df_comb$Readcount*100

colnames(bwa.counts_sr)=c("Run","Catalogue")
bwa.counts_sr=merge(bwa.counts_sr,readcounts,by="Run")
bwa.counts_sr$`Assembly Approach`="Single run"
bwa.counts_sr$percmags=bwa.counts_sr$Catalogue/bwa.counts_sr$Readcount*100

colnames(bwa.counts_coas)=c("Run","Catalogue")
bwa.counts_coas=merge(bwa.counts_coas,readcounts,by="Run")
bwa.counts_coas$`Assembly Approach`="Coassembly"
bwa.counts_coas$percmags=bwa.counts_coas$Catalogue/bwa.counts_coas$Readcount*100

bwa.counts=rbind.data.frame(bwa.counts_sr,bwa.counts_coas)
merged=merge(df_comb,bwa.counts, by="Run")

ggplot(merged, aes(x=percmags, y=percassemb, color=`Assembly Approach`)) + geom_point() +xlab("Reads mapping to MAGs (%)") + ylab("Reads mapping to assembly (%)")+theme_classic()+scale_color_manual(breaks=c("Coassembly","Single run"), values=c("#a8ddb5","#c994c7"))+xlim(0,100)+ylim(0,100)
ggsave("data/figures/perassemb_perref.png",width=5,height=5)
