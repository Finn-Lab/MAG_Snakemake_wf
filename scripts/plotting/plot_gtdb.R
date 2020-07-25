library(RColorBrewer)
library(tidyr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
gtdb.taxa=read.delim(args[1],header=TRUE,stringsAsFactor=FALSE)
ranks = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
gtdb.taxa = separate(data = gtdb.taxa, col = classification, sep = ";", into = ranks)
gtdb.bac = gtdb.taxa[which(gtdb.taxa$Domain == "d__Bacteria"),]

colors=c("#CC79A7","#0072B2","#56B4E9","#009E73")

# calculate bacteria proportions
if(exists("gtdb.fi.bac")){
  rm(gtdb.fi.bac)
}
for (r in ranks[2:(length(ranks)-1)]){
  total = nrow(gtdb.bac)   #total number of genomes
  rank.present = table(gtdb.bac[!grepl("__$", gtdb.bac[,r]),r]) 
  total.rank = sum(rank.present)
  numb_ranks=nrow(rank.present)
  numb_ranks_display=min(3,round(2/3*numb_ranks))
  gtdb.ranks = data.frame(sort(rank.present)[(length(rank.present)-numb_ranks_display):length(rank.present)])
  gtdb.ranks$Prop = gtdb.ranks$Freq/total*100
  gtdb.ranks$Rank = r
  gtdb.ranks$Colour = colors[1:(nrow(gtdb.ranks))]
  ranks.fi=gtdb.ranks

  other.name = paste(tolower(substr(r, 1, 1)), "__Other", sep="")
  novel.name = paste(tolower(substr(r, 1, 1)), "__Novel", sep="")

  other.freq = total.rank-sum(gtdb.ranks$Freq)
  novel.freq=total-total.rank  #novel genomes

  other.prop = other.freq/total*100
  novel.prop = novel.freq/total*100

  other.ranks = data.frame(Var1=other.name, Freq=other.freq, Prop=other.prop)
  novel.ranks = data.frame(Var1=novel.name, Freq=novel.freq, Prop=novel.prop)

  other.ranks$Colour="#DDDDDD"
  novel.ranks$Colour="#A9A9A9"

  other.ranks$Rank = r
  novel.ranks$Rank = r
  

  if (other.ranks$Freq!=0){ranks.fi = rbind(other.ranks, ranks.fi)}
  if (novel.ranks$Freq!=0){ranks.fi = rbind(novel.ranks, ranks.fi)}


  if(exists("gtdb.fi.bac")){
    gtdb.fi.bac = rbind(gtdb.fi.bac, ranks.fi)
  } else {
    gtdb.fi.bac = ranks.fi
  }
}

gtdb.fi.bac$Rank=factor(gtdb.fi.bac$Rank,levels=c("Phylum","Class","Order","Family","Genus"))


# plot stacked plot taxa counts
print(ggplot(gtdb.fi.bac, aes(x=Rank, y=Prop, fill=Var1)) 
      + geom_bar(stat="identity", colour="darkgrey", alpha=0.7, size=0.2, width=0.7)
      + theme_bw()
      + ylab("Proportion of species (%)")
      + coord_flip()
      + scale_fill_manual(values=as.vector(gtdb.fi.bac$Colour), name="Taxa")
      + scale_x_discrete(limits=ranks[2:(length(ranks)-1)])
      + theme(legend.position="right")
      + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      + theme(axis.title.x = element_text(size=14))
      + theme(axis.text.y = element_text(size=12))
      + theme(axis.title.y = element_blank())
      + theme(axis.text.x = element_text(size=12)))

ggsave("data/figures/gtdb_bacteria.png")
