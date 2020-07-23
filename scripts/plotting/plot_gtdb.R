library(RColorBrewer)
library(tidyr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

gtdb.taxa=read.delim(args[1],header=TRUE,stringsAsFactor=FALSE)

ranks = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
gtdb.taxa = separate(data = gtdb.taxa, col = classification, sep = ";", into = ranks)
gtdb.bac = gtdb.taxa[which(gtdb.taxa$Domain == "d__Bacteria"),]


# calculate bacteria proportions
if(exists("gtdb.fi.bac")){
  rm(gtdb.fi.bac)
}
for (r in ranks[2:(length(ranks)-1)]){
  total = nrow(gtdb.bac)
  rank.present = table(gtdb.bac[!grepl("__$", gtdb.bac[,r]),r])
  total.rank = sum(rank.present)
  gtdb.ranks = data.frame(sort(rank.present)[(length(rank.present)-6):length(rank.present)])
  gtdb.ranks$Prop = gtdb.ranks$Freq/total*100
  other.name = paste(tolower(substr(r, 1, 1)), "__Other", sep="")
  other.freq = total.rank-sum(gtdb.ranks$Freq)
  other.prop = other.freq/total*100
  other.ranks = data.frame(Var1=other.name, Freq=other.freq, Prop=other.prop)
  ranks.fi = rbind(other.ranks, gtdb.ranks)
  ranks.fi$Rank = r
  ranks.fi$Colour = c("grey", rev(brewer.pal(7,"Set3")))
  if(exists("gtdb.fi.bac")){
    gtdb.fi.bac = rbind(gtdb.fi.bac, ranks.fi)
  } else {
    gtdb.fi.bac = ranks.fi
  }
}
# plot stacked plot taxa counts
print(ggplot(gtdb.fi.bac, aes(x=Rank, y=Prop, fill=Var1)) 
      + geom_bar(stat="identity", colour="darkgrey", alpha=0.7, size=0.2, width=0.7)
      + theme_bw()
      + ylab("Proportion of species (%)")
      + coord_flip()
      + scale_fill_manual(values=as.vector(gtdb.fi.bac$Colour), name="Taxa")
      #+ guides(fill=FALSE)
      + scale_x_discrete(limits=ranks[2:(length(ranks)-1)])
      + theme(legend.position="right")
      + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      + theme(axis.title.x = element_text(size=14))
      + theme(axis.text.y = element_text(size=12))
      + theme(axis.title.y = element_blank())
      + theme(axis.text.x = element_text(size=12)))

ggsave("figures/gtdb_bacteria.png")
