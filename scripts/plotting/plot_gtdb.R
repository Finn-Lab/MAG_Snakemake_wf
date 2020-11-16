# This file is part of MAG Snakemake workflow.
#
# MAG Snakemake workflow is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# MAG Snakemake workflow is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with MAG Snakemake workflow.  If not, see <https://www.gnu.org/licenses/>.

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
  total = nrow(gtdb.bac)   #total number of genomes
  rank.present = table(gtdb.bac[!grepl("__$", gtdb.bac[,r]),r])   #summary of genus
  total.rank = sum(rank.present) #sum of genus counts 
  gtdb.ranks = data.frame(sort(rank.present)[(length(rank.present)-6):length(rank.present)])
  gtdb.ranks$Prop = gtdb.ranks$Freq/total*100
  other.name = paste(tolower(substr(r, 1, 1)), "__Other", sep="")
  novel.name = paste(tolower(substr(r, 1, 1)), "__Novel", sep="")

  other.freq = total.rank-sum(gtdb.ranks$Freq)
  novel.freq=total-total.rank  #novel genomes

  other.prop = other.freq/total*100
  novel.prop = novel.freq/total*100

  other.ranks = data.frame(Var1=other.name, Freq=other.freq, Prop=other.prop)
  novel.ranks = data.frame(Var1=novel.name, Freq=novel.freq, Prop=novel.prop)

  ranks.fi = rbind(other.ranks, gtdb.ranks)
  ranks.fi = rbind(novel.ranks, ranks.fi)

  ranks.fi$Rank = r
  ranks.fi$Colour=c("#999999","#CAB2D6","#A6CEE3","#FF7F00","#FB9A99","#E31A1C","#B2DF8A","#33A02C","#1F78B4")
  
  if(exists("gtdb.fi.bac")){
    gtdb.fi.bac = rbind(gtdb.fi.bac, ranks.fi)
  } else {
    gtdb.fi.bac = ranks.fi
  }
}

gtdb.fi.bac$Rank=factor(gtdb.fi.bac$Rank,levels=c("Phylum","Class","Order","Family","Genus"))
# plot stacked plot taxa counts
p<-print(ggplot(gtdb.fi.bac, aes(x=Rank, y=Prop, fill=Var1)) 
      + geom_bar(stat="identity", colour="darkgrey", alpha=0.5, size=0.2, width=0.7)
      + theme_bw()
      + ylab("Proportion of species (%)")
      + coord_flip()
      + scale_fill_manual(values=as.vector(gtdb.fi.bac$Colour), name="Taxa")
      #+ guides(fill=FALSE)
      + scale_x_discrete(limits=ranks[(length(ranks)-1):2])
      + theme(legend.position="right")
      + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      + theme(axis.title.x = element_text(size=14))
      + theme(axis.text.y = element_text(size=12))
      + theme(axis.title.y = element_blank())
      + theme(axis.text.x = element_text(size=12)))



ggsave("data/figures/gtdb_bacteria.png",p,width=15)
