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


library(ggplot2)
library(gridExtra)

args <- commandArgs(trailingOnly = TRUE)

checkm_sr=read.csv(args[1],header=TRUE,stringsAsFactor=FALSE)
checkm_coas=read.csv(args[2],header=TRUE,stringsAsFactor=FALSE)
dnadiff=read.table(args[3],header=TRUE,stringsAsFactor=FALSE)


colnames(dnadiff)=c("Ref","genome","Ref_length", "Refcovered","Query_length","Queryaligned", "ANI")
dnadiff$genome= gsub(".*/","",dnadiff$genome)
dnadiff$genome=gsub(".fa","",dnadiff$genome)

checkm_sr$`Assembly approach`="Single run"
checkm_coas$`Assembly approach`="Co-assembly"
checkm=rbind.data.frame(checkm_sr,checkm_coas)

checkm$Quality="Medium quality"
checkm$Quality[checkm$completeness>=90&checkm$contamination<=5]="High quality"
checkm$genome=gsub(".fa","",checkm$genome)

comb=merge(dnadiff,checkm,by="genome")
left<-ggplot(comb, aes(x = Queryaligned, y = Refcovered, color=`Assembly approach`, shape=Quality)) +
         geom_point(stroke = 1,size=1)+xlab("MAG Aligned (%)")+ylab("Reference Aligned (%)") + theme_classic() +scale_color_manual(breaks=c("Single run","Co-assembly"), values=c("#BBBBBB","#4477AA"))+xlim(0,100)+ ylim(0,100)+scale_shape_manual(values = c(1,2))+guides(color=guide_legend(title="Approach"))
right<-
ggplot(comb, aes(x=ANI,color=`Assembly approach`)) + geom_density()+theme_classic()+xlab("ANI") +xlim(0,100)+ylab("Density")+scale_color_manual(breaks=c("Single run","Co-assembly"), values=c("#BBBBBB","#4477AA"))+guides(color=guide_legend(title="Approach"))

p<-grid.arrange(left, right, nrow = 1)
ggsave("data/figures/dnadiff.png",p, width=10,height=5)



