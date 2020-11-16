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
