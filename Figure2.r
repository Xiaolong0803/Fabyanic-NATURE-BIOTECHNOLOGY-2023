#!/usr/bin/Rscript
# or open it with jupyter notebook
# Title: "NATURE BIOTECHNOLOGY README"
# Author: "Peng"

# load previous seurat object
# load clustering results
library("Seurat")
library(dplyr)
library(tidyr)
library(cowplot)
library(ggplot2)

# Figure2d
pbmc.snhmcALL <- readRDS(file="20210922_5hmCG_100kb_bin_925nuclei.rds")

options(repr.plot.width=8,repr.plot.height=4)

DimPlot(pbmc.snhmcALL, reduction = "umap",group.by = "FACS", pt.size = .2,split.by = "tech",label = T)
#ggsave("Fig2_snhmCseq_group_by_FACS_channel.pdf",width=8,height=4)

DimPlot(pbmc.snhmcALL, reduction = "umap",group.by = "res.name", pt.size = .2,split.by = "tech",label = T)
#ggsave("Fig2_snhmCseq_group_by_cluster.pdf",width=8,height=4)


# Figure2e
library(cowplot)
library(RColorBrewer)

table(pbmc.snhmcALL@meta.data$res.name,pbmc.snhmcALL@meta.data$FACS,pbmc.snhmcALL@meta.data$tech)
108/(108+20)
156/(156+30)

datComp <- data.frame(type = rep(c("FACS","non-split","split"),each=2),
                      cellNumber = c(7371,1488,108,20,156,30),
                      cluster = rep(c("Ex","Inh"),3))
datComp

options(repr.plot.width =4,repr.plot.height=4)
p6 <- ggplot(datComp,aes(type,cellNumber,fill=cluster)) + 
 geom_bar(position="fill", stat="identity")+
 scale_fill_manual("",values=brewer.pal(n = 3 , name = "Set1")) +
 theme_cowplot() + 
 geom_text(aes(label=cellNumber),position = position_fill(vjust = 0.5),color="white") +
 xlab("") + ylab("Prop. of nuclei") 
p6