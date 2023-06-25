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

#==========
# Figure 4a
#==========
pbmc.snmc <- readRDS(file="snmCseq2_split_CpH/20210917_pbmc.snmc.rds")
pbmc.snmc

DefaultAssay(pbmc.snmc) <- "gbch"

pbmc.snmc <- ScaleData(pbmc.snmc)

options(repr.plot.width = 5, repr.plot.height = 4)
DimPlot(pbmc.snmc,pt.size = .5,group.by = "res.name.new",label = TRUE)

Idents(pbmc.snmc) <- "res.name.new"

pbmc.snmc2 <- subset(pbmc.snmc,idents = "In",invert=T)

options(repr.plot.width = 5, repr.plot.height = 4)
DimPlot(pbmc.snmc2,pt.size = .5,group.by = "res.name.new",label = TRUE)


#==========
# Figure 4b & 4c
#==========

plot_box_ser2 <- function(features=c("Pvalb","Sst","Ndnf","Sox6","Reln","Cacna2d2","Lhx6","Gria1","Vip"),
                         ser.obj=pbmc.snmc,
                         assay="gbch",
                         clusterID="res.name.new",
                         nRow = 1
                         ) {
    
    metadata_df <- ser.obj@meta.data[,c("cell",clusterID)]
    mat_dat <- ser.obj[[assay]]@data
    #dim(mat_dat)
    genes.in <- intersect(features,rownames(mat_dat))
    #print(genes.in)
    mat_dat <- mat_dat[genes.in,row.names(metadata_df)]
    mat_dat2 <- as.matrix(mat_dat)
    data.use <- as.data.frame(t(mat_dat2))
    
    df.plot <- cbind(metadata_df,data.use)
    df.plot2 <- reshape2::melt(df.plot,by=c("cell",clusterID))
    df.plot2$cluster <- df.plot2[,clusterID]
    
    df.plot2 %>% select(variable,value) %>% 
             #group_by(cluster,variable) %>% 
             group_by(variable) %>% 
             arrange(desc(value)) %>% dplyr::slice(1) %>% dplyr::rename(max = value) -> df.max.value
    
    df.plot3 <- df.plot2 %>% 
                #group_by(cluster,variable) %>% 
                group_by(variable) %>% 
                left_join(df.max.value,by=c("variable")) %>% ungroup %>%
                dplyr::mutate(value2 = value/max)
    
    #cluster.order.our <- c("L2/3","L4","L4/5","L6","DL","In","Vip/Ndnf","Pv","Sst")
    cluster.order.our <- c("L2/3","L4","L4/5","L6","DL","Vip/Ndnf","Pv","Sst")
    df.plot3$cluster <- factor(df.plot3$cluster,levels = cluster.order.our)
    
    p <- ggplot(df.plot3,aes(cluster,value2)) + 
         geom_boxplot(width=.6,outlier.shape = NA) + 
         facet_wrap(.~variable,scales = "free",nrow = nRow) +
         stat_summary(fun="mean", geom="point", shape=23, size=2,col="red") +
         theme_cowplot()+
         theme(axis.text.x = element_text(angle=45,hjust=1)) +
         ylab("Scale max to 1") + xlab("")
    
    return(p)
}

options(repr.plot.width = 12 ,repr.plot.height = 6)

gene.in.sub <- c("Arpp21","Cux2","Rorb","Tle4","Slc6a1","Erbb4","Adarb2","Cacna2d2")


plot_box_ser2(features=gene.in.sub,
             ser.obj=pbmc.snmc2,
             assay="gbch",
             clusterID="res.name.new",
             nRow = 2)
			 
options(repr.plot.width=12,repr.plot.height=5)
FeaturePlot(pbmc.snmc2,features = gene.in.sub,
            slot = "scale.data",
            ncol = 4,
            pt.size = .1,
            cols = c("red", "grey")
            #min.cutoff = "q10", 
            #max.cutoff = "q90"
           )& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
		   
#==========
# Figure 4d
#==========

pbmc.JE <- readRDS(file="20210917_CpH_JE_clustering.rds")
table(pbmc.JE$Neuron.type2)
Idents(pbmc.JE) <- "Neuron.type2"

pbmc.JE.removeIn <- subset(pbmc.JE,idents = c("In"),invert=T)
pbmc.JE.removeIn@meta.data <- droplevels(pbmc.JE.removeIn@meta.data)
table(pbmc.JE.removeIn$Neuron.type2)

options(repr.plot.width = 6,repr.plot.height = 5)
#DimPlot(pbmc.JE.removeIn,group.by = "Neuron.type2")
DimPlot(pbmc.JE.removeIn,pt.size = .5,group.by = "Neuron.type2",label = TRUE)

# define the reference
ad.ref <- pbmc.JE.removeIn
# define the query
ad.query <- pbmc.snmc2
# predictions

ad.anchors <- FindTransferAnchors(reference = ad.ref,
                                  query = ad.query,
                                  features = unique(c(VariableFeatures(object = ad.ref), VariableFeatures(object = ad.query))),
                                  reduction = "cca",
                                  dims = 1:30)

predictions <- TransferData(anchorset = ad.anchors, 
                            refdata = ad.ref$Neuron.type2,
                            weight.reduction = ad.query[["pca"]],
                            dims = 1:30)

ad.query <- AddMetaData(ad.query, metadata = predictions)

datScore <- ad.query@meta.data
dim(datScore)
all.equal(rownames(datScore),rownames(pbmc.snmc@meta.data))

colnames(datScore)









