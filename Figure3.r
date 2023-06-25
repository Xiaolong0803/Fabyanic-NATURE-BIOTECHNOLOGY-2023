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
# Figure3a
#==========
CG_comb_rate <- readRDS("snmCseq2_split_CpG/20210902_CG_comb_rate.rds")
CG_comb_nsite <- readRDS("20210902_CG_comb_nsite.rds")

## assign the NA when < 20 CpG sites
CG_comb_rate[CG_comb_nsite < 20] <- NA

# calculate the bin with NA
na_num_per_bin <- apply(CG_comb_rate,1,function(x){
    sum(is.na(x))
})

dim(CG_comb_rate) # 27348 bins

keep.features <- names(na_num_per_bin)[na_num_per_bin <= ncol(CG_comb_rate) * (1-0.9)] # less than 10%
length(keep.features)

CG_comb_rate2 <- CG_comb_rate[keep.features,]

dim(CG_comb_rate2) # 8649 552

# impute with mean
m <- CG_comb_rate2
## from https://stackoverflow.com/questions/6918086/replace-na-values-by-row-means
k <- which(is.na(m), arr.ind=TRUE)
m[k] <- rowMeans(m, na.rm=TRUE)[k[,1]]
CG_comb_rate3 <- m

# load 100kb bin
pbmc.snmc <- CreateSeuratObject(counts = CG_comb_rate3)

pbmc.snmc@meta.data$cell <- rownames(pbmc.snmc@meta.data)

require(dplyr)
require(tidyr)

# add the cellinformation from table provided by Hao
# load the cellinfo
cellinfo <- read.table("/mnt/data1/peng/workbase2/scACEseq/update_from_20210909/All_snmCseq2_removePercent.txt",h=T,sep="\t")

cellinfo$cell <- paste("CpG_mm10",
                       cellinfo$FACS,
                      cellinfo$Exp_date,
                      cellinfo$Seq_run_id,
                      cellinfo$Assay_type,
                      cellinfo$pool_id,
                      cellinfo$barcode_id,
                      sep = "_")

pbmc.snmc@meta.data$cell <- rownames(pbmc.snmc@meta.data)

cellinfo$cell <- gsub("2021","21",cellinfo$cell)
pbmc.snmc@meta.data$cell <- gsub("2021","21",pbmc.snmc@meta.data$cell)

length(intersect(cellinfo$cell,pbmc.snmc@meta.data$cell))

rownames(cellinfo) <- cellinfo$cell

df <- pbmc.snmc@meta.data
df <- df %>% separate(cell, c("type", "genome", "input", "date","sequencer","method","barcode1","barcode2"),sep = "[_]",remove = FALSE)

#head(df,n=2)
#head(pbmc.snmc@meta.data,n=2)

df2 <- df %>% group_by(cell) %>% left_join(cellinfo,by="cell") %>% ungroup %>% as.data.frame

rownames(df2) <- rownames(pbmc.snmc@meta.data)

pbmc.snmc@meta.data <- as.data.frame(df2)

pbmc.snmc <- pbmc.snmc %>% 
             NormalizeData(verbose = F) %>%
             FindVariableFeatures(verbose = F) %>%
             ScaleData(verbose = F) %>%
             RunPCA(npcs = 30,verbose = F) %>%
             FindNeighbors(dims = 1:30,verbose = F) %>%
             FindClusters(resolution = c(.2,0.4,.6,.8,1,2),verbose = F) %>%
             RunUMAP(dims = 1:30,verbose = F) %>%
             identity()

library(cowplot)

options(repr.plot.width=16,repr.plot.height=3)

p1 <- DimPlot(pbmc.snmc, reduction = "umap",group.by = "input", pt.size = .3)
p2 <- DimPlot(pbmc.snmc, reduction = "umap",group.by = "RNA_snn_res.0.6", pt.size = .3,label = T)
p3 <- DimPlot(pbmc.snmc, reduction = "umap",group.by = "RNA_snn_res.1", pt.size = .3,label = T)
p4 <- DimPlot(pbmc.snmc, reduction = "umap",group.by = "RNA_snn_res.2", pt.size = .3,label = T)

plot_grid(p1,p2,p3,p4,nrow=1,align = "vh")

pbmc.snmc$input2 <- plyr::mapvalues(pbmc.snmc$input,
                                          from = c("Ex","Inh","NeuNpos","NeuNneg"),
                                          to = c('NeuNpos','NeuNpos','NeuNpos','NeuNneg'))

table(pbmc.snmc$input2)

colnames(pbmc.snmc@meta.data)

pbmc.snmc@meta.data$res.name <- plyr::mapvalues(pbmc.snmc@meta.data$RNA_snn_res.2,
                                from = c(0,2,4,5,1,3,7,9,6,8,10),
                                to = c("Ex","Ex","Ex","Ex","Inh","Inh","Inh","Striatum","Oligo","Astro","MG"))

table(pbmc.snmc@meta.data$input)

pbmc.snmc$input3 <- plyr::mapvalues(pbmc.snmc$input,
                                          from = c("Ex","Inh","NeuNpos","NeuNneg"),
                                          to = c('NeuN+/Neurod6+','NeuN+/Neurod6-','NeuN+','NeuN-'))

options(repr.plot.width=12,repr.plot.height=4)

p1 <- DimPlot(pbmc.snmc, reduction = "umap",group.by = "input3", pt.size = .3)
p2 <- DimPlot(pbmc.snmc, reduction = "umap",group.by = "res.name", pt.size = .3,label = T)
#p3 <- DimPlot(pbmc.snmc, reduction = "umap",group.by = "RNA_snn_res.1", pt.size = .3,label = T)
#p4 <- DimPlot(pbmc.snmc, reduction = "umap",group.by = "RNA_snn_res.2", pt.size = .3,label = T)

plot_grid(p1,p2,nrow=1,align = "vh")

#==========
# Figure3b
#==========
# load function
plot_box4_value_nonNeu <- function(genes.in=c("Tle4","Gad2","Gad1","Rorb"),mat1=CG_genebody,yLab="GBmCG",nRow=4){
    mat2 <- mat1[intersect(genes.in,rownames(mat1)),]
    mat.df <- as.data.frame(cbind(t(mat2),as.character(pbmc.snmc@meta.data$res.name)))
    names(mat.df)[ncol(mat.df)] <- "cluster"
    df.plot <- reshape2::melt(mat.df,id=c("cluster"))
    df.plot$value <- as.numeric(df.plot$value)
    df.plot %>% ungroup %>% select(variable,value) %>% group_by(variable) %>% 
      arrange(desc(value)) %>% dplyr::slice(1) %>% dplyr::rename(max = value) -> df.max.value
    cluster.order <- c('Ex','Inh','Striatum','Astro','Oligo','MG')
    df.plot2 <- df.plot %>% ungroup %>% group_by(variable) %>% 
                  left_join(df.max.value,by=c("variable")) %>% ungroup %>%
                  dplyr::mutate(value2 = value/max) %>% dplyr::filter(cluster %in% cluster.order) %>% droplevels
    df.plot2$cluster <- factor(df.plot2$cluster,levels=cluster.order)
    library(ggplot2)
    library(cowplot)
    p1 <- ggplot(df.plot2,aes(cluster,value)) + geom_boxplot(width=.6,outlier.shape = NA) + 
          ggrastr::geom_point_rast(position = "jitter",col="grey",width = 0.25,size=3) + 
          facet_wrap(~variable,nrow = nRow) +
          stat_summary(fun.y="mean", geom="point", shape=23, size=2,col="red") +
          theme_cowplot() + theme(axis.text.x = element_text(angle=45,hjust=1)) + ylab("Meth rate")#+ ggtitle(yLab)
    return(p1)
}


# load genebody methlyation
CG_comb_rate <- readRDS(file="/mnt/data1/peng/workbase2/scACEseq/update_from_20210902/20210915_552nuclei_genebody_CG_comb_rate.rds")
CG_comb_nsite <- readRDS(file="/mnt/data1/peng/workbase2/scACEseq/update_from_20210902/20210915_552nuclei_genebody_CG_comb_nsite.rds")

filter_by_nSite_percent <- function(comb_nSite,comb_rate,nSite = 5, percent = 0.9){
    comb_rate[comb_nSite < nSite] <- NA
    na_num_per_bin <- apply(comb_rate,1,function(x){
        sum(is.na(x))
    })
    print(table(na_num_per_bin <= ncol(comb_rate) * (1-percent)))
}

filter_impute <- function(comb_nSite,comb_rate,nSite = 5, percent = 0.9){
    comb_rate[comb_nSite < nSite] <- NA
    na_num_per_bin <- apply(comb_rate,1,function(x){
        sum(is.na(x))
    })
    keep.features <- names(na_num_per_bin)[na_num_per_bin <= ncol(comb_rate) * (1-percent)] # less than 45%
    comb_rate2 <- comb_rate[keep.features,]
    # impute with mean
    m <- comb_rate2
    ## from https://stackoverflow.com/questions/6918086/replace-na-values-by-row-means
    k <- which(is.na(m), arr.ind=TRUE)
    m[k] <- rowMeans(m, na.rm=TRUE)[k[,1]]
    comb_rate3 <- m
    return(comb_rate3)
}

# gene body CpG, 1 nSites, 80%
CG_genebody <- filter_impute(comb_nSite=CG_comb_nsite,
                             comb_rate=CG_comb_rate,
                             nSite=1,
                             percent=0.8)
							 
							 
genes.in.final <- c("Slc17a7","Lingo1","Sv2b","Neurod6","Kcnma1", # Ex marker
                    "Gad1","Gad2","Abat","Slc6a1","Stat5b","Erbb4", # Inh marker
                    "Ppp1r1b", "Meis2", # Str marker
                    "Slc1a2","Trim9","Paqr8","Cpe","Slc1a3", # Astro marker
                    "Mobp","Mbp","Mog","Pde4b","Klh12", # Oligo marker
                    "Csf1r","P2ry13","Fcrls","Zfp710","Laptm5" # MG marker
                   )


genes.in.final <- c("Slc17a7","Gad2","Ppp1r1b","Paqr8","Mobp","Csf1r",
                   "Lingo1","Slc6a1", "Meis2","Slc1a3","Mog","Fcrls")
				   
				   
options(repr.plot.width = 12,repr.plot.height = 6)
plot_box4_value_nonNeu(genes.in.final,CG_genebody,"GBmCG",nRow=2)
table(pbmc.snmc@meta.data$res.name)

#==========
# Figure3d
#==========
dim(hmC_matrix)
Idents(snRNA3) <- "res.name.6groups"
RNA.6g.mak <- FindAllMarkers(snRNA3,only.pos = T,verbose = F)
cluster.order.6g <- c("Ex","Inh","Striatum","Astro","Oligo","MG")
RNA.6g.mak$cluster <- factor(RNA.6g.mak$cluster,levels = cluster.order.6g)
RNA.6g.mak.pval01 <- subset(RNA.6g.mak, p_val_adj < 0.1)
RNA.6g.mak2 <- RNA.6g.mak.pval01 %>% arrange(cluster,p_val_adj) 
gene.6g.order <- RNA.6g.mak2$gene
length(gene.6g.order)

gene.6g.order2 <- intersect(gene.6g.order,rownames(pbmc.snmc3@assays$tmCG@data))
length(gene.6g.order2)

gene.order.by.rna.6g <- gene.6g.order2

gene.order <- gene.order.by.rna.6g
cluster.order <- cluster.order.6g

mat_hmCG_celltype_6g_raw <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name.6g == x]
    return(rowMeans(hmC_matrix[gene.order,cells.in]))
})

mat_tmC_celltype_6g_raw <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name.6g == x]
    return(rowMeans(tmC_matrix[gene.order,cells.in]))
})

mat_mCG_celltype_6g_raw <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name.6g == x]
    return(rowMeans(mC_matrix[gene.order,cells.in]))
})

mat_mCH_celltype_6g_raw <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name.6g == x]
    return(rowMeans(mCH_matrix[gene.order,cells.in]))
})

mat_ratio_celltype_6g_raw <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name.6g == x]
    return(rowMeans(ratio_matrix[gene.order,cells.in]))
})

mat_rna_celltype_6g_raw <- sapply(cluster.order,function(x){
    cells.in <- rownames(snRNA3@meta.data)[snRNA3@meta.data$res.name.6g == x]
    mat <- expm1(snRNA3@assays$RNA@data)
    return(Matrix::rowMeans(mat[gene.order,cells.in]))
})

dim(mat_mCG_celltype_6g_raw)
mat_CG_celltype_6g_raw <- 1 - mat_mCG_celltype_6g_raw
dim(mat_CG_celltype_6g_raw)

options(repr.plot.width=3,repr.plot.height=3)

h0 <- pheatmap::pheatmap(mat_CG_celltype_6g_raw,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="CG",color=rev(colorRampPalette(brewer.pal(n = 9,name="BuPu"))(50)))

h1 <- pheatmap::pheatmap(mat_tmC_celltype_6g_raw,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="tmCG",color=rev(colorRampPalette(brewer.pal(n = 9,name="BuPu"))(50)))

h2 <- pheatmap::pheatmap(mat_mCG_celltype_6g_raw,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="BSmCG",color=rev(colorRampPalette(brewer.pal(n = 9,name="BuPu"))(50)))

h3 <- pheatmap::pheatmap(mat_hmCG_celltype_6g_raw,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="hmCG",color=rev(colorRampPalette(brewer.pal(n = 9,name="BuPu"))(50)))

h4 <- pheatmap::pheatmap(mat_ratio_celltype_6g_raw,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="hmCG/BSmCG",color=colorRampPalette(brewer.pal(n = 9,name="PuRd"))(50))

h5 <- pheatmap::pheatmap(mat_mCH_celltype_6g_raw,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="BSmCH",color=rev(colorRampPalette(brewer.pal(n = 9,name="BuPu"))(50)))

h6 <- pheatmap::pheatmap(mat_rna_celltype_6g_raw,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="RNA",colours = rev(colorRampPalette(brewer.pal(n = 11, name = "RdBu"))(50)))

# combine all the plots
plot_list=list()
plot_list[["CG"]] <- h0[[4]]
plot_list[["BSmCG"]] <- h2[[4]]
plot_list[["tmCG"]] <- h1[[4]]
plot_list[["hmCG"]] <- h3[[4]]
#plot_list[["hmCG/BSmCG"]] <- h4[[4]]
plot_list[["BSmCH"]] <- h5[[4]]
plot_list[["RNA"]] <- h6[[4]]

library(gridExtra)

options(repr.plot.width = 12,repr.plot.height = 3.5)
#grid.arrange(grobs=plot_list,ncol = 4)
g2 <- grid.arrange(grobs=plot_list,ncol = 6)


#==========
# Figure3e
#==========
pbmc.snmc3
plot_vln2 <- function(genes.in=c("Tle4","Gad2","Gad1","Rorb"),mat1=as.matrix(snRNA4@assays$RNA@data),yLab="RNA"){
    mat2 <- mat1[intersect(genes.in,rownames(mat1)),]
    mat.df <- as.data.frame(cbind(t(mat2),as.character(snRNA3@meta.data$res.name.6groups)))
    names(mat.df)[ncol(mat.df)] <- "cluster"
    df.plot <- reshape2::melt(mat.df,id=c("cluster"))
    df.plot$value <- as.numeric(df.plot$value)
    df.plot %>% ungroup %>% select(variable,value) %>% group_by(variable) %>% 
      arrange(desc(value)) %>% dplyr::slice(1) %>% dplyr::rename(max = value) -> df.max.value
    df.plot2 <- df.plot %>% ungroup %>% group_by(variable) %>% 
                  left_join(df.max.value,by=c("variable")) %>% ungroup %>%
                  dplyr::mutate(value2 = value/max)
    cluster.order <- c('Ex','Inh','Striatum','Astro','Oligo','MG')
    df.plot2$cluster <- factor(df.plot2$cluster,levels=cluster.order)
    library(ggplot2)
    library(cowplot)
    p1 <- ggplot(df.plot2,aes(cluster,value2)) + 
#          geom_boxplot(width=.6,outlier.shape = NA) + 
          geom_violin(scale="width",trim=T,alpha=0.8,adjust=2) + 
          ggrastr::geom_point_rast(position = "jitter",col="grey",width = 0.25,size=0.002) +
          #geom_jitter(size = .1,col="grey") + 
          facet_grid(~variable,scales = "free") +
#          facet_grid(~variable,nrow = nRow) +
          stat_summary(fun.y="mean", geom="point", shape=23, size=2,col="red") +
          theme_cowplot() + theme(axis.text.x = element_text(angle=45,hjust=1)) + ylab("Scale Max to 1") + ggtitle(yLab)
    return(p1)
}

genes.in2 <- c("Gad2","Slc1a3","Mobp","Csf1r")

pp1 <- plot_box4_value_nonNeu(genes.in=genes.in2,mat1=hmC_matrix,yLab="gene body hmCG", nRow=1)
pp2 <- plot_box4_value_nonNeu(genes.in=genes.in2,mat1=tmC_matrix,yLab="gene body True mCG", nRow=1)
pp3 <- plot_vln2(genes.in=genes.in2,mat1=as.matrix(expm1(snRNA3@assays$RNA@data)),yLab="RNA")

options(repr.plot.width = 12,repr.plot.height = 12)
plot_grid(pp1,pp2,pp3,nrow=3)


#==========
# Figure3f
#==========
nrow(tmC_rate_post_impute)

#Using True 5mCG to align RNA
dim(tmC_rate_post_impute)
pbmc.snmc2
table(pbmc.snmc2@meta.data$res.name.6g)

# using True 5mCG to align RNA
mat <- tmC_rate_post_impute #[,rownames(pbmc.snmc3@meta.data)]
activity.matrix <- log(1/(mat + 0.01))

snmc <- CreateSeuratObject(counts = activity.matrix)
snmc <- NormalizeData(snmc,verbose = F)
snmc <- FindVariableFeatures(snmc,verbose = F)
snmc <- ScaleData(snmc,verbose = F)
snmc <- RunPCA(snmc,verbose = F)
snmc <- RunUMAP(snmc,verbose = F,dims = 1:30)
dff.1 <- pbmc.snmc2@meta.data
rownames(dff.1) <- dff.1$cellID
snmc@meta.data <- dff.1[rownames(snmc@meta.data),]

# define the reference
ad.ref <- snRNA3
# define the query
ad.query <- snmc
# predictions
ad.anchors <- FindTransferAnchors(reference = ad.ref,
                                  query = ad.query,
                                  features = unique(c(VariableFeatures(object = ad.ref), VariableFeatures(object = ad.query))),
                                  #features = VariableFeatures(object = ad.ref),
                                  reduction = "cca",
                                  dims = 1:30)

predictions <- TransferData(anchorset = ad.anchors, 
                            refdata = ad.ref$res.name.6groups,
                            weight.reduction = ad.query[["pca"]],
                            dims = 1:30)
ad.query <- AddMetaData(ad.query, metadata = predictions)
datScore <- ad.query@meta.data

head(datScore,n=3)

plot_different_res_snRNA3.1 <- function(datScore,keyWord){

    nLen <- length(grep("prediction",colnames(datScore),value=T))
    keyPrediction <- grep("prediction",colnames(datScore),value=T)[1:nLen-1]
    datScore2 <- datScore
    datScore2$combined <- datScore[[keyWord]]
    datDF <- reshape2::melt(datScore2[,c(keyPrediction,"combined")],id="combined")
    datDF %>% group_by(combined,variable) %>% dplyr::summarise(score = median(value)) -> datDF.median
    datDF.median$variable2 <- gsub("prediction.score.","",datDF.median$variable)
    

    cluster.oder.our <- c("Ex","Inh","Striatum","Astro","Oligo","MG")
    cluster.oder.snRNA <- c("Ex","Inh","Striatum","Astro","Oligo","MG")
    datDF.median <- droplevels(datDF.median[datDF.median$combined %in% cluster.oder.our,])
    datDF.median <- droplevels(datDF.median[datDF.median$variable2 %in% cluster.oder.our,])
    datDF.median$combined <- factor(datDF.median$combined,levels=cluster.oder.our)
    datDF.median$variable2 <- factor(datDF.median$variable2,levels=cluster.oder.snRNA)
    
    # plot
    p <- ggplot(data = datDF.median,aes(x = combined, y = variable2, fill = score)) + 
         geom_tile(col="white") +
         theme_cowplot() +
         geom_text(aes(label=round(score,2))) +
         #geom_tile(width=0.975, height=0.975) + theme_bw() + coord_equal() + 
         theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
         scale_fill_gradient("Median score",low = "#568dba", high = "#f42222") +
         xlab("") + ylab("") + ggtitle(keyWord)

    return(p)
}

options(repr.plot.width = 5.2, repr.plot.height = 4)
px9.1 <- plot_different_res_snRNA3.1(datScore,keyWord = "res.name.6g") + ggtitle("Using tmCG signal")
px9.1


mat <- mC_rate_post_impute
activity.matrix <- log(1/(mat + 0.01))

class(activity.matrix)

snmc <- CreateSeuratObject(counts = activity.matrix)
snmc <- NormalizeData(snmc,verbose = F)
snmc <- FindVariableFeatures(snmc,verbose = F)
snmc <- ScaleData(snmc,verbose = F)
snmc <- RunPCA(snmc,verbose = F)
snmc <- RunUMAP(snmc,verbose = F,dims = 1:30)
dff.1 <- pbmc.snmc2@meta.data
rownames(dff.1) <- dff.1$cellID
snmc@meta.data <- dff.1[rownames(snmc@meta.data),]

# define the reference
ad.ref <- snRNA3
# define the query
ad.query <- snmc
# predictions
ad.anchors <- FindTransferAnchors(reference = ad.ref,
                                  query = ad.query,
                                  features = unique(c(VariableFeatures(object = ad.ref), VariableFeatures(object = ad.query))),
                                  #features = VariableFeatures(object = ad.ref),
                                  reduction = "cca",
                                  dims = 1:30)

predictions <- TransferData(anchorset = ad.anchors, 
                            refdata = ad.ref$res.name.6groups,
                            weight.reduction = ad.query[["pca"]],
                            dims = 1:30)
ad.query <- AddMetaData(ad.query, metadata = predictions)
datScore <- ad.query@meta.data

head(datScore,n=3)
options(repr.plot.width = 5.2, repr.plot.height = 4)
px9.2 <- plot_different_res_snRNA3.1(datScore,keyWord = "res.name.6g") + ggtitle("Using BSmCG signal")
px9.2


mCH_rate_post_impute2 <- readRDS(file="20211016_552_genebody_CpH_matrix_post_imputation.rds")
colnames(mCH_rate_post_impute2)[1:3]
colnames(pbmc.snmc2)[1:2]

mCH_rate_post_impute2.2 <- mCH_rate_post_impute2[,rownames(pbmc.snmc2@meta.data)]
dim(mCH_rate_post_impute2.2)

mat <- mCH_rate_post_impute2.2
activity.matrix <- log(1/(mat + 0.01))

class(activity.matrix)

snmc <- CreateSeuratObject(counts = activity.matrix)
snmc <- NormalizeData(snmc,verbose = F)
snmc <- FindVariableFeatures(snmc,verbose = F)
snmc <- ScaleData(snmc,verbose = F)
snmc <- RunPCA(snmc,verbose = F)
snmc <- RunUMAP(snmc,verbose = F,dims = 1:30)

snmc@meta.data <- pbmc.snmc2@meta.data

# define the reference
ad.ref <- snRNA3
# define the query
ad.query <- snmc
# predictions
ad.anchors <- FindTransferAnchors(reference = ad.ref,
                                  query = ad.query,
                                  features = unique(c(VariableFeatures(object = ad.ref), VariableFeatures(object = ad.query))),
                                  #features = VariableFeatures(object = ad.ref),
                                  reduction = "cca",
                                  dims = 1:30)

predictions <- TransferData(anchorset = ad.anchors, 
                            refdata = ad.ref$res.name.6groups,
                            weight.reduction = ad.query[["pca"]],
                            dims = 1:30)
ad.query <- AddMetaData(ad.query, metadata = predictions)
datScore <- ad.query@meta.data

head(datScore,n=3)
options(repr.plot.width = 5.2, repr.plot.height = 4)
px9.3 <- plot_different_res_snRNA3.1(datScore,keyWord = "res.name.6g") + ggtitle("Using BSmCG signal")
px9.3

options(repr.plot.width = 15, repr.plot.height = 3.5)

px9.1 + px9.2 + px9.3




