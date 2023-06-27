#!/usr/bin/Rscript
# or open it with jupyter notebook
# Title: "NATURE BIOTECHNOLOGY README"
# Author: "Peng"


# load the CpG clsutering results
pbmc.snmc <- readRDS(file="20210921_with_final_annotation_pbmc.snmc.rds")

library("Seurat")

library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(pheatmap)

require(dplyr)
require(tidyr)

pbmc.snmc

head(pbmc.snmc@meta.data)

table(pbmc.snmc@meta.data$res.name2)

pbmc.snmc@meta.data$res.name3 <- plyr::mapvalues(pbmc.snmc@meta.data$res.name2,
                from = c("L2/3","L4","L4/5","L6","DL","Vip/Ndnf","Pv","Sst","Striatum","In"),
                to = c("Ex-up","Ex-mid","Ex-mid","Ex-deep","Ex-deep","Inh-Vip/Ndnf","Inh-Pv","Inh-Sst","Striatum","In"))

Idents(pbmc.snmc) <- "res.name3"

table(pbmc.snmc@meta.data$res.name3)

pbmc.snmc@meta.data$res.name.6g <- plyr::mapvalues(pbmc.snmc@meta.data$res.name2,
                from = c("L2/3","L4","L4/5","L6","DL","Vip/Ndnf","Pv","Sst","Striatum","In"),
                to = c("Ex","Ex","Ex","Ex","Ex","Inh","Inh","Inh","Striatum","In"))

Idents(pbmc.snmc) <- "res.name.6g"

table(pbmc.snmc@meta.data$res.name.6g)

options(repr.plot.width = 5,repr.plot.height = 4)
DimPlot(pbmc.snmc, reduction = "tsne",group.by = "res.name2", pt.size = .3,label = T,repel = T)

options(repr.plot.width = 5,repr.plot.height = 4)
DimPlot(pbmc.snmc, reduction = "umap",group.by = "res.name2", pt.size = .3,label = T,repel = T)

options(repr.plot.width = 5,repr.plot.height = 4)
DimPlot(pbmc.snmc, reduction = "umap",group.by = "res.name3", pt.size = .3,label = T,repel = T)

getwd()

mC_nsite_rename <- readRDS(file="20210930_GeneBody_mC_nsite_rename.rds")
hmC_rate_rename <- readRDS(file="20210930_GeneBody_hmC_rate_rename.rds")
hmC_nsite_rename <- readRDS(file="20210930_GeneBody_hmC_nsite_rename.rds")

class(mC_rate_rename)
dim(mC_rate_rename)
summary(unique(rownames(mC_rate_rename)))
head(mC_rate_rename)
head(hmC_rate_rename)

filter_by_nSite_percent <- function(comb_nSite,comb_rate,nSite = 5, percent = 0.9){
    comb_rate[comb_nSite < nSite] <- NA
    na_num_per_bin <- apply(comb_rate,1,function(x){
        sum(is.na(x))
    })
    print(table(na_num_per_bin <= ncol(comb_rate) * (1-percent)))
}

filter_by_nsite_neg_percent_impute <- function(mC_nsite,hmC_nsite,mC_rate,hmC_rate,percent,nsite_cutoff){
    
    # define <0 as NA
    tmC_rate <- mC_rate - hmC_rate
#    tmC_rate[tmC_rate < 0] <- NA  # Peng's old analysis 
     tmC_rate[tmC_rate < 0] <- 0
    
    # define both nsite < nsite-cutoff as NA
    f1 <- (mC_nsite < nsite_cutoff)
    tmC_rate[f1] <- NA
    
    f2 <- (hmC_nsite < nsite_cutoff)
    tmC_rate[f2] <- NA
    
    na_num_per_bin <- apply(tmC_rate,1,function(x){
        sum(is.na(x))
    })
    print(dim(tmC_rate))
    print(table(na_num_per_bin <= ncol(tmC_rate) * (1-percent)))
    keep.features <- names(na_num_per_bin)[na_num_per_bin <= ncol(tmC_rate) * (1-percent)] # less than 5%
    comb_rate2 <- tmC_rate[keep.features,]
    # impute with mean
    m <- comb_rate2
    ## from https://stackoverflow.com/questions/6918086/replace-na-values-by-row-means
    k <- which(is.na(m), arr.ind=TRUE)
    m[k] <- rowMeans(m, na.rm=TRUE)[k[,1]]
    comb_rate3 <- m
    return(comb_rate3)
}

tmC_rate_post_impute2 <- filter_by_nsite_neg_percent_impute(mC_nsite=mC_nsite_rename,
                                                           hmC_nsite=hmC_nsite_rename,
                                                           mC_rate=mC_rate_rename,
                                                           hmC_rate=hmC_rate_rename,
                                                           percent=.8,
                                                           nsite_cutoff=5)
#if nsite_cutoff=5:
#[1] 53656   545
#FALSE  TRUE 
#45087  8569 

# Peng's old analysis
#[1] 53656   545
# FALSE  TRUE 
# 48389  5267 

tmC_rate_post_impute2 <- filter_by_nsite_neg_percent_impute(mC_nsite=mC_nsite_rename,
                                                           hmC_nsite=hmC_nsite_rename,
                                                           mC_rate=mC_rate_rename,
                                                           hmC_rate=hmC_rate_rename,
                                                           percent=.8,
                                                           nsite_cutoff=2)
#if nsite_cutoff=5:
#[1] 53656   545
#FALSE  TRUE 
#45087  8569 

# Peng's old analysis
#[1] 53656   545
# FALSE  TRUE 
# 48389  5267 

tmC_rate_post_impute <- filter_by_nsite_neg_percent_impute(mC_nsite=mC_nsite_rename,
                                                           hmC_nsite=hmC_nsite_rename,
                                                           mC_rate=mC_rate_rename,
                                                           hmC_rate=hmC_rate_rename,
                                                           percent=.9,
                                                           nsite_cutoff=5)
# Peng's old analysis
#[1] 53656   545

# FALSE  TRUE 
# 48389  5267 

dim(tmC_rate_post_impute2) # 12583 genes, 545 nuclei /// 2 CpG sites, 80% cells
head(tmC_rate_post_impute2)

dim(tmC_rate_post_impute2) # 5 CpG sites, 80% cells
head(tmC_rate_post_impute2)

dim(tmC_rate_post_impute) # 5 CpG sites, 90% cells

# impute here
filter_impute3 <- function(comb_nSite,comb_rate,nSite,keep.features){
    comb_rate[comb_nSite < nSite] <- NA
    comb_rate2 <- comb_rate[keep.features,]
    # impute with mean
    m <- comb_rate2
    ## from https://stackoverflow.com/questions/6918086/replace-na-values-by-row-means
    k <- which(is.na(m), arr.ind=TRUE)
    m[k] <- rowMeans(m, na.rm=TRUE)[k[,1]]
    comb_rate3 <- m
    return(comb_rate3)
}

mC_rate_post_impute <- filter_impute3(mC_nsite_rename,
                                       mC_rate_rename,
                                       nSite = 2,
                                       rownames(tmC_rate_post_impute2))
dim(mC_rate_post_impute) # 5267 545

hmC_rate_post_impute <- filter_impute3(hmC_nsite_rename,
                                       hmC_rate_rename,
                                       nSite = 2,
                                       rownames(tmC_rate_post_impute2))
dim(hmC_rate_post_impute) # 5267 545

dim(hmC_rate_post_impute)
head(hmC_rate_post_impute)

## load gene body CpH Here
GB_mCH_rate <- readRDS("20210915_552nuclei_genebody_CH_comb_rate.rds")
GB_mCH_nsite <- readRDS("20210915_552nuclei_genebody_CH_comb_nsite.rds")

mCH_rate_post_impute <- filter_impute3(GB_mCH_nsite, 
                                       GB_mCH_rate, 
                                       nSite =2, 
                                       rownames(tmC_rate_post_impute2))

dim(mCH_rate_post_impute) # 5267 552 here

colnames(mCH_rate_post_impute)[1:5]

colnames(mCH_rate_post_impute) <- gsub("CpH","CpG",colnames(mCH_rate_post_impute))

colnames(mCH_rate_post_impute)[1:3]

dff <- readRDS(file="sc-bACE-seq/20210921_metacell_snhmCseq.rds")
dff2 <- dff[,c('cell_snhmCseq','res.name2','cell_snmCseq2','cell_in_snmCseq2_matrix','cell_in_snhmCseq_matrix')]
dff2$cellID <- paste("Paired",1:545,sep="")
dff$cellID <- dff2$cellID
table(dff2$res.name)
table(dff2$res.name2)

rownames(pbmc.snmc@meta.data)[1:3]
dim(dff)
head(dff)
head(dff2)

# get the common cell id here
pbmc.snmc2 <- subset(pbmc.snmc,cells = dff2$cell_snmCseq2)

pbmc.snmc # 552 nuclei
pbmc.snmc2 # after align to paired seq, 545 nuclei retrieved

table(pbmc.snmc2@meta.data$input)
table(pbmc.snmc2@meta.data$res.name.6g)


options(repr.plot.width = 15,repr.plot.height = 4)
p1 <- DimPlot(pbmc.snmc2, reduction = "umap",group.by = "res.name2", pt.size = .5,label = T,repel = T)
p2 <- DimPlot(pbmc.snmc2, reduction = "umap",group.by = "res.name3", pt.size = .5,label = T) + scale_colour_manual(values = rev(brewer.pal(n = 11, name = "Set3")))
p3 <- DimPlot(pbmc.snmc2, reduction = "umap",group.by = "input", pt.size = .5,label = T) + scale_colour_manual(values = brewer.pal(n = 4, name = "Set1"))

plot_grid(p1,p2,p3, nrow=1)

options(repr.plot.width = 15,repr.plot.height = 6)
#p1 <- DimPlot(pbmc.snmc2, reduction = "umap",group.by = "res.name2", pt.size = .5,label = T,repel = T)
p2 <- DimPlot(pbmc.snmc2, reduction = "umap",group.by = "res.name.6g", pt.size = .5,label = T) + scale_colour_manual(values = rev(brewer.pal(n = 11, name = "Set3")))
p3 <- DimPlot(pbmc.snmc2, reduction = "umap",group.by = "input", pt.size = .5,label = T) + scale_colour_manual(values = brewer.pal(n = 4, name = "Set1"))

plot_grid(p3,p2, nrow=1)

Idents(pbmc.snmc2) <- 'res.name2'

pbmc.snmc3 <- subset(pbmc.snmc2,idents = c("In"),invert=T)
pbmc.snmc3@meta.data <- droplevels(pbmc.snmc3@meta.data)

table(pbmc.snmc3$res.name2)
table(pbmc.snmc3$res.name.6g)

options(repr.plot.width = 15,repr.plot.height = 4)
p1 <- DimPlot(pbmc.snmc3, reduction = "umap",group.by = "res.name2", pt.size = .5,label = T,repel = T)
p2 <- DimPlot(pbmc.snmc3, reduction = "umap",group.by = "res.name3", pt.size = .5,label = T) + scale_colour_manual(values = rev(brewer.pal(n = 11, name = "Set3")))
p3 <- DimPlot(pbmc.snmc3, reduction = "umap",group.by = "input", pt.size = .5,label = T) + scale_colour_manual(values = brewer.pal(n = 4, name = "Set1"))

plot_grid(p1,p2,p3, nrow=1)

# https://satijalab.org/seurat/reference/dimplot

options(repr.plot.width = 10,repr.plot.height = 4)
#p1 <- DimPlot(pbmc.snmc3, reduction = "umap",group.by = "res.name2", pt.size = .5,label = T,repel = T)
p.6g <- DimPlot(pbmc.snmc3, reduction = "umap",group.by = "res.name.6g", pt.size = 3,label = T, raster = T) # + scale_colour_manual(values = brewer.pal(n = 11, name = "Set2"))
p.input <- DimPlot(pbmc.snmc3, reduction = "umap",group.by = "input", pt.size = 3,label = T, raster = T) + scale_colour_manual(values = brewer.pal(n = 4, name = "Spectral")) # Spectral Set3

plot_grid(p.input,p.6g, nrow=1)

ggsave("/mnt/data1/haowu/nfs2/sc-bACE-seq/Fig3_all.nuclei.UMAP_raster_01.04.2022.pdf",width=10,height=4)

pbmc.snmc3 # 512 nuclei after exculding 33 In nuclei

replace_cell_name_vector <- dff2$cell_snmCseq2
names(replace_cell_name_vector) <- dff2$cellID
replace_cell_name_vector[1:3]

hmC_matrix <- hmC_rate_post_impute
mC_matrix <- mC_rate_post_impute
tmC_matrix <- tmC_rate_post_impute2
mCH_matrix <- mCH_rate_post_impute
dim(hmC_matrix)
dim(mC_matrix)
dim(tmC_matrix)
dim(mCH_matrix)

colnames(hmC_matrix) <- replace_cell_name_vector[colnames(hmC_matrix)]
colnames(mC_matrix) <- replace_cell_name_vector[colnames(mC_matrix)]
colnames(tmC_matrix) <- replace_cell_name_vector[colnames(tmC_matrix)]
#colnames(mCH_matrix) <- replace_cell_name_vector[colnames(mCH_matrix)]

hmC_matrix <- hmC_matrix[,rownames(pbmc.snmc3@meta.data)]
mC_matrix <- mC_matrix[,rownames(pbmc.snmc3@meta.data)]
tmC_matrix <- tmC_matrix[,rownames(pbmc.snmc3@meta.data)]
mCH_matrix <- mCH_matrix[,rownames(pbmc.snmc3@meta.data)]

dim(hmC_matrix)
dim(mC_matrix)
dim(tmC_matrix)
dim(mCH_matrix)

colnames(mC_matrix)[1:3]
colnames(mCH_matrix)[1:3]
head(mC_matrix)

ratio_matrix <- hmC_matrix/mC_matrix

sum( hmC_matrix >= mC_matrix)
sum( hmC_matrix < mC_matrix)

ratio_matrix[ hmC_matrix >= mC_matrix ] <- 1
#ratio_matrix[ mC_matrix = 0] <- 0

min(ratio_matrix)
max(ratio_matrix)

pbmc.snmc3[["hmCG"]] <- CreateAssayObject(counts = hmC_matrix)
pbmc.snmc3 <- NormalizeData(pbmc.snmc3,assay = "hmCG")
pbmc.snmc3 <- ScaleData(pbmc.snmc3,assay = "hmCG")

pbmc.snmc3[["tmCG"]] <- CreateAssayObject(counts = tmC_matrix)
pbmc.snmc3 <- NormalizeData(pbmc.snmc3,assay = "tmCG")
pbmc.snmc3 <- ScaleData(pbmc.snmc3,assay = "tmCG")

pbmc.snmc3[["BSmCG"]] <- CreateAssayObject(counts = mC_matrix)
pbmc.snmc3 <- NormalizeData(pbmc.snmc3,assay = "BSmCG")
pbmc.snmc3 <- ScaleData(pbmc.snmc3,assay = "BSmCG")

pbmc.snmc3[["BSmCH"]] <- CreateAssayObject(counts = mCH_matrix)
pbmc.snmc3 <- NormalizeData(pbmc.snmc3,assay = "BSmCH")
pbmc.snmc3 <- ScaleData(pbmc.snmc3,assay = "BSmCH")

dim(hmC_matrix)
head(hmC_matrix)

dim(pbmc.snmc3@assays$hmCG@data)
#head(pbmc.snmc3@assays$hmCG@data)

pbmc.snmc3[["ratio"]] <- CreateAssayObject(counts = ratio_matrix)
pbmc.snmc3 <- NormalizeData(pbmc.snmc3,assay = "ratio")
pbmc.snmc3 <- ScaleData(pbmc.snmc3,assay = "ratio")

colnames(pbmc.snmc3@meta.data)

?FeaturePlot

options(repr.plot.width=20,repr.plot.height=8)
FeaturePlot(pbmc.snmc3,c("mm10_dedup_filter_CG",
                    "mm10_dedup_filter_CH",
                    "mm10_dedup_filter_covered_CG_num",
                    "mm10_dedup_filter_covered_CH_num",
                    "lambda_dedup_filter_CG",
                    "lambda_dedup_filter_CH",
                    "mm10_dedup_filter_read_num"
                    ),
            min.cutoff = 0,
            max.cutoff = 4,
            keep.scale = "all",
            ncol=4)

options(repr.plot.width=20,repr.plot.height=8)
FeaturePlot(pbmc.snmc3,c("mm10_dedup_filter_CG",
                    "mm10_dedup_filter_CH",
                    "mm10_dedup_filter_covered_CG_num",
                    "mm10_dedup_filter_covered_CH_num",
                    "lambda_dedup_filter_CG",
                    "lambda_dedup_filter_CH",
                    "mm10_dedup_filter_read_num"),ncol=4)

ggsave("HW-Neuronal.NonNeu.UMAP.w.In_01.05.2022.pdf",width=20,height=8)

genes.in <- c("Satb2","Slc17a7","Gad1","Gad2", ## 
             "Pvalb","Sst","Reln","Ndnf","Vip","Rorb","Foxp2","Tle4",
             "Cux1","Cux2","Slc6a1","Sulf1","Prox1","Grik3","Slc1a2","Csf1r",
             "Mbp","Cspg4","Nxn","Lcp2","Mef2c","Pdgfra","Ppp1r1b", "Meis2")

genes.in.final <- c("Slc17a7","Lingo1","Sv2b","Neurod6","Kcnma1", # Ex marker
                    "Gad1","Gad2","Abat","Slc6a1","Stat5b","Erbb4", # Inh marker
                    "Ppp1r1b", "Meis2", # Str marker
                    "Slc1a2","Trim9","Paqr8","Cpe","Slc1a3", # Astro marker
                    "Mobp","Mbp","Mog","Pde4b","Klh12", # Oligo marker
                    "Csf1r","P2ry13","Fcrls","Zfp710","Laptm5" # MG marker
                   ) 

genes.in2 <- intersect(genes.in.final,rownames(tmC_matrix))

DefaultAssay(pbmc.snmc3) <- "BSmCG"
#pbmc.snmc3
options(repr.plot.width = 12,repr.plot.height = 12)
FeaturePlot(pbmc.snmc3,features = genes.in2,raster = TRUE,ncol = 4,pt.size = 4,slot = "scale.data")& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

FeaturePlot(pbmc.snmc3,features = genes.in2,
            min.cutoff = -2,
            max.cutoff = 2,
            raster = TRUE,ncol = 4,pt.size = 4,slot = "scale.data")& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
#FeaturePlot(pbmc.snmc2,features = genes.in2,raster = TRUE,ncol = 4,pt.size = 4,slot = "counts",min.cutoff = 0)
#ggsave("Fig4_removeIn_BSmCG.featureplot_umap.pdf",width=12,height=20)

DefaultAssay(pbmc.snmc3) <- "BSmCG"
#pbmc.snmc3
options(repr.plot.width = 12,repr.plot.height = 6)
FeaturePlot(pbmc.snmc3,features = c("Gad1","Gad2"),
#            min.cutoff = -2,
#            max.cutoff = 2,
            raster = TRUE,ncol = 2,pt.size = 4,slot = "scale.data")& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
#FeaturePlot(pbmc.snmc2,features = genes.in2,raster = TRUE,ncol = 4,pt.size = 4,slot = "counts",min.cutoff = 0)
#ggsave("Fig4_removeIn_BSmCG.featureplot_umap.pdf",width=12,height=20)

DefaultAssay(pbmc.snmc3) <- "tmCG"
#pbmc.snmc3
options(repr.plot.width = 12,repr.plot.height = 12)
FeaturePlot(pbmc.snmc3,features = genes.in2,raster = TRUE,ncol = 4,pt.size = 4,slot = "scale.data")& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
#FeaturePlot(pbmc.snmc2,features = genes.in2,raster = TRUE,ncol = 4,pt.size = 4,slot = "counts",min.cutoff = 0)
#ggsave("Fig4_removeIn_BSmCG.featureplot_umap.pdf",width=12,height=20)

DefaultAssay(pbmc.snmc3) <- "ratio"
#pbmc.snmc3
options(repr.plot.width = 12,repr.plot.height = 12)
FeaturePlot(pbmc.snmc3,features = genes.in2,raster = TRUE,ncol = 4,pt.size = 4,slot = "scale.data")& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
#FeaturePlot(pbmc.snmc2,features = genes.in2,raster = TRUE,ncol = 4,pt.size = 4,slot = "counts",min.cutoff = 0)
#ggsave("Fig4_removeIn_BSmCG.featureplot_umap.pdf",width=12,height=20)

DefaultAssay(pbmc.snmc3) <- "hmCG"
#pbmc.snmc3
options(repr.plot.width = 12,repr.plot.height = 12)
FeaturePlot(pbmc.snmc3,features = genes.in2,raster = TRUE,ncol = 4,pt.size = 4,slot = "scale.data")& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
#FeaturePlot(pbmc.snmc2,features = genes.in2,raster = TRUE,ncol = 4,pt.size = 4,slot = "counts",min.cutoff = 0)
#ggsave("Fig4_removeIn_BSmCG.featureplot_umap.pdf",width=12,height=20)

genes.in.final <- c("Slc17a7","Lingo1", # "Sv2b","Neurod6","Kcnma1", # Ex marker
                    "Slc6a1","Gad2", #"Gad2","Abat","Stat5b","Erbb4", # Inh marker
                    "Ppp1r1b", "Meis2", # Str marker
                    "Slc1a3", "Paqr8", # "Cpe","Slc1a2","Trim9", # Astro marker
                    "Mobp","Mog", # "Mbp","Pde4b","Klh12", # Oligo marker
                    "Csf1r","Fcrls"  # "P2ry13","Zfp710","Laptm5" # MG marker
                   ) 

genes.in2 <- intersect(genes.in.final,rownames(tmC_matrix))

DefaultAssay(pbmc.snmc3) <- "tmCG"
#pbmc.snmc3
options(repr.plot.width = 18,repr.plot.height = 8)
#FeaturePlot(pbmc.snmc3,features = genes.in2,raster = TRUE,ncol = 4,pt.size = 4,slot = "scale.data")& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

FeaturePlot(pbmc.snmc3,features = genes.in2,
            min.cutoff = -2,
            max.cutoff = 2,
            raster = TRUE,ncol = 4,pt.size = 4,slot = "scale.data")& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))


DefaultAssay(pbmc.snmc3) <- "hmCG"
#pbmc.snmc3
options(repr.plot.width = 18,repr.plot.height = 8)
#FeaturePlot(pbmc.snmc3,features = genes.in2,raster = TRUE,ncol = 4,pt.size = 4,slot = "scale.data")& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

FeaturePlot(pbmc.snmc3,features = genes.in2,
            min.cutoff = -2,
            max.cutoff = 2,
            raster = TRUE,ncol = 4,pt.size = 4,slot = "scale.data")& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))


DefaultAssay(pbmc.snmc3) <- "ratio"
#pbmc.snmc3
options(repr.plot.width = 18,repr.plot.height = 8)
#FeaturePlot(pbmc.snmc3,features = genes.in2,raster = TRUE,ncol = 4,pt.size = 4,slot = "scale.data")& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

FeaturePlot(pbmc.snmc3,features = genes.in2,
#            min.cutoff = -2,
#            max.cutoff = 2,
            raster = TRUE,ncol = 4,pt.size = 4,slot = "scale.data")& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))



getwd()

saveRDS(pbmc.snmc3,file="2022.06.13_pbmc.snmc3_12583genes.rds") # 80% cells, 2 CpG

saveRDS(pbmc.snmc3,file="2022.03.17_pbmc.snmc3_14590genes.rds") # 80% cells, 1 CpG

saveRDS(pbmc.snmc3,file="2022.03.09_pbmc.snmc3_8569genes.rds") # 80% cells, 5 CpGs

saveRDS(pbmc.snmc3,file="2022.03.09_pbmc.snmc3_11994genes.rds") # 75% cells, 3 CpGs

# load the snRNA matrix
snRNA3 <- readRDS(file="20211016-rename-snRNA3-match-JE.rds")
snRNA3

colnames(snRNA3@meta.data)

table(snRNA3@meta.data$res.name.JE2)
summary(snRNA3@meta.data$res.name.JE)

head(snRNA3@meta.data)

sort(unique(pbmc.snmc3$res.name2))

snRNA3@meta.data$res.name.JE3 <- plyr::mapvalues(snRNA3@meta.data$res.name.JE2,
                                    from = c("DL_1-2","L2_3","L4_5","Ndnf_Vip","Str"),
                                    to = c("DL","L2/3","L4/5","Vip/Ndnf","Striatum"))

sort(unique(snRNA3@meta.data$res.name.JE3))

cluster.oder <- c("L2/3","L4","L4/5","L6","DL","Vip/Ndnf","Pv","Sst","Striatum","Astro","Oligo","MG")

options(repr.plot.width = 15,repr.plot.height = 6)
#DimPlot(snRNA3, reduction = "umap",group.by = "res.name.JE3", pt.size = .1,label = T,repel = T)

pt1 <- DimPlot(snRNA3, reduction = "umap",group.by = "orig.ident", pt.size = .2,label = T) + ggtitle("FANS channel") + scale_colour_manual(values = brewer.pal(n = 4, name = "Set1"))
pt2 <- DimPlot(snRNA3, reduction = "umap",group.by = "res.name.JE3", pt.size = .2,label = T) + ggtitle("Cell types") + scale_colour_manual(values = brewer.pal(n = 12, name = "Set3"))

plot_grid(pt1,pt2,nrow=1,align = "vh")

ggsave("HW-sNucDropseq_Neuronal.NonNeu.UMAP_06.15.2022.pdf",width=15,height=6)

pbmc.snmc3

table(Idents(pbmc.snmc3))

install.packages('BiocManager')
BiocManager::install('limma')

tmC.gene <- FindAllMarkers(pbmc.snmc3,logfc.threshold = 0,assay = "tmCG")

head(tmC.gene,n=5)
dim(tmC.gene) # 90% cells: 19119 rows; 80% cells, 2 CpGs: 27006 genes; 
summary(unique(tmC.gene$gene)) # 90% cells: 5539 genes; 80% cells, 2 CpGs; 9346 genes; 

head(tmC.gene,n=2) # old version

require(dplyr)
require(tidyr)

tmC.gene2 <- tmC.gene %>% dplyr::filter(avg_log2FC < 0,p_val_adj < .1) %>% droplevels

cluster.order <- c("L2/3","L4","L4/5","L6","DL","Vip/Ndnf","Pv","Sst","Striatum","Astro","Oligo","MG")
tmC.gene2$cluster <- factor(tmC.gene2$cluster,levels = cluster.order)

tmC.gene3 <- tmC.gene2 %>% arrange(cluster,p_val_adj) 

gene.order <- tmC.gene3$gene

length(gene.order) # p_val_adj<.05: 3263 genes; p_val_adj<.1: 3586 genes (90% cells), 3931 genes (80% cells)

gene.order.by.tmCG <- intersect(gene.order,rownames(snRNA3@assays$RNA@data))
length(gene.order.by.tmCG) # p_val_adj<.05: 1777 genes; p_val_adj<.1: 1917 genes (90% cells), 1920 genes (80% cells)

cluster.order

table(pbmc.snmc3@meta.data$res.name2)
table(pbmc.snmc3@meta.data$res.name3)

dim(pbmc.snmc3@meta.data)
head(rownames(pbmc.snmc3@meta.data))

table(Idents(pbmc.snmc3))

hmC.gene <- FindAllMarkers(pbmc.snmc3,logfc.threshold = 0,assay = "hmCG")

dim(hmC.gene)
summary(unique(hmC.gene$gene))
head(hmC.gene, n=5)

hmC.gene2 <- hmC.gene %>% dplyr::filter(avg_log2FC >0, p_val_adj < .1) %>% droplevels

cluster.order <- c("L2/3","L4","L4/5","L6","DL","Vip/Ndnf","Pv","Sst","Striatum","Astro","Oligo","MG")
hmC.gene2$cluster <- factor(hmC.gene2$cluster,levels = cluster.order)

hmC.gene3 <- hmC.gene2 %>% arrange(cluster,p_val_adj) 

gene.order <- hmC.gene3$gene
length(gene.order)

table(hmC.gene3$cluster)
hmC.gene3 %>% filter(cluster=="L2/3") 

hmC.gene3 %>% filter(cluster=="L6") 

hmC.gene3 %>% filter(cluster=="Pv")

hmC.gene3 %>% filter(cluster=="MG")

gene.order.by.hmCG <- intersect(gene.order,rownames(snRNA3@assays$RNA@data))
length(gene.order.by.hmCG)
head(gene.order.by.hmCG)

table(Idents(pbmc.snmc3))
pbmc.snmc3.neu <- subset(pbmc.snmc3,idents = c("L2/3","L4","L4/5","L6","DL","Pv","Sst","Vip/Ndnf"))

table(Idents(pbmc.snmc3.neu))

hmC.gene.neu <- FindAllMarkers(pbmc.snmc3.neu,logfc.threshold = 0,assay = "hmCG")

cluster.order.neu <- c("L2/3","L4","L4/5","L6","DL","Vip/Ndnf","Pv","Sst")
hmC.gene.neu$cluster <- factor(hmC.gene.neu$cluster,levels = cluster.order.neu)

dim(hmC.gene.neu)
#summary(unique(hmC.gene$gene))
head(hmC.gene.neu, n=5)

hmC.gene2.neu <- hmC.gene.neu %>% dplyr::filter(avg_log2FC >0, p_val_adj < .1) %>% droplevels
hmC.gene3.neu <- hmC.gene2.neu %>% arrange(cluster,p_val_adj) 
length(hmC.gene3.neu$gene)
table(hmC.gene3.neu$cluster)

hmC.gene3.neu %>% filter(cluster=="L6") 

hmC.gene3.neu %>% filter(cluster=="Pv") 

gene.order.by.hmCG.neu <- intersect(hmC.gene3.neu$gene,rownames(snRNA3@assays$RNA@data))
length(gene.order.by.hmCG.neu)

dim(hmC_matrix)

head(hmC_matrix)

# new analysis (80% cells, 2 CpGs):
gene.order.all <- rownames(hmC_matrix)
length(gene.order.all)
length(rownames(snRNA3@assays$RNA@data))
gene.order.all.rna <- intersect(gene.order.all,rownames(snRNA3@assays$RNA@data))
length(gene.order.all.rna)
unique(pbmc.snmc3@meta.data$res.name2)

# old analysis (80% cells, 5 CpGs):
gene.order.all <- rownames(hmC_matrix)
length(gene.order.all)
length(rownames(snRNA3@assays$RNA@data))
gene.order.all.rna <- intersect(gene.order.all,rownames(snRNA3@assays$RNA@data))
length(gene.order.all.rna)
unique(pbmc.snmc3@meta.data$res.name2)

cluster.order <- c("L2/3","L4","L4/5","L6","DL","Vip/Ndnf","Pv","Sst","Striatum","Astro","Oligo","MG")

### how to use sapply: https://r-coder.com/sapply-function-r/
mat_hmCG_celltype_all.genes <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name2 == x]
    return(rowMeans(hmC_matrix[gene.order.all.rna,cells.in]))
})

head(mat_hmCG_celltype_all.genes)

mat_hmCG_celltype_all.genes <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name2 == x]
    return(rowMeans(hmC_matrix[gene.order.all.rna,cells.in]))
})

mat_tmCG_celltype_all.genes <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name2 == x]
    return(rowMeans(tmC_matrix[gene.order.all.rna,cells.in]))
})

mat_mCG_celltype_all.genes <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name2 == x]
    return(rowMeans(mC_matrix[gene.order.all.rna,cells.in]))
})

mat_mCH_celltype_all.genes <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name2 == x]
    return(rowMeans(mCH_matrix[gene.order.all.rna,cells.in]))
})

mat_ratio_celltype_all.genes <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name2 == x]
    return(rowMeans(ratio_matrix[gene.order.all.rna,cells.in]))
})

mat_rna_celltype_all.genes <- sapply(cluster.order,function(x){
    cells.in <- rownames(snRNA3@meta.data)[snRNA3@meta.data$res.name.JE3 == x]
    mat <- expm1(snRNA3@assays$RNA@data)
    return(Matrix::rowMeans(mat[gene.order.all.rna,cells.in]))
})


dim(mat_hmCG_celltype_all.genes)
dim(mat_rna_celltype_all.genes)
head(mat_tmCG_celltype_all.genes)

class(mat_hmCG_celltype_all.genes)
class(hmC_matrix)

df_mCH_celltype_all.genes<-as.data.frame(mat_mCH_celltype_all.genes)
head(df_mCH_celltype_all.genes)
dim(df_mCH_celltype_all.genes)

mat_CG_celltype_all.genes_v2 <- 1 - mat_mCG_celltype_all.genes
sum(mat_CG_celltype_all.genes_v2 >= 0)
sum(mat_CG_celltype_all.genes_v2 < 0)
head(mat_CG_celltype_all.genes_v2)
head(1 - mat_tmCG_celltype_all.genes - mat_hmCG_celltype_all.genes)

mat_CG_celltype_all.genes <- 1 - mat_tmCG_celltype_all.genes - mat_hmCG_celltype_all.genes
#head(mat_CG_celltype_all.genes)
#mat_CG_celltype_all.genes[,"L2/3"]
#mat_CG_celltype_all.genes["Erbb4",]

df_CG_celltype_all.genes<-as.data.frame(mat_CG_celltype_all.genes)
colnames(df_CG_celltype_all.genes) <- paste(colnames(df_CG_celltype_all.genes),"CG",sep="_")

df_tmCG_celltype_all.genes<-as.data.frame(mat_tmCG_celltype_all.genes)
colnames(df_tmCG_celltype_all.genes) <- paste(colnames(df_tmCG_celltype_all.genes),"tmCG",sep="_")

df_hmCG_celltype_all.genes<-as.data.frame(mat_hmCG_celltype_all.genes)
colnames(df_hmCG_celltype_all.genes) <- paste(colnames(df_hmCG_celltype_all.genes),"hmCG",sep="_")

df_mCH_celltype_all.genes<-as.data.frame(mat_mCH_celltype_all.genes)
colnames(df_mCH_celltype_all.genes) <- paste(colnames(df_mCH_celltype_all.genes),"mCH",sep="_")

df_rna_celltype_all.genes<-as.data.frame(mat_rna_celltype_all.genes)
colnames(df_rna_celltype_all.genes) <- paste(colnames(df_rna_celltype_all.genes),"rna",sep="_")

head(df_CG_celltype_all.genes)
head(df_tmCG_celltype_all.genes)
head(df_hmCG_celltype_all.genes)
head(df_mCH_celltype_all.genes)
head(df_rna_celltype_all.genes)

df_5modalities_celltype_all.genes<-cbind(df_CG_celltype_all.genes, 
                                         df_tmCG_celltype_all.genes, 
                                         df_hmCG_celltype_all.genes, 
                                         df_mCH_celltype_all.genes,
                                         df_rna_celltype_all.genes)

colnames(df_5modalities_celltype_all.genes) = gsub("/", ".", colnames(df_5modalities_celltype_all.genes))

head(df_5modalities_celltype_all.genes)


df_4modalities_celltype_all.genes<-cbind(df_CG_celltype_all.genes, 
                                         df_tmCG_celltype_all.genes, 
                                         df_hmCG_celltype_all.genes, 
                                         df_rna_celltype_all.genes)

colnames(df_4modalities_celltype_all.genes) = gsub("/", ".", colnames(df_4modalities_celltype_all.genes))

head(df_4modalities_celltype_all.genes)

#sum( hmC_matrix >= mC_matrix)
#sum( hmC_matrix < mC_matrix)

#ratio_matrix[ hmC_matrix >= mC_matrix ] <- 1
#ratio_matrix[ mC_matrix = 0] <- 0

df_rna_celltype_all.genes[c("Arpp21","Enpp2","Cux2","Rorb","Tle4","Slc6a1","Erbb4","Reln","Adarb2","Cacna2d2"), ]

L2.3_CG.low_hmC.hi_rna.hi_genes <- df_4modalities_celltype_all.genes %>% filter(L2.3_CG<0.2 & L2.3_hmCG>0.3 & L2.3_rna>2)
L2.3_CG.low_hmC.hi_rna.hi_genes

L2.3_CG.hi_hmC.low_rna.hi_genes <- df_4modalities_celltype_all.genes %>% filter(L2.3_CG>0.5 & L2.3_hmCG<0.3 & L2.3_rna>2)
L2.3_CG.hi_hmC.low_rna.hi_genes

df_4modalities_celltype_all.genes["Tle4", ]
L6_CG.low_hmC.hi_rna.hi_genes <- df_4modalities_celltype_all.genes %>% filter(L6_CG<0.2 & L6_hmCG>0.3 & L6_rna>2)
dim(L6_CG.low_hmC.hi_rna.hi_genes)

L6_CG.hi_hmC.low_rna.hi_genes <- df_4modalities_celltype_all.genes %>% filter(L6_CG>0.5 & L6_hmCG<0.3 & L6_rna>2)
dim(L6_CG.hi_hmC.low_rna.hi_genes)

library(ggtern)

sessionInfo()  # 10-29-2022

sessionInfo()  # 09-13-2022

colnames(df_5modalities_celltype_all.genes)

#pdf("./Figure5/ggtern/ggtern_ATAC.VIP_density.rgbw.pdf")
ggtern(data=df_4modalities_celltype_all.genes, aes(x=L2.3_tmCG, y=L2.3_hmCG, z=L2.3_CG)) +

#  geom_point(alpha = 0.01, size=1/100) +
# stat_density_tern(geom='polygon', aes(fill=..level.., alpha=..level..)) + 
  stat_density_tern(geom='polygon', aes(fill=..level..), base='identity', bins=100) +

  scale_fill_gradient(low='blue', high='red') +
#  scale_fill_gradient2(high='red') +

#  geom_segment(data=lines, aes(x, y, z, xend = xend, yend = yend, zend = zend), color = 'gray20', size = 0.5) +

#  theme_bw() +
  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +

  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +

#  guides(fill = guide_colourbar(order=3)) +
#  labs(fill = '5hmC signals') +
  ggtitle('Ex-L2.3')

library(gridExtra)

options(repr.plot.width = 15,repr.plot.height = 4)


L2.3<-ggtern(data=df_4modalities_celltype_all.genes, aes(x=L2.3_tmCG, y=L2.3_hmCG, z=L2.3_CG)) +
  stat_density_tern(geom='polygon', aes(fill=..level..), base='identity', bins=100) +
  scale_fill_gradient(low='blue', high='red') +
   theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +
  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +
  ggtitle('Ex-L2.3')

L4<-ggtern(data=df_4modalities_celltype_all.genes, aes(x=L4_tmCG, y=L4_hmCG, z=L4_CG)) +
  stat_density_tern(geom='polygon', aes(fill=..level..), base='identity', bins=100) +
  scale_fill_gradient(low='blue', high='red') +
   theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +
  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +
  ggtitle('Ex-L4')

L5<-ggtern(data=df_4modalities_celltype_all.genes, aes(x=L4.5_tmCG, y=L4.5_hmCG, z=L4.5_CG)) +
  stat_density_tern(geom='polygon', aes(fill=..level..), base='identity', bins=100) +
  scale_fill_gradient(low='blue', high='red') +
   theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +
  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +
  ggtitle('Ex-L5')

L6<-ggtern(data=df_4modalities_celltype_all.genes, aes(x=L6_tmCG, y=L6_hmCG, z=L6_CG)) +
  stat_density_tern(geom='polygon', aes(fill=..level..), base='identity', bins=100) +
  scale_fill_gradient(low='blue', high='red') +
   theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +
  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +
  ggtitle('Ex-L6')

DL<-ggtern(data=df_4modalities_celltype_all.genes, aes(x=DL_tmCG, y=DL_hmCG, z=DL_CG)) +
  stat_density_tern(geom='polygon', aes(fill=..level..), base='identity', bins=100) +
  scale_fill_gradient(low='blue', high='red') +
   theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +
  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +
  ggtitle('Ex-DL')

Sst<-ggtern(data=df_4modalities_celltype_all.genes, aes(x=Sst_tmCG, y=Sst_hmCG, z=Sst_CG)) +
  stat_density_tern(geom='polygon', aes(fill=..level..), base='identity', bins=100) +
  scale_fill_gradient(low='blue', high='red') +
   theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +
  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +
  ggtitle('Inh-Sst')

Pv<-ggtern(data=df_4modalities_celltype_all.genes, aes(x=Pv_tmCG, y=Pv_hmCG, z=Pv_CG)) +
  stat_density_tern(geom='polygon', aes(fill=..level..), base='identity', bins=100) +
  scale_fill_gradient(low='blue', high='red') +
   theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +
  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +
  ggtitle('Inh-Pv')

Vip.Ndnf<-ggtern(data=df_4modalities_celltype_all.genes, aes(x=Vip.Ndnf_tmCG, y=Vip.Ndnf_hmCG, z=Vip.Ndnf_CG)) +
  stat_density_tern(geom='polygon', aes(fill=..level..), base='identity', bins=100) +
  scale_fill_gradient(low='blue', high='red') +
   theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +
  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +
  ggtitle('Inh-Vip.Ndnf')

Striatum<-ggtern(data=df_4modalities_celltype_all.genes, aes(x=Striatum_tmCG, y=Striatum_hmCG, z=Striatum_CG)) +
  stat_density_tern(geom='polygon', aes(fill=..level..), base='identity', bins=100) +
  scale_fill_gradient(low='blue', high='red') +
   theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +
  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +
  ggtitle('Striatum')

Oligo<-ggtern(data=df_4modalities_celltype_all.genes, aes(x=Oligo_tmCG, y=Oligo_hmCG, z=Oligo_CG)) +
  stat_density_tern(geom='polygon', aes(fill=..level..), base='identity', bins=100) +
  scale_fill_gradient(low='blue', high='red') +
   theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +
  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +
  ggtitle('Oligo')

MG<-ggtern(data=df_4modalities_celltype_all.genes, aes(x=MG_tmCG, y=MG_hmCG, z=MG_CG)) +
  stat_density_tern(geom='polygon', aes(fill=..level..), base='identity', bins=100) +
  scale_fill_gradient(low='blue', high='red') +
   theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +
  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +
  ggtitle('MG')

Astro<-ggtern(data=df_4modalities_celltype_all.genes, aes(x=Astro_tmCG, y=Astro_hmCG, z=Astro_CG)) +
  stat_density_tern(geom='polygon', aes(fill=..level..), base='identity', bins=100) +
  scale_fill_gradient(low='blue', high='red') +
   theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +
  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +
  ggtitle('Astro')

###
plot_list=list()
plot_list[["L2.3"]] <- L2.3
plot_list[["L4"]] <- L4
plot_list[["L5"]] <- L5
plot_list[["L6"]] <- L6
plot_list[["DL"]] <- DL
plot_list[["Pv"]] <- Pv
plot_list[["Sst"]] <- Sst
plot_list[["Vip.Ndnf"]] <- Vip.Ndnf
plot_list[["Striatum"]] <- Striatum
plot_list[["Astro"]] <- Astro
plot_list[["Oligo"]] <- Oligo
plot_list[["MG"]] <- MG

#plot_list

#library(gridExtra)

options(repr.plot.width = 15,repr.plot.height = 9)

ternary_plot <- ggtern::grid.arrange(grobs=plot_list,nrow = 3)

ggsave("Fig6_12ternary.plots_09.13.2022.pdf",ternary_plot,width=15,height=9)

options(repr.plot.width = 15,repr.plot.height = 4)

myPalette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(-6, 6))


L2.3<-ggtern(data=df_4modalities_celltype_all.genes, aes(x=L2.3_tmCG, y=L2.3_hmCG, z=L2.3_CG)) +
  geom_mask() +
  geom_point(aes(colour=log2(L2.3_rna)), size=0.2, alpha=0.5) +
  sc + 
   theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +
  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +
  ggtitle('Ex-L2.3')

L4<-ggtern(data=df_4modalities_celltype_all.genes, aes(x=L4_tmCG, y=L4_hmCG, z=L4_CG)) +
  geom_mask() +
  geom_point(aes(colour=log2(L4_rna)), size=0.2, alpha=0.5) +
  sc + 
   theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +
  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +
  ggtitle('Ex-L4')

L5<-ggtern(data=df_4modalities_celltype_all.genes, aes(x=L4.5_tmCG, y=L4.5_hmCG, z=L4.5_CG)) +
  geom_mask() +
  geom_point(aes(colour=log2(L4.5_rna)), size=0.2, alpha=0.5) +
  sc + 
   theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +
  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +
  ggtitle('Ex-L5')

L6<-ggtern(data=df_4modalities_celltype_all.genes, aes(x=L6_tmCG, y=L6_hmCG, z=L6_CG)) +
  geom_mask() +
  geom_point(aes(colour=log2(L6_rna)), size=0.2, alpha=0.5) +
  sc + 
   theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +
  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +
  ggtitle('Ex-L6')

DL<-ggtern(data=df_4modalities_celltype_all.genes, aes(x=DL_tmCG, y=DL_hmCG, z=DL_CG)) +
  geom_mask() +
  geom_point(aes(colour=log2(DL_rna)), size=0.2, alpha=0.5) +
  sc + 
   theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +
  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +
  ggtitle('Ex-DL')

Sst<-ggtern(data=df_4modalities_celltype_all.genes, aes(x=Sst_tmCG, y=Sst_hmCG, z=Sst_CG)) +
  geom_mask() +
  geom_point(aes(colour=log2(Sst_rna)), size=0.2, alpha=0.5) +
  sc + 
   theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +
  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +
  ggtitle('Inh-Sst')

Pv<-ggtern(data=df_4modalities_celltype_all.genes, aes(x=Pv_tmCG, y=Pv_hmCG, z=Pv_CG)) +
  geom_mask() +
  geom_point(aes(colour=log2(Pv_rna)), size=0.2, alpha=0.5) +
  sc + 
   theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +
  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +
  ggtitle('Inh-Pv')

Vip.Ndnf<-ggtern(data=df_4modalities_celltype_all.genes, aes(x=Vip.Ndnf_tmCG, y=Vip.Ndnf_hmCG, z=Vip.Ndnf_CG)) +
  geom_mask() +
  geom_point(aes(colour=log2(Vip.Ndnf_rna)), size=0.2, alpha=0.5) +
  sc + 
   theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +
  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +
  ggtitle('Inh-Vip.Ndnf')

Striatum<-ggtern(data=df_4modalities_celltype_all.genes, aes(x=Striatum_tmCG, y=Striatum_hmCG, z=Striatum_CG)) +
  geom_mask() +
  geom_point(aes(colour=log2(Striatum_rna)), size=0.2, alpha=0.5) +
  sc + 
   theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +
  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +
  ggtitle('Striatum')

Oligo<-ggtern(data=df_4modalities_celltype_all.genes, aes(x=Oligo_tmCG, y=Oligo_hmCG, z=Oligo_CG)) +
  geom_mask() +
  geom_point(aes(colour=log2(Oligo_rna)), size=0.2, alpha=0.5) +
  sc + 
   theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +
  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +
  ggtitle('Oligo')

MG<-ggtern(data=df_4modalities_celltype_all.genes, aes(x=MG_tmCG, y=MG_hmCG, z=MG_CG)) +
  geom_mask() +
  geom_point(aes(colour=log2(MG_rna)), size=0.2, alpha=0.5) +
  sc + 
   theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +
  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +
  ggtitle('MG')

Astro<-ggtern(data=df_4modalities_celltype_all.genes, aes(x=Astro_tmCG, y=Astro_hmCG, z=Astro_CG)) +
  geom_mask() +
  geom_point(aes(colour=log2(Astro_rna)), size=0.2, alpha=0.5) +
  sc + 
   theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +
  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +
  ggtitle('Astro')

###
plot_list=list()
plot_list[["L2.3"]] <- L2.3
plot_list[["L4"]] <- L4
plot_list[["L5"]] <- L5
plot_list[["L6"]] <- L6
plot_list[["DL"]] <- DL
plot_list[["Pv"]] <- Pv
plot_list[["Sst"]] <- Sst
plot_list[["Vip.Ndnf"]] <- Vip.Ndnf
plot_list[["Striatum"]] <- Striatum
plot_list[["Astro"]] <- Astro
plot_list[["Oligo"]] <- Oligo
plot_list[["MG"]] <- MG

library(gridExtra)

options(repr.plot.width = 15,repr.plot.height = 9)
#g <- grid.arrange(grobs=plot_list,nrow = 1)
ternary_plot <- ggtern::grid.arrange(grobs=plot_list,nrow = 3)

ggsave("Fig6_12ternary.plots.RNA_log2RNA.-6to6_09.13.2022.pdf",ternary_plot,width=15,height=9)

options(repr.plot.width = 15,repr.plot.height = 4)

myPalette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(-8, 8))


L2.3<-ggtern(data=df_4modalities_celltype_all.genes, aes(x=L2.3_tmCG, y=L2.3_hmCG, z=L2.3_CG)) +
  geom_mask() +
  geom_point(aes(colour=log2(L2.3_rna)), size=0.2, alpha=0.5) +
  sc + 
   theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +
  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +
  ggtitle('Ex-L2.3')

L4<-ggtern(data=df_4modalities_celltype_all.genes, aes(x=L4_tmCG, y=L4_hmCG, z=L4_CG)) +
  geom_mask() +
  geom_point(aes(colour=log2(L4_rna)), size=0.2, alpha=0.5) +
  sc + 
   theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +
  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +
  ggtitle('Ex-L4')

L5<-ggtern(data=df_4modalities_celltype_all.genes, aes(x=L4.5_tmCG, y=L4.5_hmCG, z=L4.5_CG)) +
  geom_mask() +
  geom_point(aes(colour=log2(L4.5_rna)), size=0.2, alpha=0.5) +
  sc + 
   theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +
  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +
  ggtitle('Ex-L5')

L6<-ggtern(data=df_4modalities_celltype_all.genes, aes(x=L6_tmCG, y=L6_hmCG, z=L6_CG)) +
  geom_mask() +
  geom_point(aes(colour=log2(L6_rna)), size=0.2, alpha=0.5) +
  sc + 
   theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +
  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +
  ggtitle('Ex-L6')

DL<-ggtern(data=df_4modalities_celltype_all.genes, aes(x=DL_tmCG, y=DL_hmCG, z=DL_CG)) +
  geom_mask() +
  geom_point(aes(colour=log2(DL_rna)), size=0.2, alpha=0.5) +
  sc + 
   theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +
  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +
  ggtitle('Ex-DL')

Sst<-ggtern(data=df_4modalities_celltype_all.genes, aes(x=Sst_tmCG, y=Sst_hmCG, z=Sst_CG)) +
  geom_mask() +
  geom_point(aes(colour=log2(Sst_rna)), size=0.2, alpha=0.5) +
  sc + 
   theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +
  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +
  ggtitle('Inh-Sst')

Pv<-ggtern(data=df_4modalities_celltype_all.genes, aes(x=Pv_tmCG, y=Pv_hmCG, z=Pv_CG)) +
  geom_mask() +
  geom_point(aes(colour=log2(Pv_rna)), size=0.2, alpha=0.5) +
  sc + 
   theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +
  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +
  ggtitle('Inh-Pv')

Vip.Ndnf<-ggtern(data=df_4modalities_celltype_all.genes, aes(x=Vip.Ndnf_tmCG, y=Vip.Ndnf_hmCG, z=Vip.Ndnf_CG)) +
  geom_mask() +
  geom_point(aes(colour=log2(Vip.Ndnf_rna)), size=0.2, alpha=0.5) +
  sc + 
   theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +
  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +
  ggtitle('Inh-Vip.Ndnf')

Striatum<-ggtern(data=df_4modalities_celltype_all.genes, aes(x=Striatum_tmCG, y=Striatum_hmCG, z=Striatum_CG)) +
  geom_mask() +
  geom_point(aes(colour=log2(Striatum_rna)), size=0.2, alpha=0.5) +
  sc + 
   theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +
  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +
  ggtitle('Striatum')

Oligo<-ggtern(data=df_4modalities_celltype_all.genes, aes(x=Oligo_tmCG, y=Oligo_hmCG, z=Oligo_CG)) +
  geom_mask() +
  geom_point(aes(colour=log2(Oligo_rna)), size=0.2, alpha=0.5) +
  sc + 
   theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +
  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +
  ggtitle('Oligo')

MG<-ggtern(data=df_4modalities_celltype_all.genes, aes(x=MG_tmCG, y=MG_hmCG, z=MG_CG)) +
  geom_mask() +
  geom_point(aes(colour=log2(MG_rna)), size=0.2, alpha=0.5) +
  sc + 
   theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +
  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +
  ggtitle('MG')

Astro<-ggtern(data=df_4modalities_celltype_all.genes, aes(x=Astro_tmCG, y=Astro_hmCG, z=Astro_CG)) +
  geom_mask() +
  geom_point(aes(colour=log2(Astro_rna)), size=0.2, alpha=0.5) +
  sc + 
   theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +
  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +
  ggtitle('Astro')

###
plot_list=list()
plot_list[["L2.3"]] <- L2.3
plot_list[["L4"]] <- L4
plot_list[["L5"]] <- L5
plot_list[["L6"]] <- L6
plot_list[["DL"]] <- DL
plot_list[["Pv"]] <- Pv
plot_list[["Sst"]] <- Sst
plot_list[["Vip.Ndnf"]] <- Vip.Ndnf
plot_list[["Striatum"]] <- Striatum
plot_list[["Astro"]] <- Astro
plot_list[["Oligo"]] <- Oligo
plot_list[["MG"]] <- MG

library(gridExtra)

options(repr.plot.width = 15,repr.plot.height = 9)
#g <- grid.arrange(grobs=plot_list,nrow = 1)
ternary_plot <- grid.arrange(grobs=plot_list,nrow = 3)

ggsave("Fig4_12ternary.plots_log2RNA.-8to8_12.28.2021.pdf",ternary_plot,width=15,height=9)

library(gghighlight)

# https://bjnnowak.netlify.app/2021/07/26/r-plotting-soil-textures-example-of-water-storage-capacity/
df_L6_CG.low_hmC.hi_rna.hi_genes <- df_4modalities_celltype_all.genes %>% filter(L6_CG<0.2 & L6_hmCG>0.3 & L6_rna>2)

df_L6_CG.hi_hmC.low_rna.hi_genes <- df_4modalities_celltype_all.genes %>% filter(L6_CG>0.5 & L6_hmCG<0.3 & L6_rna>2)

myPalette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(-8, 8))


ggtern(data=df_4modalities_celltype_all.genes, aes(x=L6_tmCG, y=L6_hmCG, z=L6_CG)) +
#ggtern(data=df_4modalities_celltype_all.genes %>% filter(L6_rna>10), aes(x=L6_tmCG, y=L6_hmCG, z=L6_CG)) +

  geom_mask() +
  geom_point(aes(colour=log2(L6_rna)), size=0.5, alpha=0.8) +
  sc + 
#  scale_color_continuous(
#    low="dodgerblue",
#    high="tomato"
#  )+

#  gghighlight(L6_rna>4) + 
#gghighlight(L6_CG >= 0.5) +
#gghighlight(L6_hmCG/L6_tmCG > 1) +
#gghighlight(L6_CG<0.2 & L6_hmCG>0.3 & L6_rna>2) + 

geom_point(data=df_4modalities_celltype_all.genes[c("Tle4","Arpp21","Erbb4"), ], colour="red", size=2) +

#geom_point(data=df_L6_CG.low_hmC.hi_rna.hi_genes, 
#             aes(x=lifeExp,y=gdpPercap), 
#             color='red',
#             size=2) +

#  geom_point(size=3, aes(shape = DNA_methylation,
#                         fill = Germline_imprints)) +
#scales
#  scale_shape_manual(values = c(21, 24)) +
#  scale_fill_manual(values = c("blue","yellow","yellow")) + # Two 'NA's were not filled with color.
#  scale_size_continuous(range = c(2.5, 7.5)) +
  
#  stat_density_tern(geom='polygon', aes(fill=..level.., alpha=..level..)) + 
#  stat_density_tern(geom='polygon', aes(fill=..level..), base='identity', bins=100) +

#  geom_segment(data=lines, aes(x, y, z, xend = xend, yend = yend, zend = zend), color = 'gray20', size = 0.5) +

  theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +

  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +

#  guides(shape = guide_legend(override.aes=list(size=3))) +

#  geom_Rline(Rintercept=.5, colour='gray') +

#  labs(shape = 'Parental origin', 
#       fill = 'Germline imprints') +

  ggtitle('Ex-L6')



# https://bjnnowak.netlify.app/2021/07/26/r-plotting-soil-textures-example-of-water-storage-capacity/
myPalette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(0, 6))


L6_mCH <- ggtern(data=df_5modalities_celltype_all.genes, aes(x=L6_tmCG, y=L6_hmCG, z=L6_CG)) +
#ggtern(data=df_4modalities_celltype_all.genes %>% filter(L6_rna>10), aes(x=L6_tmCG, y=L6_hmCG, z=L6_CG)) +

  geom_mask() +
  geom_point(aes(colour=L6_mCH*100), size=0.5, alpha=0.8) +
  sc + 
#  scale_color_continuous(
#    low="dodgerblue",
#    high="tomato"
#  )+

#  gghighlight(L6_rna>4) + 
#gghighlight(L6_CG >= 0.5) +
#gghighlight(L6_hmCG/L6_tmCG > 1) +
#gghighlight(L6_CG<0.2 & L6_hmCG>0.3 & L6_rna>2) + 

geom_point(data=df_4modalities_celltype_all.genes[c("Tle4","Arpp21","Erbb4"), ], colour="red", size=2) +

#  geom_segment(data=lines, aes(x, y, z, xend = xend, yend = yend, zend = zend), color = 'gray20', size = 0.5) +

  theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +

  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +

#  guides(shape = guide_legend(override.aes=list(size=3))) +

#  geom_Rline(Rintercept=.5, colour='gray') +

#  labs(shape = 'Parental origin', 
#       fill = 'Germline imprints') +

  ggtitle('Ex-L6')

Pv_mCH <- ggtern(data=df_5modalities_celltype_all.genes, aes(x=Pv_tmCG, y=Pv_hmCG, z=Pv_CG)) +
#ggtern(data=df_4modalities_celltype_all.genes %>% filter(L6_rna>10), aes(x=L6_tmCG, y=L6_hmCG, z=L6_CG)) +

  geom_mask() +
  geom_point(aes(colour=Pv_mCH*100), size=0.5, alpha=0.8) +
  sc + 
#  scale_color_continuous(
#    low="dodgerblue",
#    high="tomato"
#  )+

#  gghighlight(L6_rna>4) + 
#gghighlight(L6_CG >= 0.5) +
#gghighlight(L6_hmCG/L6_tmCG > 1) +
#gghighlight(L6_CG<0.2 & L6_hmCG>0.3 & L6_rna>2) + 

geom_point(data=df_4modalities_celltype_all.genes[c("Tle4","Arpp21","Erbb4"), ], colour="red", size=2) +

#  geom_segment(data=lines, aes(x, y, z, xend = xend, yend = yend, zend = zend), color = 'gray20', size = 0.5) +

  theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +

  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +

#  guides(shape = guide_legend(override.aes=list(size=3))) +

#  geom_Rline(Rintercept=.5, colour='gray') +

#  labs(shape = 'Parental origin', 
#       fill = 'Germline imprints') +

  ggtitle('Inh-Pv')

plot_list=list()
#plot_list[["L2.3"]] <- L2.3
#plot_list[["L4"]] <- L4
#plot_list[["L5"]] <- L5
plot_list[["L6"]] <- L6_mCH
#plot_list[["DL"]] <- DL
plot_list[["Pv"]] <- Pv_mCH
#plot_list[["Sst"]] <- Sst
#plot_list[["Vip.Ndnf"]] <- Vip.Ndnf
#plot_list[["Striatum"]] <- Striatum
#plot_list[["Astro"]] <- Astro
#plot_list[["Oligo"]] <- Oligo
#plot_list[["MG"]] <- MG

library(gridExtra)

options(repr.plot.width = 15,repr.plot.height = 9)


#g <- grid.arrange(grobs=plot_list,nrow = 1)
ternary_plot <- ggtern::grid.arrange(grobs=plot_list,nrow = 1)

# https://bjnnowak.netlify.app/2021/07/26/r-plotting-soil-textures-example-of-water-storage-capacity/

myPalette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
#myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(-6, 6))

L6 <- ggtern(data=df_4modalities_celltype_all.genes, aes(x=L6_tmCG, y=L6_hmCG, z=L6_CG)) +
#ggtern(data=df_4modalities_celltype_all.genes %>% filter(L6_rna>10), aes(x=L6_tmCG, y=L6_hmCG, z=L6_CG)) +

  geom_mask() +
  geom_point(aes(colour=log2(L6_rna)), size=0.5, alpha=0.5) + sc +
#  scale_color_continuous(
#    low="dodgerblue",
#    high="tomato"
#  )+

  gghighlight(L6_rna>1) + 
#gghighlight(L6_CG >= 0.5) +
#gghighlight(L6_hmCG/L6_tmCG > 1) +
#gghighlight(L6_CG<0.2 & L6_hmCG>0.3 & L6_rna>2) + 

geom_point(data=df_4modalities_celltype_all.genes[c("Tle4"), ], colour="red", size=2) +
geom_point(data=df_4modalities_celltype_all.genes[c("Adarb2"), ], colour="green", size=2) +

  theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +

  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +

  ggtitle('Ex-L6')


Vip.Ndnf <- ggtern(data=df_4modalities_celltype_all.genes, aes(x=Vip.Ndnf_tmCG, y=Vip.Ndnf_hmCG, z=Vip.Ndnf_CG)) +
#ggtern(data=df_4modalities_celltype_all.genes %>% filter(L6_rna>10), aes(x=L6_tmCG, y=L6_hmCG, z=L6_CG)) +

  geom_mask() +
  geom_point(aes(colour=log2(Vip.Ndnf_rna)), size=0.5, alpha=0.5) + sc +
#  scale_color_continuous(
#    low="dodgerblue",
#    high="tomato"
#  )+

  gghighlight(Vip.Ndnf_rna>1) + 
#gghighlight(L6_CG >= 0.5) +
#gghighlight(L6_hmCG/L6_tmCG > 1) +
#gghighlight(L6_CG<0.2 & L6_hmCG>0.3 & L6_rna>2) + 

geom_point(data=df_4modalities_celltype_all.genes[c("Tle4"), ], colour="red", size=2) +
geom_point(data=df_4modalities_celltype_all.genes[c("Adarb2"), ], colour="green", size=2) +

  theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +

  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +

  ggtitle('Inh-Vip.Ndnf')

plot_list=list()
#plot_list[["L2.3"]] <- L2.3
#plot_list[["L4"]] <- L4
#plot_list[["L5"]] <- L5
plot_list[["L6"]] <- L6
#plot_list[["DL"]] <- DL
#plot_list[["Pv"]] <- Pv
#plot_list[["Sst"]] <- Sst
plot_list[["Vip.Ndnf"]] <- Vip.Ndnf
#plot_list[["Striatum"]] <- Striatum
#plot_list[["Astro"]] <- Astro
#plot_list[["Oligo"]] <- Oligo
#plot_list[["MG"]] <- MG

library(gridExtra)

options(repr.plot.width = 15,repr.plot.height = 9)
#g <- grid.arrange(grobs=plot_list,nrow = 1)

#plot_grid(L6, Vip.Ndnf, nrow=1)
ggtern::grid.arrange(grobs=plot_list,nrow = 1)




df_4modalities_celltype_all.genes %>% arrange(desc(L6_rna)) %>% filter(log2(L6_rna)>4) # sort L6_rna

cluster.order
dim(hmC_matrix)
head(hmC_matrix)
#table(pbmc.snmc3@meta.data$res.name)
table(pbmc.snmc3@meta.data$res.name2)
table(pbmc.snmc3@meta.data$res.name3)
head(pbmc.snmc3@meta.data)

mat_hmCG_celltype <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name2 == x]
    return(rowMeans(hmC_matrix[gene.order.by.tmCG,cells.in]))
})

mat_tmC_celltype <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name2 == x]
    return(rowMeans(tmC_matrix[gene.order.by.tmCG,cells.in]))
})

mat_mCG_celltype <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name2 == x]
    return(rowMeans(mC_matrix[gene.order.by.tmCG,cells.in]))
})

mat_mCH_celltype <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name2 == x]
    return(rowMeans(mCH_matrix[gene.order.by.tmCG,cells.in]))
})

mat_ratio_celltype <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name2 == x]
    return(rowMeans(ratio_matrix[gene.order.by.tmCG,cells.in]))
})

mat_rna_celltype <- sapply(cluster.order,function(x){
    cells.in <- rownames(snRNA3@meta.data)[snRNA3@meta.data$res.name.JE3 == x]
    mat <- expm1(snRNA3@assays$RNA@data)
    return(Matrix::rowMeans(mat[gene.order.by.tmCG,cells.in]))
})



mat_CG_celltype <- 1 - mat_mCG_celltype

dim(mat_CG_celltype)
head(mat_CG_celltype)

head(mat_hmCG_celltype)
dim(mat_hmCG_celltype)

library(RColorBrewer)

#install.packages("pheatmap")
library(pheatmap)

options(repr.plot.width=4,repr.plot.height=4)


h1 <- pheatmap::pheatmap(mat_CG_celltype,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="CG",color=rev(colorRampPalette(brewer.pal(n = 9,name="BuPu"))(50)))

h2 <- pheatmap::pheatmap(mat_tmC_celltype,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="tmCG",color=rev(colorRampPalette(brewer.pal(n = 9,name="BuPu"))(50)))
#h2 <- pheatmap::pheatmap(mat_mCG_celltype,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
#                   scale="row",clustering_method = "ward.D2",main="BSmCG",color=rev(colorRampPalette(brewer.pal(n = 9,name="BuPu"))(50)))

h3 <- pheatmap::pheatmap(mat_hmCG_celltype,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="hmCG",color=rev(colorRampPalette(brewer.pal(n = 9,name="BuPu"))(50)))

h4 <- pheatmap::pheatmap(mat_ratio_celltype,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="hmCG/BSmCG",color=colorRampPalette(brewer.pal(n = 9,name="BuPu"))(50))

h5 <- pheatmap::pheatmap(mat_mCH_celltype,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="BSmCH",color=rev(colorRampPalette(brewer.pal(n = 9,name="BuPu"))(50)))

h6 <- pheatmap::pheatmap(mat_rna_celltype,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="RNA")

# combine all the plots; 80% cells
plot_list=list()
plot_list[["CG"]] <- h1[[4]]
plot_list[["tmCG"]] <- h2[[4]]
plot_list[["hmCG"]] <- h3[[4]]
#plot_list[["hmCG/BSmCG"]] <- h4[[4]]
plot_list[["BSmCH"]] <- h5[[4]]
plot_list[["RNA"]] <- h6[[4]]

library(gridExtra)

options(repr.plot.width = 16,repr.plot.height = 4)
#grid.arrange(grobs=plot_list,ncol = 4)
g2.by.tmC <- grid.arrange(grobs=plot_list,ncol = 6)

# Peng's analysis:
# combine all the plots
plot_list=list()
plot_list[["tmCG"]] <- h1[[4]]
plot_list[["BSmCG"]] <- h2[[4]]
plot_list[["hmCG"]] <- h3[[4]]
plot_list[["hmCG/BSmCG"]] <- h4[[4]]
plot_list[["BSmCH"]] <- h5[[4]]
plot_list[["RNA"]] <- h6[[4]]

library(gridExtra)

options(repr.plot.width = 12,repr.plot.height = 3)
#grid.arrange(grobs=plot_list,ncol = 4)
g2 <- grid.arrange(grobs=plot_list,ncol = 6)

ggsave("Fig6_all.cell.types_6heatmaps_order_by_tmCG_09.16.2022_1920genes.pdf",g2.by.tmC,width=16,height=4)

gene.order <- gene.order.by.hmCG

mat_hmCG_celltype <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name2 == x]
    return(rowMeans(hmC_matrix[gene.order,cells.in]))
})

mat_tmC_celltype <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name2 == x]
    return(rowMeans(tmC_matrix[gene.order,cells.in]))
})

mat_mCG_celltype <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name2 == x]
    return(rowMeans(mC_matrix[gene.order,cells.in]))
})

mat_mCH_celltype <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name2 == x]
    return(rowMeans(mCH_matrix[gene.order,cells.in]))
})

mat_ratio_celltype <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name2 == x]
    return(rowMeans(ratio_matrix[gene.order,cells.in]))
})

mat_rna_celltype <- sapply(cluster.order,function(x){
    cells.in <- rownames(snRNA3@meta.data)[snRNA3@meta.data$res.name.JE3 == x]
    mat <- expm1(snRNA3@assays$RNA@data)
    return(Matrix::rowMeans(mat[gene.order,cells.in]))
})

mat_CG_celltype <- 1 - mat_mCG_celltype

dim(mat_CG_celltype)
head(mat_CG_celltype)

options(repr.plot.width=4,repr.plot.height=4)


h1 <- pheatmap::pheatmap(mat_CG_celltype,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="CG",color=rev(colorRampPalette(brewer.pal(n = 9,name="BuPu"))(50)))

h2 <- pheatmap::pheatmap(mat_tmC_celltype,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="tmCG",color=rev(colorRampPalette(brewer.pal(n = 9,name="BuPu"))(50)))
#h2 <- pheatmap::pheatmap(mat_mCG_celltype,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
#                   scale="row",clustering_method = "ward.D2",main="BSmCG",color=rev(colorRampPalette(brewer.pal(n = 9,name="BuPu"))(50)))

h3 <- pheatmap::pheatmap(mat_hmCG_celltype,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="hmCG",color=rev(colorRampPalette(brewer.pal(n = 9,name="BuPu"))(50)))

h4 <- pheatmap::pheatmap(mat_ratio_celltype,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="hmCG/BSmCG",color=colorRampPalette(brewer.pal(n = 9,name="BuPu"))(50))

h5 <- pheatmap::pheatmap(mat_mCH_celltype,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="BSmCH",color=rev(colorRampPalette(brewer.pal(n = 9,name="BuPu"))(50)))

h6 <- pheatmap::pheatmap(mat_rna_celltype,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="RNA")

# combine all the plots; 80% cells
plot_list=list()
plot_list[["CG"]] <- h1[[4]]
plot_list[["tmCG"]] <- h2[[4]]
plot_list[["hmCG"]] <- h3[[4]]
#plot_list[["hmCG/BSmCG"]] <- h4[[4]]
plot_list[["BSmCH"]] <- h5[[4]]
plot_list[["RNA"]] <- h6[[4]]

library(gridExtra)

options(repr.plot.width = 16,repr.plot.height = 4)
#grid.arrange(grobs=plot_list,ncol = 4)
g2.by.hmC <- grid.arrange(grobs=plot_list,ncol = 6)

ggsave("Fig6_all.cell.types_6heatmaps_order_by_hmCG_09.16.2022_1087genes.pdf",g2.by.hmC,width=16,height=4)

Idents(snRNA3) <- "res.name.JE3"

table(Idents(snRNA3))

rna.gene <- FindAllMarkers(snRNA3,only.pos = T)

length(unique(rna.gene$gene)) # Peng: 2771; Hao: 3317 genes

cluster.order <- c("L2/3","L4","L4/5","L6","DL","Vip/Ndnf","Pv","Sst","Striatum","Astro","Oligo","MG")
rna.gene$cluster <- factor(rna.gene$cluster,levels = cluster.order)
rna.gene2 <- rna.gene %>% arrange(cluster,p_val_adj) 
#rna.gene2 <- rna.gene %>% arrange(cluster,avg_log2FC) 
#gene.order <- rna.gene2$gene

dim(rna.gene)
head(rna.gene)
tail(rna.gene)
rna.gene.pval01 <- subset(rna.gene, p_val_adj < 0.1)
rna.gene2.pval01 <- rna.gene.pval01 %>% arrange(cluster,p_val_adj) 
dim(rna.gene2.pval01)
#head(rna.gene2.pval01)
gene.order <- rna.gene2.pval01$gene

head(subset(rna.gene2, cluster=="MG"))
dim(subset(rna.gene2, cluster=="MG"))
dim(subset(rna.gene2, cluster=="MG" & p_val_adj < 0.1))
dim(subset(rna.gene2.pval01, cluster=="MG"))

length(rownames(pbmc.snmc3@assays$tmCG@data))

gene.order.by.rna <- intersect(gene.order,rownames(pbmc.snmc3@assays$tmCG@data))
length(gene.order.by.rna) 
# Peng: 1544; 
# Hao (90% cells, 5 CpGs): 1501 (pval<0.1); 1295 (pval<0.05); 
# Hao (80% cells, 5 CpGs): 1653 (pval<0.1)
# Hao (80% cells, 1 CpGs): 1873 (pval<0.1)
# Hao (80% cells, 2 CpGs): 1825 (pval<0.1)

# Hao (75% cells, 3 CpGs): 1805 (pval<0.1)

# old analysis: 
gene.order.by.rna <- intersect(gene.order,rownames(pbmc.snmc3@assays$tmCG@data))
length(gene.order.by.rna) # Peng: 1544; Hao (90% cells): 1295 (pval<0.05); 1501 (pval<0.1); (80% cells): 1653 (pval<0.1)

head(snRNA3@meta.data)
head(gene.order.by.rna)
head(snRNA3@meta.data$res.name.JE3)
#snRNA3@assays$RNA@data

gene.order <- gene.order.by.rna

mat_hmCG_celltype <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name2 == x]
    return(rowMeans(hmC_matrix[gene.order,cells.in]))
})

mat_tmC_celltype <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name2 == x]
    return(rowMeans(tmC_matrix[gene.order,cells.in]))
})

mat_mCG_celltype <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name2 == x]
    return(rowMeans(mC_matrix[gene.order,cells.in]))
})

mat_mCH_celltype <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name2 == x]
    return(rowMeans(mCH_matrix[gene.order,cells.in]))
})

mat_ratio_celltype <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name2 == x]
    return(rowMeans(ratio_matrix[gene.order,cells.in]))
})

mat_rna_celltype <- sapply(cluster.order,function(x){
    cells.in <- rownames(snRNA3@meta.data)[snRNA3@meta.data$res.name.JE3 == x]
    mat <- expm1(snRNA3@assays$RNA@data)
    return(Matrix::rowMeans(mat[gene.order,cells.in]))
})


mat_CG_celltype <- 1 - mat_mCG_celltype

dim(mat_CG_celltype)
dim(mat_tmC_celltype)
head(mat_CG_celltype)

#scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))
head(mat_tmC_celltype, n=10)

options(repr.plot.width=3,repr.plot.height=3)

h1 <- pheatmap::pheatmap(mat_CG_celltype,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="CG",color=rev(colorRampPalette(brewer.pal(n = 9,name="BuPu"))(50)))

h2 <- pheatmap::pheatmap(mat_tmC_celltype,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="tmCG",color=rev(colorRampPalette(brewer.pal(n = 9,name="BuPu"))(50)))
#h2 <- pheatmap::pheatmap(mat_mCG_celltype,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
#                   scale="row",clustering_method = "ward.D2",main="BSmCG",color=rev(colorRampPalette(brewer.pal(n = 9,name="BuPu"))(50)))

h3 <- pheatmap::pheatmap(mat_hmCG_celltype,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="hmCG",color=rev(colorRampPalette(brewer.pal(n = 9,name="BuPu"))(50)))

h4 <- pheatmap::pheatmap(mat_ratio_celltype,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="hmCG/BSmCG",color=colorRampPalette(brewer.pal(n = 9,name="BuPu"))(50))

h5 <- pheatmap::pheatmap(mat_mCH_celltype,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="BSmCH",color=rev(colorRampPalette(brewer.pal(n = 9,name="BuPu"))(50)))

h6 <- pheatmap::pheatmap(mat_rna_celltype,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="RNA",colours = rev(colorRampPalette(brewer.pal(n = 11, name = "RdBu"))(50)))

h6 <- pheatmap::pheatmap(mat_rna_celltype,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="RNA",colours = colorRampPalette(brewer.pal(9,"RdBu"))(50))

# Hao's new analysis -  1825 genes 
# combine all the plots
plot_list=list()
plot_list[["CG"]] <- h1[[4]]
plot_list[["tmCG"]] <- h2[[4]]
plot_list[["hmCG"]] <- h3[[4]]
#plot_list[["hmCG/BSmCG"]] <- h4[[4]]
plot_list[["BSmCH"]] <- h5[[4]]
plot_list[["RNA"]] <- h6[[4]]

library(gridExtra)

options(repr.plot.width = 16,repr.plot.height = 4)
#grid.arrange(grobs=plot_list,ncol = 4)
g2.by.rna <- grid.arrange(grobs=plot_list,ncol = 6)

# Hao's old analysis - 1653 genes 
# combine all the plots
plot_list=list()
plot_list[["CG"]] <- h1[[4]]
plot_list[["tmCG"]] <- h2[[4]]
plot_list[["hmCG"]] <- h3[[4]]
plot_list[["hmCG/BSmCG"]] <- h4[[4]]
plot_list[["BSmCH"]] <- h5[[4]]
plot_list[["RNA"]] <- h6[[4]]

library(gridExtra)

options(repr.plot.width = 12,repr.plot.height = 4)
#grid.arrange(grobs=plot_list,ncol = 4)
g2 <- grid.arrange(grobs=plot_list,ncol = 6)

# Peng's old plots:
# combine all the plots
plot_list=list()
plot_list[["tmCG"]] <- h1[[4]]
plot_list[["BSmCG"]] <- h2[[4]]
plot_list[["hmCG"]] <- h3[[4]]
plot_list[["hmCG/BSmCG"]] <- h4[[4]]
plot_list[["BSmCH"]] <- h5[[4]]
plot_list[["RNA"]] <- h6[[4]]

library(gridExtra)

options(repr.plot.width = 12,repr.plot.height = 3)
#grid.arrange(grobs=plot_list,ncol = 4)
g2 <- grid.arrange(grobs=plot_list,ncol = 6)

ggsave("Fig4_6heatmaps_order_by_RNA_03.17.2022_1873genes.pdf",g2,width=12,height=4)

ggsave("Fig6_all.cell.types_6heatmaps_order_by_RNA_09.16.2022_1825genes.pdf",g2.by.rna,width=16,height=4)

cluster.order <- c("L2/3","L4","L4/5","L6","DL","Vip/Ndnf","Pv","Sst","Striatum","Astro","Oligo","MG")
length(gene.order.by.hmCG) # 1544

gene.order2 <- gene.order.by.hmCG

mat_hmCG_celltype_norm <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name2 == x]
    m <- hmC_matrix
    mat <- sweep(m,2,colSums(m),`/`)
    return(rowMeans(mat[gene.order2,cells.in]))
})

mat_tmC_celltype_norm <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name2 == x]
    m <- tmC_matrix
    mat <- sweep(m,2,colSums(m),`/`)
    return(rowMeans(mat[gene.order2,cells.in]))
    #return(rowMeans(tmC_matrix[gene.order2.rna,cells.in]))
})

mat_mCG_celltype_norm <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name2 == x]
    m <- mC_matrix
    mat <- sweep(m,2,colSums(m),`/`)
    return(rowMeans(mat[gene.order2,cells.in]))
    #return(rowMeans(mC_matrix[gene.order2.rna,cells.in]))
})

mat_mCH_celltype_norm <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name2 == x]
    m <- mCH_matrix
    mat <- sweep(m,2,colSums(m),`/`)
    return(rowMeans(mat[gene.order2,cells.in]))
    #return(rowMeans(mCH_matrix[gene.order2.rna,cells.in]))
})

mat_ratio_celltype_norm <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name2 == x]
    m <- ratio_matrix
    mat <- sweep(m,2,colSums(m),`/`)
    return(rowMeans(mat[gene.order2,cells.in]))
    #return(rowMeans(ratio_matrix[gene.order2.rna,cells.in]))
})

mat_rna_celltype_norm <- sapply(cluster.order,function(x){
    cells.in <- rownames(snRNA3@meta.data)[snRNA3@meta.data$res.name.JE3 == x]
    mat <- expm1(snRNA3@assays$RNA@data)
    return(Matrix::rowMeans(mat[gene.order2,cells.in]))
})

options(repr.plot.width=3,repr.plot.height=3)

h1 <- pheatmap::pheatmap(mat_tmC_celltype_norm,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="tmCG",color=rev(colorRampPalette(brewer.pal(n = 9,name="BuPu"))(50)))

h2 <- pheatmap::pheatmap(mat_mCG_celltype_norm,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="BSmCG",color=rev(colorRampPalette(brewer.pal(n = 9,name="BuPu"))(50)))

h3 <- pheatmap::pheatmap(mat_hmCG_celltype_norm,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="hmCG",color=rev(colorRampPalette(brewer.pal(n = 9,name="BuPu"))(50)))

h4 <- pheatmap::pheatmap(mat_ratio_celltype_norm,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="hmCG/BSmCG",color=colorRampPalette(brewer.pal(n = 9,name="BuPu"))(50))

h5 <- pheatmap::pheatmap(mat_mCH_celltype_norm,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="BSmCH",color=rev(colorRampPalette(brewer.pal(n = 9,name="BuPu"))(50)))

h6 <- pheatmap::pheatmap(mat_rna_celltype_norm,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="RNA",colours = rev(colorRampPalette(brewer.pal(n = 11, name = "RdBu"))(50)))

h6 <- pheatmap::pheatmap(mat_rna_celltype_norm,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="RNA",colours = colorRampPalette(c("blue","white","red"))(50))

# combine all the plots
plot_list=list()
plot_list[["tmCG"]] <- h1[[4]]
plot_list[["BSmCG"]] <- h2[[4]]
plot_list[["hmCG"]] <- h3[[4]]
#plot_list[["hmCG/BSmCG"]] <- h4[[4]]
plot_list[["BSmCH"]] <- h5[[4]]
plot_list[["RNA"]] <- h6[[4]]

library(gridExtra)

options(repr.plot.width = 16,repr.plot.height = 4)
#grid.arrange(grobs=plot_list,ncol = 4)
g2.by.hmC.col.scaled <- grid.arrange(grobs=plot_list,ncol = 6)

ggsave("Fig6_all.cell.types_6heatmaps_order_by_hmCG.col.scaled_09.16.2022_1087genes.pdf",g2.by.hmC.col.scaled,width=16,height=4)

plot_box3_value_rmIn_nonNeu <- function(genes.in=c("Tle4","Gad2","Gad1","Rorb"),mat1=CG_genebody,yLab="GBmCG",nRow=4){
    mat2 <- mat1[intersect(genes.in,rownames(mat1)),]
    mat.df <- as.data.frame(cbind(t(mat2),as.character(pbmc.snmc3$res.name2)))
    names(mat.df)[ncol(mat.df)] <- "cluster"
    df.plot <- reshape2::melt(mat.df,id=c("cluster"))
    df.plot$value <- as.numeric(df.plot$value)
    df.plot %>% ungroup %>% select(variable,value) %>% group_by(variable) %>% 
      arrange(desc(value)) %>% dplyr::slice(1) %>% dplyr::rename(max = value) -> df.max.value
    cluster.order <- c('L2/3','L4','L4/5','L6','DL','Vip/Ndnf','Pv','Sst','Striatum','Oligo','Astro','MG')
    df.plot2 <- df.plot %>% ungroup %>% group_by(variable) %>% 
                  left_join(df.max.value,by=c("variable")) %>% ungroup %>%
                  dplyr::mutate(value2 = value/max) %>% dplyr::filter(cluster %in% cluster.order) %>% droplevels
    df.plot2$cluster <- factor(df.plot2$cluster,levels=cluster.order)
    library(ggplot2)
    library(cowplot)
    p1 <- ggplot(df.plot2,aes(cluster,value)) + geom_boxplot(width=.6,outlier.shape = NA) + 
          ggrastr::geom_point_rast(position = "jitter",col="grey",width = 0.25,size=0.25) + 
          facet_wrap(~variable,nrow = nRow) +
          stat_summary(fun.y="mean", geom="point", shape=23, size=2,col="red") +
          theme_cowplot() + theme(axis.text.x = element_text(angle=45,hjust=1)) + ylab("Meth rate") + ggtitle(yLab)
    return(p1)
}

plot_vln2 <- function(genes.in=c("Tle4","Gad2","Gad1","Rorb"),mat1=as.matrix(snRNA3@assays$RNA@data),yLab="RNA"){
    mat2 <- mat1[intersect(genes.in,rownames(mat1)),]
    mat.df <- as.data.frame(cbind(t(mat2),as.character(snRNA3@meta.data$res.name.JE3)))
    names(mat.df)[ncol(mat.df)] <- "cluster"
    df.plot <- reshape2::melt(mat.df,id=c("cluster"))
    df.plot$value <- as.numeric(df.plot$value)
    df.plot %>% ungroup %>% select(variable,value) %>% group_by(variable) %>% 
      arrange(desc(value)) %>% dplyr::slice(1) %>% dplyr::rename(max = value) -> df.max.value
    df.plot2 <- df.plot %>% ungroup %>% group_by(variable) %>% 
                  left_join(df.max.value,by=c("variable")) %>% ungroup %>%
                  dplyr::mutate(value2 = value/max)
    cluster.order <- c('L2/3','L4','L4/5','L6','DL','Vip/Ndnf','Pv','Sst','Striatum','Oligo','Astro','MG')
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


genes.in <- c("Satb2","Slc17a7","Gad1","Gad2", ## 
             "Pvalb","Sst","Reln","Ndnf","Vip","Rorb","Foxp2","Tle4",
             "Cux1","Cux2","Slc6a1","Sulf1","Prox1","Grik3","Slc1a2","Csf1r",
             "Mbp","Cspg4","Nxn","Lcp2","Mef2c","Pdgfra","Ppp1r1b", "Meis2")

pp2 <- plot_box3_value_rmIn_nonNeu(genes.in=gene.list,mat1=hmC_matrix,yLab="hmCG", nRow=4)
pp3 <- plot_box3_value_rmIn_nonNeu(genes.in=gene.list,mat1=tmC_matrix,yLab="tmCG", nRow=4)
#pp8 <- plot_vln2(genes.in=gene.list.final,mat1=as.matrix(expm1(snRNA3@assays$RNA@data)),yLab="RNA")
#pp4 <- plot_box4_value_nonNeu(genes.in=gene.list,mat1=ratio_matrix,yLab="hmCG/BSmCG", nRow=4)
#pp5 <- plot_box4_value_nonNeu(genes.in=gene.list,mat1=as.matrix(pbmc.snmc2@assays$gbch@counts),yLab="GBmCH")
#pp6 <- plot_box4_value_nonNeu(genes.in=gene.list,mat1=C_matrix,yLab="CG", nRow=4)
#pp7 <- plot_box4_value_nonNeu(genes.in=gene.list,mat1=ratio2_matrix,yLab="hmCG/tmCG")

options(repr.plot.width = 18,repr.plot.height = 24)
plot_grid(pp2,pp3,nrow=2)

C_matrix <- 1-mC_matrix

gene.list.final <- c(
                    "Sv2b", "Arpp21", # Ex marker, "Lingo1",
                    "Gad2",  # Inh marker
                    "Meis2", # Str marker "Ppp1r1b",
                    "Slc1a3", # Astro marker
                    "Mobp", # Oligo marker
                    "Csf1r","Serinc3" # MG marker
                   ) 

#gene.list <- c("Arpp21","Cux2","Rorb","Tle4","Slc6a1","Erbb4","Adarb2","Cacna2d2")
#gene.list2 <- c("Arpp21","Rorb","Tle4","Slc6a1","Erbb4","Adarb2","Gabrg3","Cntnap5c","Il1rapl2")
#pp1 <- plot_box4_value_nonNeu(genes.in=gene.list.final,mat1=mC_matrix,yLab="BSmCG")
pp2 <- plot_box3_value_rmIn_nonNeu(genes.in=gene.list.final,mat1=hmC_matrix,yLab="hmCG", nRow=1)
pp3 <- plot_box3_value_rmIn_nonNeu(genes.in=gene.list.final,mat1=tmC_matrix,yLab="tmCG", nRow=1)
pp4 <- plot_box3_value_rmIn_nonNeu(genes.in=gene.list.final,mat1=C_matrix,yLab="CG", nRow=1)
#pp5 <- plot_box3_value_rmIn_nonNeu(genes.in=gene.list.final,mat1=as.matrix(pbmc.snmc3@assays$gbch@counts),yLab="GBmCH")
#pp6 <- plot_box4_value_nonNeu(genes.in=gene.list,mat1=C_matrix,yLab="CG")
#pp7 <- plot_box4_value_nonNeu(genes.in=gene.list,mat1=ratio2_matrix,yLab="hmCG/tmCG")

pp8 <- plot_vln2(genes.in=gene.list.final,mat1=as.matrix(expm1(snRNA3@assays$RNA@data)),yLab="RNA")

options(repr.plot.width = 16,repr.plot.height = 16)
plot_grid(pp3,pp2,pp4,pp8,nrow=4)

table(Idents(snRNA3))

table(snRNA3@meta.data$res.name.JE3)

snRNA3@meta.data$res.name.6groups <- plyr::mapvalues(snRNA3@meta.data$res.name.JE3,
                from = c("L2/3","L4","L4/5","L6","DL","Pv","Sst","Vip/Ndnf"),
                to = c("Ex","Ex","Ex","Ex","Ex","Inh","Inh","Inh"))

table(snRNA3@meta.data$res.name.6groups)

Idents(snRNA3) <- "res.name.6groups"

RNA.6g.mak <- FindAllMarkers(snRNA3,only.pos = T)

cluster.order.6g <- c("Ex","Inh","Striatum","Astro","Oligo","MG")

dim(RNA.6g.mak)

RNA.6g.mak$cluster <- factor(RNA.6g.mak$cluster,levels = cluster.order.6g)

RNA.6g.mak.pval01 <- subset(RNA.6g.mak, p_val_adj < 0.1)
RNA.6g.mak2 <- RNA.6g.mak.pval01 %>% arrange(cluster,p_val_adj) 

gene.6g.order <- RNA.6g.mak2$gene

length(gene.6g.order) # Peng: 2857 genes; no p_val_adj filtering: 3569 genes; p_val_adj<0.1: 2281 genes

gene.6g.order2 <- intersect(gene.6g.order,rownames(pbmc.snmc3@assays$tmCG@data))

length(gene.6g.order2)
# Peng: 1371 genes; 
# Hao: no p_val_adj filtering: 1867 genes; 
#      p_val_adj<0.1: 1383 genes (90% cells, 5 CGs); 
#                     1525 genes (80% cells, 5 CGs); 
#                     1668 genes (75% cells, 3 CGs);
#                     1732 genes (80% cells, 1 CG);
#                     1685 genes (80% cells, 2 CG)

# Hao's old analysis:
length(gene.6g.order) # Peng: 2857 genes; no p_val_adj filtering: 3569 genes; p_val_adj<0.1: 2281 genes

gene.6g.order2 <- intersect(gene.6g.order,rownames(pbmc.snmc3@assays$tmCG@data))

length(gene.6g.order2) 
# Peng: 1371 genes; Hao: no p_val_adj filtering: 1867 genes; p_val_adj<0.1: 1383 genes (90% cells); 1525 genes (80% cells)

table(pbmc.snmc3@meta.data$res.name2)

pbmc.snmc3@meta.data$res.name.6groups <- plyr::mapvalues(pbmc.snmc3@meta.data$res.name2,
                from = c("L2/3","L4","L4/5","L6","DL","Pv","Sst","Vip/Ndnf"),
                to = c("Ex","Ex","Ex","Ex","Ex","Inh","Inh","Inh"))

table(pbmc.snmc3@meta.data$res.name.6groups)
table(pbmc.snmc3@meta.data$res.name)

gene.order.by.rna.6g <- gene.6g.order2

length(gene.order.by.rna.6g) # 1371; 1867; Hao: 1383 genes (90% cells); 1685 genes (80% cells, 2 CpGs)

gene.order <- gene.order.by.rna.6g
cluster.order <- cluster.order.6g

mat_hmCG_celltype_6g_raw <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name.6groups == x]
    return(rowMeans(hmC_matrix[gene.order,cells.in]))
})

mat_tmC_celltype_6g_raw <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name.6groups == x]
    return(rowMeans(tmC_matrix[gene.order,cells.in]))
})

mat_mCG_celltype_6g_raw <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name.6groups == x]
    return(rowMeans(mC_matrix[gene.order,cells.in]))
})

mat_mCH_celltype_6g_raw <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name.6groups == x]
    return(rowMeans(mCH_matrix[gene.order,cells.in]))
})

mat_ratio_celltype_6g_raw <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name.6groups == x]
    return(rowMeans(ratio_matrix[gene.order,cells.in]))
})

mat_rna_celltype_6g_raw <- sapply(cluster.order,function(x){
    cells.in <- rownames(snRNA3@meta.data)[snRNA3@meta.data$res.name.6groups == x]
    mat <- expm1(snRNA3@assays$RNA@data)
    return(Matrix::rowMeans(mat[gene.order,cells.in]))
})

dim(mat_mCG_celltype_6g_raw)
head(mat_mCG_celltype_6g_raw)

mat_CG_celltype_6g_raw <- 1 - mat_mCG_celltype_6g_raw
dim(mat_CG_celltype_6g_raw)
head(mat_CG_celltype_6g_raw)

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

#h6 <- pheatmap::pheatmap(mat_rna_celltype_6g_raw,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
#                   scale="row",clustering_method = "ward.D2",main="RNA",colours = colorRampPalette(c("blue","white","red"))(50))

# 1685 genes (80% cells, 2 CpGs), 6-16-2022 (final)
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

ggsave("HW-Fig4_6heatmaps_6groups_06.16.2022_1685genes.pdf",g2,width=12,height=3.5)

# 1732 genes 
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

options(repr.plot.width = 12,repr.plot.height = 3)
#grid.arrange(grobs=plot_list,ncol = 4)
g2 <- grid.arrange(grobs=plot_list,ncol = 6)

# 1668 genes 
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

options(repr.plot.width = 12,repr.plot.height = 3)
#grid.arrange(grobs=plot_list,ncol = 4)
g2 <- grid.arrange(grobs=plot_list,ncol = 6)

# 1525 genes
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

options(repr.plot.width = 12,repr.plot.height = 4)
#grid.arrange(grobs=plot_list,ncol = 4)
g2 <- grid.arrange(grobs=plot_list,ncol = 6)

# Peng's old analysis
# combine all the plots
plot_list=list()
plot_list[["tmCG"]] <- h1[[4]]
plot_list[["BSmCG"]] <- h2[[4]]
plot_list[["hmCG"]] <- h3[[4]]
plot_list[["hmCG/BSmCG"]] <- h4[[4]]
plot_list[["BSmCH"]] <- h5[[4]]
plot_list[["RNA"]] <- h6[[4]]

library(gridExtra)

options(repr.plot.width = 12,repr.plot.height = 3)
#grid.arrange(grobs=plot_list,ncol = 4)
g2 <- grid.arrange(grobs=plot_list,ncol = 6)

ggsave("HW-Fig4_5heatmaps_6groups_03.17.2022_1732genes.pdf",g2,width=12,height=3)

ggsave("sc-bACE-seq/HW-Fig4_5heatmaps_6groups_03.12.2022_1668genes.pdf",g2,width=12,height=3)

ggsave("HW-Fig4_5heatmaps_6groups_01.04.2022.pdf",g2,width=12,height=4)

df_CG_celltype_6g_raw<-as.data.frame(mat_CG_celltype_6g_raw)
colnames(df_CG_celltype_6g_raw) <- paste(colnames(df_CG_celltype_6g_raw),"CG",sep="_")

df_tmCG_celltype_6g_raw<-as.data.frame(mat_tmC_celltype_6g_raw)
colnames(df_tmCG_celltype_6g_raw) <- paste(colnames(df_tmCG_celltype_6g_raw),"tmCG",sep="_")

df_hmCG_celltype_6g_raw<-as.data.frame(mat_hmCG_celltype_6g_raw)
colnames(df_hmCG_celltype_6g_raw) <- paste(colnames(df_hmCG_celltype_6g_raw),"hmCG",sep="_")

df_BSmCG_celltype_6g_raw<-as.data.frame(mat_mCG_celltype_6g_raw)
colnames(df_BSmCG_celltype_6g_raw) <- paste(colnames(df_BSmCG_celltype_6g_raw),"BSmCG",sep="_")

df_mCH_celltype_6g_raw<-as.data.frame(mat_mCH_celltype_6g_raw)
colnames(df_mCH_celltype_6g_raw) <- paste(colnames(df_mCH_celltype_6g_raw),"mCH",sep="_")

df_rna_celltype_6g_raw<-as.data.frame(mat_rna_celltype_6g_raw)
colnames(df_rna_celltype_6g_raw) <- paste(colnames(df_rna_celltype_6g_raw),"rna",sep="_")

head(df_CG_celltype_6g_raw)
head(df_tmCG_celltype_6g_raw)
head(df_hmCG_celltype_6g_raw)
head(df_BSmCG_celltype_6g_raw)
head(df_mCH_celltype_6g_raw)
head(df_rna_celltype_6g_raw)

df_6modalities_celltype_markers.6g_raw<-cbind(df_CG_celltype_6g_raw, 
                                         df_tmCG_celltype_6g_raw, 
                                         df_hmCG_celltype_6g_raw, 
                                         df_BSmCG_celltype_6g_raw,
                                         df_mCH_celltype_6g_raw,
                                         df_rna_celltype_6g_raw)

#colnames(df_5modalities_celltype_6g_raw) = gsub("/", ".", colnames(df_5modalities_celltype_6g_raw))
dim(df_6modalities_celltype_markers.6g_raw)
head(df_6modalities_celltype_markers.6g_raw)

head(RNA.6g.mak2)
#head(pbmc.snmc3@assays$tmCG@data)
Ex_rna.specific.gene<-intersect(RNA.6g.mak2 %>% filter(cluster=="Ex") %>% pull(gene),rownames(df_6modalities_celltype_markers.6g_raw))
length(Ex_rna.specific.gene)

Inh_rna.specific.gene<-intersect(RNA.6g.mak2 %>% filter(cluster=="Inh") %>% pull(gene),rownames(df_6modalities_celltype_markers.6g_raw))
length(Inh_rna.specific.gene)

Striatum_rna.specific.gene<-intersect(RNA.6g.mak2 %>% filter(cluster=="Striatum") %>% pull(gene),rownames(df_6modalities_celltype_markers.6g_raw))
length(Striatum_rna.specific.gene)

Astro_rna.specific.gene<-intersect(RNA.6g.mak2 %>% filter(cluster=="Astro") %>% pull(gene),rownames(df_6modalities_celltype_markers.6g_raw))
length(Astro_rna.specific.gene)

Oligo_rna.specific.gene<-intersect(RNA.6g.mak2 %>% filter(cluster=="Oligo") %>% pull(gene),rownames(df_6modalities_celltype_markers.6g_raw))
length(Oligo_rna.specific.gene)

MG_rna.specific.gene<-intersect(RNA.6g.mak2 %>% filter(cluster=="MG") %>% pull(gene),rownames(df_6modalities_celltype_markers.6g_raw))
length(MG_rna.specific.gene)


# Hao's old analysis: 
head(RNA.6g.mak2)
#head(pbmc.snmc3@assays$tmCG@data)
Ex_rna.specific.gene<-intersect(RNA.6g.mak2 %>% filter(cluster=="Ex") %>% pull(gene),rownames(df_6modalities_celltype_markers.6g_raw))
length(Ex_rna.specific.gene)

Inh_rna.specific.gene<-intersect(RNA.6g.mak2 %>% filter(cluster=="Inh") %>% pull(gene),rownames(df_6modalities_celltype_markers.6g_raw))
length(Inh_rna.specific.gene)

Striatum_rna.specific.gene<-intersect(RNA.6g.mak2 %>% filter(cluster=="Striatum") %>% pull(gene),rownames(df_6modalities_celltype_markers.6g_raw))
length(Striatum_rna.specific.gene)

Astro_rna.specific.gene<-intersect(RNA.6g.mak2 %>% filter(cluster=="Astro") %>% pull(gene),rownames(df_6modalities_celltype_markers.6g_raw))
length(Astro_rna.specific.gene)

Oligo_rna.specific.gene<-intersect(RNA.6g.mak2 %>% filter(cluster=="Oligo") %>% pull(gene),rownames(df_6modalities_celltype_markers.6g_raw))
length(Oligo_rna.specific.gene)

MG_rna.specific.gene<-intersect(RNA.6g.mak2 %>% filter(cluster=="MG") %>% pull(gene),rownames(df_6modalities_celltype_markers.6g_raw))
length(MG_rna.specific.gene)


df_6modalities_celltype_markers.6g_raw[c("Gad2","Slc1a3","Mobp","Csf1r","Serinc3"), ]
df_6modalities_celltype_markers.6g_raw["Csf1r", ] %>% select(MG_CG, MG_tmCG, MG_hmCG, Ex_rna, Inh_rna, Striatum_rna, Astro_rna, Oligo_rna,MG_rna)
df.6g_MG <- df_6modalities_celltype_markers.6g_raw[MG_rna.specific.gene, ] %>%  mutate(MG_hmCG.tmCG.ratio = MG_hmCG / MG_tmCG)%>% arrange(desc(MG_hmCG)) %>% filter(log2(MG_rna)>0) %>% select(MG_CG, MG_tmCG, MG_hmCG, MG_mCH, MG_hmCG.tmCG.ratio, Ex_rna, Inh_rna, Striatum_rna, Astro_rna, Oligo_rna,MG_rna) # sort Vip.Ndnf_rna

dim(df.6g_MG)
df.6g_MG

df_CG_celltype_6g_raw["Csf1r", ]
df_tmCG_celltype_6g_raw["Csf1r", ]
df_hmCG_celltype_6g_raw["Csf1r", ]
df_mCH_celltype_6g_raw["Csf1r", ]
df_rna_celltype_6g_raw["Csf1r", ]


Ex.markers <- df_6modalities_celltype_markers.6g_raw[Ex_rna.specific.gene, ]
Inh.markers <- df_6modalities_celltype_markers.6g_raw[Inh_rna.specific.gene, ]
Striatum.markers <- df_6modalities_celltype_markers.6g_raw[Striatum_rna.specific.gene, ]
Oligo.markers <- df_6modalities_celltype_markers.6g_raw[Oligo_rna.specific.gene, ]
Astro.markers <- df_6modalities_celltype_markers.6g_raw[Astro_rna.specific.gene, ]
MG.markers <- df_6modalities_celltype_markers.6g_raw[MG_rna.specific.gene, ]

dim(Ex.markers)
dim(Inh.markers)
dim(Striatum.markers)
dim(Oligo.markers)
dim(Astro.markers)
dim(MG.markers)

df.stat <- data.frame(cell_type=c("Ex","Inh","Striatum","Astro","Oligo","MG"), pvalue=1:6, cor=1:6)
rownames(df.stat) <- df.stat$cell_type

ts.Ex <- cor.test(Ex.markers$Ex_tmCG, Ex.markers$Ex_hmCG)
ts.Inh <- cor.test(Inh.markers$Inh_tmCG, Inh.markers$Inh_hmCG)
ts.Striatum <- cor.test(Striatum.markers$Striatum_tmCG, Striatum.markers$Striatum_hmCG)
ts.Oligo <- cor.test(Oligo.markers$Oligo_tmCG, Oligo.markers$Oligo_hmCG)
ts.Astro <- cor.test(Astro.markers$Astro_tmCG, Astro.markers$Astro_hmCG)
ts.MG <- cor.test(MG.markers$MG_tmCG, MG.markers$MG_hmCG)

df.stat["Ex",c("pvalue","cor")] <- t(c(as.numeric(ts.Ex$p.value), as.numeric(ts.Ex$estimate)))
df.stat["Inh",2:3] <- t(c(as.numeric(ts.Inh$p.value), as.numeric(ts.Inh$estimate)))
df.stat["Striatum",2:3] <- t(c(as.numeric(ts.Striatum$p.value), as.numeric(ts.Striatum$estimate)))
df.stat["Oligo",2:3] <- t(c(as.numeric(ts.Oligo$p.value), as.numeric(ts.Oligo$estimate)))
df.stat["Astro",2:3] <- t(c(as.numeric(ts.Astro$p.value), as.numeric(ts.Astro$estimate)))
df.stat["MG",2:3] <- t(c(as.numeric(ts.MG$p.value), as.numeric(ts.MG$estimate)))

df.stat



options(repr.plot.width = 12,repr.plot.height = 12)

#myPalette <- colorRampPalette(rev(brewer.pal(11, "RdBu"))) 
myPalette <- colorRampPalette(brewer.pal(11, "YlOrRd")) 
#myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(-1,8))
#sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(0.5,4))

Ex.markers_hmCG.tmCG <- ggplot(Ex.markers,aes(x=Ex_tmCG*100, y=Ex_hmCG*100)) +
    geom_point(aes(colour=log2(Ex_rna)), size=2, alpha=0.6) + sc +
    theme_bw() + 
    theme(aspect.ratio=1) +
    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
    scale_y_continuous(breaks=seq(0,50,10), labels=seq(0,50,10), limits=c(0,50)) + 
    geom_text(data = df.stat["Ex", ],
                 aes(x = 80, y=45, label = formatC(pvalue,digits = 2)),color = "red",hjust = 0) +
    geom_text(data = df.stat["Ex", ],
                 aes(x = 80, y=40, label = round(cor,3)),color = "red",hjust = 0) +
    geom_smooth(method='lm', formula= y~x)

Inh.markers_hmCG.tmCG <- ggplot(Inh.markers,aes(x=Inh_tmCG*100, y=Inh_hmCG*100)) +
    geom_point(aes(colour=log2(Inh_rna)), size=2, alpha=0.6) + sc +
    theme_bw() + 
    theme(aspect.ratio=1) +
    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
    scale_y_continuous(breaks=seq(0,50,10), labels=seq(0,50,10), limits=c(0,50)) + 
    geom_text(data = df.stat["Inh", ],
                 aes(x = 80, y=45, label = formatC(pvalue,digits = 2)),color = "red",hjust = 0) +
    geom_text(data = df.stat["Inh", ],
                 aes(x = 80, y=40, label = round(cor,3)),color = "red",hjust = 0) +
    geom_smooth(method='lm', formula= y~x)

Striatum.markers_hmCG.tmCG <- ggplot(Striatum.markers,aes(x=Striatum_tmCG*100, y=Striatum_hmCG*100)) +
    geom_point(aes(colour=log2(Striatum_rna)), size=2, alpha=0.6) + sc +
    theme_bw() + 
    theme(aspect.ratio=1) +
    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
    scale_y_continuous(breaks=seq(0,50,10), labels=seq(0,50,10), limits=c(0,50)) + 
    geom_text(data = df.stat["Striatum", ],
                 aes(x = 80, y=45, label = formatC(pvalue,digits = 2)),color = "red",hjust = 0) +
    geom_text(data = df.stat["Striatum", ],
                 aes(x = 80, y=40, label = round(cor,3)),color = "red",hjust = 0) +
    geom_smooth(method='lm', formula= y~x)

Oligo.markers_hmCG.tmCG <- ggplot(Oligo.markers,aes(x=Oligo_tmCG*100, y=Oligo_hmCG*100)) +
    geom_point(aes(colour=log2(Oligo_rna)), size=2, alpha=0.6) + sc +
    theme_bw() + 
    theme(aspect.ratio=1) +
    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
    scale_y_continuous(breaks=seq(0,50,10), labels=seq(0,50,10), limits=c(0,50)) + 
    geom_text(data = df.stat["Oligo", ],
                 aes(x = 80, y=45, label = formatC(pvalue,digits = 2)),color = "red",hjust = 0) +
    geom_text(data = df.stat["Oligo", ],
                 aes(x = 80, y=40, label = round(cor,3)),color = "red",hjust = 0) +
    geom_smooth(method='lm', formula= y~x)

Astro.markers_hmCG.tmCG <- ggplot(Astro.markers,aes(x=Astro_tmCG*100, y=Astro_hmCG*100)) +
    geom_point(aes(colour=log2(Astro_rna)), size=2, alpha=0.6) + sc +
    theme_bw() + 
    theme(aspect.ratio=1) +
    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
    scale_y_continuous(breaks=seq(0,50,10), labels=seq(0,50,10), limits=c(0,50)) + 
    geom_text(data = df.stat["Astro", ],
                 aes(x = 80, y=45, label = formatC(pvalue,digits = 2)),color = "red",hjust = 0) +
    geom_text(data = df.stat["Astro", ],
                 aes(x = 80, y=40, label = round(cor,3)),color = "red",hjust = 0) +
    geom_smooth(method='lm', formula= y~x)

MG.markers_hmCG.tmCG <- ggplot(MG.markers,aes(x=MG_tmCG*100, y=MG_hmCG*100)) +
    geom_point(aes(colour=log2(MG_rna)), size=2, alpha=0.6) + sc +
    theme_bw() + 
    theme(aspect.ratio=1) +
    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
    scale_y_continuous(breaks=seq(0,50,10), labels=seq(0,50,10), limits=c(0,50)) +  
    geom_text(data = df.stat["MG", ],
                 aes(x = 80, y=45, label = formatC(pvalue,digits = 2)),color = "red",hjust = 0) +
    geom_text(data = df.stat["MG", ],
                 aes(x = 80, y=40, label = round(cor,3)),color = "red",hjust = 0) +
    geom_smooth(method='lm', formula= y~x)


plot_grid(
    Ex.markers_hmCG.tmCG,   
    Astro.markers_hmCG.tmCG,
    Inh.markers_hmCG.tmCG,
    Oligo.markers_hmCG.tmCG,
    Striatum.markers_hmCG.tmCG,
    MG.markers_hmCG.tmCG,
    nrow=3,align = "vh")

options(repr.plot.width = 15,repr.plot.height = 6)

Ex.markers_hmCG.rna <- ggplot(Ex.markers,aes(x=Ex_hmCG*100, y=log2(Ex_rna))) +
    geom_point(aes(colour=log2(Ex_rna)), size=1, alpha=0.6) + sc +
    theme_bw() + 
    theme(aspect.ratio=1) +
    geom_smooth(method='lm', formula= y~x)

Ex.markers_tmCG.rna <- ggplot(Ex.markers,aes(x=Ex_tmCG*100, y=log2(Ex_rna))) +
    geom_point(aes(colour=log2(Ex_rna)), size=1, alpha=0.6) + sc +
    theme_bw() + 
    theme(aspect.ratio=1) +
    geom_smooth(method='lm', formula= y~x)

MG.markers_hmCG.rna <- ggplot(MG.markers,aes(x=MG_hmCG*100, y=log2(MG_rna))) +
    geom_point(aes(colour=log2(MG_rna)), size=1, alpha=0.6) + sc +
    theme_bw() + 
    theme(aspect.ratio=1) +
    geom_smooth(method='lm', formula= y~x)

MG.markers_tmCG.rna <- ggplot(MG.markers,aes(x=MG_tmCG*100, y=log2(MG_rna))) +
    geom_point(aes(colour=log2(MG_rna)), size=1, alpha=0.6) + sc +
    theme_bw() + 
    theme(aspect.ratio=1) +
    geom_smooth(method='lm', formula= y~x)

plot_grid(
    Ex.markers_hmCG.rna,   
    Ex.markers_tmCG.rna,
    Ex.markers_hmCG.tmCG,
    MG.markers_hmCG.rna,   
    MG.markers_tmCG.rna,
    MG.markers_hmCG.tmCG,
    nrow=2,align = "vh")

options(repr.plot.width = 12,repr.plot.height = 12)

#myPalette <- colorRampPalette(rev(brewer.pal(11, "RdBu"))) 
myPalette <- colorRampPalette(brewer.pal(11, "YlOrRd")) 
#myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(-1,8))
#sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(0.5,4))

Ex.markers_hmCG.tmCG <- ggplot(Ex.markers,aes(x=Ex_tmCG*100, y=Ex_hmCG*100)) +
    geom_point(aes(colour=log2(Ex_rna)), size=2, alpha=0.6) + sc +
    theme_bw() + 
    theme(aspect.ratio=1) +
    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
#    scale_y_continuous(breaks=seq(0,50,10), labels=seq(0,50,10), limits=c(0,50)) + 
    geom_text(data = df.stat["Ex", ],
                 aes(x = 80, y=15, label = formatC(pvalue,digits = 2)),color = "red",hjust = 0) +
    geom_text(data = df.stat["Ex", ],
                 aes(x = 80, y=10, label = round(cor,3)),color = "red",hjust = 0) +
    geom_smooth(method='lm', formula= y~x)
  #   stat_smooth(method="lm", se=TRUE, formula=y ~ poly(x, 3, raw=FALSE))



Inh.markers_hmCG.tmCG <- ggplot(Inh.markers,aes(x=Inh_tmCG*100, y=Inh_hmCG*100)) +
    geom_point(aes(colour=log2(Inh_rna)), size=2, alpha=0.6) + sc +
    theme_bw() + 
    theme(aspect.ratio=1) +
    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
#    scale_y_continuous(breaks=seq(0,50,10), labels=seq(0,50,10), limits=c(0,50)) + 
    geom_text(data = df.stat["Inh", ],
                 aes(x = 80, y=15, label = formatC(pvalue,digits = 2)),color = "red",hjust = 0) +
    geom_text(data = df.stat["Inh", ],
                 aes(x = 80, y=10, label = round(cor,3)),color = "red",hjust = 0) +
    geom_smooth(method='lm', formula= y~x)

Striatum.markers_hmCG.tmCG <- ggplot(Striatum.markers,aes(x=Striatum_tmCG*100, y=Striatum_hmCG*100)) +
    geom_point(aes(colour=log2(Striatum_rna)), size=2, alpha=0.6) + sc +
    theme_bw() + 
    theme(aspect.ratio=1) +
    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
#    scale_y_continuous(breaks=seq(0,50,10), labels=seq(0,50,10), limits=c(0,50)) + 
    geom_text(data = df.stat["Striatum", ],
                 aes(x = 80, y=15, label = formatC(pvalue,digits = 2)),color = "red",hjust = 0) +
    geom_text(data = df.stat["Striatum", ],
                 aes(x = 80, y=10, label = round(cor,3)),color = "red",hjust = 0) +
    geom_smooth(method='lm', formula= y~x)

Oligo.markers_hmCG.tmCG <- ggplot(Oligo.markers,aes(x=Oligo_tmCG*100, y=Oligo_hmCG*100)) +
    geom_point(aes(colour=log2(Oligo_rna)), size=2, alpha=0.6) + sc +
    theme_bw() + 
    theme(aspect.ratio=1) +
    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
#    scale_y_continuous(breaks=seq(0,50,10), labels=seq(0,50,10), limits=c(0,50)) + 
    geom_text(data = df.stat["Oligo", ],
                 aes(x = 80, y=15, label = formatC(pvalue,digits = 2)),color = "red",hjust = 0) +
    geom_text(data = df.stat["Oligo", ],
                 aes(x = 80, y=10, label = round(cor,3)),color = "red",hjust = 0) +
    geom_smooth(method='lm', formula= y~x)

Astro.markers_hmCG.tmCG <- ggplot(Astro.markers,aes(x=Astro_tmCG*100, y=Astro_hmCG*100)) +
    geom_point(aes(colour=log2(Astro_rna)), size=2, alpha=0.6) + sc +
    theme_bw() + 
    theme(aspect.ratio=1) +
    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
#    scale_y_continuous(breaks=seq(0,50,10), labels=seq(0,50,10), limits=c(0,50)) + 
    geom_text(data = df.stat["Astro", ],
                 aes(x = 80, y=15, label = formatC(pvalue,digits = 2)),color = "red",hjust = 0) +
    geom_text(data = df.stat["Astro", ],
                 aes(x = 80, y=10, label = round(cor,3)),color = "red",hjust = 0) +
    geom_smooth(method='lm', formula= y~x)

MG.markers_hmCG.tmCG <- ggplot(MG.markers,aes(x=MG_tmCG*100, y=MG_hmCG*100)) +
    geom_point(aes(colour=log2(MG_rna)), size=2, alpha=0.6) + sc +
    theme_bw() + 
    theme(aspect.ratio=1) +
    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
#    scale_y_continuous(breaks=seq(0,50,10), labels=seq(0,50,10), limits=c(0,50)) +  
    geom_text(data = df.stat["MG", ],
                 aes(x = 80, y=15, label = formatC(pvalue,digits = 2)),color = "red",hjust = 0) +
    geom_text(data = df.stat["MG", ],
                 aes(x = 80, y=10, label = round(cor,3)),color = "red",hjust = 0) +
    geom_smooth(method='lm', formula= y~x)


plot_grid(
    Ex.markers_hmCG.tmCG,   
    Astro.markers_hmCG.tmCG,
    Inh.markers_hmCG.tmCG,
    Oligo.markers_hmCG.tmCG,
    Striatum.markers_hmCG.tmCG,
    MG.markers_hmCG.tmCG,
    nrow=3,align = "vh")

ggsave("Fig3_celltype.6g.removeIn_genebody.tmCG.vs.hmCG_06-16-2022_scatterplot.pdf",width = 20 ,height = 20)

options(repr.plot.width = 12,repr.plot.height = 12)

#myPalette <- colorRampPalette(rev(brewer.pal(11, "RdBu"))) 
myPalette <- colorRampPalette(brewer.pal(11, "YlOrRd")) 
#myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(-1,8))
#sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(0.5,4))

Ex.markers_hmCG.tmCG <- ggplot(Ex.markers,aes(x=Ex_tmCG*100, y=Ex_hmCG*100)) +
    geom_point(aes(colour=log2(Ex_rna)), size=2, alpha=0.6) + sc +
    theme_bw() + 
    theme(aspect.ratio=1) +
    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
#    scale_y_continuous(breaks=seq(0,50,10), labels=seq(0,50,10), limits=c(0,50)) + 
    geom_text(data = df.stat["Ex", ],
                 aes(x = 80, y=15, label = formatC(pvalue,digits = 2)),color = "red",hjust = 0) +
    geom_text(data = df.stat["Ex", ],
                 aes(x = 80, y=10, label = round(cor,3)),color = "red",hjust = 0) +
  #  geom_smooth(method='lm', formula= y~x)
     stat_smooth(method="lm", se=TRUE, formula=y ~ poly(x, 3, raw=FALSE))



Inh.markers_hmCG.tmCG <- ggplot(Inh.markers,aes(x=Inh_tmCG*100, y=Inh_hmCG*100)) +
    geom_point(aes(colour=log2(Inh_rna)), size=2, alpha=0.6) + sc +
    theme_bw() + 
    theme(aspect.ratio=1) +
    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
#    scale_y_continuous(breaks=seq(0,50,10), labels=seq(0,50,10), limits=c(0,50)) + 
    geom_text(data = df.stat["Inh", ],
                 aes(x = 80, y=15, label = formatC(pvalue,digits = 2)),color = "red",hjust = 0) +
    geom_text(data = df.stat["Inh", ],
                 aes(x = 80, y=10, label = round(cor,3)),color = "red",hjust = 0) +
   #  geom_smooth(method='lm', formula= y~x)
     stat_smooth(method="lm", se=TRUE, formula=y ~ poly(x, 3, raw=FALSE))

Striatum.markers_hmCG.tmCG <- ggplot(Striatum.markers,aes(x=Striatum_tmCG*100, y=Striatum_hmCG*100)) +
    geom_point(aes(colour=log2(Striatum_rna)), size=2, alpha=0.6) + sc +
    theme_bw() + 
    theme(aspect.ratio=1) +
    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
#    scale_y_continuous(breaks=seq(0,50,10), labels=seq(0,50,10), limits=c(0,50)) + 
    geom_text(data = df.stat["Striatum", ],
                 aes(x = 80, y=15, label = formatC(pvalue,digits = 2)),color = "red",hjust = 0) +
    geom_text(data = df.stat["Striatum", ],
                 aes(x = 80, y=10, label = round(cor,3)),color = "red",hjust = 0) +
  #  geom_smooth(method='lm', formula= y~x)
     stat_smooth(method="lm", se=TRUE, formula=y ~ poly(x, 3, raw=FALSE))

Oligo.markers_hmCG.tmCG <- ggplot(Oligo.markers,aes(x=Oligo_tmCG*100, y=Oligo_hmCG*100)) +
    geom_point(aes(colour=log2(Oligo_rna)), size=2, alpha=0.6) + sc +
    theme_bw() + 
    theme(aspect.ratio=1) +
    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
#    scale_y_continuous(breaks=seq(0,50,10), labels=seq(0,50,10), limits=c(0,50)) + 
    geom_text(data = df.stat["Oligo", ],
                 aes(x = 80, y=15, label = formatC(pvalue,digits = 2)),color = "red",hjust = 0) +
    geom_text(data = df.stat["Oligo", ],
                 aes(x = 80, y=10, label = round(cor,3)),color = "red",hjust = 0) +
   #  geom_smooth(method='lm', formula= y~x)
     stat_smooth(method="lm", se=TRUE, formula=y ~ poly(x, 3, raw=FALSE))

Astro.markers_hmCG.tmCG <- ggplot(Astro.markers,aes(x=Astro_tmCG*100, y=Astro_hmCG*100)) +
    geom_point(aes(colour=log2(Astro_rna)), size=2, alpha=0.6) + sc +
    theme_bw() + 
    theme(aspect.ratio=1) +
    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
#    scale_y_continuous(breaks=seq(0,50,10), labels=seq(0,50,10), limits=c(0,50)) + 
    geom_text(data = df.stat["Astro", ],
                 aes(x = 80, y=15, label = formatC(pvalue,digits = 2)),color = "red",hjust = 0) +
    geom_text(data = df.stat["Astro", ],
                 aes(x = 80, y=10, label = round(cor,3)),color = "red",hjust = 0) +
   #  geom_smooth(method='lm', formula= y~x)
     stat_smooth(method="lm", se=TRUE, formula=y ~ poly(x, 3, raw=FALSE))

MG.markers_hmCG.tmCG <- ggplot(MG.markers,aes(x=MG_tmCG*100, y=MG_hmCG*100)) +
    geom_point(aes(colour=log2(MG_rna)), size=2, alpha=0.6) + sc +
    theme_bw() + 
    theme(aspect.ratio=1) +
    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
#    scale_y_continuous(breaks=seq(0,50,10), labels=seq(0,50,10), limits=c(0,50)) +  
    geom_text(data = df.stat["MG", ],
                 aes(x = 80, y=15, label = formatC(pvalue,digits = 2)),color = "red",hjust = 0) +
    geom_text(data = df.stat["MG", ],
                 aes(x = 80, y=10, label = round(cor,3)),color = "red",hjust = 0) +
  #  geom_smooth(method='lm', formula= y~x)
     stat_smooth(method="lm", se=TRUE, formula=y ~ poly(x, 3, raw=FALSE))


plot_grid(
    Ex.markers_hmCG.tmCG,   
    Astro.markers_hmCG.tmCG,
    Inh.markers_hmCG.tmCG,
    Oligo.markers_hmCG.tmCG,
    Striatum.markers_hmCG.tmCG,
    MG.markers_hmCG.tmCG,
    nrow=3,align = "vh")

options(repr.plot.width = 12,repr.plot.height = 12)

#myPalette <- colorRampPalette(rev(brewer.pal(11, "RdBu"))) 
#myPalette <- colorRampPalette(brewer.pal(11, "YlOrRd")) 
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(-1,6))
#sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(0.5,4))

Ex.markers <- df_6modalities_celltype_markers.6g_raw[Ex_rna.specific.gene, ] 
Inh.markers <- df_6modalities_celltype_markers.6g_raw[Inh_rna.specific.gene, ]
Striatum.markers <- df_6modalities_celltype_markers.6g_raw[Striatum_rna.specific.gene, ]
Oligo.markers <- df_6modalities_celltype_markers.6g_raw[Oligo_rna.specific.gene, ]
Astro.markers <- df_6modalities_celltype_markers.6g_raw[Astro_rna.specific.gene, ]
MG.markers <- df_6modalities_celltype_markers.6g_raw[MG_rna.specific.gene, ]

Ex.markers_hmCG.CG <- ggplot(Ex.markers,aes(x=Ex_CG*100, y=Ex_hmCG*100)) +
    geom_point(aes(colour=log2(Ex_rna)), size=2, alpha=0.6) + sc +
    theme_bw() + 
    theme(aspect.ratio=1) +
    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
#  geom_smooth(method='lm', formula= y~x)
     stat_smooth(method="lm", se=TRUE, formula=y ~ poly(x, 3, raw=FALSE))

Inh.markers_hmCG.CG <- ggplot(Inh.markers,aes(x=Inh_CG*100, y=Inh_hmCG*100)) +
    geom_point(aes(colour=log2(Inh_rna)), size=2, alpha=0.6) + sc +
    theme_bw() + 
    theme(aspect.ratio=1) +
    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
 #  geom_smooth(method='lm', formula= y~x)
     stat_smooth(method="lm", se=TRUE, formula=y ~ poly(x, 3, raw=FALSE))

Striatum.markers_hmCG.CG <- ggplot(Striatum.markers,aes(x=Striatum_CG*100, y=Striatum_hmCG*100)) +
    geom_point(aes(colour=log2(Striatum_rna)), size=2, alpha=0.6) + sc +
    theme_bw() + 
    theme(aspect.ratio=1) +
    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
 #  geom_smooth(method='lm', formula= y~x)
     stat_smooth(method="lm", se=TRUE, formula=y ~ poly(x, 3, raw=FALSE))

Oligo.markers_hmCG.CG <- ggplot(Oligo.markers,aes(x=Oligo_CG*100, y=Oligo_hmCG*100)) +
    geom_point(aes(colour=log2(Oligo_rna)), size=2, alpha=0.6) + sc +
    theme_bw() + 
    theme(aspect.ratio=1) +
    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
#  geom_smooth(method='lm', formula= y~x)
     stat_smooth(method="lm", se=TRUE, formula=y ~ poly(x, 3, raw=FALSE))

Astro.markers_hmCG.CG <- ggplot(Astro.markers,aes(x=Astro_CG*100, y=Astro_hmCG*100)) +
    geom_point(aes(colour=log2(Astro_rna)), size=2, alpha=0.6) + sc +
    theme_bw() + 
    theme(aspect.ratio=1) +
    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
 #  geom_smooth(method='lm', formula= y~x)
     stat_smooth(method="lm", se=TRUE, formula=y ~ poly(x, 3, raw=FALSE))

MG.markers_hmCG.CG <- ggplot(MG.markers,aes(x=MG_CG*100, y=MG_hmCG*100)) +
    geom_point(aes(colour=log2(MG_rna)), size=2, alpha=0.6) + sc +
    theme_bw() + 
    theme(aspect.ratio=1) +
    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
#    scale_y_continuous(breaks=seq(0,4,1), labels=seq(0,4,1), limits=c(0,4)) + 
#  geom_smooth(method='lm', formula= y~x)
     stat_smooth(method="lm", se=TRUE, formula=y ~ poly(x, 3, raw=FALSE))


plot_grid(
    Ex.markers_hmCG.CG,   
    Astro.markers_hmCG.CG,
    Inh.markers_hmCG.CG,
    Oligo.markers_hmCG.CG,
    Striatum.markers_hmCG.CG,
    MG.markers_hmCG.CG,
      #    L6.markers_ratio.CG,
#    Vip.Ndnf.markers_ratio.CG,
#    L6.markers_mCH.CG,
#    Vip.Ndnf.markers_mCH.CG,
    nrow=3,align = "vh")

ggsave("Fig3_celltype.6g.removeIn_genebody.hmCG.vs.CG_poly_scatterplot.pdf", width = 12 ,height = 12)  


Ex.markers.log2rna1 <- df_6modalities_celltype_markers.6g_raw[Ex_rna.specific.gene, ] %>% filter(log2(Ex_rna)>1)
Inh.markers.log2rna1 <- df_6modalities_celltype_markers.6g_raw[Inh_rna.specific.gene, ] %>% filter(log2(Inh_rna)>1)
Striatum.markers.log2rna1 <- df_6modalities_celltype_markers.6g_raw[Striatum_rna.specific.gene, ] %>% filter(log2(Striatum_rna)>1)
Oligo.markers.log2rna1 <- df_6modalities_celltype_markers.6g_raw[Oligo_rna.specific.gene, ] %>% filter(log2(Oligo_rna)>1)
Astro.markers.log2rna1 <- df_6modalities_celltype_markers.6g_raw[Astro_rna.specific.gene, ] %>% filter(log2(Astro_rna)>1)
MG.markers.log2rna1 <- df_6modalities_celltype_markers.6g_raw[MG_rna.specific.gene, ] %>% filter(log2(MG_rna)>1)

dim(Ex.markers.log2rna1)
dim(Inh.markers.log2rna1)
dim(Striatum.markers.log2rna1)
dim(Oligo.markers.log2rna1)
dim(Astro.markers.log2rna1)
dim(MG.markers.log2rna1)

df.stat <- data.frame(cell_type=c("Ex","Inh","Striatum","Astro","Oligo","MG"), pvalue=1:6, cor=1:6)
rownames(df.stat) <- df.stat$cell_type

ts.Ex <- cor.test(Ex.markers.log2rna1$Ex_CG, Ex.markers.log2rna1$Ex_hmCG)
ts.Inh <- cor.test(Inh.markers.log2rna1$Inh_CG, Inh.markers.log2rna1$Inh_hmCG)
ts.Striatum <- cor.test(Striatum.markers.log2rna1$Striatum_CG, Striatum.markers.log2rna1$Striatum_hmCG)
ts.Oligo <- cor.test(Oligo.markers.log2rna1$Oligo_CG, Oligo.markers.log2rna1$Oligo_hmCG)
ts.Astro <- cor.test(Astro.markers.log2rna1$Astro_CG, Astro.markers.log2rna1$Astro_hmCG)
ts.MG <- cor.test(MG.markers.log2rna1$MG_CG, MG.markers.log2rna1$MG_hmCG)

df.stat["Ex",c("pvalue","cor")] <- t(c(as.numeric(ts.Ex$p.value), as.numeric(ts.Ex$estimate)))
df.stat["Inh",2:3] <- t(c(as.numeric(ts.Inh$p.value), as.numeric(ts.Inh$estimate)))
df.stat["Striatum",2:3] <- t(c(as.numeric(ts.Striatum$p.value), as.numeric(ts.Striatum$estimate)))
df.stat["Oligo",2:3] <- t(c(as.numeric(ts.Oligo$p.value), as.numeric(ts.Oligo$estimate)))
df.stat["Astro",2:3] <- t(c(as.numeric(ts.Astro$p.value), as.numeric(ts.Astro$estimate)))
df.stat["MG",2:3] <- t(c(as.numeric(ts.MG$p.value), as.numeric(ts.MG$estimate)))

df.stat


options(repr.plot.width = 12,repr.plot.height = 12)

myPalette <- colorRampPalette(rev(brewer.pal(11, "RdBu"))) 
#myPalette <- colorRampPalette(brewer.pal(11, "YlOrRd")) 
#myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(-1,6))
#sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(0.5,4))

Ex.markers_hmCG.CG <- ggplot(Ex.markers.log2rna1,aes(x=Ex_CG*100, y=Ex_hmCG*100)) +
    geom_point(aes(colour=log2(Ex_rna)), size=2, alpha=0.6) + sc +
    theme_bw() + 
    theme(aspect.ratio=1) +
    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
#  geom_smooth(method='lm', formula= y~x)
     stat_smooth(method="lm", se=TRUE, formula=y ~ poly(x, 3, raw=FALSE))

Inh.markers_hmCG.CG <- ggplot(Inh.markers.log2rna1,aes(x=Inh_CG*100, y=Inh_hmCG*100)) +
    geom_point(aes(colour=log2(Inh_rna)), size=2, alpha=0.6) + sc +
    theme_bw() + 
    theme(aspect.ratio=1) +
    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
 #  geom_smooth(method='lm', formula= y~x)
     stat_smooth(method="lm", se=TRUE, formula=y ~ poly(x, 3, raw=FALSE))

Striatum.markers_hmCG.CG <- ggplot(Striatum.markers.log2rna1,aes(x=Striatum_CG*100, y=Striatum_hmCG*100)) +
    geom_point(aes(colour=log2(Striatum_rna)), size=2, alpha=0.6) + sc +
    theme_bw() + 
    theme(aspect.ratio=1) +
    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
 #  geom_smooth(method='lm', formula= y~x)
     stat_smooth(method="lm", se=TRUE, formula=y ~ poly(x, 3, raw=FALSE))

Oligo.markers_hmCG.CG <- ggplot(Oligo.markers.log2rna1,aes(x=Oligo_CG*100, y=Oligo_hmCG*100)) +
    geom_point(aes(colour=log2(Oligo_rna)), size=2, alpha=0.6) + sc +
    theme_bw() + 
    theme(aspect.ratio=1) +
    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
#  geom_smooth(method='lm', formula= y~x)
     stat_smooth(method="lm", se=TRUE, formula=y ~ poly(x, 3, raw=FALSE))

Astro.markers_hmCG.CG <- ggplot(Astro.markers.log2rna1,aes(x=Astro_CG*100, y=Astro_hmCG*100)) +
    geom_point(aes(colour=log2(Astro_rna)), size=2, alpha=0.6) + sc +
    theme_bw() + 
    theme(aspect.ratio=1) +
    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
 #  geom_smooth(method='lm', formula= y~x)
     stat_smooth(method="lm", se=TRUE, formula=y ~ poly(x, 3, raw=FALSE))

MG.markers_hmCG.CG <- ggplot(MG.markers.log2rna1,aes(x=MG_CG*100, y=MG_hmCG*100)) +
    geom_point(aes(colour=log2(MG_rna)), size=2, alpha=0.6) + sc +
    theme_bw() + 
    theme(aspect.ratio=1) +
    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
#    scale_y_continuous(breaks=seq(0,4,1), labels=seq(0,4,1), limits=c(0,4)) + 
#  geom_smooth(method='lm', formula= y~x)
     stat_smooth(method="lm", se=TRUE, formula=y ~ poly(x, 3, raw=FALSE))


plot_grid(
    Ex.markers_hmCG.CG,   
    Astro.markers_hmCG.CG,
    Inh.markers_hmCG.CG,
    Oligo.markers_hmCG.CG,
    Striatum.markers_hmCG.CG,
    MG.markers_hmCG.CG,
      #    L6.markers_ratio.CG,
#    Vip.Ndnf.markers_ratio.CG,
#    L6.markers_mCH.CG,
#    Vip.Ndnf.markers_mCH.CG,
    nrow=3,align = "vh")

ggsave("Fig3_celltype.6g.removeIn_genebody.hmCG.vs.CG_markers.log2rna1_scatterplot.pdf",width = 12 ,height = 12) 

options(repr.plot.width = 20,repr.plot.height = 20)

myPalette <- colorRampPalette(rev(brewer.pal(11, "RdBu"))) 
#myPalette <- colorRampPalette(brewer.pal(11, "YlOrRd")) 
#myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(-1,8))
#sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(0.5,4))


Ex.markers.log2rna1_hmCG.CG <- ggplot(Ex.markers.log2rna1,aes(x=Ex_CG*100, y=Ex_hmCG*100)) +
    geom_point(aes(colour=log2(Ex_rna)), size=2, alpha=0.6) + sc +
    theme_bw() + 
    theme(aspect.ratio=1) +
    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
    scale_y_continuous(breaks=seq(0,50,10), labels=seq(0,50,10), limits=c(0,50)) + 
    geom_text(data = df.stat["Ex", ],
                 aes(x = 80, y=45, label = formatC(pvalue,digits = 2)),color = "red",hjust = 0) +
    geom_text(data = df.stat["Ex", ],
                 aes(x = 80, y=40, label = round(cor,3)),color = "red",hjust = 0) +
#  geom_smooth(method='lm', formula= y~x)
     stat_smooth(method="lm", se=TRUE, formula=y ~ poly(x, 3, raw=FALSE))

Inh.markers.log2rna1_hmCG.CG <- ggplot(Inh.markers.log2rna1,aes(x=Inh_CG*100, y=Inh_hmCG*100)) +
    geom_point(aes(colour=log2(Inh_rna)), size=2, alpha=0.6) + sc +
    theme_bw() + 
    theme(aspect.ratio=1) +
    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
    scale_y_continuous(breaks=seq(0,50,10), labels=seq(0,50,10), limits=c(0,50)) + 
    geom_text(data = df.stat["Inh", ],
                 aes(x = 80, y=45, label = formatC(pvalue,digits = 2)),color = "red",hjust = 0) +
    geom_text(data = df.stat["Inh", ],
                 aes(x = 80, y=40, label = round(cor,3)),color = "red",hjust = 0) +
   #  geom_smooth(method='lm', formula= y~x)
     stat_smooth(method="lm", se=TRUE, formula=y ~ poly(x, 3, raw=FALSE))

Striatum.markers.log2rna1_hmCG.CG <- ggplot(Striatum.markers.log2rna1,aes(x=Striatum_CG*100, y=Striatum_hmCG*100)) +
    geom_point(aes(colour=log2(Striatum_rna)), size=2, alpha=0.6) + sc +
    theme_bw() + 
    theme(aspect.ratio=1) +
    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
    scale_y_continuous(breaks=seq(0,50,10), labels=seq(0,50,10), limits=c(0,50)) + 
    geom_text(data = df.stat["Striatum", ],
                 aes(x = 80, y=45, label = formatC(pvalue,digits = 2)),color = "red",hjust = 0) +
    geom_text(data = df.stat["Striatum", ],
                 aes(x = 80, y=40, label = round(cor,3)),color = "red",hjust = 0) +
  #  geom_smooth(method='lm', formula= y~x)
     stat_smooth(method="lm", se=TRUE, formula=y ~ poly(x, 3, raw=FALSE))

Oligo.markers.log2rna1_hmCG.CG <- ggplot(Oligo.markers.log2rna1,aes(x=Oligo_CG*100, y=Oligo_hmCG*100)) +
    geom_point(aes(colour=log2(Oligo_rna)), size=2, alpha=0.6) + sc +
    theme_bw() + 
    theme(aspect.ratio=1) +
    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
    scale_y_continuous(breaks=seq(0,50,10), labels=seq(0,50,10), limits=c(0,50)) + 
    geom_text(data = df.stat["Oligo", ],
                 aes(x = 80, y=45, label = formatC(pvalue,digits = 2)),color = "red",hjust = 0) +
    geom_text(data = df.stat["Oligo", ],
                 aes(x = 80, y=40, label = round(cor,3)),color = "red",hjust = 0) +
   #  geom_smooth(method='lm', formula= y~x)
     stat_smooth(method="lm", se=TRUE, formula=y ~ poly(x, 3, raw=FALSE))

Astro.markers.log2rna1_hmCG.CG <- ggplot(Astro.markers.log2rna1,aes(x=Astro_CG*100, y=Astro_hmCG*100)) +
    geom_point(aes(colour=log2(Astro_rna)), size=2, alpha=0.6) + sc +
    theme_bw() + 
    theme(aspect.ratio=1) +
    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
    scale_y_continuous(breaks=seq(0,50,10), labels=seq(0,50,10), limits=c(0,50)) + 
    geom_text(data = df.stat["Astro", ],
                 aes(x = 80, y=45, label = formatC(pvalue,digits = 2)),color = "red",hjust = 0) +
    geom_text(data = df.stat["Astro", ],
                 aes(x = 80, y=40, label = round(cor,3)),color = "red",hjust = 0) +
    #  geom_smooth(method='lm', formula= y~x)
     stat_smooth(method="lm", se=TRUE, formula=y ~ poly(x, 3, raw=FALSE))

#L6.markers_hmCG.CG <- ggplot(df_L6,aes(x=L6_CG, y=L6_hmCG)) + 
#L6.markers_hmCG.rna <- ggplot(df_L6,aes(x=log2(L6_rna), y=L6_hmCG)) +
#L6.markers_ratio.rna <- ggplot(df_L6,aes(x=log2(L6_rna), y=L6_hmCG.tmCG.ratio)) +
MG.markers.log2rna1_hmCG.CG <- ggplot(MG.markers.log2rna1,aes(x=MG_CG*100, y=MG_hmCG*100)) +
#ggplot(df_L6,aes(x=L6_hmCG.tmCG.ratio, y=L6_mCH*100)) +
#   geom_mask() +
#    geom_point(aes(size=log10GeneLength.kb, colour=log2(L6_rna)), alpha=0.8) + sc +
#    geom_point(aes(size=L6_mCH*100, colour=log2(L6_rna)), alpha=0.5) + sc +
    geom_point(aes(colour=log2(MG_rna)), size=2, alpha=0.6) + sc +
#    gghighlight(L6_rna>2) + 
    theme_bw() + 
    theme(aspect.ratio=1) +
#    scale_size_continuous(range = c(0.5,4)) +
#    scale_size_manual(c(0.5,4)) +
    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
    scale_y_continuous(breaks=seq(0,50,10), labels=seq(0,50,10), limits=c(0,50)) + 
    geom_text(data = df.stat["MG", ],
                 aes(x = 80, y=45, label = formatC(pvalue,digits = 2)),color = "red",hjust = 0) +
    geom_text(data = df.stat["MG", ],
                 aes(x = 80, y=40, label = round(cor,3)),color = "red",hjust = 0) +
  #  geom_smooth(method='lm', formula= y~x)
     stat_smooth(method="lm", se=TRUE, formula=y ~ poly(x, 3, raw=FALSE))



plot_grid(

    Ex.markers.log2rna1_hmCG.CG,   
    Astro.markers.log2rna1_hmCG.CG,
    Inh.markers.log2rna1_hmCG.CG,
    Oligo.markers.log2rna1_hmCG.CG,
    Striatum.markers.log2rna1_hmCG.CG,
    MG.markers.log2rna1_hmCG.CG,
#    L6.markers_ratio.CG,
#    Vip.Ndnf.markers_ratio.CG,
#    L6.markers_mCH.CG,
#    Vip.Ndnf.markers_mCH.CG,
    nrow=3,align = "vh")

ggsave("Fig3_celltype.6g.removeIn_genebody.hmCG.vs.CG_v2_scatterplot.pdf",width = 20 ,height = 20) 

head(pbmc.snmc3@meta.data)

plot_box4_value_nonNeu <- function(genes.in=c("Tle4","Gad2","Gad1","Rorb"),mat1=CG_genebody,yLab="GBmCG",nRow=1){
    mat2 <- mat1[intersect(genes.in,rownames(mat1)),]
    mat.df <- as.data.frame(cbind(t(mat2),as.character(pbmc.snmc3@meta.data$res.name.6groups)))
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
    p1 <- ggplot(df.plot2,aes(cluster,value*100)) + geom_boxplot(width=.6,outlier.shape = NA) + 
          ggrastr::geom_point_rast(position = "jitter",col="grey",width = 0.25,size=0.25) + 
          facet_wrap(~variable,nrow = nRow) +
          stat_summary(fun.y="mean", geom="point", shape=23, size=2,col="red") +
          theme_cowplot() + theme(axis.text.x = element_text(angle=45,hjust=1)) + ylab("Meth rate") + ggtitle(yLab)
    return(p1)
}

table(pbmc.snmc3@meta.data$res.name.6groups)
table(pbmc.snmc3$res.name2)

dim(mC_matrix)
head(mC_matrix)
#dim(C_matrix)


# old analysis:
table(pbmc.snmc3@meta.data$res.name.6groups)

dim(mC_matrix)
head(mC_matrix)

#options(repr.plot.width = 12,repr.plot.height = 3)
gene.list <- c("Arpp21","Lingo1","Sv2b","Neurod6","Kcnma1", # Ex marker
                    "Gad1","Gad2","Abat","Slc6a1","Stat5b","Erbb4", # Inh marker
                    "Ppp1r1b", "Meis2", # Str marker
                    "Slc1a2","Trim9","Paqr8","Cpe","Slc1a3", # Astro marker
                    "Mobp","Mbp","Mog","Pde4b","Klh12", # Oligo marker
                    "Csf1r","P2ry13","Fcrls","Zfp710","Laptm5","Serinc3","Itgb5" # MG marker
                   ) 

#gene.list <- c("Arpp21","Cux2","Rorb","Tle4","Slc6a1","Erbb4","Adarb2","Cacna2d2")
#gene.list2 <- c("Arpp21","Rorb","Tle4","Slc6a1","Erbb4","Adarb2","Gabrg3","Cntnap5c","Il1rapl2")
#pp1 <- plot_box4_value_nonNeu(genes.in=gene.list,mat1=mC_matrix,yLab="BSmCG")
pp2 <- plot_box4_value_nonNeu(genes.in=gene.list,mat1=hmC_matrix,yLab="hmCG", nRow=4)
pp3 <- plot_box4_value_nonNeu(genes.in=gene.list,mat1=tmC_matrix,yLab="tmCG", nRow=4)
#pp4 <- plot_box4_value_nonNeu(genes.in=gene.list,mat1=ratio_matrix,yLab="hmCG/BSmCG", nRow=4)
#pp5 <- plot_box4_value_nonNeu(genes.in=gene.list,mat1=as.matrix(pbmc.snmc2@assays$gbch@counts),yLab="GBmCH")
#pp6 <- plot_box4_value_nonNeu(genes.in=gene.list,mat1=C_matrix,yLab="CG", nRow=4)
#pp7 <- plot_box4_value_nonNeu(genes.in=gene.list,mat1=ratio2_matrix,yLab="hmCG/tmCG")

options(repr.plot.width = 18,repr.plot.height = 24)
plot_grid(pp2,pp3,nrow=2)

# old analysis: 
#options(repr.plot.width = 12,repr.plot.height = 3)
gene.list <- c("Slc17a7","Lingo1","Sv2b","Neurod6","Kcnma1", # Ex marker
                    "Gad1","Gad2","Abat","Slc6a1","Stat5b","Erbb4", # Inh marker
                    "Ppp1r1b", "Meis2", # Str marker
                    "Slc1a2","Trim9","Paqr8","Cpe","Slc1a3", # Astro marker
                    "Mobp","Mbp","Mog","Pde4b","Klh12", # Oligo marker
                    "Csf1r","P2ry13","Fcrls","Zfp710","Laptm5" # MG marker
                   ) 

#gene.list <- c("Arpp21","Cux2","Rorb","Tle4","Slc6a1","Erbb4","Adarb2","Cacna2d2")
#gene.list2 <- c("Arpp21","Rorb","Tle4","Slc6a1","Erbb4","Adarb2","Gabrg3","Cntnap5c","Il1rapl2")
#pp1 <- plot_box4_value_nonNeu(genes.in=gene.list,mat1=mC_matrix,yLab="BSmCG")
pp2 <- plot_box4_value_nonNeu(genes.in=gene.list,mat1=hmC_matrix,yLab="hmCG", nRow=4)
pp3 <- plot_box4_value_nonNeu(genes.in=gene.list,mat1=tmC_matrix,yLab="tmCG", nRow=4)
#pp4 <- plot_box4_value_nonNeu(genes.in=gene.list,mat1=ratio_matrix,yLab="hmCG/BSmCG")
#pp5 <- plot_box4_value_nonNeu(genes.in=gene.list,mat1=as.matrix(pbmc.snmc2@assays$gbch@counts),yLab="GBmCH")
#pp6 <- plot_box4_value_nonNeu(genes.in=gene.list,mat1=C_matrix,yLab="CG")
#pp7 <- plot_box4_value_nonNeu(genes.in=gene.list,mat1=ratio2_matrix,yLab="hmCG/tmCG")

options(repr.plot.width = 9,repr.plot.height = 18)
plot_grid(pp2,pp3,nrow=2)


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



gene.list.final <- c(
                    "Sv2b", # Ex marker, "Lingo1",
                    "Gad2", # Inh marker
                    #"Ppp1r1b","Meis2", # Str marker
                    "Slc1a3", # Astro marker
                    "Mobp", # Oligo marker
                    "Csf1r","Serinc3" # MG marker
                   ) 

#gene.list <- c("Arpp21","Cux2","Rorb","Tle4","Slc6a1","Erbb4","Adarb2","Cacna2d2")
#gene.list2 <- c("Arpp21","Rorb","Tle4","Slc6a1","Erbb4","Adarb2","Gabrg3","Cntnap5c","Il1rapl2")
#pp1 <- plot_box4_value_nonNeu(genes.in=gene.list.final,mat1=mC_matrix,yLab="BSmCG")
pp2 <- plot_box4_value_nonNeu(genes.in=gene.list.final,mat1=hmC_matrix,yLab="hmCG", nRow=1)
pp3 <- plot_box4_value_nonNeu(genes.in=gene.list.final,mat1=tmC_matrix,yLab="tmCG", nRow=1)
pp4 <- plot_box4_value_nonNeu(genes.in=gene.list.final,mat1=ratio_matrix,yLab="hmCG/BSmCG", nRow=1)
pp5 <- plot_box4_value_nonNeu(genes.in=gene.list.final,mat1=as.matrix(pbmc.snmc2@assays$gbch@counts),yLab="GBmCH")
#pp6 <- plot_box4_value_nonNeu(genes.in=gene.list,mat1=C_matrix,yLab="CG")
#pp7 <- plot_box4_value_nonNeu(genes.in=gene.list,mat1=ratio2_matrix,yLab="hmCG/tmCG")

pp8 <- plot_vln2(genes.in=gene.list.final,mat1=as.matrix(expm1(snRNA3@assays$RNA@data)),yLab="RNA")

options(repr.plot.width = 12,repr.plot.height = 16)
plot_grid(pp3,pp2,pp5,pp8, nrow=4)




ggsave("Fig3_celltype.6g.removeIn_rawMeth_marker_genes_boxplot_v2.pdf",width = 7 ,height = 9)

dim(hmC_matrix) # Peng's analysis: 6343 genes x 512 nuclei

gene.order2.rna <- gene.order.by.rna.6g

mat_hmCG_celltype_6g_norm <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name.6groups == x]
    m <- hmC_matrix
    mat <- sweep(m,2,colSums(m),`/`)
    return(rowMeans(mat[gene.order2.rna,cells.in]))
})

mat_tmC_celltype_6g_norm <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name.6groups == x]
    m <- tmC_matrix
    mat <- sweep(m,2,colSums(m),`/`)
    return(rowMeans(mat[gene.order2.rna,cells.in]))
    #return(rowMeans(tmC_matrix[gene.order2.rna,cells.in]))
})

mat_mCG_celltype_6g_norm <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name.6groups == x]
    m <- mC_matrix
    mat <- sweep(m,2,colSums(m),`/`)
    return(rowMeans(mat[gene.order2.rna,cells.in]))
    #return(rowMeans(mC_matrix[gene.order2.rna,cells.in]))
})

mat_mCH_celltype_6g_norm <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name.6groups == x]
    m <- mCH_matrix
    mat <- sweep(m,2,colSums(m),`/`)
    return(rowMeans(mat[gene.order2.rna,cells.in]))
    #return(rowMeans(mCH_matrix[gene.order2.rna,cells.in]))
})

mat_ratio_celltype_6g_norm <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name.6groups == x]
    m <- ratio_matrix
    mat <- sweep(m,2,colSums(m),`/`)
    return(rowMeans(mat[gene.order2.rna,cells.in]))
    #return(rowMeans(ratio_matrix[gene.order2.rna,cells.in]))
})

mat_rna_celltype_6g_norm <- sapply(cluster.order,function(x){
    cells.in <- rownames(snRNA3@meta.data)[snRNA3@meta.data$res.name.6groups == x]
    mat <- expm1(snRNA3@assays$RNA@data)
    return(Matrix::rowMeans(mat[gene.order2.rna,cells.in]))
})

options(repr.plot.width=3,repr.plot.height=3)

h1 <- pheatmap::pheatmap(mat_tmC_celltype_6g_norm,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="tmCG",color=rev(colorRampPalette(brewer.pal(n = 9,name="BuPu"))(50)))

h2 <- pheatmap::pheatmap(mat_mCG_celltype_6g_norm,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="BSmCG",color=rev(colorRampPalette(brewer.pal(n = 9,name="BuPu"))(50)))

h3 <- pheatmap::pheatmap(mat_hmCG_celltype_6g_norm,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="hmCG",color=rev(colorRampPalette(brewer.pal(n = 9,name="BuPu"))(50)))

h4 <- pheatmap::pheatmap(mat_ratio_celltype_6g_norm,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="hmCG/BSmCG",color=colorRampPalette(brewer.pal(n = 9,name="PuRd"))(50))

h5 <- pheatmap::pheatmap(mat_mCH_celltype_6g_norm,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="BSmCH",color=rev(colorRampPalette(brewer.pal(n = 9,name="BuPu"))(50)))

h6 <- pheatmap::pheatmap(mat_rna_celltype_6g_norm,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="RNA",colours = rev(colorRampPalette(brewer.pal(n = 11, name = "RdBu"))(50)))

h6 <- pheatmap::pheatmap(mat_rna_celltype_6g_norm,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="RNA",colours = colorRampPalette(c("blue","white","red"))(50))

# combine all the plots
plot_list=list()
plot_list[["tmCG"]] <- h1[[4]]
plot_list[["BSmCG"]] <- h2[[4]]
plot_list[["hmCG"]] <- h3[[4]]
plot_list[["hmCG/BSmCG"]] <- h4[[4]]
plot_list[["BSmCH"]] <- h5[[4]]
plot_list[["RNA"]] <- h6[[4]]

#library(gridExtra)

options(repr.plot.width = 12,repr.plot.height = 3)
#grid.arrange(grobs=plot_list,ncol = 4)
g2 <- grid.arrange(grobs=plot_list,ncol = 6)

# Peng's analysis
options(repr.plot.width = 12,repr.plot.height = 3)
#grid.arrange(grobs=plot_list,ncol = 4)
g2 <- grid.arrange(grobs=plot_list,ncol = 6)

ggsave("Fig4_6heatmaps_6groups_normalized_1029.pdf",g2,width=12,height=3)

table(Idents(snRNA3))

table(snRNA3@meta.data$res.name.JE3)

Idents(snRNA3) <- "res.name.JE3"

snRNA3.neu <- subset(snRNA3,idents = c("L2/3","L4","L4/5","L6","DL","Pv","Sst","Vip/Ndnf"))

neu.rna.mak <- FindAllMarkers(snRNA3.neu,only.pos = T)

cluster.order.neu <- c("L2/3","L4","L4/5","L6","DL","Vip/Ndnf","Pv","Sst")

neu.rna.mak$cluster <- factor(neu.rna.mak$cluster,levels = cluster.order.neu)
neu.rna.mak2 <- neu.rna.mak %>% filter(p_val_adj < 0.1) %>% arrange(cluster,p_val_adj) 
gene.neu.order <- neu.rna.mak2$gene
length(gene.neu.order) # Peng: 2315 genes; Hao: 3033 genes (w/o pval flitering), 1988 genes (pval.adj<0.1)
head(neu.rna.mak2)

gene.order.by.rna.neu <- intersect(gene.neu.order,rownames(pbmc.snmc3@assays$tmCG@data))

length(gene.order.by.rna.neu) 
# Peng: 989 genes; 
# Hao: 1548 genes (w/o pval flitering), 
#      1101 genes (pval.adj<0.1, 5 CpGs, 80%); 
#      1172 genes (pval.adj<0.1, 3 CpGs, 75%);
#      1202 genes (pval.adj<0.1, 1 CpG, 80%);
#      1179 genes (pval.adj<0.1, 2 CpG, 80%)



length(rownames(pbmc.snmc3@assays$tmCG@data))
table(neu.rna.mak2$cluster) # 1988 genes

pbmc.snmc3@assays

cluster.order <- cluster.order.neu
gene.order <- gene.order.by.rna.neu

mat_hmCG_celltype_neu_raw <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name2 == x]
    return(rowMeans(hmC_matrix[gene.order,cells.in]))
})

mat_tmC_celltype_neu_raw <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name2 == x]
    return(rowMeans(tmC_matrix[gene.order,cells.in]))
})

mat_mCG_celltype_neu_raw <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name2 == x]
    return(rowMeans(mC_matrix[gene.order,cells.in]))
})

mat_mCH_celltype_neu_raw <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name2 == x]
    return(rowMeans(mCH_matrix[gene.order,cells.in]))
})

mat_ratio_celltype_neu_raw <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name2 == x]
    return(rowMeans(ratio_matrix[gene.order,cells.in]))
})

mat_rna_celltype_neu_raw <- sapply(cluster.order,function(x){
    cells.in <- rownames(snRNA3.neu@meta.data)[snRNA3.neu@meta.data$res.name.JE3 == x]
    mat <- expm1(snRNA3.neu@assays$RNA@data)
    return(Matrix::rowMeans(mat[gene.order,cells.in]))
})

dim(mat_mCG_celltype_neu_raw)
head(mat_mCG_celltype_neu_raw)

mat_CG_celltype_neu_raw <- 1 - mat_mCG_celltype_neu_raw
dim(mat_CG_celltype_neu_raw)
head(mat_CG_celltype_neu_raw)

options(repr.plot.width=3,repr.plot.height=3)

h0 <- pheatmap::pheatmap(mat_CG_celltype_neu_raw,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="CG",color=rev(colorRampPalette(brewer.pal(n = 9,name="BuPu"))(50)))

h1 <- pheatmap::pheatmap(mat_tmC_celltype_neu_raw,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="tmCG",color=rev(colorRampPalette(brewer.pal(n = 9,name="BuPu"))(50)))

h2 <- pheatmap::pheatmap(mat_mCG_celltype_neu_raw,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="BSmCG",color=rev(colorRampPalette(brewer.pal(n = 9,name="BuPu"))(50)))

h3 <- pheatmap::pheatmap(mat_hmCG_celltype_neu_raw,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="hmCG",color=rev(colorRampPalette(brewer.pal(n = 9,name="BuPu"))(50)))

h4 <- pheatmap::pheatmap(mat_ratio_celltype_neu_raw,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="hmCG/BSmCG",color=colorRampPalette(brewer.pal(n = 9,name="PuRd"))(50))

h5 <- pheatmap::pheatmap(mat_mCH_celltype_neu_raw,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="BSmCH",color=rev(colorRampPalette(brewer.pal(n = 9,name="BuPu"))(50)))

h6 <- pheatmap::pheatmap(mat_rna_celltype_neu_raw,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="RNA",colours = rev(colorRampPalette(brewer.pal(n = 11, name = "RdBu"))(50)))

h7 <- pheatmap::pheatmap(mat_rna_celltype_neu_raw,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="RNA",colours = colorRampPalette(c("blue","white","red"))(50))

# 6-16-2022 final analysis:  1179 genes 
# combine all the plots
plot_list=list()
plot_list[["CG"]] <- h0[[4]]
plot_list[["tmCG"]] <- h1[[4]]
#plot_list[["BSmCG"]] <- h2[[4]]
plot_list[["hmCG"]] <- h3[[4]]
#plot_list[["hmCG/BSmCG"]] <- h4[[4]]
plot_list[["BSmCH"]] <- h5[[4]]
plot_list[["RNA"]] <- h6[[4]]

library(gridExtra)

options(repr.plot.width = 12,repr.plot.height = 4)
#grid.arrange(grobs=plot_list,ncol = 4)
g2 <- grid.arrange(grobs=plot_list,ncol = 5)

# old analysis:  1101 genes 
# combine all the plots
plot_list=list()
plot_list[["CG"]] <- h0[[4]]
plot_list[["tmCG"]] <- h1[[4]]
#plot_list[["BSmCG"]] <- h2[[4]]
plot_list[["hmCG"]] <- h3[[4]]
#plot_list[["hmCG/BSmCG"]] <- h4[[4]]
plot_list[["BSmCH"]] <- h5[[4]]
plot_list[["RNA"]] <- h6[[4]]

library(gridExtra)

options(repr.plot.width = 12,repr.plot.height = 4)
#grid.arrange(grobs=plot_list,ncol = 4)
g2 <- grid.arrange(grobs=plot_list,ncol = 5)

ggsave("HW-Fig5_6heatmaps_NeuClusters_06.16.2022_1179genes.pdf",g2,width=12,height=4)

ggsave("HW-Fig4_6heatmaps_NeuClusters_03.17.2022.pdf",g2,width=12,height=4)

ggsave("HW-Fig4_6heatmaps_NeuClusters_01.02.2022.pdf",g2,width=12,height=4)

# neu.rna.mak2 <- neu.rna.mak %>% filter(p_val_adj < 0.1) %>% arrange(cluster,p_val_adj) 
#gene.neu.order <- neu.rna.mak2$gene
#length(gene.neu.order) # Peng: 2315 genes; Hao: 3033 genes (w/o pval flitering), 1988 genes (pval.adj<0.1)
dim(neu.rna.mak2)
dim(df_4modalities_celltype_all.genes)
table(neu.rna.mak2$cluster)
head(neu.rna.mak2)
head(df_4modalities_celltype_all.genes)

# old analysis:
# neu.rna.mak2 <- neu.rna.mak %>% filter(p_val_adj < 0.1) %>% arrange(cluster,p_val_adj) 
#gene.neu.order <- neu.rna.mak2$gene
#length(gene.neu.order) # Peng: 2315 genes; Hao: 3033 genes (w/o pval flitering), 1988 genes (pval.adj<0.1)
dim(neu.rna.mak2)
dim(df_4modalities_celltype_all.genes)
table(neu.rna.mak2$cluster)
head(neu.rna.mak2)
head(df_4modalities_celltype_all.genes)

length(unique(neu.rna.mak2$gene)) # 1988 genes; unique 1378 genes
length(intersect(neu.rna.mak2$gene,rownames(df_4modalities_celltype_all.genes)))
#neu.rna.mak2$gene
#neu.rna.mak2 %>% filter(cluster=="L6") %>% pull(gene)

L2.3_rna.specific.gene<-intersect(neu.rna.mak2 %>% filter(cluster=="L2/3") %>% pull(gene),rownames(df_4modalities_celltype_all.genes))
length(L2.3_rna.specific.gene)
L6_rna.specific.gene<-intersect(neu.rna.mak2 %>% filter(cluster=="L6") %>% pull(gene),rownames(df_4modalities_celltype_all.genes))
length(L6_rna.specific.gene)
Vip.Ndnf_rna.specific.gene<-intersect(neu.rna.mak2 %>% filter(cluster=="Vip/Ndnf") %>% pull(gene),rownames(df_4modalities_celltype_all.genes))
length(Vip.Ndnf_rna.specific.gene)
Sst_rna.specific.gene<-intersect(neu.rna.mak2 %>% filter(cluster=="Sst") %>% pull(gene),rownames(df_4modalities_celltype_all.genes))
length(Sst_rna.specific.gene)
Pv_rna.specific.gene<-intersect(neu.rna.mak2 %>% filter(cluster=="Pv") %>% pull(gene),rownames(df_4modalities_celltype_all.genes))
length(Pv_rna.specific.gene)

L2.3_rna.specific.gene
Pv_rna.specific.gene
intersect(L2.3_rna.specific.gene, Pv_rna.specific.gene)

summary(df_4modalities_celltype_all.genes[L6_rna.specific.gene, c("L6_rna","Pv_rna","Sst_rna","Vip.Ndnf_rna")])
summary(df_4modalities_celltype_all.genes[Pv_rna.specific.gene, c("L6_rna","Pv_rna","Sst_rna","Vip.Ndnf_rna")])
summary(df_4modalities_celltype_all.genes[Sst_rna.specific.gene, c("L6_rna","Pv_rna","Sst_rna","Vip.Ndnf_rna")])
summary(df_4modalities_celltype_all.genes[Vip.Ndnf_rna.specific.gene, c("L6_rna","Pv_rna","Sst_rna","Vip.Ndnf_rna")])

length(gene.order.by.hmCG.neu)
gene.order.by.hmCG.neu

cluster.order <- cluster.order.neu
gene.order <- gene.order.by.hmCG.neu

mat_hmCG_celltype_neu_raw <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name2 == x]
    return(rowMeans(hmC_matrix[gene.order,cells.in]))
})

mat_tmC_celltype_neu_raw <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name2 == x]
    return(rowMeans(tmC_matrix[gene.order,cells.in]))
})

mat_mCG_celltype_neu_raw <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name2 == x]
    return(rowMeans(mC_matrix[gene.order,cells.in]))
})

mat_mCH_celltype_neu_raw <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name2 == x]
    return(rowMeans(mCH_matrix[gene.order,cells.in]))
})

mat_ratio_celltype_neu_raw <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name2 == x]
    return(rowMeans(ratio_matrix[gene.order,cells.in]))
})

mat_rna_celltype_neu_raw <- sapply(cluster.order,function(x){
    cells.in <- rownames(snRNA3.neu@meta.data)[snRNA3.neu@meta.data$res.name.JE3 == x]
    mat <- expm1(snRNA3.neu@assays$RNA@data)
    return(Matrix::rowMeans(mat[gene.order,cells.in]))
})

mat_CG_celltype_neu_raw <- 1 - mat_mCG_celltype_neu_raw

options(repr.plot.width=3,repr.plot.height=3)

h0 <- pheatmap::pheatmap(mat_CG_celltype_neu_raw,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="CG",color=rev(colorRampPalette(brewer.pal(n = 9,name="BuPu"))(50)))

h1 <- pheatmap::pheatmap(mat_tmC_celltype_neu_raw,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="tmCG",color=rev(colorRampPalette(brewer.pal(n = 9,name="BuPu"))(50)))

h2 <- pheatmap::pheatmap(mat_mCG_celltype_neu_raw,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="BSmCG",color=rev(colorRampPalette(brewer.pal(n = 9,name="BuPu"))(50)))

h3 <- pheatmap::pheatmap(mat_hmCG_celltype_neu_raw,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="hmCG",color=rev(colorRampPalette(brewer.pal(n = 9,name="BuPu"))(50)))

h4 <- pheatmap::pheatmap(mat_ratio_celltype_neu_raw,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="hmCG/BSmCG",color=colorRampPalette(brewer.pal(n = 9,name="PuRd"))(50))

h5 <- pheatmap::pheatmap(mat_mCH_celltype_neu_raw,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="BSmCH",color=rev(colorRampPalette(brewer.pal(n = 9,name="BuPu"))(50)))

h6 <- pheatmap::pheatmap(mat_rna_celltype_neu_raw,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="RNA",colours = rev(colorRampPalette(brewer.pal(n = 11, name = "RdBu"))(50)))


Genebody.info<-read.table("regions.all_genes.bed",sep="\t",quote="\"",fill=TRUE, col.names=c("chr","start","end","strand","id1","id2","GeneName","Transcript_type"))
head(Genebody.info)
Genebody.info %>% filter(GeneName=="Zc3h11a")
length(Genebody.info$GeneName)
length(unique(Genebody.info$GeneName))

Genebody.info.uniq <- Genebody.info %>% distinct(GeneName, .keep_all = TRUE) # remove all non-unique entries, e.g. Zc3h11a
dim(Genebody.info.uniq)

df_5modalities_celltype_all.genes$GeneName<-rownames(df_5modalities_celltype_all.genes)
dim(df_5modalities_celltype_all.genes)
head(df_5modalities_celltype_all.genes)

df_5modalities_celltype_all.genes_genebody <- inner_join(df_5modalities_celltype_all.genes, Genebody.info.uniq, by=c("GeneName"="GeneName"))
dim(df_5modalities_celltype_all.genes)
dim(df_5modalities_celltype_all.genes_genebody)
rownames(df_5modalities_celltype_all.genes_genebody)<-df_5modalities_celltype_all.genes_genebody$GeneName
df_5modalities_celltype_all.genes_genelength <- df_5modalities_celltype_all.genes_genebody %>% mutate(GeneLength = end - start) %>% mutate(log10GeneLength.kb = log10(GeneLength/1000))
head(df_5modalities_celltype_all.genes_genelength)
summary(df_5modalities_celltype_all.genes_genelength$log10GeneLength.kb)

Genebody.info.uniq_genelength <- Genebody.info.uniq %>% mutate(GeneLength = end - start) %>% mutate(log10GeneLength.kb = log10(GeneLength/1000))
table(Genebody.info.uniq_genelength$Transcript_type)

Genebody.info.uniq_genelength %>% arrange(desc(log10GeneLength.kb)) %>% filter(Transcript_type == "protein_coding")

# Histogram with density plot
# protein_coding
ggplot(Genebody.info.uniq_genelength %>% filter(Transcript_type == "protein_coding"), aes(x=GeneLength)) + 
 geom_histogram(aes(y=..density..), colour="black", fill="white")+
 geom_density(alpha=.2, fill="#FF6666") 

ggplot(Genebody.info.uniq_genelength %>% filter(Transcript_type == "protein_coding"), aes(x=log10GeneLength.kb)) + 
 geom_histogram(aes(y=..density..), colour="black", fill="white")+
 geom_density(alpha=.2, fill="#FF6666") 

# protein_coding + lincRNA
ggplot(Genebody.info.uniq_genelength %>% filter(Transcript_type %in% c("protein_coding", "lincRNA")), aes(x=log10GeneLength.kb, color=Transcript_type, fill=Transcript_type)) +
geom_histogram(aes(y=..density..), position="identity", alpha=0.5)+
geom_density(alpha=0.6)+
#geom_vline(data=Genebody.info.uniq_genelength, aes(xintercept=mean(log10GeneLength.kb), color=Transcript_type),linetype="dashed")+
scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
#labs(title="Weight histogram plot",x="Weight(kg)", y = "Density")+
theme_classic()

# Histogram with density plot
ggplot(df_5modalities_celltype_all.genes_genelength, aes(x=GeneLength)) + 
 geom_histogram(aes(y=..density..), colour="black", fill="white")+
 geom_density(alpha=.2, fill="#FF6666") 

ggplot(df_5modalities_celltype_all.genes_genelength, aes(x=log10GeneLength.kb)) + 
 geom_histogram(aes(y=..density..), colour="black", fill="white")+
 geom_density(alpha=.2, fill="#FF6666") 

df_5modalities_celltype_all.genes_genebody <- inner_join(df_5modalities_celltype_all.genes, Genebody.info.uniq, by=c("GeneName"="GeneName"))
dim(df_5modalities_celltype_all.genes)
dim(df_5modalities_celltype_all.genes_genebody)
rownames(df_5modalities_celltype_all.genes_genebody)<-df_5modalities_celltype_all.genes_genebody$GeneName
df_5modalities_celltype_all.genes_genelength <- df_5modalities_celltype_all.genes_genebody %>% mutate(GeneLength = end - start) %>% mutate(log10GeneLength.kb = log10(GeneLength/1000))
head(df_5modalities_celltype_all.genes_genelength)
summary(df_5modalities_celltype_all.genes_genelength$log10GeneLength.kb)

# Old analysis:  
df_5modalities_celltype_all.genes_genebody <- inner_join(df_5modalities_celltype_all.genes, Genebody.info.uniq, by=c("GeneName"="GeneName"))
dim(df_5modalities_celltype_all.genes)
dim(df_5modalities_celltype_all.genes_genebody)
rownames(df_5modalities_celltype_all.genes_genebody)<-df_5modalities_celltype_all.genes_genebody$GeneName
df_5modalities_celltype_all.genes_genelength <- df_5modalities_celltype_all.genes_genebody %>% mutate(GeneLength = end - start) %>% mutate(log10GeneLength.kb = log10(GeneLength/1000))
head(df_5modalities_celltype_all.genes_genelength)
summary(df_5modalities_celltype_all.genes_genelength$log10GeneLength.kb)

df_5modalities_celltype_all.genes_genelength["Arpp21", ] %>% select(L2.3_CG, L2.3_tmCG, L2.3_hmCG, L2.3_rna)
df_L2.3 <- df_5modalities_celltype_all.genes_genelength[L2.3_rna.specific.gene, ] %>%  mutate(L2.3_hmCG.tmCG.ratio = L2.3_hmCG / L2.3_tmCG) %>% arrange(desc(L2.3_CG)) %>% filter(log2(L2.3_rna)>0) %>% select(GeneLength, log10GeneLength.kb, L2.3_CG, L2.3_tmCG, L2.3_hmCG, L2.3_mCH, L2.3_hmCG.tmCG.ratio, L2.3_rna,L4_rna,L4.5_rna,L6_rna,DL_rna,Pv_rna,Vip.Ndnf_rna, Sst_rna) # sort L2.3_rna
length(L2.3_rna.specific.gene)
dim(df_L2.3)
df_L2.3[c("Cux2","Cux1"), ]
df_L2.3

df_5modalities_celltype_all.genes_genelength["Arpp21", ] %>% select(L6_CG, L6_tmCG, L6_hmCG, L6_rna)
df_L6 <- df_5modalities_celltype_all.genes_genelength[L6_rna.specific.gene, ] %>%  mutate(L6_hmCG.tmCG.ratio = L6_hmCG / L6_tmCG) %>% arrange(desc(L6_CG)) %>% filter(log2(L6_rna)>0) %>% select(GeneLength, log10GeneLength.kb, L6_CG, L6_tmCG, L6_hmCG, L6_mCH, L6_hmCG.tmCG.ratio, L2.3_rna,L4_rna,L4.5_rna,L6_rna,DL_rna,Pv_rna,Vip.Ndnf_rna, Sst_rna) # sort L6_rna
length(L6_rna.specific.gene)
dim(df_L6)
df_L6["Tle4", ]
df_L6

df_5modalities_celltype_all.genes_genelength["Slc6a1", ] %>% select(Vip.Ndnf_CG, Vip.Ndnf_tmCG, Vip.Ndnf_hmCG, Vip.Ndnf_rna)
df_Vip.Ndnf <- df_5modalities_celltype_all.genes_genelength[Vip.Ndnf_rna.specific.gene, ] %>%  mutate(Vip.Ndnf_hmCG.tmCG.ratio = Vip.Ndnf_hmCG / Vip.Ndnf_tmCG)%>% arrange(desc(Vip.Ndnf_CG)) %>% filter(log2(Vip.Ndnf_rna)>0) %>% select(GeneLength, log10GeneLength.kb, Vip.Ndnf_CG, Vip.Ndnf_tmCG, Vip.Ndnf_hmCG, Vip.Ndnf_mCH, Vip.Ndnf_hmCG.tmCG.ratio, L2.3_rna,L4_rna,L4.5_rna,L6_rna,DL_rna,Pv_rna,Vip.Ndnf_rna, Sst_rna) # sort Vip.Ndnf_rna
length(Vip.Ndnf_rna.specific.gene)
dim(df_Vip.Ndnf)
df_Vip.Ndnf["Adarb2", ]
df_Vip.Ndnf

df_5modalities_celltype_all.genes_genelength["Slc6a1", ] %>% select(Pv_CG, Pv_tmCG, Pv_hmCG, Pv_rna)
df_Pv <- df_5modalities_celltype_all.genes_genelength[Pv_rna.specific.gene, ] %>%  mutate(Pv_hmCG.tmCG.ratio = Pv_hmCG / Pv_tmCG)%>% arrange(desc(Pv_CG)) %>% filter(log2(Pv_rna)>0) %>% select(GeneLength, log10GeneLength.kb, Pv_CG, Pv_tmCG, Pv_hmCG, Pv_mCH, Pv_hmCG.tmCG.ratio, L2.3_rna,L4_rna,L4.5_rna,L6_rna,DL_rna,Pv_rna,Vip.Ndnf_rna, Sst_rna) # sort Pv_rna
length(Pv_rna.specific.gene)
dim(df_Pv)
df_Pv["Erbb4", ]
head(df_Pv)

options(repr.plot.width = 20,repr.plot.height = 10)

myPalette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
#myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(-8, 8))

#L6.all_rna.CG <- ggplot(df_4modalities_celltype_all.genes,aes(x=log1p(L6_rna), y=L6_CG)) + 
#L6.all_rna.CG <- ggplot(df_4modalities_celltype_all.genes,aes(x=log2(L6_rna), y=L6_CG)) + 
L6.all_rna.hmCG <- ggplot(df_5modalities_celltype_all.genes,aes(x=log2(L6_rna), y=L6_hmCG*100)) + 
#   geom_mask() +
    geom_point(aes(colour=log2(L6_rna)), size=0.5, alpha=0.5) + sc +
#    gghighlight(L6_rna>2) + 
    theme_bw() + 
    theme(aspect.ratio=1) +
    geom_smooth(method='lm', formula= y~x)

L6.all_rna.CG <- ggplot(df_5modalities_celltype_all.genes,aes(x=log2(L6_rna), y=L6_CG*100)) + 
#L6.all_rna.hmCG <- ggplot(df_4modalities_celltype_all.genes,aes(x=log2(L6_rna), y=L6_hmCG*100)) + 
#   geom_mask() +
    geom_point(aes(colour=log2(L6_rna)), size=0.5, alpha=0.5) + sc +
#    gghighlight(L6_rna>2) + 
    theme_bw() + 
    theme(aspect.ratio=1) +
    geom_smooth(method='lm', formula= y~x)

L6.all_rna.tmCG <- ggplot(df_5modalities_celltype_all.genes,aes(x=log2(L6_rna), y=L6_tmCG*100)) + 
#L6.all_rna.hmCG <- ggplot(df_4modalities_celltype_all.genes,aes(x=log2(L6_rna), y=L6_hmCG*100)) + 
#   geom_mask() +
    geom_point(aes(colour=log2(L6_rna)), size=0.5, alpha=0.5) + sc +
#    gghighlight(L6_rna>2) + 
    theme_bw() + 
    theme(aspect.ratio=1) +
    geom_smooth(method='lm', formula= y~x)

L6.all_rna.mCH <- ggplot(df_5modalities_celltype_all.genes,aes(x=log2(L6_rna), y=L6_mCH)) +
#L6.all_rna.ratio <- ggplot(df_5modalities_celltype_all.genes,aes(x=log2(L6_rna), y=L6_hmCG/L6_tmCG)) + 
#   geom_mask() +
    geom_point(aes(colour=log2(L6_rna)), size=0.5, alpha=0.5) + sc +
#    gghighlight(L6_rna>2) + 
    theme_bw() + 
    theme(aspect.ratio=1) +
    geom_smooth(method='lm', formula= y~x)

L6.all_rna.hmCG <- ggplot(df_5modalities_celltype_all.genes,aes(x=log2(L6_rna), y=L6_hmCG*100)) + 
#   geom_mask() +
    geom_point(aes(colour=log2(L6_rna)), size=0.5, alpha=0.5) + sc +
#    gghighlight(L6_rna>2) + 
    theme_bw() + 
    theme(aspect.ratio=1) +
    geom_smooth(method='lm', formula= y~x)

L6.all_rna.CG <- ggplot(df_5modalities_celltype_all.genes,aes(x=log2(L6_rna), y=L6_CG*100)) + 
#L6.all_rna.hmCG <- ggplot(df_4modalities_celltype_all.genes,aes(x=log2(L6_rna), y=L6_hmCG*100)) + 
#   geom_mask() +
    geom_point(aes(colour=log2(L6_rna)), size=0.5, alpha=0.5) + sc +
#    gghighlight(L6_rna>2) + 
    theme_bw() + 
    theme(aspect.ratio=1) +
    geom_smooth(method='lm', formula= y~x)

L6.all_rna.tmCG <- ggplot(df_5modalities_celltype_all.genes,aes(x=log2(L6_rna), y=L6_tmCG*100)) + 
#L6.all_rna.hmCG <- ggplot(df_4modalities_celltype_all.genes,aes(x=log2(L6_rna), y=L6_hmCG*100)) + 
#   geom_mask() +
    geom_point(aes(colour=log2(L6_rna)), size=0.5, alpha=0.5) + sc +
#    gghighlight(L6_rna>2) + 
    theme_bw() + 
    theme(aspect.ratio=1) +
    geom_smooth(method='lm', formula= y~x)

L6.all_rna.mCH <- ggplot(df_5modalities_celltype_all.genes,aes(x=log2(L6_rna), y=L6_mCH)) +
#L6.all_rna.ratio <- ggplot(df_5modalities_celltype_all.genes,aes(x=log2(L6_rna), y=L6_hmCG/L6_tmCG)) + 
#   geom_mask() +
    geom_point(aes(colour=log2(L6_rna)), size=0.5, alpha=0.5) + sc +
#    gghighlight(L6_rna>2) + 
    theme_bw() + 
    theme(aspect.ratio=1) +
    geom_smooth(method='lm', formula= y~x)

#L6.markers_hmCG.CG <- ggplot(df_L6,aes(x=L6_CG, y=L6_hmCG)) + 
#L6.markers_hmCG.rna <- ggplot(df_L6,aes(x=log2(L6_rna), y=L6_hmCG)) +
#L6.markers_ratio.rna <- ggplot(df_L6,aes(x=log2(L6_rna), y=L6_hmCG.tmCG.ratio)) +
L6.markers_ratio.CG <- ggplot(df_L6,aes(x=L6_CG*100, y=L6_hmCG.tmCG.ratio)) +
#   geom_mask() +
    geom_point(aes(colour=L6_mCH*100), size=0.5, alpha=0.5) + sc +
#    gghighlight(L6_rna>2) + 
    theme_bw() + 
    theme(aspect.ratio=1) +
    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
    geom_smooth(method='lm', formula= y~x)

#L6.markers_ratio2.CG <- ggplot(df_L6,aes(x=1/(L6_hmCG+L6_tmCG), y=L6_CG)) + 
L6.markers_mCH.tmCG <- ggplot(df_L6,aes(x=L6_mCH, y=L6_tmCG)) + 
#   geom_mask() +
    geom_point(aes(colour=log2(L6_rna)), size=0.5, alpha=0.5) + sc +
#    gghighlight(L6_rna>2) + 
    theme_bw() + 
    theme(aspect.ratio=1) +
    geom_smooth(method='lm', formula= y~x)


L6.markers_hmCG.rna <- ggplot(df_L6,aes(x=log2(L6_rna), y=L6_hmCG)) +
#L6.markers_tmCG.rna <- ggplot(df_L6,aes(x=log2(L6_rna), y=L6_tmCG)) +
#L6.markers_CG.rna <- ggplot(df_L6,aes(x=log2(L6_rna), y=L6_CG)) +
#L6.markers_ratio.rna <- ggplot(df_L6,aes(x=log2(L6_rna), y=L6_hmCG.tmCG.ratio)) +
#   geom_mask() +
    geom_point(aes(colour=log2(L6_rna)), size=0.5, alpha=0.5) + sc +
#    gghighlight(L6_rna>2) + 
    theme_bw() + 
    theme(aspect.ratio=1) +
#    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
    geom_smooth(method='lm', formula= y~x)

#L6.markers_hmCG.rna <- ggplot(df_L6,aes(x=log2(L6_rna), y=L6_hmCG)) +
L6.markers_tmCG.rna <- ggplot(df_L6,aes(x=log2(L6_rna), y=L6_tmCG)) +
#L6.markers_CG.rna <- ggplot(df_L6,aes(x=log2(L6_rna), y=L6_CG)) +
#L6.markers_ratio.rna <- ggplot(df_L6,aes(x=log2(L6_rna), y=L6_hmCG.tmCG.ratio)) +
#   geom_mask() +
    geom_point(aes(colour=log2(L6_rna)), size=0.5, alpha=0.5) + sc +
#    gghighlight(L6_rna>2) + 
    theme_bw() + 
    theme(aspect.ratio=1) +
#    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
    geom_smooth(method='lm', formula= y~x)

#L6.markers_hmCG.rna <- ggplot(df_L6,aes(x=log2(L6_rna), y=L6_hmCG)) +
#L6.markers_tmCG.rna <- ggplot(df_L6,aes(x=log2(L6_rna), y=L6_tmCG)) +
L6.markers_CG.rna <- ggplot(df_L6,aes(x=log2(L6_rna), y=L6_CG)) +
#L6.markers_ratio.rna <- ggplot(df_L6,aes(x=log2(L6_rna), y=L6_hmCG.tmCG.ratio)) +
#   geom_mask() +
    geom_point(aes(colour=log2(L6_rna)), size=0.5, alpha=0.5) + sc +
#    gghighlight(L6_rna>2) + 
    theme_bw() + 
    theme(aspect.ratio=1) +
#    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
    geom_smooth(method='lm', formula= y~x)

#L6.markers_hmCG.rna <- ggplot(df_L6,aes(x=log2(L6_rna), y=L6_hmCG)) +
#L6.markers_tmCG.rna <- ggplot(df_L6,aes(x=log2(L6_rna), y=L6_tmCG)) +
L6.markers_CG.rna <- ggplot(df_L6,aes(x=log2(L6_rna), y=L6_CG)) +
#L6.markers_ratio.rna <- ggplot(df_L6,aes(x=log2(L6_rna), y=L6_hmCG.tmCG.ratio)) +
#   geom_mask() +
    geom_point(aes(colour=log2(L6_rna)), size=0.5, alpha=0.5) + sc +
#    gghighlight(L6_rna>2) + 
    theme_bw() + 
    theme(aspect.ratio=1) +
#    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
    geom_smooth(method='lm', formula= y~x)

#L6.markers_hmCG.rna <- ggplot(df_L6,aes(x=log2(L6_rna), y=L6_hmCG)) +
#L6.markers_tmCG.rna <- ggplot(df_L6,aes(x=log2(L6_rna), y=L6_tmCG)) +
#L6.markers_CG.rna <- ggplot(df_L6,aes(x=log2(L6_rna), y=L6_CG)) +
#L6.markers_ratio.rna <- ggplot(df_L6,aes(x=log2(L6_rna), y=L6_hmCG.tmCG.ratio)) +
L6.markers_mCH.rna <- ggplot(df_L6,aes(x=log2(L6_rna), y=L6_mCH)) +
#   geom_mask() +
    geom_point(aes(colour=log2(L6_rna)), size=0.5, alpha=0.5) + sc +
#    gghighlight(L6_rna>2) + 
    theme_bw() + 
    theme(aspect.ratio=1) +
#    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
    geom_smooth(method='lm', formula= y~x)

Vip.Ndnf.markers_ratio.CG <- ggplot(df_Vip.Ndnf,aes(x=Vip.Ndnf_hmCG.tmCG.ratio, y=Vip.Ndnf_CG)) + 
    geom_point() + 
    gghighlight(Vip.Ndnf_rna>2) + 
    theme_bw() + 
    theme(aspect.ratio=1) +
    geom_smooth(method='lm', formula= y~x)

plot_grid(
#    L6.all_rna.CG, 
     L6.all_rna.hmCG,
     L6.all_rna.tmCG,
     L6.all_rna.CG,
     L6.all_rna.mCH,
#     L6.all_rna.ratio,
#    L6.markers_ratio.rna,
#    L6.markers_hmCG.rna,
#    L6.markers_mCH.tmCG,
#    L6.markers_ratio.CG,
#    L6.markers_ratio2.CG, 
    L6.markers_hmCG.rna,
    L6.markers_tmCG.rna,
    L6.markers_CG.rna,
    L6.markers_mCH.rna,
#    L6.markers_ratio.rna,
#    Vip.Ndnf.markers_ratio.CG,
    nrow=2,align = "vh")

dim(df_L6)
head(df_L6)

options(repr.plot.width = 15,repr.plot.height = 12)

myPalette <- colorRampPalette(rev(brewer.pal(11, "RdBu"))) 
#myPalette <- colorRampPalette(brewer.pal(11, "YlOrRd")) 
#myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(1,3.75))
#sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(0.5,4))

#L6.markers_hmCG.CG <- ggplot(df_L6 %>% filter(log2(L6_rna)>1 & L6_CG*100 < 30),aes(x=L6_CG*100, y=L6_hmCG*100)) +
L6.markers_hmCG.CG <- ggplot(df_L6 %>% filter(log2(L6_rna)>1),aes(x=L6_CG*100, y=L6_hmCG*100)) + 
#L6.markers_hmCG.rna <- ggplot(df_L6,aes(x=log2(L6_rna), y=L6_hmCG)) +
#L6.markers_ratio.rna <- ggplot(df_L6,aes(x=log2(L6_rna), y=L6_hmCG.tmCG.ratio)) +
#L6.markers_ratio.CG <- ggplot(df_L6 %>% filter(log2(L6_rna)>1),aes(x=L6_CG*100, y=L6_hmCG.tmCG.ratio)) +
#ggplot(df_L6,aes(x=L6_hmCG.tmCG.ratio, y=L6_mCH*100)) +
#   geom_mask() +
    geom_point(aes(size=log2(L6_rna), colour=log10GeneLength.kb), alpha=0.8) + sc +
#    geom_point(aes(size=L6_mCH*100, colour=log2(L6_rna)), alpha=0.5) + sc +
#    geom_point(aes(colour=log10GeneLength.kb), size=2, alpha=0.8) + sc +
#    gghighlight(L6_rna>2) + 
    theme_bw() + 
    theme(aspect.ratio=1) +
    scale_size_continuous(range = c(0.5,4)) +
#    scale_size_manual(c(0.5,4)) +
    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
    scale_y_continuous(breaks=seq(0,50,10), labels=seq(0,50,10), limits=c(0,50)) + 
    stat_smooth(method="lm", se=TRUE, formula=y ~ poly(x, 3, raw=FALSE))




#L6.markers_hmCG.CG <- ggplot(df_L6,aes(x=L6_CG, y=L6_hmCG)) + 
#L6.markers_hmCG.rna <- ggplot(df_L6,aes(x=log2(L6_rna), y=L6_hmCG)) +
#L6.markers_ratio.rna <- ggplot(df_L6,aes(x=log2(L6_rna), y=L6_hmCG.tmCG.ratio)) +
L6.markers_tmCG.CG <- ggplot(df_L6 %>% filter(log2(L6_rna)>1),aes(x=L6_CG*100, y=L6_tmCG*100)) +
#ggplot(df_L6,aes(x=L6_hmCG.tmCG.ratio, y=L6_mCH*100)) +
#   geom_mask() +
#    geom_point(aes(size=log10GeneLength.kb, colour=log2(L6_rna)), alpha=0.8) + sc +
    geom_point(aes(size=log2(L6_rna), colour=log10GeneLength.kb), alpha=0.8) + sc +
#    geom_point(aes(size=L6_mCH*100, colour=log2(L6_rna)), alpha=0.5) + sc +
#    geom_point(aes(colour=log10GeneLength.kb), size=2, alpha=0.8) + sc +
#    gghighlight(L6_rna>2) + 
    theme_bw() + 
    theme(aspect.ratio=1) +
    scale_size_continuous(range = c(0.5,4)) +
#    scale_size_manual(c(0.5,4)) +
    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
#    scale_y_continuous(breaks=seq(0,4,1), labels=seq(0,4,1), limits=c(0,4)) + 
    geom_smooth(method='lm', formula= y~x)

#L6.markers_hmCG.CG <- ggplot(df_L6,aes(x=L6_CG, y=L6_hmCG)) + 
#L6.markers_hmCG.rna <- ggplot(df_L6,aes(x=log2(L6_rna), y=L6_hmCG)) +
#L6.markers_ratio.rna <- ggplot(df_L6,aes(x=log2(L6_rna), y=L6_hmCG.tmCG.ratio)) +
L6.markers_mCH.CG <- ggplot(df_L6 %>% filter(log2(L6_rna)>1),aes(x=L6_CG*100, y=L6_mCH*100)) +
#ggplot(df_L6,aes(x=L6_hmCG.tmCG.ratio, y=L6_mCH*100)) +
#   geom_mask() +
#    geom_point(aes(size=log10GeneLength.kb, colour=log2(L6_rna)), alpha=0.8) + sc +
    geom_point(aes(size=log2(L6_rna), colour=log10GeneLength.kb), alpha=0.8) + sc +
#    geom_point(aes(size=L6_mCH*100, colour=log2(L6_rna)), alpha=0.5) + sc +
#    geom_point(aes(colour=log10GeneLength.kb), size=2, alpha=0.8) + sc +
#    gghighlight(L6_rna>2) + 
    theme_bw() + 
    theme(aspect.ratio=1) +
    scale_size_continuous(range = c(0.5,4)) +
#    scale_size_manual(c(0.5,4)) +
    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
#    scale_y_continuous(breaks=seq(0,4,1), labels=seq(0,4,1), limits=c(0,4)) + 
    geom_smooth(method='lm', formula= y~x)


Pv.markers_hmCG.CG <- ggplot(df_Pv %>% filter(log2(Pv_rna)>1),aes(x=Pv_CG*100, y=Pv_hmCG*100)) + 
#Vip.Ndnf.markers_hmCG.CG <- ggplot(df_Vip.Ndnf %>% filter(log2(Vip.Ndnf_rna)>1),aes(x=Vip.Ndnf_CG*100, y=Vip.Ndnf_hmCG*100)) + 
#    gghighlight(Vip.Ndnf_rna>2) + 
#   geom_mask() +
    geom_point(aes(size=log2(Pv_rna), colour=log10GeneLength.kb), alpha=0.8) + sc +
#    geom_point(aes(size=L6_mCH*100, colour=log2(L6_rna)), alpha=0.5) + sc +
#    geom_point(aes(colour=log10GeneLength.kb), size=2, alpha=0.8) + sc +
#    gghighlight(L6_rna>2) + 
    theme_bw() + 
    theme(aspect.ratio=1) +
    scale_size_continuous(range = c(0.5,4)) +
#    scale_size_manual(c(0.5,4)) +
    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
#    scale_y_continuous(breaks=seq(0,4,1), labels=seq(0,4,1), limits=c(0,4)) + 
     stat_smooth(method="lm", se=TRUE, formula=y ~ poly(x, 3, raw=FALSE))



Pv.markers_tmCG.CG <- ggplot(df_Pv %>% filter(log2(Pv_rna)>1),aes(x=Pv_CG*100, y=Pv_tmCG*100)) + 
#Vip.Ndnf.markers_tmCG.CG <- ggplot(df_Vip.Ndnf %>% filter(log2(Vip.Ndnf_rna)>1),aes(x=Vip.Ndnf_CG*100, y=Vip.Ndnf_tmCG*100)) + 
#    gghighlight(Vip.Ndnf_rna>2) + 
#   geom_mask() +
#    geom_point(aes(size=log10GeneLength.kb, colour=log2(Pv_rna)), alpha=0.8) + sc +
    geom_point(aes(size=log2(Pv_rna), colour=log10GeneLength.kb), alpha=0.8) + sc +
#    geom_point(aes(size=L6_mCH*100, colour=log2(L6_rna)), alpha=0.5) + sc +
#    geom_point(aes(colour=log10GeneLength.kb), size=2, alpha=0.8) + sc +
#    gghighlight(L6_rna>2) + 
    theme_bw() + 
    theme(aspect.ratio=1) +
    scale_size_continuous(range = c(0.5,4)) +
#    scale_size_manual(c(0.5,4)) +
    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
#    scale_y_continuous(breaks=seq(0,4,1), labels=seq(0,4,1), limits=c(0,4)) + 
    geom_smooth(method='lm', formula= y~x)

Pv.markers_mCH.CG <- ggplot(df_Pv %>% filter(log2(Pv_rna)>1),aes(x=Pv_CG*100, y=Pv_mCH*100)) + 
#Vip.Ndnf.markers_mCH.CG <- ggplot(df_Vip.Ndnf %>% filter(log2(Vip.Ndnf_rna)>1),aes(x=Vip.Ndnf_CG*100, y=Vip.Ndnf_mCH*100)) + 
#    gghighlight(Vip.Ndnf_rna>2) + 
#   geom_mask() +
#   geom_point(aes(size=log10GeneLength.kb, colour=log2(Pv_rna)), alpha=0.8) + sc +
    geom_point(aes(size=log2(Pv_rna), colour=log10GeneLength.kb), alpha=0.8) + sc +
#    geom_point(aes(size=L6_mCH*100, colour=log2(L6_rna)), alpha=0.5) + sc +
#    geom_point(aes(colour=log10GeneLength.kb), size=2, alpha=0.8) + sc +
#    gghighlight(L6_rna>2) + 
    theme_bw() + 
    theme(aspect.ratio=1) +
    scale_size_continuous(range = c(0.5,4)) +
#    scale_size_manual(c(0.5,4)) +
    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
#    scale_y_continuous(breaks=seq(0,4,1), labels=seq(0,4,1), limits=c(0,4)) + 
    geom_smooth(method='lm', formula= y~x)



#==========
# Figure6e
#==========
plot_grid(
    L6.markers_hmCG.CG,
    Pv.markers_hmCG.CG,
    L6.markers_tmCG.CG,
    Pv.markers_tmCG.CG,
    L6.markers_mCH.CG,
    Pv.markers_mCH.CG,
    nrow=3,align = "vh")

ggsave("Fig6_L6.vs.Pv.markers_hmCG.vs.CG_11.05.2022.pdf",  width=15,height=12)   

candidates<-df_5modalities_celltype_all.genes_genelength[c("Arpp21","Tle4","Khdrbs3","Atp2b1","Agbl4","Mgat4c","Thsd7b","Cdh18","Slc6a1","Erbb4","Cacna2d2","Gabrg3","Cntnap5c","Il1rapl2"), ]
candidates %>% select(L6_CG, L6_hmCG, L6_tmCG, Pv_CG, Pv_hmCG, Pv_tmCG, L6_rna, Pv_rna)
rownames(candidates)
dim(df_L6 %>% filter(log2(L6_rna)>1))

# Add selective labels to a subset of data points
# https://stackoverflow.com/questions/15015356/how-to-do-selective-labeling-with-ggplot-geom-point

myPalette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
#myPalette <- colorRampPalette(brewer.pal(11, "YlOrRd"))
#myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(-1,7))

L6_geneList <- c("Arpp21","Tle4","Khdrbs3","Atp2b1","Agbl4","Mgat4c","Thsd7b","Cdh18")
Pv_geneList <- c("Slc6a1","Erbb4","Cacna2d2","Gabrg3","Cntnap5c","Il1rapl2")

L6_geneList
Pv_geneList

#rep.genes <- c("Tle4","Arpp21","Cdh18","Erbb4","Slc6a1","Cntnap5c")

#L6.markers_hmCG.CG <- ggplot(df_L6,aes(x=L6_CG, y=L6_hmCG)) + 
#L6.markers_hmCG.rna <- ggplot(df_L6,aes(x=log2(L6_rna), y=L6_hmCG)) +
#L6.markers_ratio.rna <- ggplot(df_L6,aes(x=log2(L6_rna), y=L6_hmCG.tmCG.ratio)) +
L6.markers_genelength.hmCG <- ggplot(df_L6 %>% filter(log2(L6_rna)>1),aes(x=L6_hmCG*100, y=log10GeneLength.kb)) +
#ggplot(df_L6,aes(x=L6_hmCG.tmCG.ratio, y=L6_mCH*100)) +
#   geom_mask() +
    geom_point(aes(colour=log2(L6_rna)), size=2, alpha=0.8) + sc +
    geom_point(data=df_5modalities_celltype_all.genes_genelength[L6_geneList, ], colour="red", size=3, shape=1) +
    geom_point(data=df_5modalities_celltype_all.genes_genelength[Pv_geneList, ], colour="blue", size=3, shape=1) +
##    gghighlight(L6_rna>2) + 
    theme_bw() + 
    geom_text(data=candidates, aes(x=L6_hmCG*100, y=log10GeneLength.kb, label = rownames(candidates)), hjust=-0.1,vjust=-0.1,size=5) +
#   geom_text(aes(label = rownames(df_L6 %>% filter(log2(L6_rna)>1))), hjust=-0.1,vjust=-0.1,size=3) +
#   geom_text(hjust=-0.1,vjust=-0.1,size=3) +
    theme(aspect.ratio=1) +
#    scale_size_continuous(range = c(1,5)) +
#    scale_size_manual(c(2,6)) +
    scale_x_continuous(breaks=seq(10,50,10), labels=seq(10,50,10), limits=c(10,50)) + 
    scale_y_continuous(breaks=seq(0.5,3.5,1), labels=seq(0.5,3.5,1), limits=c(0.5,3.5)) +  
    geom_smooth(method='lm', formula= y~x)

L6.markers_genelength.CG <- ggplot(df_L6 %>% filter(log2(L6_rna)>1),aes(x=L6_CG*100, y=log10GeneLength.kb)) +
#ggplot(df_L6,aes(x=L6_hmCG.tmCG.ratio, y=L6_mCH*100)) +
#   geom_mask() +
    geom_point(aes(colour=log2(L6_rna)), size=2, alpha=0.8) + sc +
    geom_point(data=df_5modalities_celltype_all.genes_genelength[L6_geneList, ], colour="red", size=3, shape=1) +
    geom_point(data=df_5modalities_celltype_all.genes_genelength[Pv_geneList, ], colour="blue", size=3, shape=1) +
##    gghighlight(L6_rna>2) + 
    geom_text(data=candidates, aes(x=L6_CG*100, y=log10GeneLength.kb, label = rownames(candidates)), hjust=-0.1,vjust=-0.1,size=5) +
    theme_bw() + 
    theme(aspect.ratio=1) +
#    scale_size_continuous(range = c(1,5)) +
#    scale_size_manual(c(2,6)) +
    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
    scale_y_continuous(breaks=seq(0.5,3.5,1), labels=seq(0.5,3.5,1), limits=c(0.5,3.5)) +  
    geom_smooth(method='lm', formula= y~x)


Vip.Ndnf.markers_genelength.CG <- ggplot(df_Vip.Ndnf %>% filter(log2(Vip.Ndnf_rna)>1),aes(x=Vip.Ndnf_CG*100, y=log10GeneLength.kb)) + 
    geom_point(aes(colour=log2(Vip.Ndnf_rna)), size=2, alpha=0.8) + sc +
    geom_point(data=df_5modalities_celltype_all.genes_genelength[c("Tle4","Arpp21","Cdh18"), ], colour="red", size=2, shape=1) +
    geom_point(data=df_5modalities_celltype_all.genes_genelength[c("Erbb4","Slc6a1","Cntnap5c","Cacna2d2"), ], colour="blue", size=2, shape=1) +
#    geom_point(data=df_5modalities_celltype_all.genes_genelength[c("Erbb4","Adarb2","Slc6a1"), ], colour="red", size=2, shape=1) +
#    gghighlight(Vip.Ndnf_rna>2) + 
    theme_bw() + 
    theme(aspect.ratio=1) +
    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
#    scale_y_continuous(breaks=seq(0,4,1), labels=seq(0,4,1), limits=c(0,4)) + 
    geom_smooth(method='lm', formula= y~x)

Pv.markers_genelength.hmCG <- ggplot(df_Pv %>% filter(log2(Pv_rna)>1),aes(x=Pv_hmCG*100, y=log10GeneLength.kb)) +
#ggplot(df_L6,aes(x=L6_hmCG.tmCG.ratio, y=L6_mCH*100)) +
#   geom_mask() +
    geom_point(aes(colour=log2(Pv_rna)), size=2, alpha=0.8) + sc +
    geom_point(data=df_5modalities_celltype_all.genes_genelength[L6_geneList, ], colour="red", size=3, shape=1) +
    geom_point(data=df_5modalities_celltype_all.genes_genelength[Pv_geneList, ], colour="blue", size=3, shape=1) +
##    gghighlight(L6_rna>2) + 
    theme_bw() + 
    geom_text(data=candidates, aes(x=Pv_hmCG*100, y=log10GeneLength.kb, label = rownames(candidates)), hjust=-0.1,vjust=-0.1,size=5) +
    theme(aspect.ratio=1) +
#    scale_size_continuous(range = c(1,5)) +
#    scale_size_manual(c(2,6)) +
    scale_x_continuous(breaks=seq(10,50,10), labels=seq(10,50,10), limits=c(10,50)) + 
    scale_y_continuous(breaks=seq(0.5,3.5,1), labels=seq(0.5,3.5,1), limits=c(0.5,3.5)) +  
    geom_smooth(method='lm', formula= y~x)

Pv.markers_genelength.CG <- ggplot(df_Pv %>% filter(log2(Pv_rna)>1),aes(x=Pv_CG*100, y=log10GeneLength.kb)) + 
    geom_point(aes(colour=log2(Pv_rna)), size=2, alpha=0.8) + sc +
    geom_point(data=df_5modalities_celltype_all.genes_genelength[L6_geneList, ], colour="red", size=3, shape=1) +
    geom_point(data=df_5modalities_celltype_all.genes_genelength[Pv_geneList, ], colour="blue", size=3, shape=1) +
##    gghighlight(L6_rna>2) + 
    theme_bw() + 
    geom_text(data=candidates, aes(x=Pv_CG*100, y=log10GeneLength.kb, label = rownames(candidates)), hjust=-0.1,vjust=-0.1,size=5) +
    theme(aspect.ratio=1) +
    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
    scale_y_continuous(breaks=seq(0.5,3.5,1), labels=seq(0.5,3.5,1), limits=c(0.5,3.5)) + 
    geom_smooth(method='lm', formula= y~x)




L6.markers_genelength.tmCG <- ggplot(df_L6 %>% filter(log2(L6_rna)>1),aes(x=L6_tmCG*100, y=log10GeneLength.kb)) + 
    geom_point(aes(colour=log2(L6_rna)), size=2, alpha=0.8) + sc +
    geom_point(data=df_5modalities_celltype_all.genes_genelength[L6_geneList, ], colour="red", size=3, shape=1) +
    geom_point(data=df_5modalities_celltype_all.genes_genelength[Pv_geneList, ], colour="blue", size=3, shape=1) +
##    gghighlight(L6_rna>2) + 
    theme_bw() + 
    geom_text(data=candidates, aes(x=L6_tmCG*100, y=log10GeneLength.kb, label = rownames(candidates)), hjust=-0.1,vjust=-0.1,size=5) +
    theme(aspect.ratio=1) +
    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
    scale_y_continuous(breaks=seq(0.5,3.5,1), labels=seq(0.5,3.5,1), limits=c(0.5,3.5)) + 
    geom_smooth(method='lm', formula= y~x)


Pv.markers_genelength.tmCG <- ggplot(df_Pv %>% filter(log2(Pv_rna)>1),aes(x=Pv_tmCG*100, y=log10GeneLength.kb)) + 
    geom_point(aes(colour=log2(Pv_rna)), size=2, alpha=0.8) + sc +
    geom_point(data=df_5modalities_celltype_all.genes_genelength[L6_geneList, ], colour="red", size=3, shape=1) +
    geom_point(data=df_5modalities_celltype_all.genes_genelength[Pv_geneList, ], colour="blue", size=3, shape=1) +
##    gghighlight(L6_rna>2) + 
    theme_bw() + 
    geom_text(data=candidates, aes(x=Pv_tmCG*100, y=log10GeneLength.kb, label = rownames(candidates)), hjust=-0.1,vjust=-0.1,size=5) +
    theme(aspect.ratio=1) +
    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
    scale_y_continuous(breaks=seq(0.5,3.5,1), labels=seq(0.5,3.5,1), limits=c(0.5,3.5)) + 
    geom_smooth(method='lm', formula= y~x)

options(repr.plot.width = 15,repr.plot.height = 18)

plot_grid(
    L6.markers_genelength.hmCG,
    #Vip.Ndnf.markers_genelength.CG,
    Pv.markers_genelength.hmCG,
    L6.markers_genelength.CG,
    Pv.markers_genelength.CG,
    L6.markers_genelength.tmCG,
    Pv.markers_genelength.tmCG,
    nrow=3,align = "vh")

ggsave("Fig6_L6.vs.Pv.markers_geneLength.vs.hmCG.CG.tmCG_RdBu_11.07.2022.pdf",  width=15,height=18)   


myPalette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
#myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(-1,7))

#L6.markers_hmCG.CG <- ggplot(df_L6,aes(x=L6_CG, y=L6_hmCG)) + 
#L6.markers_hmCG.rna <- ggplot(df_L6,aes(x=log2(L6_rna), y=L6_hmCG)) +
#L6.markers_ratio.rna <- ggplot(df_L6,aes(x=log2(L6_rna), y=L6_hmCG.tmCG.ratio)) +
L6.markers_genelength.CG <- ggplot(df_5modalities_celltype_all.genes_genelength %>% filter(log2(L6_rna)>2),aes(x=L6_CG*100, y=log10GeneLength.kb)) +
#ggplot(df_L6,aes(x=L6_hmCG.tmCG.ratio, y=L6_mCH*100)) +
#   geom_mask() +
    geom_point(aes(colour=log2(L6_rna)), size=2, alpha=0.8) + sc +
#    gghighlight(L6_rna>2) + 
    theme_bw() + 
    theme(aspect.ratio=1) +
#    scale_size_continuous(range = c(1,5)) +
#    scale_size_manual(c(2,6)) +
#    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
#    scale_y_continuous(breaks=seq(0,4,1), labels=seq(0,4,1), limits=c(0,4)) + 
    geom_smooth(method='lm', formula= y~x)

L2.3.markers_genelength.CG <- ggplot(df_5modalities_celltype_all.genes_genelength %>% filter(log2(L2.3_rna)>2),aes(x=L2.3_CG*100, y=log10GeneLength.kb)) + 
    geom_point(aes(colour=log2(L2.3_rna)), size=2, alpha=0.8) + sc +
#    gghighlight(Vip.Ndnf_rna>2) + 
    theme_bw() + 
    theme(aspect.ratio=1) +
#    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
#    scale_y_continuous(breaks=seq(0,4,1), labels=seq(0,4,1), limits=c(0,4)) + 
    geom_smooth(method='lm', formula= y~x)

Vip.Ndnf.markers_genelength.CG <- ggplot(df_5modalities_celltype_all.genes_genelength %>% filter(log2(Vip.Ndnf_rna)>2),aes(x=Vip.Ndnf_CG*100, y=log10GeneLength.kb)) + 
    geom_point(aes(colour=log2(Vip.Ndnf_rna)), size=2, alpha=0.8) + sc +
#    gghighlight(Vip.Ndnf_rna>2) + 
    theme_bw() + 
    theme(aspect.ratio=1) +
#    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
#    scale_y_continuous(breaks=seq(0,4,1), labels=seq(0,4,1), limits=c(0,4)) + 
    geom_smooth(method='lm', formula= y~x)


Pv.markers_genelength.CG <- ggplot(df_5modalities_celltype_all.genes_genelength %>% filter(log2(Pv_rna)>2),aes(x=Pv_CG*100, y=log10GeneLength.kb)) + 
    geom_point(aes(colour=log2(Pv_rna)), size=2, alpha=0.8) + sc +
#    gghighlight(Vip.Ndnf_rna>2) + 
    theme_bw() + 
    theme(aspect.ratio=1) +
#    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
#    scale_y_continuous(breaks=seq(0,4,1), labels=seq(0,4,1), limits=c(0,4)) + 
    geom_smooth(method='lm', formula= y~x)


options(repr.plot.width = 15,repr.plot.height = 12)

plot_grid(
    L6.markers_genelength.CG,
    Pv.markers_genelength.CG,
    L2.3.markers_genelength.CG,
    Vip.Ndnf.markers_genelength.CG,
    nrow=2,align = "vh")

ggsave("Fig6_L6.L23.Pv.VipNdnf.markers_geneLength.vs.CG_11.07.2022.pdf",  width=15,height=12)   

options(repr.plot.width = 15,repr.plot.height = 12)

myPalette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
#myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(0,4))

#L6.all_rna.CG <- ggplot(df_4modalities_celltype_all.genes,aes(x=log1p(L6_rna), y=L6_CG)) + 
#L6.all_rna.CG <- ggplot(df_4modalities_celltype_all.genes,aes(x=log2(L6_rna), y=L6_CG)) + 
#L6.all_rna.mCH <- 
L6.all_rna.mCH <- ggplot(df_5modalities_celltype_all.genes_genelength,aes(x=log2(L6_rna), y=L6_mCH*100)) + 
#   geom_mask() +
    geom_point(aes(colour=log10GeneLength.kb), size=0.5, alpha=0.5) + sc +
#    gghighlight(L6_rna>2) + 
    theme_bw() + 
    theme(aspect.ratio=1) +
    scale_y_continuous(breaks=seq(0,6,1), labels=seq(0,6,1), limits=c(0,6)) + 
    geom_smooth(method='lm', formula= y~x)


#L6.all_rna.CG <- ggplot(df_4modalities_celltype_all.genes,aes(x=log1p(L6_rna), y=L6_CG)) + 
#L6.all_rna.CG <- ggplot(df_4modalities_celltype_all.genes,aes(x=log2(L6_rna), y=L6_CG)) + 
L6.all_rna.hmCG <- ggplot(df_5modalities_celltype_all.genes_genelength,aes(x=log2(L6_rna), y=L6_hmCG*100)) + 
#   geom_mask() +
    geom_point(aes(colour=log10GeneLength.kb), size=0.5, alpha=0.5) + sc +
#    gghighlight(L6_rna>2) + 
    theme_bw() + 
    theme(aspect.ratio=1) +
    geom_smooth(method='lm', formula= y~x)

L6.all_rna.tmCG <- ggplot(df_5modalities_celltype_all.genes_genelength,aes(x=log2(L6_rna), y=L6_tmCG*100)) + 
#   geom_mask() +
    geom_point(aes(colour=log10GeneLength.kb), size=0.5, alpha=0.5) + sc +
#    gghighlight(L6_rna>2) + 
    theme_bw() + 
    theme(aspect.ratio=1) +
#    scale_y_continuous(breaks=seq(0,6,1), labels=seq(0,6,1), limits=c(0,6)) + 
    geom_smooth(method='lm', formula= y~x)


#L6.all_rna.CG <- ggplot(df_5modalities_celltype_all.genes_genelength,aes(x=log2(L6_rna), y=L6_CG*100)) + 
L6.all_rna.hmCG.CG <- ggplot(df_5modalities_celltype_all.genes_genelength,aes(x=log2(L6_rna), y=(L6_hmCG+L6_CG)*100)) + 
#   geom_mask() +
    geom_point(aes(colour=log10GeneLength.kb), size=0.5, alpha=0.5) + sc +
#    gghighlight(L6_rna>2) + 
    theme_bw() + 
    theme(aspect.ratio=1) +
    geom_smooth(method='lm', formula= y~x)

plot_grid(
    L6.all_rna.mCH,
    L6.all_rna.tmCG,
    L6.all_rna.hmCG,
#    L6.all_rna.CG,
    L6.all_rna.hmCG.CG,
    nrow=2,align = "vh")

ggsave("Fig6_L6_all.genes_11.06.2022.pdf", width=15,height=18)

options(repr.plot.width = 15,repr.plot.height = 12)

myPalette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
#myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(-8,8))

#L6.all_rna.CG <- ggplot(df_4modalities_celltype_all.genes,aes(x=log1p(L6_rna), y=L6_CG)) + 
#L6.all_rna.CG <- ggplot(df_4modalities_celltype_all.genes,aes(x=log2(L6_rna), y=L6_CG)) + 
#L6.all_rna.mCH <- 
L6.all_CG.mCH <- ggplot(df_5modalities_celltype_all.genes_genelength,aes(x=L6_CG*100, y=L6_mCH*100)) + 
#   geom_mask() +
    geom_point(aes(colour=log2(L6_rna)), size=0.8, alpha=0.5) + sc +
#    gghighlight(L6_rna>2) + 
    theme_bw() + 
    theme(aspect.ratio=1) +
    scale_y_continuous(breaks=seq(0,6,1), labels=seq(0,6,1), limits=c(0,6)) + 
    geom_smooth(method='lm', formula= y~x)


#L6.all_rna.CG <- ggplot(df_4modalities_celltype_all.genes,aes(x=log1p(L6_rna), y=L6_CG)) + 
#L6.all_rna.CG <- ggplot(df_4modalities_celltype_all.genes,aes(x=log2(L6_rna), y=L6_CG)) + 
L6.all_CG.hmCG <- ggplot(df_5modalities_celltype_all.genes_genelength,aes(x=L6_CG*100, y=L6_hmCG*100)) + 
#   geom_mask() +
    geom_point(aes(colour=log2(L6_rna)), size=0.8, alpha=0.5) + sc +
#    gghighlight(L6_rna>2) + 
    theme_bw() + 
    theme(aspect.ratio=1) +
    geom_smooth(method='lm', formula= y~x)

L6.all_CG.tmCG <- ggplot(df_5modalities_celltype_all.genes_genelength,aes(x=L6_CG*100, y=L6_tmCG*100)) + 
#   geom_mask() +
    geom_point(aes(colour=log2(L6_rna)), size=0.8, alpha=0.5) + sc +
#    gghighlight(L6_rna>2) + 
    theme_bw() + 
    theme(aspect.ratio=1) +
#    scale_y_continuous(breaks=seq(0,6,1), labels=seq(0,6,1), limits=c(0,6)) + 
    geom_smooth(method='lm', formula= y~x)


   L6.all_tmCG.mCH <- ggplot(df_5modalities_celltype_all.genes_genelength,aes(x=L6_tmCG*100, y=L6_mCH*100)) + 
#   geom_mask() +
    geom_point(aes(colour=log2(L6_rna)), size=0.8, alpha=0.5) + sc +
#    gghighlight(L6_rna>2) + 
    theme_bw() + 
    theme(aspect.ratio=1) +
    scale_y_continuous(breaks=seq(0,6,1), labels=seq(0,6,1), limits=c(0,6)) + 
    geom_smooth(method='lm', formula= y~x)


#L6.all_rna.CG <- ggplot(df_4modalities_celltype_all.genes,aes(x=log1p(L6_rna), y=L6_CG)) + 
#L6.all_rna.CG <- ggplot(df_4modalities_celltype_all.genes,aes(x=log2(L6_rna), y=L6_CG)) + 
L6.all_tmCG.CG <- ggplot(df_5modalities_celltype_all.genes_genelength,aes(x=L6_tmCG*100, y=L6_CG*100)) + 
#   geom_mask() +
    geom_point(aes(colour=log2(L6_rna)), size=0.8, alpha=0.5) + sc +
#    gghighlight(L6_rna>2) + 
    theme_bw() + 
    theme(aspect.ratio=1) +
    geom_smooth(method='lm', formula= y~x)

L6.all_tmCG.hmCG <- ggplot(df_5modalities_celltype_all.genes_genelength,aes(x=L6_tmCG*100, y=L6_hmCG*100)) + 
#   geom_mask() +
    geom_point(aes(colour=log2(L6_rna)), size=0.8, alpha=0.5) + sc +
#    gghighlight(L6_rna>2) + 
    theme_bw() + 
    theme(aspect.ratio=1) +
#    scale_y_continuous(breaks=seq(0,6,1), labels=seq(0,6,1), limits=c(0,6)) + 
    geom_smooth(method='lm', formula= y~x)



plot_grid(
    L6.all_CG.mCH,
    L6.all_CG.tmCG,
    L6.all_CG.hmCG,
    L6.all_tmCG.mCH,
    L6.all_tmCG.CG,
    L6.all_tmCG.hmCG,
#    L6.all_rna.CG,
#    L6.all_genelength.CG,
    nrow=2,align = "vh")

options(repr.plot.width = 15,repr.plot.height = 12)

myPalette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
#myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(-8,8))

#L6.all_rna.CG <- ggplot(df_4modalities_celltype_all.genes,aes(x=log1p(L6_rna), y=L6_CG)) + 
#L6.all_rna.CG <- ggplot(df_4modalities_celltype_all.genes,aes(x=log2(L6_rna), y=L6_CG)) + 
#L6.all_rna.mCH <- 
L6.all_genelength.mCH <- ggplot(df_5modalities_celltype_all.genes_genelength,aes(x=log10GeneLength.kb, y=L6_mCH*100)) + 
#   geom_mask() +
    geom_point(aes(colour=log2(L6_rna)), size=0.8, alpha=0.5) + sc +
#    gghighlight(L6_rna>2) + 
    theme_bw() + 
    theme(aspect.ratio=1) +
    scale_y_continuous(breaks=seq(0,6,1), labels=seq(0,6,1), limits=c(0,6)) + 
    geom_smooth(method='lm', formula= y~x)


#L6.all_rna.CG <- ggplot(df_4modalities_celltype_all.genes,aes(x=log1p(L6_rna), y=L6_CG)) + 
#L6.all_rna.CG <- ggplot(df_4modalities_celltype_all.genes,aes(x=log2(L6_rna), y=L6_CG)) + 
L6.all_genelength.hmCG <- ggplot(df_5modalities_celltype_all.genes_genelength,aes(x=log10GeneLength.kb, y=L6_hmCG*100)) + 
#   geom_mask() +
    geom_point(aes(colour=log2(L6_rna)), size=0.8, alpha=0.5) + sc +
#    gghighlight(L6_rna>2) + 
    theme_bw() + 
    theme(aspect.ratio=1) +
    geom_smooth(method='lm', formula= y~x)

L6.all_genelength.tmCG <- ggplot(df_5modalities_celltype_all.genes_genelength,aes(x=log10GeneLength.kb, y=L6_tmCG*100)) + 
#   geom_mask() +
    geom_point(aes(colour=log2(L6_rna)), size=0.8, alpha=0.5) + sc +
#    gghighlight(L6_rna>2) + 
    theme_bw() + 
    theme(aspect.ratio=1) +
#    scale_y_continuous(breaks=seq(0,6,1), labels=seq(0,6,1), limits=c(0,6)) + 
    geom_smooth(method='lm', formula= y~x)


#L6.all_rna.CG <- ggplot(df_5modalities_celltype_all.genes_genelength,aes(x=log2(L6_rna), y=L6_CG*100)) + 
L6.all_genelength.CG <- ggplot(df_5modalities_celltype_all.genes_genelength,aes(x=log10GeneLength.kb, y=L6_CG*100)) + 
#   geom_mask() +
    geom_point(aes(colour=log2(L6_rna)), size=0.8, alpha=0.5) + sc +
#    gghighlight(L6_rna>2) + 
    theme_bw() + 
    theme(aspect.ratio=1) +
    geom_smooth(method='lm', formula= y~x)

plot_grid(
    L6.all_genelength.mCH,
    L6.all_genelength.tmCG,
    L6.all_genelength.hmCG,
#    L6.all_rna.CG,
    L6.all_genelength.CG,
    nrow=2,align = "vh")


myPalette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
#myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(0,4))

Pv.all_rna.hmCG <- ggplot(df_5modalities_celltype_all.genes_genelength,aes(x=log2(Pv_rna), y=Pv_hmCG*100)) + 
#   geom_mask() +
    geom_point(aes(colour=log10GeneLength.kb), size=0.5, alpha=0.5) + sc +
#    gghighlight(Pv_rna>2) + 
#    theme_bw() + theme(legend.position="top") +
    theme_bw() + theme(legend.position="none") +
    theme(aspect.ratio=1) +
    geom_smooth(method='lm', formula= y~x)

Pv.all_rna.CG <- ggplot(df_5modalities_celltype_all.genes_genelength,aes(x=log2(Pv_rna), y=Pv_CG*100)) + 
#Pv.all_rna.hmCG <- ggplot(df_4modalities_celltype_all.genes,aes(x=log2(Pv_rna), y=Pv_hmCG*100)) + 
#   geom_mask() +
    geom_point(aes(colour=log10GeneLength.kb), size=0.5, alpha=0.5) + sc +
#    gghighlight(Pv_rna>2) + 
    theme_bw() + theme(legend.position="none") +
    theme(aspect.ratio=1) +
    geom_smooth(method='lm', formula= y~x)

Pv.all_rna.tmCG <- ggplot(df_5modalities_celltype_all.genes_genelength,aes(x=log2(Pv_rna), y=Pv_tmCG*100)) + 
#Pv.all_rna.hmCG <- ggplot(df_4modalities_celltype_all.genes,aes(x=log2(Pv_rna), y=Pv_hmCG*100)) + 
#   geom_mask() +
    geom_point(aes(colour=log10GeneLength.kb), size=0.5, alpha=0.5) + sc +
#    gghighlight(Pv_rna>2) + 
    theme_bw() + theme(legend.position="none") +
    theme(aspect.ratio=1) +
    geom_smooth(method='lm', formula= y~x)

Pv.all_rna.mCH <- ggplot(df_5modalities_celltype_all.genes_genelength,aes(x=log2(Pv_rna), y=Pv_mCH)) +
#Pv.all_rna.ratio <- ggplot(df_5modalities_celltype_all.genes,aes(x=log2(Pv_rna), y=Pv_hmCG/Pv_tmCG)) + 
#   geom_mask() +
    geom_point(aes(colour=log10GeneLength.kb), size=0.5, alpha=0.5) + sc +
#    gghighlight(Pv_rna>2) + 
    theme_bw() + theme(legend.position="none") +
    theme(aspect.ratio=1) +
    geom_smooth(method='lm', formula= y~x)



#L6.all_rna.CG <- ggplot(df_4modalities_celltype_all.genes,aes(x=log1p(L6_rna), y=L6_CG)) + 
#L6.all_rna.CG <- ggplot(df_4modalities_celltype_all.genes,aes(x=log2(L6_rna), y=L6_CG)) + 
L6.all_rna.hmCG <- ggplot(df_5modalities_celltype_all.genes_genelength,aes(x=log2(L6_rna), y=L6_hmCG*100)) + 
#   geom_mask() +
    geom_point(aes(colour=log10GeneLength.kb), size=0.5, alpha=0.5) + sc +
#    gghighlight(L6_rna>2) + 
    theme_bw() + theme(legend.position="none") +
    theme(aspect.ratio=1) +
    geom_smooth(method='lm', formula= y~x)

L6.all_rna.CG <- ggplot(df_5modalities_celltype_all.genes_genelength,aes(x=log2(L6_rna), y=L6_CG*100)) + 
#L6.all_rna.hmCG <- ggplot(df_4modalities_celltype_all.genes,aes(x=log2(L6_rna), y=L6_hmCG*100)) + 
#   geom_mask() +
    geom_point(aes(colour=log10GeneLength.kb), size=0.5, alpha=0.5) + sc +
#    gghighlight(L6_rna>2) + 
    theme_bw() + theme(legend.position="none") +
    theme(aspect.ratio=1) +
    geom_smooth(method='lm', formula= y~x)

L6.all_rna.tmCG <- ggplot(df_5modalities_celltype_all.genes_genelength,aes(x=log2(L6_rna), y=L6_tmCG*100)) + 
#L6.all_rna.hmCG <- ggplot(df_4modalities_celltype_all.genes,aes(x=log2(L6_rna), y=L6_hmCG*100)) + 
#   geom_mask() +
    geom_point(aes(colour=log10GeneLength.kb), size=0.5, alpha=0.5) + sc +
#    gghighlight(L6_rna>2) + 
    theme_bw() + theme(legend.position="none") +
    theme(aspect.ratio=1) +
    geom_smooth(method='lm', formula= y~x)

L6.all_rna.mCH <- ggplot(df_5modalities_celltype_all.genes_genelength,aes(x=log2(L6_rna), y=L6_mCH)) +
#L6.all_rna.ratio <- ggplot(df_5modalities_celltype_all.genes,aes(x=log2(L6_rna), y=L6_hmCG/L6_tmCG)) + 
#   geom_mask() +
    geom_point(aes(colour=log10GeneLength.kb), size=0.5, alpha=0.5) + sc +
#    gghighlight(L6_rna>2) + 
    theme_bw() + theme(legend.position="none") +
    theme(aspect.ratio=1) +
    geom_smooth(method='lm', formula= y~x)


L2.3.all_rna.hmCG <- ggplot(df_5modalities_celltype_all.genes_genelength,aes(x=log2(L2.3_rna), y=L2.3_hmCG*100)) + 
#   geom_mask() +
    geom_point(aes(colour=log10GeneLength.kb), size=0.5, alpha=0.5) + sc +
#    gghighlight(L2.3_rna>2) + 
    theme_bw() + theme(legend.position="none") +
    theme(aspect.ratio=1) +
    geom_smooth(method='lm', formula= y~x)

L2.3.all_rna.CG <- ggplot(df_5modalities_celltype_all.genes_genelength,aes(x=log2(L2.3_rna), y=L2.3_CG*100)) + 
#L2.3.all_rna.hmCG <- ggplot(df_4modalities_celltype_all.genes,aes(x=log2(L2.3_rna), y=L2.3_hmCG*100)) + 
#   geom_mask() +
    geom_point(aes(colour=log10GeneLength.kb), size=0.5, alpha=0.5) + sc +
#    gghighlight(L2.3_rna>2) + 
    theme_bw() + theme(legend.position="none") +
    theme(aspect.ratio=1) +
    geom_smooth(method='lm', formula= y~x)

L2.3.all_rna.tmCG <- ggplot(df_5modalities_celltype_all.genes_genelength,aes(x=log2(L2.3_rna), y=L2.3_tmCG*100)) + 
#L2.3.all_rna.hmCG <- ggplot(df_4modalities_celltype_all.genes,aes(x=log2(L2.3_rna), y=L2.3_hmCG*100)) + 
#   geom_mask() +
    geom_point(aes(colour=log10GeneLength.kb), size=0.5, alpha=0.5) + sc +
#    gghighlight(L2.3_rna>2) + 
    theme_bw() + theme(legend.position="none") +
    theme(aspect.ratio=1) +
    geom_smooth(method='lm', formula= y~x)

L2.3.all_rna.mCH <- ggplot(df_5modalities_celltype_all.genes_genelength,aes(x=log2(L2.3_rna), y=L2.3_mCH)) +
#L2.3.all_rna.ratio <- ggplot(df_5modalities_celltype_all.genes,aes(x=log2(L2.3_rna), y=L2.3_hmCG/L2.3_tmCG)) + 
#   geom_mask() +
    geom_point(aes(colour=log10GeneLength.kb), size=0.5, alpha=0.5) + sc +
#    gghighlight(L2.3_rna>2) + 
    theme_bw() + theme(legend.position="none") +
    theme(aspect.ratio=1) +
    geom_smooth(method='lm', formula= y~x)


Vip.Ndnf.all_rna.hmCG <- ggplot(df_5modalities_celltype_all.genes_genelength,aes(x=log2(Vip.Ndnf_rna), y=Vip.Ndnf_hmCG*100)) + 
#   geom_mask() +
    geom_point(aes(colour=log10GeneLength.kb), size=0.5, alpha=0.5) + sc +
#    gghighlight(Vip.Ndnf_rna>2) + 
    theme_bw() + theme(legend.position="none") +
    theme(aspect.ratio=1) +
    geom_smooth(method='lm', formula= y~x)

Vip.Ndnf.all_rna.CG <- ggplot(df_5modalities_celltype_all.genes_genelength,aes(x=log2(Vip.Ndnf_rna), y=Vip.Ndnf_CG*100)) + 
#Vip.Ndnf.all_rna.hmCG <- ggplot(df_4modalities_celltype_all.genes,aes(x=log2(Vip.Ndnf_rna), y=Vip.Ndnf_hmCG*100)) + 
#   geom_mask() +
    geom_point(aes(colour=log10GeneLength.kb), size=0.5, alpha=0.5) + sc +
#    gghighlight(Vip.Ndnf_rna>2) + 
    theme_bw() + theme(legend.position="none") +
    theme(aspect.ratio=1) +
    geom_smooth(method='lm', formula= y~x)

Vip.Ndnf.all_rna.tmCG <- ggplot(df_5modalities_celltype_all.genes_genelength,aes(x=log2(Vip.Ndnf_rna), y=Vip.Ndnf_tmCG*100)) + 
#Vip.Ndnf.all_rna.hmCG <- ggplot(df_4modalities_celltype_all.genes,aes(x=log2(Vip.Ndnf_rna), y=Vip.Ndnf_hmCG*100)) + 
#   geom_mask() +
    geom_point(aes(colour=log10GeneLength.kb), size=0.5, alpha=0.5) + sc +
#    gghighlight(Vip.Ndnf_rna>2) + 
    theme_bw() + theme(legend.position="none") +
    theme(aspect.ratio=1) +
    geom_smooth(method='lm', formula= y~x)

Vip.Ndnf.all_rna.mCH <- ggplot(df_5modalities_celltype_all.genes_genelength,aes(x=log2(Vip.Ndnf_rna), y=Vip.Ndnf_mCH)) +
#Vip.Ndnf.all_rna.ratio <- ggplot(df_5modalities_celltype_all.genes,aes(x=log2(Vip.Ndnf_rna), y=Vip.Ndnf_hmCG/Vip.Ndnf_tmCG)) + 
#   geom_mask() +
    geom_point(aes(colour=log10GeneLength.kb), size=0.5, alpha=0.5) + sc +
#    gghighlight(Vip.Ndnf_rna>2) + 
    theme_bw() + theme(legend.position="none") +
    theme(aspect.ratio=1) +
    geom_smooth(method='lm', formula= y~x)

options(repr.plot.width = 20,repr.plot.height = 20)

plot_grid(
#    L6.all_rna.CG, 
     L2.3.all_rna.hmCG,
     L2.3.all_rna.tmCG,
     L2.3.all_rna.CG,
     L2.3.all_rna.mCH,
     L6.all_rna.hmCG,
     L6.all_rna.tmCG,
     L6.all_rna.CG,
     L6.all_rna.mCH,
     Pv.all_rna.hmCG,
     Pv.all_rna.tmCG,
     Pv.all_rna.CG,
     Pv.all_rna.mCH,
     Vip.Ndnf.all_rna.hmCG,
     Vip.Ndnf.all_rna.tmCG,
     Vip.Ndnf.all_rna.CG,
     Vip.Ndnf.all_rna.mCH,
    
#     L6.all_rna.ratio,
#    L6.markers_ratio.rna,
#    L6.markers_hmCG.rna,
#    L6.markers_mCH.tmCG,
#    L6.markers_ratio.CG,
#    L6.markers_ratio2.CG, 
#    L6.markers_hmCG.rna,
#    L6.markers_tmCG.rna,
#    L6.markers_CG.rna,
#    L6.markers_mCH.rna,
#    L6.markers_ratio.rna,
#    Vip.Ndnf.markers_ratio.CG,
    nrow=4,align = "vh")


ggsave("Fig6_L23.L6.Pv.VipNdnf_all.genes_11.06.2022.pdf", width=15,height=18)

myPalette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
#myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(-8, 8))

L6_rna<-ggtern(data=df_4modalities_celltype_all.genes, aes(x=L6_tmCG, y=L6_hmCG, z=L6_CG)) +
  geom_mask() +
  geom_point(aes(colour=log2(L6_rna)), size=1, alpha=0.5) +
  sc + 
   theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +
  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +
  L2.3.all_rna.hmCG <- ggplot(df_5modalities_celltype_all.genes_genelength,aes(x=log2(L2.3_rna), y=L2.3_hmCG*100)) + 
#   geom_mask() +
    geom_point(aes(colour=log10GeneLength.kb), size=0.5, alpha=0.5) + sc +
#    gghighlight(L2.3_rna>2) + 
    theme_bw() + theme(legend.position="top") +
    theme(aspect.ratio=1) +
    geom_smooth(method='lm', formula= y~x)

L2.3.all_rna.CG <- ggplot(df_5modalities_celltype_all.genes_genelength,aes(x=log2(L2.3_rna), y=L2.3_CG*100)) + 
#L2.3.all_rna.hmCG <- ggplot(df_4modalities_celltype_all.genes,aes(x=log2(L2.3_rna), y=L2.3_hmCG*100)) + 
#   geom_mask() +
    geom_point(aes(colour=log10GeneLength.kb), size=0.5, alpha=0.5) + sc +
#    gghighlight(L2.3_rna>2) + 
    theme_bw() + theme(legend.position="top") +
    theme(aspect.ratio=1) +
    geom_smooth(method='lm', formula= y~x)

L2.3.all_rna.tmCG <- ggplot(df_5modalities_celltype_all.genes_genelength,aes(x=log2(L2.3_rna), y=L2.3_tmCG*100)) + 
#L2.3.all_rna.hmCG <- ggplot(df_4modalities_celltype_all.genes,aes(x=log2(L2.3_rna), y=L2.3_hmCG*100)) + 
#   geom_mask() +
    geom_point(aes(colour=log10GeneLength.kb), size=0.5, alpha=0.5) + sc +
#    gghighlight(L2.3_rna>2) + 
    theme_bw() + theme(legend.position="top") +
    theme(aspect.ratio=1) +
    geom_smooth(method='lm', formula= y~x)

L2.3.all_rna.mCH <- ggplot(df_5modalities_celltype_all.genes_genelength,aes(x=log2(L2.3_rna), y=L2.3_mCH)) +
#L2.3.all_rna.ratio <- ggplot(df_5modalities_celltype_all.genes,aes(x=log2(L2.3_rna), y=L2.3_hmCG/L2.3_tmCG)) + 
#   geom_mask() +
    geom_point(aes(colour=log10GeneLength.kb), size=0.5, alpha=0.5) + sc +
#    gghighlight(L2.3_rna>2) + 
    theme_bw() + theme(legend.position="top") +
    theme(aspect.ratio=1) +
    geom_smooth(method='lm', formula= y~x)

ggtitle('Ex-L6')

Pv_rna<-ggtern(data=df_4modalities_celltype_all.genes, aes(x=Pv_tmCG, y=Pv_hmCG, z=Pv_CG)) +
  geom_mask() +
  geom_point(aes(colour=log2(Pv_rna)), size=1, alpha=0.5) +
  sc + 
   theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +
  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +
  ggtitle('Inh-Pv')


plot_list=list()
plot_list[["L6_rna"]] <- L6_rna
plot_list[["Pv_rna"]] <- Pv_rna

#==========
# Figure6a
#==========
options(repr.plot.width = 15,repr.plot.height = 6)
#g <- grid.arrange(grobs=plot_list,nrow = 1)
ternary_plot <- ggtern::grid.arrange(grobs=plot_list,nrow = 1)

# https://bjnnowak.netlify.app/2021/07/26/r-plotting-soil-textures-example-of-water-storage-capacity/

#myPalette <- colorRampPalette(rev(brewer.pal(11, "RdBu"))) 
#myPalette <- colorRampPalette(brewer.pal(11, "YlOrRd")) 
#myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
#sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(0,7))
#sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(0.5,4))

L6_genelength <- ggtern(data=df_5modalities_celltype_all.genes_genelength, aes(x=L6_tmCG, y=L6_hmCG, z=L6_CG)) +
#ggtern(data=df_4modalities_celltype_all.genes %>% filter(L6_rna>10), aes(x=L6_tmCG, y=L6_hmCG, z=L6_CG)) +

  geom_mask() +
  geom_point(aes(colour=log10GeneLength.kb), size=0.6, alpha=0.8) + 
#  sc +
  scale_color_continuous(
    low="dodgerblue",
    high="tomato"
  )+

  gghighlight(L6_rna>1) + 
#gghighlight(L6_CG >= 0.5) +
#gghighlight(L6_hmCG/L6_tmCG > 1) +
#gghighlight(L6_CG<0.2 & L6_hmCG>0.3 & L6_rna>2) + 

geom_point(data=df_5modalities_celltype_all.genes_genelength[c("Tle4","Arpp21"), ], colour="red", size=2, shape=1) +
geom_point(data=df_5modalities_celltype_all.genes_genelength[c("Erbb4","Slc6a1"), ], colour="blue", size=2, shape=1) + # ,"Adarb2"
geom_point(data=df_5modalities_celltype_all.genes_genelength[c("Cdh18"), ], colour="green", size=2, shape=1) +
geom_point(data=df_5modalities_celltype_all.genes_genelength[c("Cntnap5c"), ], colour="black", size=2, shape=1) + # ,"Adarb2"


#geom_point(data=df_4modalities_celltype_all.genes, colour="gray", size=0.5, alpha=0.25) +
#geom_point(data=df_4modalities_celltype_all.genes[L6_rna.specific.gene, ] %>% filter(log2(L6_rna)>1), colour="tomato", size=0.5, alpha=0.25) +
#geom_point(data=df_4modalities_celltype_all.genes[L2.3_rna.specific.gene, ], colour="green", size=0.5, alpha=0.25) +
#geom_point(data=df_4modalities_celltype_all.genes[Pv_rna.specific.gene, ] %>% filter(log2(Pv_rna)>1), colour="dodgerblue", size=0.5, alpha=0.25) +

#geom_point(data=df_4modalities_celltype_all.genes[c("Tle4"), ], colour="red", size=2) +
#geom_point(data=df_4modalities_celltype_all.genes[c("Adarb2"), ], colour="green", size=2) +

  theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +

  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +

  ggtitle('Ex-L6')


Pv_genelength <- ggtern(data=df_5modalities_celltype_all.genes_genelength, aes(x=Pv_tmCG, y=Pv_hmCG, z=Pv_CG)) +
#ggtern(data=df_4modalities_celltype_all.genes %>% filter(L6_rna>10), aes(x=L6_tmCG, y=L6_hmCG, z=L6_CG)) +

  geom_mask() +
  geom_point(aes(colour=log10GeneLength.kb), size=0.6, alpha=0.8) + 
# sc +
  scale_color_continuous(
    low="dodgerblue",
    high="tomato"
  )+

  gghighlight(Pv_rna>1) + 
#gghighlight(L6_CG >= 0.5) +
#gghighlight(L6_hmCG/L6_tmCG > 1) +
#gghighlight(L6_CG<0.2 & L6_hmCG>0.3 & L6_rna>2) + 

geom_point(data=df_5modalities_celltype_all.genes_genelength[c("Tle4","Arpp21"), ], colour="red", size=2, shape=1) +
geom_point(data=df_5modalities_celltype_all.genes_genelength[c("Erbb4","Slc6a1"), ], colour="blue", size=2, shape=1) + # ,"Adarb2"
geom_point(data=df_5modalities_celltype_all.genes_genelength[c("Cdh18"), ], colour="green", size=2, shape=1) +
geom_point(data=df_5modalities_celltype_all.genes_genelength[c("Cntnap5c"), ], colour="black", size=2, shape=1) + # ,"Adarb2"


#geom_point(data=df_4modalities_celltype_all.genes, colour="gray", size=0.5, alpha=0.25) +
#geom_point(data=df_4modalities_celltype_all.genes[L6_rna.specific.gene, ] %>% filter(log2(L6_rna)>1), colour="tomato", size=0.5, alpha=0.25) +
#geom_point(data=df_4modalities_celltype_all.genes[L2.3_rna.specific.gene, ], colour="green", size=0.5, alpha=0.25) +
#geom_point(data=df_4modalities_celltype_all.genes[Pv_rna.specific.gene, ] %>% filter(log2(Pv_rna)>1), colour="dodgerblue", size=0.5, alpha=0.25) +

#geom_point(data=df_4modalities_celltype_all.genes[c("Tle4"), ], colour="red", size=2) +
#geom_point(data=df_4modalities_celltype_all.genes[c("Adarb2"), ], colour="green", size=2) +


  theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +

  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +

  ggtitle('Inh-Pv')



plot_list=list()
#plot_list[["L2.3"]] <- L2.3
#plot_list[["L4"]] <- L4
#plot_list[["L5"]] <- L5
plot_list[["L6"]] <- L6_genelength
#plot_list[["DL"]] <- DL
plot_list[["Pv"]] <- Pv_genelength
#plot_list[["Sst"]] <- Sst
#plot_list[["Vip.Ndnf"]] <- Vip.Ndnf
#plot_list[["Striatum"]] <- Striatum
#plot_list[["Astro"]] <- Astro
#plot_list[["Oligo"]] <- Oligo
#plot_list[["MG"]] <- MG

library(gridExtra)

#==========
# Figure6b
#==========
options(repr.plot.width = 15,repr.plot.height = 9)
#g <- grid.arrange(grobs=plot_list,nrow = 1)
ternary_plot <- ggtern::grid.arrange(grobs=plot_list,nrow = 1)

# https://bjnnowak.netlify.app/2021/07/26/r-plotting-soil-textures-example-of-water-storage-capacity/

#myPalette <- colorRampPalette(rev(brewer.pal(11, "RdBu"))) 
#myPalette <- colorRampPalette(brewer.pal(11, "YlOrRd")) 
#myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
#sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(0,7))
#sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(0.5,4))

L6_genelength <- ggtern(data=df_5modalities_celltype_all.genes_genelength, aes(x=L6_tmCG, y=L6_hmCG, z=L6_CG)) +
#ggtern(data=df_4modalities_celltype_all.genes %>% filter(L6_rna>10), aes(x=L6_tmCG, y=L6_hmCG, z=L6_CG)) +

  geom_mask() +
  geom_point(aes(colour=log10GeneLength.kb), size=0.6, alpha=0.8) + 
#  sc +
  scale_color_continuous(
    low="dodgerblue",
    high="tomato"
  )+

#  gghighlight(L6_rna>1) + 
#gghighlight(L6_CG >= 0.5) +
#gghighlight(L6_hmCG/L6_tmCG > 1) +
#gghighlight(L6_CG<0.2 & L6_hmCG>0.3 & L6_rna>2) + 

#geom_point(data=df_5modalities_celltype_all.genes_genelength[c("Tle4","Arpp21"), ], colour="red", size=2, shape=1) +
#geom_point(data=df_5modalities_celltype_all.genes_genelength[c("Erbb4","Slc6a1"), ], colour="blue", size=2, shape=1) + # ,"Adarb2"
#geom_point(data=df_5modalities_celltype_all.genes_genelength[c("Cdh18"), ], colour="green", size=2, shape=1) +
#geom_point(data=df_5modalities_celltype_all.genes_genelength[c("Cntnap5c"), ], colour="black", size=2, shape=1) + # ,"Adarb2"


#geom_point(data=df_4modalities_celltype_all.genes, colour="gray", size=0.5, alpha=0.25) +
#geom_point(data=df_4modalities_celltype_all.genes[L6_rna.specific.gene, ] %>% filter(log2(L6_rna)>1), colour="tomato", size=0.5, alpha=0.25) +
#geom_point(data=df_4modalities_celltype_all.genes[L2.3_rna.specific.gene, ], colour="green", size=0.5, alpha=0.25) +
#geom_point(data=df_4modalities_celltype_all.genes[Pv_rna.specific.gene, ] %>% filter(log2(Pv_rna)>1), colour="dodgerblue", size=0.5, alpha=0.25) +

#geom_point(data=df_4modalities_celltype_all.genes[c("Tle4"), ], colour="red", size=2) +
#geom_point(data=df_4modalities_celltype_all.genes[c("Adarb2"), ], colour="green", size=2) +

  theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +

  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +

  ggtitle('Ex-L6')


Pv_genelength <- ggtern(data=df_5modalities_celltype_all.genes_genelength, aes(x=Pv_tmCG, y=Pv_hmCG, z=Pv_CG)) +
#ggtern(data=df_4modalities_celltype_all.genes %>% filter(L6_rna>10), aes(x=L6_tmCG, y=L6_hmCG, z=L6_CG)) +

  geom_mask() +
  geom_point(aes(colour=log10GeneLength.kb), size=0.6, alpha=0.8) + 
# sc +
  scale_color_continuous(
    low="dodgerblue",
    high="tomato"
  )+

#  gghighlight(Pv_rna>1) + 
#gghighlight(L6_CG >= 0.5) +
#gghighlight(L6_hmCG/L6_tmCG > 1) +
#gghighlight(L6_CG<0.2 & L6_hmCG>0.3 & L6_rna>2) + 

#geom_point(data=df_5modalities_celltype_all.genes_genelength[c("Tle4","Arpp21"), ], colour="red", size=2, shape=1) +
#geom_point(data=df_5modalities_celltype_all.genes_genelength[c("Erbb4","Slc6a1"), ], colour="blue", size=2, shape=1) + # ,"Adarb2"
#geom_point(data=df_5modalities_celltype_all.genes_genelength[c("Cdh18"), ], colour="green", size=2, shape=1) +
#geom_point(data=df_5modalities_celltype_all.genes_genelength[c("Cntnap5c"), ], colour="black", size=2, shape=1) + # ,"Adarb2"


#geom_point(data=df_4modalities_celltype_all.genes, colour="gray", size=0.5, alpha=0.25) +
#geom_point(data=df_4modalities_celltype_all.genes[L6_rna.specific.gene, ] %>% filter(log2(L6_rna)>1), colour="tomato", size=0.5, alpha=0.25) +
#geom_point(data=df_4modalities_celltype_all.genes[L2.3_rna.specific.gene, ], colour="green", size=0.5, alpha=0.25) +
#geom_point(data=df_4modalities_celltype_all.genes[Pv_rna.specific.gene, ] %>% filter(log2(Pv_rna)>1), colour="dodgerblue", size=0.5, alpha=0.25) +

#geom_point(data=df_4modalities_celltype_all.genes[c("Tle4"), ], colour="red", size=2) +
#geom_point(data=df_4modalities_celltype_all.genes[c("Adarb2"), ], colour="green", size=2) +


  theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +

  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +

  ggtitle('Inh-Pv')



plot_list=list()
#plot_list[["L2.3"]] <- L2.3
#plot_list[["L4"]] <- L4
#plot_list[["L5"]] <- L5
plot_list[["L6"]] <- L6_genelength
#plot_list[["DL"]] <- DL
plot_list[["Pv"]] <- Pv_genelength
#plot_list[["Sst"]] <- Sst
#plot_list[["Vip.Ndnf"]] <- Vip.Ndnf
#plot_list[["Striatum"]] <- Striatum
#plot_list[["Astro"]] <- Astro
#plot_list[["Oligo"]] <- Oligo
#plot_list[["MG"]] <- MG

library(gridExtra)

options(repr.plot.width = 15,repr.plot.height = 9)

#g <- grid.arrange(grobs=plot_list,nrow = 1)
ternary_plot <- ggtern::grid.arrange(grobs=plot_list,nrow = 1)

# in R [(R_4.1.1) haowu@methylome2:]:  install.packages("directlabels")
library("directlabels")

# https://stackoverflow.com/questions/29267679/avoid-overlapping-labels-with-ggplot2-and-directlabels-using-ggtern
options(repr.plot.width = 8,repr.plot.height = 6)
data(Feldspar)
base <- ggtern(Feldspar,
               aes(x=Ab,y=An,z=Or,group=Feldspar,colour=Feldspar)) + 
        geom_point() +
        geom_mask() 


.approvedgeom = c(ggtern:::.approvedgeom,c('GeomDl'))
assignInNamespace('.approvedgeom',
                   .approvedgeom,
                   envir=as.environment('package:ggtern'))
direct.label(base)


# https://stackoverflow.com/questions/64555039/add-variable-function-label-to-ternary-plot-in-r
options(repr.plot.width = 8,repr.plot.height = 6)
df <- data.frame(x=c(10,20,30), y=c(15,25,35), z=c(75,55,35), VALUE=c(1,2,3), lab=c("One", "", "Three"))
df 

ggtern(data=df, aes(x=x, y=y, z=z, value = VALUE)) +
  geom_point(aes(fill = VALUE), size = 2, stroke = 0, shape = 21) + 
  scale_fill_gradient(low = "red",high = "yellow", guide = F) + 
  scale_color_gradient(low = "red",high = "yellow", guide = F) +
  geom_text(aes(label = lab), vjust=1)

candidates<-df_5modalities_celltype_all.genes_genelength[c("Tle4","Arpp21","Cdh18","Erbb4","Slc6a1","Cntnap5c"), ]

candidates

#options(repr.plot.width = 8,repr.plot.height = 6)

L6_markers <- ggtern(data=candidates, aes(x=L6_tmCG, y=L6_hmCG, z=L6_CG)) +
    geom_point(data=df_5modalities_celltype_all.genes_genelength[c("Tle4","Arpp21","Cdh18"), ], colour="red", size=2, shape=1) +
    geom_point(data=df_5modalities_celltype_all.genes_genelength[c("Erbb4","Slc6a1","Cntnap5c"), ], colour="blue", size=2, shape=1) + # ,"Adarb2"
    geom_text(aes(label = rownames(candidates)), hjust=-0.2,vjust=-0.2,size=3) +
    theme_bw() +
    theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") 

# "Tle4","Arpp21","Cdh18","Thsd7b","Gabrg3","Erbb4","Slc6a1","Cntnap5c","Cacna2d2"
Pv_markers <- ggtern(data=candidates, aes(x=Pv_tmCG, y=Pv_hmCG, z=Pv_CG)) +
    geom_point(data=df_5modalities_celltype_all.genes_genelength[c("Tle4","Arpp21","Cdh18","Thsd7b","Gabrg3"), ], colour="red", size=2, shape=1) +
    geom_point(data=df_5modalities_celltype_all.genes_genelength[c("Erbb4","Cacna2d2","Slc6a1","Cntnap5c"), ], colour="blue", size=2, shape=1) + # ,"Adarb2"
    geom_text(aes(label = rownames(candidates)), hjust=-0.1,vjust=-0.1,size=3) +
    theme_bw() +
    theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") 

plot_list=list()

plot_list[["L6_markers"]] <- L6_markers
#plot_list[["DL"]] <- DL
plot_list[["Pv_markers"]] <- Pv_markers



options(repr.plot.width = 15,repr.plot.height = 6)
#g <- grid.arrange(grobs=plot_list,nrow = 1)
ternary_plot_rep.genes <- ggtern::grid.arrange(grobs=plot_list,nrow = 1)



ggsave("/mnt/data1/haowu/nfs2/sc-bACE-seq/Fig6_L6.vs.Pv.markers_ternary.plot_rep.genes_10.31.2022.pdf", ternary_plot_rep.genes, width=15,height=6)

# https://bjnnowak.netlify.app/2021/07/26/r-plotting-soil-textures-example-of-water-storage-capacity/

myPalette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
#myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(-6, 6))



L6_markers <- ggtern(data=df_4modalities_celltype_all.genes, aes(x=L6_tmCG, y=L6_hmCG, z=L6_CG)) +
#ggtern(data=df_4modalities_celltype_all.genes %>% filter(L6_rna>10), aes(x=L6_tmCG, y=L6_hmCG, z=L6_CG)) +

  geom_mask() +
#  geom_point(aes(colour=log2(L6_rna)), size=0.5, alpha=0.5) + sc +
#  scale_color_continuous(
#    low="dodgerblue",
#    high="tomato"
#  )+

#  gghighlight(L6_rna>2) + 
#gghighlight(L6_CG >= 0.5) +
#gghighlight(L6_hmCG/L6_tmCG > 1) +
#gghighlight(L6_CG<0.2 & L6_hmCG>0.3 & L6_rna>2) + 

geom_point(data=df_4modalities_celltype_all.genes, colour="gray", size=1, alpha=0.25) +
geom_point(data=df_4modalities_celltype_all.genes[L6_rna.specific.gene, ] %>% filter(log2(L6_rna)>0), colour="#66C2A5", size=1, alpha=0.6) +
#geom_point(data=df_4modalities_celltype_all.genes[L2.3_rna.specific.gene, ], colour="green", size=0.5, alpha=0.25) +
geom_point(data=df_4modalities_celltype_all.genes[Pv_rna.specific.gene, ] %>% filter(log2(Pv_rna)>0), colour="#E78AC3", size=1, alpha=0.6) +
geom_point(data=df_5modalities_celltype_all.genes_genelength[c("Tle4","Arpp21","Cdh18"), ], colour="red", size=2, shape=1) +
geom_point(data=df_5modalities_celltype_all.genes_genelength[c("Erbb4","Slc6a1","Cntnap5c"), ], colour="blue", size=2, shape=1) + # ,"Adarb2"

#geom_text(aes(label=rownames(candidates)), size=3, fontface="bold", color="navyblue", alpha=.5,hjust=-0.1, position = "jitter")

#geom_point(data=df_4modalities_celltype_all.genes[c("Tle4"), ], colour="red", size=2) +
#geom_point(data=df_4modalities_celltype_all.genes[c("Adarb2"), ], colour="green", size=2) +

  theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +

  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +

  ggtitle('Ex-L6')




Pv_markers <- ggtern(data=df_4modalities_celltype_all.genes, aes(x=Pv_tmCG, y=Pv_hmCG, z=Pv_CG)) +
#ggtern(data=df_4modalities_celltype_all.genes %>% filter(L6_rna>10), aes(x=L6_tmCG, y=L6_hmCG, z=L6_CG)) +

  geom_mask() +
#  geom_point(aes(colour=log2(Sst_rna)), size=0.5, alpha=0.5) + sc +
#  scale_color_continuous(
#    low="dodgerblue",
#    high="tomato"
#  )+

#  gghighlight(L6_rna>2) + 
#gghighlight(L6_CG >= 0.5) +
#gghighlight(L6_hmCG/L6_tmCG > 1) +
#gghighlight(L6_CG<0.2 & L6_hmCG>0.3 & L6_rna>2) + 

geom_point(data=df_4modalities_celltype_all.genes, colour="gray", size=1, alpha=0.25) +
geom_point(data=df_4modalities_celltype_all.genes[L6_rna.specific.gene, ] %>% filter(log2(L6_rna)>0), colour="#66C2A5", size=1, alpha=0.6) +
#geom_point(data=df_4modalities_celltype_all.genes[L2.3_rna.specific.gene, ], colour="green", size=0.5, alpha=0.25) +
geom_point(data=df_4modalities_celltype_all.genes[Pv_rna.specific.gene, ] %>% filter(log2(Pv_rna)>0), colour="#E78AC3", size=1, alpha=0.6) +
geom_point(data=df_5modalities_celltype_all.genes_genelength[c("Tle4","Arpp21","Cdh18"), ], colour="red", size=2, shape=1) +
geom_point(data=df_5modalities_celltype_all.genes_genelength[c("Erbb4","Slc6a1","Cntnap5c"), ], colour="blue", size=2, shape=1) + # ,"Adarb2"

#geom_point(data=df_4modalities_celltype_all.genes[c("Tle4"), ], colour="red", size=2) +
#geom_point(data=df_4modalities_celltype_all.genes[c("Adarb2"), ], colour="green", size=2) +


  theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +

  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +

  ggtitle('Inh-Pv')



L2.3 <- ggtern(data=df_4modalities_celltype_all.genes, aes(x=L2.3_tmCG, y=L2.3_hmCG, z=L2.3_CG)) +
#ggtern(data=df_4modalities_celltype_all.genes %>% filter(L6_rna>10), aes(x=L6_tmCG, y=L6_hmCG, z=L6_CG)) +

  geom_mask() +
#  geom_point(aes(colour=log2(Vip.Ndnf_rna)), size=0.5, alpha=0.5) + sc +
#  scale_color_continuous(
#    low="dodgerblue",
#    high="tomato"
#  )+

#  gghighlight(L6_rna>2) + 
#gghighlight(L6_CG >= 0.5) +
#gghighlight(L6_hmCG/L6_tmCG > 1) +
#gghighlight(L6_CG<0.2 & L6_hmCG>0.3 & L6_rna>2) + 

geom_point(data=df_4modalities_celltype_all.genes, colour="gray", size=1, alpha=0.25) +
geom_point(data=df_4modalities_celltype_all.genes[L6_rna.specific.gene, ] %>% filter(log2(L6_rna)>0), colour="tomato", size=1, alpha=0.6) +
geom_point(data=df_4modalities_celltype_all.genes[L2.3_rna.specific.gene, ] %>% filter(log2(L2.3_rna)>0), colour="green", size=1, alpha=0.6) +
#geom_point(data=df_4modalities_celltype_all.genes[Pv_rna.specific.gene, ], colour="dodgerblue", size=0.5, alpha=0.25) +

#geom_point(data=df_4modalities_celltype_all.genes[c("Tle4"), ], colour="red", size=2) +
#geom_point(data=df_4modalities_celltype_all.genes[c("Adarb2"), ], colour="green", size=2) +


  theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +

  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +

  ggtitle('Ex-L2.3')

plot_list=list()
#plot_list[["L2.3"]] <- L2.3
#plot_list[["L4"]] <- L4
#plot_list[["L5"]] <- L5
plot_list[["L6_rna"]] <- L6_rna
plot_list[["Pv_rna"]] <- Pv_rna
plot_list[["L6_markers"]] <- L6_markers
#plot_list[["DL"]] <- DL
plot_list[["Pv_markers"]] <- Pv_markers
plot_list[["L6_genelength"]] <- L6_genelength
plot_list[["Pv_genelength"]] <- Pv_genelength
#plot_list[["Sst"]] <- Sst
#plot_list[["Vip.Ndnf"]] <- Vip.Ndnf
#plot_list[["Striatum"]] <- Striatum
#plot_list[["Astro"]] <- Astro
#plot_list[["Oligo"]] <- Oligo
#plot_list[["MG"]] <- MG



options(repr.plot.width = 15,repr.plot.height = 18)
#g <- grid.arrange(grobs=plot_list,nrow = 1)
ternary_plot <- ggtern::grid.arrange(grobs=plot_list,nrow = 3)

ggsave("/mnt/data1/haowu/nfs2/sc-bACE-seq/Fig6_L6.vs.Pv.markers_ternary.plot_09.13.2022.pdf", ternary_plot, width=15,height=18)

ggsave("/mnt/data1/haowu/nfs2/sc-bACE-seq/Fig6_L6.vs.Pv.markers_ternary.plot_10.31.2022.pdf", ternary_plot, width=15,height=18)

# old analysis: 
# https://bjnnowak.netlify.app/2021/07/26/r-plotting-soil-textures-example-of-water-storage-capacity/

myPalette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
#myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(-6, 6))

L6_markers <- ggtern(data=df_4modalities_celltype_all.genes, aes(x=L6_tmCG, y=L6_hmCG, z=L6_CG)) +
#ggtern(data=df_4modalities_celltype_all.genes %>% filter(L6_rna>10), aes(x=L6_tmCG, y=L6_hmCG, z=L6_CG)) +

  geom_mask() +
#  geom_point(aes(colour=log2(L6_rna)), size=0.5, alpha=0.5) + sc +
#  scale_color_continuous(
#    low="dodgerblue",
#    high="tomato"
#  )+

#  gghighlight(L6_rna>2) + 
#gghighlight(L6_CG >= 0.5) +
#gghighlight(L6_hmCG/L6_tmCG > 1) +
#gghighlight(L6_CG<0.2 & L6_hmCG>0.3 & L6_rna>2) + 

geom_point(data=df_4modalities_celltype_all.genes, colour="gray", size=1, alpha=0.25) +
geom_point(data=df_4modalities_celltype_all.genes[L6_rna.specific.gene, ] %>% filter(log2(L6_rna)>0), colour="#66C2A5", size=1, alpha=0.6) +
#geom_point(data=df_4modalities_celltype_all.genes[L2.3_rna.specific.gene, ], colour="green", size=0.5, alpha=0.25) +
geom_point(data=df_4modalities_celltype_all.genes[Pv_rna.specific.gene, ] %>% filter(log2(Pv_rna)>0), colour="#E78AC3", size=1, alpha=0.6) +
geom_point(data=df_5modalities_celltype_all.genes_genelength[c("Tle4","Arpp21","Cdh18"), ], colour="red", size=2, shape=1) +
geom_point(data=df_5modalities_celltype_all.genes_genelength[c("Erbb4","Slc6a1","Cntnap5c"), ], colour="blue", size=2, shape=1) + # ,"Adarb2"

#geom_point(data=df_4modalities_celltype_all.genes[c("Tle4"), ], colour="red", size=2) +
#geom_point(data=df_4modalities_celltype_all.genes[c("Adarb2"), ], colour="green", size=2) +

  theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +

  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +

  ggtitle('Ex-L6')


Pv_markers <- ggtern(data=df_4modalities_celltype_all.genes, aes(x=Pv_tmCG, y=Pv_hmCG, z=Pv_CG)) +
#ggtern(data=df_4modalities_celltype_all.genes %>% filter(L6_rna>10), aes(x=L6_tmCG, y=L6_hmCG, z=L6_CG)) +

  geom_mask() +
#  geom_point(aes(colour=log2(Sst_rna)), size=0.5, alpha=0.5) + sc +
#  scale_color_continuous(
#    low="dodgerblue",
#    high="tomato"
#  )+

#  gghighlight(L6_rna>2) + 
#gghighlight(L6_CG >= 0.5) +
#gghighlight(L6_hmCG/L6_tmCG > 1) +
#gghighlight(L6_CG<0.2 & L6_hmCG>0.3 & L6_rna>2) + 

geom_point(data=df_4modalities_celltype_all.genes, colour="gray", size=1, alpha=0.25) +
geom_point(data=df_4modalities_celltype_all.genes[L6_rna.specific.gene, ] %>% filter(log2(L6_rna)>0), colour="#66C2A5", size=1, alpha=0.6) +
#geom_point(data=df_4modalities_celltype_all.genes[L2.3_rna.specific.gene, ], colour="green", size=0.5, alpha=0.25) +
geom_point(data=df_4modalities_celltype_all.genes[Pv_rna.specific.gene, ] %>% filter(log2(Pv_rna)>0), colour="#E78AC3", size=1, alpha=0.6) +
geom_point(data=df_5modalities_celltype_all.genes_genelength[c("Tle4","Arpp21","Cdh18"), ], colour="red", size=2, shape=1) +
geom_point(data=df_5modalities_celltype_all.genes_genelength[c("Erbb4","Slc6a1","Cntnap5c"), ], colour="blue", size=2, shape=1) + # ,"Adarb2"

#geom_point(data=df_4modalities_celltype_all.genes[c("Tle4"), ], colour="red", size=2) +
#geom_point(data=df_4modalities_celltype_all.genes[c("Adarb2"), ], colour="green", size=2) +


  theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +

  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +

  ggtitle('Inh-Pv')



L2.3 <- ggtern(data=df_4modalities_celltype_all.genes, aes(x=L2.3_tmCG, y=L2.3_hmCG, z=L2.3_CG)) +
#ggtern(data=df_4modalities_celltype_all.genes %>% filter(L6_rna>10), aes(x=L6_tmCG, y=L6_hmCG, z=L6_CG)) +

  geom_mask() +
#  geom_point(aes(colour=log2(Vip.Ndnf_rna)), size=0.5, alpha=0.5) + sc +
#  scale_color_continuous(
#    low="dodgerblue",
#    high="tomato"
#  )+

#  gghighlight(L6_rna>2) + 
#gghighlight(L6_CG >= 0.5) +
#gghighlight(L6_hmCG/L6_tmCG > 1) +
#gghighlight(L6_CG<0.2 & L6_hmCG>0.3 & L6_rna>2) + 

geom_point(data=df_4modalities_celltype_all.genes, colour="gray", size=1, alpha=0.25) +
geom_point(data=df_4modalities_celltype_all.genes[L6_rna.specific.gene, ] %>% filter(log2(L6_rna)>0), colour="tomato", size=1, alpha=0.6) +
geom_point(data=df_4modalities_celltype_all.genes[L2.3_rna.specific.gene, ] %>% filter(log2(L2.3_rna)>0), colour="green", size=1, alpha=0.6) +
#geom_point(data=df_4modalities_celltype_all.genes[Pv_rna.specific.gene, ], colour="dodgerblue", size=0.5, alpha=0.25) +

#geom_point(data=df_4modalities_celltype_all.genes[c("Tle4"), ], colour="red", size=2) +
#geom_point(data=df_4modalities_celltype_all.genes[c("Adarb2"), ], colour="green", size=2) +


  theme_bw() +
#  theme_rgbw() +
  theme_showarrows() +
    xlab("mCG") +                       #replace default axis labels
    ylab("hmCG") +
    zlab("CG") +

  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.just      = 'left') +

  ggtitle('Ex-L2.3')

plot_list=list()
#plot_list[["L2.3"]] <- L2.3
#plot_list[["L4"]] <- L4
#plot_list[["L5"]] <- L5
plot_list[["L6_rna"]] <- L6_rna
plot_list[["Pv_rna"]] <- Pv_rna
plot_list[["L6_markers"]] <- L6_markers
#plot_list[["DL"]] <- DL
plot_list[["Pv_markers"]] <- Pv_markers
plot_list[["L6_genelength"]] <- L6_genelength
plot_list[["Pv_genelength"]] <- Pv_genelength
#plot_list[["Sst"]] <- Sst
#plot_list[["Vip.Ndnf"]] <- Vip.Ndnf
#plot_list[["Striatum"]] <- Striatum
#plot_list[["Astro"]] <- Astro
#plot_list[["Oligo"]] <- Oligo
#plot_list[["MG"]] <- MG

library(gridExtra)

options(repr.plot.width = 15,repr.plot.height = 18)
#g <- grid.arrange(grobs=plot_list,nrow = 1)
ternary_plot <- grid.arrange(grobs=plot_list,nrow = 3)

df_5modalities_celltype_all.genes_genelength[L6_rna.specific.gene, ] %>% filter(log2(L6_rna)>0)

df_L6_markers <- df_5modalities_celltype_all.genes_genelength[L6_rna.specific.gene, ] %>%  mutate(L6_hmCG.tmCG.ratio = L6_hmCG / L6_tmCG) %>% arrange(desc(L6_CG)) %>% filter(log2(L6_rna)>0) %>% select(GeneName, log10GeneLength.kb, L6_CG, L6_tmCG, L6_hmCG, L6_hmCG.tmCG.ratio) 
df_L6_markers

library(reshape2)


df_L6_markers_long <- melt(df_L6_markers, id.vars=c("GeneName","log10GeneLength.kb"), variable.name="mod_type", value.name="mod_level")
df_L6_markers_long$GeneName <- factor(df_L6_markers_long$GeneName, levels = df_L6_markers$GeneName) # retain the order of original dataframe before melt
df_L6_markers_long


options(repr.plot.width = 6,repr.plot.height = 12)

ggplot(df_L6_markers_long %>% filter(mod_type != "L6_hmCG.tmCG.ratio"), aes(x=GeneName, y=mod_level, fill=mod_type)) + 
    geom_bar(position="stack", stat="identity") +
#    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
    scale_y_continuous(breaks=seq(0,1,0.1), labels=seq(0,1,0.1), limits=c(0,1)) + 
    theme(aspect.ratio=6) + 
    coord_flip() +
    scale_fill_manual(values=c("#deca0c","#fc8d62","#8DA0CB")) +  
    theme_bw()



df_L6_markers_Pv.as.ctrl <- df_5modalities_celltype_all.genes_genelength[L6_rna.specific.gene, ] %>% 
#            mutate(L6_hmCG.BSmCG.ratio = L6_hmCG / (1-L6_CG), Pv_hmCG.BSmCG.ratio = Pv_hmCG / (1-Pv_CG)) %>% 
            mutate(L6_CG.tmCG.ratio = L6_CG / L6_tmCG, Pv_CG.tmCG.ratio = Pv_CG / Pv_tmCG) %>% 
            mutate(L6_hmCG.tmCG.ratio = L6_hmCG / L6_tmCG, Pv_hmCG.tmCG.ratio = Pv_hmCG / Pv_tmCG) %>% 
            mutate(marker_class = case_when(L6_CG < 0.30 ~ 'Class1_L6', L6_CG >= 0.30 ~ 'Class2_L6')) %>% 
            arrange(desc(L6_CG)) %>% 
            filter(log2(L6_rna)>0) %>% 
            select(GeneName, log10GeneLength.kb, marker_class, L6_hmCG.tmCG.ratio, Pv_hmCG.tmCG.ratio, L6_CG.tmCG.ratio, Pv_CG.tmCG.ratio, L6_CG, L6_tmCG, L6_hmCG, L6_rna, Pv_CG, Pv_tmCG, Pv_hmCG, Pv_rna) 

table(df_L6_markers_Pv.as.ctrl$marker_class)
df_L6_markers_Pv.as.ctrl
#      Class1_L6 Class2_L6 
# 20%       43       222
# 25%       93       172
# 30%      130       135 
# 40%      208        57
# 40%      219        63 ; all markers included 



df_Pv_markers_L6.as.ctrl <- df_5modalities_celltype_all.genes_genelength[Pv_rna.specific.gene, ] %>%  
#            mutate(L6_hmCG.BSmCG.ratio = L6_hmCG / (1-L6_CG), Pv_hmCG.BSmCG.ratio = Pv_hmCG / (1-Pv_CG)) %>% 
            mutate(L6_CG.tmCG.ratio = L6_CG / L6_tmCG, Pv_CG.tmCG.ratio = Pv_CG / Pv_tmCG) %>%
            mutate(L6_hmCG.tmCG.ratio = L6_hmCG / L6_tmCG, Pv_hmCG.tmCG.ratio = Pv_hmCG / Pv_tmCG) %>% 
            mutate(marker_class = case_when(Pv_CG < 0.30 ~ 'Class1_Pv', Pv_CG >= 0.30 ~ 'Class2_Pv')) %>% 
            arrange(desc(Pv_CG)) %>% 
            filter(log2(Pv_rna)>0) %>% 
            select(GeneName, log10GeneLength.kb, marker_class, L6_hmCG.tmCG.ratio, Pv_hmCG.tmCG.ratio, L6_CG.tmCG.ratio, Pv_CG.tmCG.ratio, L6_CG, L6_tmCG, L6_hmCG, L6_rna, Pv_CG, Pv_tmCG, Pv_hmCG, Pv_rna) 

table(df_Pv_markers_L6.as.ctrl$marker_class)

df_L6_markers_Pv.as.ctrl_long <- melt(df_L6_markers_Pv.as.ctrl, id.vars=c("GeneName","log10GeneLength.kb","marker_class"), variable.name="cell.type_assay", value.name="assay.value")
df_L6_markers_Pv.as.ctrl_long$GeneName <- factor(df_L6_markers_Pv.as.ctrl_long$GeneName, levels = df_L6_markers_Pv.as.ctrl$GeneName) # retain the order of original dataframe before melt
df_L6_markers_Pv.as.ctrl_long

df_Pv_markers_L6.as.ctrl_long <- melt(df_Pv_markers_L6.as.ctrl, id.vars=c("GeneName","log10GeneLength.kb","marker_class"), variable.name="cell.type_assay", value.name="assay.value")
df_Pv_markers_L6.as.ctrl_long$GeneName <- factor(df_Pv_markers_L6.as.ctrl_long$GeneName, levels = df_Pv_markers_L6.as.ctrl$GeneName) # retain the order of original dataframe before melt
df_Pv_markers_L6.as.ctrl_long

df_L6_markers_Pv.as.ctrl_long_sep <- df_L6_markers_Pv.as.ctrl_long %>% separate(cell.type_assay, c("cell.type","assay.type"), sep="_")
head(df_L6_markers_Pv.as.ctrl_long_sep)

df_Pv_markers_L6.as.ctrl_long_sep <- df_Pv_markers_L6.as.ctrl_long %>% separate(cell.type_assay, c("cell.type","assay.type"), sep="_")
head(df_Pv_markers_L6.as.ctrl_long_sep)



#==========
# Figure6d
#==========
options(repr.plot.width = 8,repr.plot.height = 12)
L6_bargraph <- ggplot(df_L6_markers_Pv.as.ctrl_long_sep %>% filter(assay.type %in% c('CG','tmCG','hmCG')), aes(x=GeneName, y=assay.value, fill=assay.type)) + 
    geom_bar(position="stack", stat="identity") +
#    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
    scale_y_continuous(breaks=seq(0,1,0.1), labels=seq(0,1,0.1), limits=c(0,1)) + 
    theme(aspect.ratio=6) + 
    coord_flip() +
    scale_fill_manual(values=c("#deca0c","#8DA0CB","#fc8d62")) + 
    theme_bw() +
    facet_grid(. ~ cell.type)

L6_bargraph


ggsave("Fig6_L6.markers_barplot_09.13.2022.pdf", L6_bargraph, width=8,height=12)

options(repr.plot.width = 8,repr.plot.height = 12)

Pv_bargraph <- ggplot(df_Pv_markers_L6.as.ctrl_long_sep %>% filter(assay.type %in% c('CG','tmCG','hmCG')), aes(x=GeneName, y=assay.value, fill=assay.type)) + 
    geom_bar(position="stack", stat="identity") +
#    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
    scale_y_continuous(breaks=seq(0,1,0.1), labels=seq(0,1,0.1), limits=c(0,1)) + 
    theme(aspect.ratio=6) + 
    coord_flip() +
    scale_fill_manual(values=c("#deca0c","#8DA0CB","#fc8d62")) + 
    theme_bw() +
    facet_grid(. ~ cell.type)

Pv_bargraph

ggsave("Fig6_Pv.markers_barplot_09.13.2022.pdf", Pv_bargraph, width=8,height=12)

options(repr.plot.width = 6,repr.plot.height = 8)

ggplot(df_L6_markers_Pv.as.ctrl_long_sep %>% filter(assay.type %in% c('CG','tmCG','hmCG')), aes(x=cell.type, y=assay.value, fill=cell.type)) + 
        geom_boxplot(outlier.shape = NA) + geom_jitter(height = 0, width = 0.2, shape=1, size=2, alpha = 2/10) +  
#       geom_jitter(height = 0, width = 0.2, shape=21, size=1, alpha=0.6) + 
#    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
#    scale_y_continuous(breaks=seq(0,1,0.1), labels=seq(0,1,0.1), limits=c(0,1)) + 
        stat_summary(fun.y=mean, geom="point", fill="white", shape=21, size=3) + 
        theme(aspect.ratio=2) + 
        #scale_fill_brewer(palette="Set3") +
        scale_fill_manual(values=c("#ed6462", "#1f78b4")) +   # "#7F7F7F","#CCCCCC"
#    coord_flip() +
        theme_bw() +
        theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1, size=rel(0.9))) +  # adjust=4 is used to smoothen the density curve
        stat_summary(fun.y=mean, colour="white", geom="text", show.legend = FALSE, vjust=-0.7, size=4, aes(label=round(..y.., digits=2))) + 
        theme(legend.position="top") +
        facet_grid(marker_class ~ assay.type)


#==========
# Figure6f
#==========
options(repr.plot.width = 6,repr.plot.height = 8)

ggplot(df_L6_markers_Pv.as.ctrl_long_sep %>% filter(assay.type %in% c('CG','tmCG','hmCG')), aes(x=marker_class, y=assay.value, fill=marker_class)) + 
        geom_boxplot(outlier.shape = NA) + geom_jitter(height = 0, width = 0.2, shape=1, size=2, alpha = 2/10) +  
#       geom_jitter(height = 0, width = 0.2, shape=21, size=1, alpha=0.6) + 
#    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
#    scale_y_continuous(breaks=seq(0,1,0.1), labels=seq(0,1,0.1), limits=c(0,1)) + 
        stat_summary(fun.y=mean, geom="point", fill="white", shape=21, size=3) + 
        theme(aspect.ratio=2) + 
        #scale_fill_brewer(palette="Set3") +
        scale_fill_manual(values=c("#ed6462", "#1f78b4")) +   # "#7F7F7F","#CCCCCC"
#    coord_flip() +
        theme_bw() +
        theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1, size=rel(0.9))) +  # adjust=4 is used to smoothen the density curve
        stat_summary(fun.y=mean, colour="white", geom="text", show.legend = FALSE, vjust=-0.7, size=4, aes(label=round(..y.., digits=2))) + 
        theme(legend.position="top") +
        facet_grid(cell.type ~ assay.type)




ggsave("Fig6_L6.C1vsC2markers_boxplot_hmCG.CG.tmCG_RdBu_11.07.2022.pdf",  width=6,height=8)   

options(repr.plot.width = 6,repr.plot.height = 6)

ggplot(df_L6_markers_Pv.as.ctrl_long_sep %>% filter(assay.type %in% c('rna')), aes(x=cell.type, y=log2(assay.value), fill=cell.type)) + 
        geom_boxplot(outlier.shape = NA) + geom_jitter(height = 0, width = 0.2, shape=1, size=2, alpha = 2/10) + 
#        geom_boxplot(alpha=0.6,adjust=3) +  
#        geom_jitter(height = 0, width = 0.2, shape=21, size=1, alpha=0.6) + 
#    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
#    scale_y_continuous(breaks=seq(0,1,0.1), labels=seq(0,1,0.1), limits=c(0,1)) + 
        stat_summary(fun.y=mean, geom="point", fill="white", shape=21, size=3) + 
        theme(aspect.ratio=2) + 
        #scale_fill_brewer(palette="Set3") +
        scale_fill_manual(values=c("#ed6462", "#1f78b4")) +   # "#7F7F7F","#CCCCCC"
#    coord_flip() +
        theme_bw() +
        theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1, size=rel(0.9))) +  # adjust=4 is used to smoothen the density curve
        stat_summary(fun.y=mean, colour="white", geom="text", show.legend = FALSE, vjust=-0.7, size=4, aes(label=round(..y.., digits=2))) + 
        theme(legend.position="top") +
        facet_grid(. ~ marker_class)




options(repr.plot.width = 6,repr.plot.height = 6)

ggplot(df_L6_markers_Pv.as.ctrl_long_sep %>% filter(assay.type %in% c('CG.tmCG.ratio')), aes(x=cell.type, y=assay.value, fill=cell.type)) + 
        geom_boxplot(outlier.shape = NA) + geom_jitter(height = 0, width = 0.2, shape=1, size=2, alpha = 2/10) + 
#        geom_boxplot(alpha=0.6,adjust=3) +  
#        geom_jitter(height = 0, width = 0.2, shape=21, size=1, alpha=0.6) + 
#    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
#    scale_y_continuous(breaks=seq(0,1,0.1), labels=seq(0,1,0.1), limits=c(0,1)) + 
        stat_summary(fun.y=mean, geom="point", fill="white", shape=21, size=3) + 
        theme(aspect.ratio=2) + 
        #scale_fill_brewer(palette="Set3") +
        scale_fill_manual(values=c("#ed6462", "#1f78b4")) +   # "#7F7F7F","#CCCCCC"
#    coord_flip() +
        theme_bw() +
        theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1, size=rel(0.9))) +  # adjust=4 is used to smoothen the density curve
        stat_summary(fun.y=mean, colour="indianred2", geom="text", show.legend = FALSE, vjust=-0.7, size=4, aes(label=round(..y.., digits=2))) + 
        theme(legend.position="top") +
        facet_grid(. ~ marker_class)

options(repr.plot.width = 6,repr.plot.height = 6)

ggplot(df_L6_markers_Pv.as.ctrl_long_sep %>% filter(assay.type %in% c('hmCG.tmCG.ratio')), aes(x=cell.type, y=assay.value, fill=cell.type)) + 
        geom_boxplot(outlier.shape = NA) + geom_jitter(height = 0, width = 0.2, shape=1, size=2, alpha = 2/10) + 
#        geom_boxplot(alpha=0.6,adjust=3) +  
#        geom_jitter(height = 0, width = 0.2, shape=21, size=1, alpha=0.6) + 
#    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
#    scale_y_continuous(breaks=seq(0,1,0.1), labels=seq(0,1,0.1), limits=c(0,1)) + 
        stat_summary(fun.y=mean, geom="point", fill="white", shape=21, size=3) + 
        theme(aspect.ratio=2) + 
        #scale_fill_brewer(palette="Set3") +
        scale_fill_manual(values=c("#ed6462", "#1f78b4")) +   # "#7F7F7F","#CCCCCC"
#    coord_flip() +
        theme_bw() +
        theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1, size=rel(0.9))) +  # adjust=4 is used to smoothen the density curve
        stat_summary(fun.y=mean, colour="white", geom="text", show.legend = FALSE, vjust=-0.7, size=4, aes(label=round(..y.., digits=2))) + 
        theme(legend.position="top") +
        facet_grid(. ~ marker_class)

df_L6_markers_Pv.as.ctrl

options(repr.plot.width = 6,repr.plot.height = 8)

ggplot(df_L6_markers_Pv.as.ctrl, aes(x=marker_class, y=log10GeneLength.kb, fill=marker_class)) + 
        geom_boxplot(outlier.shape = NA) + geom_jitter(height = 0, width = 0.2, shape=1, size=2, alpha = 2/10) + 
#        geom_boxplot(alpha=0.6,adjust=3) +  
#        geom_jitter(height = 0, width = 0.2, shape=21, size=1, alpha=0.6) + 
#    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
#    scale_y_continuous(breaks=seq(0,1,0.1), labels=seq(0,1,0.1), limits=c(0,1)) + 
        stat_summary(fun.y=mean, geom="point", fill="white", shape=21, size=3) + 
        theme(aspect.ratio=1) + 
        scale_fill_brewer(palette="Set3") +
#        scale_fill_manual(values=c("#ed6462", "#1f78b4")) +   # "#7F7F7F","#CCCCCC"
#    coord_flip() +
        theme_bw() +
        theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1, size=rel(0.9))) +  # adjust=4 is used to smoothen the density curve
        stat_summary(fun.y=mean, colour="indianred2", geom="text", show.legend = FALSE, vjust=-0.7, size=4, aes(label=round(..y.., digits=2))) + 
#        stat_summary(fun.y=mean, colour="white", geom="text", show.legend = FALSE, vjust=-0.7, size=4, aes(label=round(..y.., digits=2))) + 
        theme(legend.position="top")
#        facet_grid(. ~ marker_class)





# https://dplyr.tidyverse.org/reference/summarise.html
mean_geneLength <- df_L6_markers_Pv.as.ctrl %>% 
    group_by(marker_class) %>%
    summarise(mean = mean(log10GeneLength.kb), n=n())

mean_geneLength



#==========
# Figure6g
#==========
options(repr.plot.width = 8, repr.plot.height = 4)
# Histogram with density plot
ggplot(df_L6_markers_Pv.as.ctrl, aes(x=log10GeneLength.kb, color=marker_class, fill=marker_class)) + 
#     geom_histogram(aes(y=..density..), position="identity", alpha=0.5)+
     geom_density(alpha=0.6) +
     theme(aspect.ratio=1) + 
#        scale_fill_brewer(palette="Set3") +
#        scale_fill_manual(values=c("#ed6462", "#1f78b4")) +   # "#7F7F7F","#CCCCCC"
#    coord_flip() +
      theme_bw() +
#      theme_classic() +
      geom_vline(data=mean_geneLength, aes(xintercept=mean, color=marker_class),linetype="dashed", size=0.75)+
      scale_color_manual(values=c("#ed6462", "#1f78b4"))+
      scale_fill_manual(values=c("#ed6462", "#1f78b4")) + 
#	  annotate("text",x=mean_geneLength[1,2, drop = TRUE],y=0,label=paste(format(round(mean_geneLength[1,2, drop = TRUE]*100,digits=1),nsmall=1), "%"),colour="#DC143C",vjust=-0.5,size=4) + 
#	  annotate("text",x=mean_geneLength[2,2, drop = TRUE],y=0,label=paste(format(round(mean_geneLength[2,2, drop = TRUE]*100,digits=1),nsmall=1), "%"),colour="#ADD8E6",vjust=-0.5,size=4)
	  annotate("text",x=mean_geneLength[1,2, drop = TRUE],y=1,label=format(round(mean_geneLength[1,2, drop = TRUE],digits=2),nsmall=1),colour="#ed6462",vjust=-0.5,size=4) + 
	  annotate("text",x=mean_geneLength[2,2, drop = TRUE],y=1,label=format(round(mean_geneLength[2,2, drop = TRUE],digits=2),nsmall=1),,colour="#1f78b4",vjust=-0.5,size=4)
#labs(title="Weight histogram plot",x="Weight(kg)", y = "Density")+




ggsave("Fig6_L6.C1vsC2markers_densityplot_11.08.2022.pdf",  width=8,height=4)   

options(repr.plot.width = 6,repr.plot.height = 8)

ggplot(df_Pv_markers_L6.as.ctrl_long_sep %>% filter(assay.type %in% c('CG','tmCG','hmCG')), aes(x=cell.type, y=assay.value, fill=cell.type)) + 
        geom_boxplot(outlier.shape = NA) + geom_jitter(height = 0, width = 0.2, shape=1, size=2, alpha = 2/10) + 
#        geom_boxplot(alpha=0.6,adjust=3) +  
#       geom_jitter(height = 0, width = 0.2, shape=21, size=1, alpha=0.6) + 
#    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
#    scale_y_continuous(breaks=seq(0,1,0.1), labels=seq(0,1,0.1), limits=c(0,1)) + 
        stat_summary(fun.y=mean, geom="point", fill="white", shape=21, size=3) + 
        theme(aspect.ratio=2) + 
        #scale_fill_brewer(palette="Set3") +
        scale_fill_manual(values=c("#ed6462", "#1f78b4")) +   # "#7F7F7F","#CCCCCC"
#    coord_flip() +
        theme_bw() +
        theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1, size=rel(0.9))) +  # adjust=4 is used to smoothen the density curve
        stat_summary(fun.y=mean, colour="white", geom="text", show.legend = FALSE, vjust=-0.7, size=4, aes(label=round(..y.., digits=2))) + 
        theme(legend.position="top") +
        facet_grid(marker_class ~ assay.type)



options(repr.plot.width = 6,repr.plot.height = 8)

ggplot(df_Pv_markers_L6.as.ctrl_long_sep %>% filter(assay.type %in% c('CG','tmCG','hmCG')), aes(x=marker_class, y=assay.value, fill=marker_class)) + 
        geom_boxplot(outlier.shape = NA) + geom_jitter(height = 0, width = 0.2, shape=1, size=2, alpha = 2/10) +  
#       geom_jitter(height = 0, width = 0.2, shape=21, size=1, alpha=0.6) + 
#    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
#    scale_y_continuous(breaks=seq(0,1,0.1), labels=seq(0,1,0.1), limits=c(0,1)) + 
        stat_summary(fun.y=mean, geom="point", fill="white", shape=21, size=3) + 
        theme(aspect.ratio=2) + 
        #scale_fill_brewer(palette="Set3") +
        scale_fill_manual(values=c("#ed6462", "#1f78b4")) +   # "#7F7F7F","#CCCCCC"
#    coord_flip() +
        theme_bw() +
        theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1, size=rel(0.9))) +  # adjust=4 is used to smoothen the density curve
        stat_summary(fun.y=mean, colour="white", geom="text", show.legend = FALSE, vjust=-0.7, size=4, aes(label=round(..y.., digits=2))) + 
        theme(legend.position="top") +
        facet_grid(cell.type ~ assay.type)





ggsave("Fig6_Pv.C1vsC2markers_boxplot_hmCG.CG.tmCG_RdBu_11.07.2022.pdf",  width=6,height=8)   

options(repr.plot.width = 6,repr.plot.height = 6)

ggplot(df_Pv_markers_L6.as.ctrl_long_sep %>% filter(assay.type %in% c('rna')), aes(x=cell.type, y=log2(assay.value), fill=cell.type)) + 
        geom_boxplot(outlier.shape = NA) + geom_jitter(height = 0, width = 0.2, shape=1, size=2, alpha = 2/10) + 
#        geom_boxplot(alpha=0.6,adjust=3) +  
#        geom_jitter(height = 0, width = 0.2, shape=21, size=1, alpha=0.6) + 
#    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
#    scale_y_continuous(breaks=seq(0,1,0.1), labels=seq(0,1,0.1), limits=c(0,1)) + 
        stat_summary(fun.y=mean, geom="point", fill="white", shape=21, size=3) + 
        theme(aspect.ratio=2) + 
        #scale_fill_brewer(palette="Set3") +
        scale_fill_manual(values=c("#ed6462", "#1f78b4")) +   # "#7F7F7F","#CCCCCC"
#    coord_flip() +
        theme_bw() +
        theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1, size=rel(0.9))) +  # adjust=4 is used to smoothen the density curve
        stat_summary(fun.y=mean, colour="white", geom="text", show.legend = FALSE, vjust=-0.7, size=4, aes(label=round(..y.., digits=2))) + 
        theme(legend.position="top") +
        facet_grid(. ~ marker_class)

options(repr.plot.width = 6,repr.plot.height = 6)

ggplot(df_Pv_markers_L6.as.ctrl_long_sep %>% filter(assay.type %in% c('CG.tmCG.ratio')), aes(x=cell.type, y=assay.value, fill=cell.type)) + 
        geom_boxplot(outlier.shape = NA) + geom_jitter(height = 0, width = 0.2, shape=1, size=2, alpha = 2/10) + 
#        geom_boxplot(alpha=0.6,adjust=3) +  
#        geom_jitter(height = 0, width = 0.2, shape=21, size=1, alpha=0.6) + 
#    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
#    scale_y_continuous(breaks=seq(0,1,0.1), labels=seq(0,1,0.1), limits=c(0,1)) + 
        stat_summary(fun.y=mean, geom="point", fill="white", shape=21, size=3) + 
        theme(aspect.ratio=2) + 
        #scale_fill_brewer(palette="Set3") +
        scale_fill_manual(values=c("#ed6462", "#1f78b4")) +   # "#7F7F7F","#CCCCCC"
#    coord_flip() +
        theme_bw() +
        theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1, size=rel(0.9))) +  # adjust=4 is used to smoothen the density curve
        stat_summary(fun.y=mean, colour="indianred2", geom="text", show.legend = FALSE, vjust=-0.7, size=4, aes(label=round(..y.., digits=2))) + 
        theme(legend.position="top") +
        facet_grid(. ~ marker_class)

options(repr.plot.width = 6,repr.plot.height = 8)

ggplot(df_Pv_markers_L6.as.ctrl, aes(x=marker_class, y=log10GeneLength.kb, fill=marker_class)) + 
        geom_boxplot(outlier.shape = NA) + geom_jitter(height = 0, width = 0.2, shape=1, size=2, alpha = 2/10) + 
#        geom_boxplot(alpha=0.6,adjust=3) +  
#        geom_jitter(height = 0, width = 0.2, shape=21, size=1, alpha=0.6) + 
#    scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20), limits=c(0,100)) + 
#    scale_y_continuous(breaks=seq(0,1,0.1), labels=seq(0,1,0.1), limits=c(0,1)) + 
        stat_summary(fun.y=mean, geom="point", fill="white", shape=21, size=3) + 
        theme(aspect.ratio=2) + 
        scale_fill_brewer(palette="Set3") +
#        scale_fill_manual(values=c("#ed6462", "#1f78b4")) +   # "#7F7F7F","#CCCCCC"
#    coord_flip() +
        theme_bw() +
        theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1, size=rel(0.9))) +  # adjust=4 is used to smoothen the density curve
        stat_summary(fun.y=mean, colour="indianred2", geom="text", show.legend = FALSE, vjust=-0.7, size=4, aes(label=round(..y.., digits=2))) + 
        theme(legend.position="top")
#        facet_grid(. ~ marker_class)



# https://dplyr.tidyverse.org/reference/summarise.html
mean_geneLength <- df_Pv_markers_L6.as.ctrl %>% 
    group_by(marker_class) %>%
    summarise(mean = mean(log10GeneLength.kb), n=n())

mean_geneLength


#==========
# Figure6g
#==========
options(repr.plot.width = 8, repr.plot.height = 4)

# Histogram with density plot
ggplot(df_Pv_markers_L6.as.ctrl, aes(x=log10GeneLength.kb, color=marker_class, fill=marker_class)) + 
#     geom_histogram(aes(y=..density..), position="identity", alpha=0.5)+
     geom_density(alpha=0.6) +
     theme(aspect.ratio=1) + 
#        scale_fill_brewer(palette="Set3") +
#        scale_fill_manual(values=c("#ed6462", "#1f78b4")) +   # "#7F7F7F","#CCCCCC"
#    coord_flip() +
      theme_bw() +
#      theme_classic() +
      geom_vline(data=mean_geneLength, aes(xintercept=mean, color=marker_class),linetype="dashed", size=0.75)+
      scale_color_manual(values=c("#ed6462", "#1f78b4"))+
      scale_fill_manual(values=c("#ed6462", "#1f78b4")) + 
#	  annotate("text",x=mean_geneLength[1,2, drop = TRUE],y=0,label=paste(format(round(mean_geneLength[1,2, drop = TRUE]*100,digits=1),nsmall=1), "%"),colour="#DC143C",vjust=-0.5,size=4) + 
#	  annotate("text",x=mean_geneLength[2,2, drop = TRUE],y=0,label=paste(format(round(mean_geneLength[2,2, drop = TRUE]*100,digits=1),nsmall=1), "%"),colour="#ADD8E6",vjust=-0.5,size=4)
	  annotate("text",x=mean_geneLength[1,2, drop = TRUE],y=1,label=format(round(mean_geneLength[1,2, drop = TRUE],digits=2),nsmall=1),colour="#ed6462",vjust=-0.5,size=4) + 
	  annotate("text",x=mean_geneLength[2,2, drop = TRUE],y=1,label=format(round(mean_geneLength[2,2, drop = TRUE],digits=2),nsmall=1),,colour="#1f78b4",vjust=-0.5,size=4)
#labs(title="Weight histogram plot",x="Weight(kg)", y = "Density")+





ggsave("Fig6_Pv.C1vsC2markers_densityplot_11.08.2022.pdf",  width=8,height=4)   

df_Pv_markers_L6.as.ctrl

head(df_5modalities_celltype_all.genes)
head(df_5modalities_celltype_all.genes_genelength)
colnames(df_5modalities_celltype_all.genes_genelength)

df_5modalities_celltype_all.genes_genelength_long <- melt(df_5modalities_celltype_all.genes_genelength, id.vars=c("GeneName",'chr','start','end','strand','id1','id2','Transcript_type','GeneLength','log10GeneLength.kb'), variable.name="mod_type", value.name="mod_level")
#df_L6_markers_long$GeneName <- factor(df_L6_markers_long$GeneName, levels = df_L6_markers$GeneName) # retain the order of original dataframe before melt
dim(df_5modalities_celltype_all.genes_genelength_long)
unique(df_5modalities_celltype_all.genes_genelength_long$mod_type)
df_5modalities_celltype_all.genes_genelength_long

df_5modalities_celltype_all.genes_genelength_long_sep <- df_5modalities_celltype_all.genes_genelength_long %>% separate(mod_type, c("cell_type","assay_type"), sep="_")
head(df_5modalities_celltype_all.genes_genelength_long_sep)
unique(df_5modalities_celltype_all.genes_genelength_long_sep$cell_type)
unique(df_5modalities_celltype_all.genes_genelength_long_sep$assay_type)

#rownames(df_5modalities_celltype_all.genes_genelength_long_sep) <- df_5modalities_celltype_all.genes_genelength_long_sep$GeneName
#df_5modalities_celltype_all.genes_genelength_long_sep[L6_rna.specific.gene, ]

# https://bioconda.github.io/recipes/bioconductor-clusterprofiler/README.html?highlight=clusterprofiler#package-package%20&#x27;bioconductor-clusterprofiler&#x27;
# conda install bioconductor-clusterprofiler
# conda install bioconductor-pathview

# http://bioconductor.org/packages/release/data/annotation/html/org.Mm.eg.db.html
# conda install bioconductor-org.mm.eg.db
library(org.Mm.eg.db)

keytypes(org.Mm.eg.db)

library(clusterProfiler)

length(rownames(df_5modalities_celltype_all.genes))

L2.3_rna.specific.gene<-intersect(neu.rna.mak2 %>% filter(cluster=="L2/3") %>% pull(gene),rownames(df_4modalities_celltype_all.genes))
length(L2.3_rna.specific.gene)
L6_rna.specific.gene<-intersect(neu.rna.mak2 %>% filter(cluster=="L6") %>% pull(gene),rownames(df_4modalities_celltype_all.genes))
length(L6_rna.specific.gene)
Vip.Ndnf_rna.specific.gene<-intersect(neu.rna.mak2 %>% filter(cluster=="Vip/Ndnf") %>% pull(gene),rownames(df_4modalities_celltype_all.genes))
length(Vip.Ndnf_rna.specific.gene)
Sst_rna.specific.gene<-intersect(neu.rna.mak2 %>% filter(cluster=="Sst") %>% pull(gene),rownames(df_4modalities_celltype_all.genes))
length(Sst_rna.specific.gene)
Pv_rna.specific.gene<-intersect(neu.rna.mak2 %>% filter(cluster=="Pv") %>% pull(gene),rownames(df_4modalities_celltype_all.genes))
length(Pv_rna.specific.gene)

class(L6_rna.specific.gene)
head(L6_rna.specific.gene)

class(rownames(df_5modalities_celltype_all.genes))

all.genes<-df_5modalities_celltype_all.genes$L6_rna
names(all.genes)<-rownames(df_5modalities_celltype_all.genes)

class(all.genes)
tail(all.genes)
length(unique(names(all.genes)))
length(all.genes)

# https://support.bioconductor.org/p/63502/
## 1. Download Gene Info directly from NCBI
### wget ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz

## 2. After downloading, set taxonomy ID for column 1 (e.g. 9606=Hs, 10090=Mm), and extract relevant columns (2=GeneID, 3=Symbol, 10=Type_of_gene). Do this using AWK/*nix.
### zcat gene_info.gz | awk 'BEGIN {FS="\t"} $1==10090 {print $2 "\t" $3 "\t" $10}' > geneInfo_Mm.txt

## 3. Then load resulting file into R-session:
GeneIDs<-read.table("geneInfo_Mm.txt",sep="\t",quote="\"",na.strings="-",fill=TRUE, col.names=c("GeneID","Symbol","TypeOfGene"))
head(GeneIDs)

df.all.genes <- as.data.frame(all.genes)
df.all.genes$gene_name <- rownames(df.all.genes)
df.all.genes <- rename(df.all.genes, L6_rna = all.genes)
head(df.all.genes)

df.all.genes_common <- inner_join(df.all.genes, GeneIDs, by=c("gene_name"="Symbol"))
head(df.all.genes_common)
dim(df.all.genes_common)

all.GeneIDs <- as.character(df.all.genes_common$GeneID)

df.L6.genes <- as.data.frame(L6_rna.specific.gene)
head(df.L6.genes)

df.L6.genes_common <- inner_join(df.L6.genes, GeneIDs, by=c("L6_rna.specific.gene"="Symbol"))
head(df.L6.genes_common)
dim(df.L6.genes_common)

L6.GeneIDs <- as.character(df.L6.genes_common$GeneID)


df.Pv.genes <- as.data.frame(Pv_rna.specific.gene)

df.Pv.genes_common <- inner_join(df.Pv.genes, GeneIDs, by=c("Pv_rna.specific.gene"="Symbol"))
dim(df.Pv.genes_common)
Pv.GeneIDs <- as.character(df.Pv.genes_common$GeneID)

# https://hbctraining.github.io/Training-modules/DGE-functional-analysis/lessons/functional_analysis_2019.html
## Run GO enrichment analysis 
L6_ego <- enrichGO(gene = L6.GeneIDs, 
                universe = all.GeneIDs,
                keyType = "ENTREZID",
                OrgDb = org.Mm.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)

L6_cluster_summary <- data.frame(L6_ego)
head(L6_cluster_summary)

Pv_ego <- enrichGO(gene = Pv.GeneIDs, 
                universe = all.GeneIDs,
#                keyType = "REFSEQ",
                keyType = "ENTREZID",
                OrgDb = org.Mm.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)

Pv_cluster_summary <- data.frame(Pv_ego)
head(Pv_cluster_summary)

options(repr.plot.width=10,repr.plot.height=10)

L6 <- dotplot(L6_ego, showCategory=30)
Pv <- dotplot(Pv_ego, showCategory=30)

L6
Pv

cluster.order

gene.order2.rna <- gene.neu.order2

mat_hmCG_celltype_neu_norm <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name.6groups == x]
    m <- hmC_matrix[gene.order2.rna,cells.in]
    mat <- sweep(m,2,colSums(m),`/`)
    return(rowMeans(mat[gene.order2.rna,cells.in]))
})

mat_tmC_celltype_neu_norm <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name.6groups == x]
    m <- tmC_matrix[gene.order2.rna,cells.in]
    mat <- sweep(m,2,colSums(m),`/`)
    return(rowMeans(mat[gene.order2.rna,cells.in]))
    #return(rowMeans(tmC_matrix[gene.order2.rna,cells.in]))
})

mat_mCG_celltype_neu_norm <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name.6groups == x]
    m <- mC_matrix[gene.order2.rna,cells.in]
    mat <- sweep(m,2,colSums(m),`/`)
    return(rowMeans(mat[gene.order2.rna,cells.in]))
    #return(rowMeans(mC_matrix[gene.order2.rna,cells.in]))
})

mat_mCH_celltype_neu_norm <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name.6groups == x]
    m <- mCH_matrix[gene.order2.rna,cells.in]
    mat <- sweep(m,2,colSums(m),`/`)
    return(rowMeans(mat[gene.order2.rna,cells.in]))
    #return(rowMeans(mCH_matrix[gene.order2.rna,cells.in]))
})

mat_ratio_celltype_neu_norm <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name.6groups == x]
    m <- ratio_matrix[gene.order2.rna,cells.in]
    mat <- sweep(m,2,colSums(m),`/`)
    return(rowMeans(mat[gene.order2.rna,cells.in]))
    #return(rowMeans(ratio_matrix[gene.order2.rna,cells.in]))
})

mat_rna_celltype_neu_norm <- sapply(cluster.order,function(x){
    cells.in <- rownames(snRNA3.neu@meta.data)[snRNA3.neu@meta.data$res.name.6groups == x]
    mat <- expm1(snRNA3.neu@assays$RNA@data)
    return(Matrix::rowMeans(mat[gene.order2.rna,cells.in]))
})

options(repr.plot.width=3,repr.plot.height=3)

h1 <- pheatmap::pheatmap(mat_tmC_celltype_neu_norm,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="tmCG",color=rev(colorRampPalette(brewer.pal(n = 9,name="BuPu"))(50)))

h2 <- pheatmap::pheatmap(mat_mCG_celltype_neu_norm,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="BSmCG",color=rev(colorRampPalette(brewer.pal(n = 9,name="BuPu"))(50)))

h3 <- pheatmap::pheatmap(mat_hmCG_celltype_neu_norm,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="hmCG",color=rev(colorRampPalette(brewer.pal(n = 9,name="BuPu"))(50)))

h4 <- pheatmap::pheatmap(mat_ratio_celltype_neu_norm,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="hmCG/BSmCG",color=colorRampPalette(brewer.pal(n = 9,name="PuRd"))(50))

h5 <- pheatmap::pheatmap(mat_mCH_celltype_neu_norm,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="BSmCH",color=rev(colorRampPalette(brewer.pal(n = 9,name="BuPu"))(50)))

h6 <- pheatmap::pheatmap(mat_rna_celltype_neu_norm,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="RNA",colours = rev(colorRampPalette(brewer.pal(n = 11, name = "RdBu"))(50)))

h6 <- pheatmap::pheatmap(mat_rna_celltype_neu_norm,cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,
                   scale="row",clustering_method = "ward.D2",main="RNA",colours = colorRampPalette(c("blue","white","red"))(50))

table(pbmc.snmc3@meta.data$res.name.6groups)
table(pbmc.snmc3@meta.data$res.name2)

table(snRNA3@meta.data$res.name.6groups)
table(snRNA3@meta.data$res.name)

cluster.order <- c("Ex","Inh","Striatum","Astro","Oligo","MG")

gene.order2.rna <- intersect(rownames(hmC_matrix),rownames(snRNA3@assays$RNA@counts))
length(gene.order2.rna) # Peng: 4255 genes; Hao: 5182 genes (90% cells); 6716 genes (80% cells)

# old analysis:
gene.order2.rna <- intersect(rownames(hmC_matrix),rownames(snRNA3@assays$RNA@counts))
length(gene.order2.rna) # Peng: 4255 genes; Hao: 5182 genes (90% cells); 6716 genes (80% cells)

mat_hmCG_celltype_6g_raw_scatter <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name.6groups == x]
    return(rowMeans(hmC_matrix[gene.order2.rna,cells.in]))
})

mat_tmCG_celltype_6g_raw_scatter <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name.6groups == x]
    return(rowMeans(tmC_matrix[gene.order2.rna,cells.in]))
})

mat_mCG_celltype_6g_raw_scatter <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name.6groups == x]
    return(rowMeans(mC_matrix[gene.order2.rna,cells.in]))
})

mat_mCH_celltype_6g_raw_scatter <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name.6groups == x]
    return(rowMeans(mCH_matrix[gene.order2.rna,cells.in]))
})

mat_ratio_celltype_6g_raw_scatter <- sapply(cluster.order,function(x){
    cells.in <- rownames(pbmc.snmc3@meta.data)[pbmc.snmc3@meta.data$res.name.6groups == x]
    return(rowMeans(ratio_matrix[gene.order2.rna,cells.in]))
})

mat_rna_celltype_6g_raw_scatter <- sapply(cluster.order,function(x){
    cells.in <- rownames(snRNA3@meta.data)[snRNA3@meta.data$res.name.6groups == x]
    mat <- expm1(snRNA3@assays$RNA@data)
    return(Matrix::rowMeans(mat[gene.order2.rna,cells.in]))
})

colnames(mat_hmCG_celltype_6g_raw_scatter) <- paste(colnames(mat_hmCG_celltype_6g_raw_scatter),"hmCG",sep="_")
colnames(mat_tmCG_celltype_6g_raw_scatter) <- paste(colnames(mat_tmCG_celltype_6g_raw_scatter),"tmCG",sep="_")
colnames(mat_mCG_celltype_6g_raw_scatter) <- paste(colnames(mat_mCG_celltype_6g_raw_scatter),"BSmCG",sep="_")
colnames(mat_mCH_celltype_6g_raw_scatter) <- paste(colnames(mat_mCH_celltype_6g_raw_scatter),"BSmCH",sep="_")
colnames(mat_ratio_celltype_6g_raw_scatter) <- paste(colnames(mat_ratio_celltype_6g_raw_scatter),"ratio",sep="_")

colnames(mat_rna_celltype_6g_raw_scatter) <- paste(colnames(mat_rna_celltype_6g_raw_scatter),"RNA",sep="_")

head(mat_rna_celltype_6g_raw_scatter)

library(ggrastr)

library(cowplot)

dat1 <- mat_hmCG_celltype_6g_raw_scatter
dat2 <- mat_mCG_celltype_6g_raw_scatter

gene.common <- intersect(rownames(dat1),rownames(dat2))

dat1.2 <- dat1[gene.common,]
dat2.2 <- dat2[gene.common,]

mat.comb <- cbind(dat1.2,dat2.2)
mat.comb2 <- reshape2::melt(mat.comb,by="row.names")
names(mat.comb2) <- c("gene","type","value")
mat.comb2$celltype <- apply(mat.comb2,1,function(x){
    unlist(strsplit(x[2],"_"))[1]
})
mat.comb2$tech <- apply(mat.comb2,1,function(x){
    unlist(strsplit(x[2],"_"))[2]
})
mat.comb3 <- reshape2::dcast(mat.comb2,gene+celltype~tech,value.var = "value")

cluster.order <- c("Ex","Inh","Striatum","Astro","Oligo","MG")

mat.comb3$celltype <- factor(mat.comb3$celltype,levels=cluster.order)

sub.stat <- sapply(cluster.order,function(i){
    tmp <- mat.comb3 %>% dplyr::filter(celltype == i) %>% droplevels
    ts <- cor.test(tmp$hmCG, tmp$BSmCG)
    return(c(as.numeric(ts$p.value),as.numeric(ts$estimate)))
})

sub.stat <- t(sub.stat)
sub.stat <- as.data.frame(sub.stat)
names(sub.stat) <- c("pvalue","cor")
sub.stat$celltype <- rownames(sub.stat)
sub.stat$celltype <- factor(sub.stat$celltype,levels=cluster.order)

options(repr.plot.width=21,repr.plot.height=3)

px1 <- ggplot(mat.comb3,aes(BSmCG,hmCG)) + 
       geom_point_rast(size = .5,alpha=.1,raster.dpi = 300) + 
       geom_smooth(method = "lm", se = F) +
       geom_text(data = sub.stat,
                 aes(x = .1,y=0.45,label = formatC(pvalue,digits = 2)),color = "red",hjust = 0) +
       geom_text(data = sub.stat,
                 aes(x = .1,y=0.35,label = round(cor,2)),color = "red",hjust = 0) +
       facet_grid(~celltype) + theme_bw() #+ theme_cowplot()
px1

dat1 <- mat_hmCG_celltype_6g_raw_scatter
dat2 <- mat_tmCG_celltype_6g_raw_scatter

gene.common <- intersect(rownames(dat1),rownames(dat2))

dat1.2 <- dat1[gene.common,]
dat2.2 <- dat2[gene.common,]

mat.comb <- cbind(dat1.2,dat2.2)
mat.comb2 <- reshape2::melt(mat.comb,by="row.names")
names(mat.comb2) <- c("gene","type","value")
mat.comb2$celltype <- apply(mat.comb2,1,function(x){
    unlist(strsplit(x[2],"_"))[1]
})
mat.comb2$tech <- apply(mat.comb2,1,function(x){
    unlist(strsplit(x[2],"_"))[2]
})
mat.comb3 <- reshape2::dcast(mat.comb2,gene+celltype~tech,value.var = "value")

cluster.order <- c("Ex","Inh","Striatum","Astro","Oligo","MG")

mat.comb3$celltype <- factor(mat.comb3$celltype,levels=cluster.order)

sub.stat <- sapply(cluster.order,function(i){
    tmp <- mat.comb3 %>% dplyr::filter(celltype == i) %>% droplevels
    ts <- cor.test(tmp$hmCG, tmp$tmCG)
    return(c(as.numeric(ts$p.value),as.numeric(ts$estimate)))
})

sub.stat <- t(sub.stat)
sub.stat <- as.data.frame(sub.stat)
names(sub.stat) <- c("pvalue","cor")
sub.stat$celltype <- rownames(sub.stat)
sub.stat$celltype <- factor(sub.stat$celltype,levels=cluster.order)

options(repr.plot.width=21,repr.plot.height=3)

px2 <- ggplot(mat.comb3,aes(tmCG,hmCG)) + 
       geom_point_rast(size = .5,alpha=.1,raster.dpi = 300) + 
       geom_smooth(method = "lm", se = FALSE) +
       geom_text(data = sub.stat,
                 aes(x = .1,y=0.45,label = formatC(pvalue,digits = 2)),color = "red",hjust = 0) +
       geom_text(data = sub.stat,
                 aes(x = .1,y=0.35,label = round(cor,2)),color = "red",hjust = 0) +
       facet_grid(~celltype) + theme_bw() #+ theme_cowplot()
px2

dat1 <- mat_mCG_celltype_6g_raw_scatter
dat2 <- mat_tmCG_celltype_6g_raw_scatter

gene.common <- intersect(rownames(dat1),rownames(dat2))

dat1.2 <- dat1[gene.common,]
dat2.2 <- dat2[gene.common,]

mat.comb <- cbind(dat1.2,dat2.2)
mat.comb2 <- reshape2::melt(mat.comb,by="row.names")
names(mat.comb2) <- c("gene","type","value")
mat.comb2$celltype <- apply(mat.comb2,1,function(x){
    unlist(strsplit(x[2],"_"))[1]
})
mat.comb2$tech <- apply(mat.comb2,1,function(x){
    unlist(strsplit(x[2],"_"))[2]
})
mat.comb3 <- reshape2::dcast(mat.comb2,gene+celltype~tech,value.var = "value")

cluster.order <- c("Ex","Inh","Striatum","Astro","Oligo","MG")

mat.comb3$celltype <- factor(mat.comb3$celltype,levels=cluster.order)

sub.stat <- sapply(cluster.order,function(i){
    tmp <- mat.comb3 %>% dplyr::filter(celltype == i) %>% droplevels
    ts <- cor.test(tmp$BSmCG, tmp$tmCG)
    return(c(as.numeric(ts$p.value),as.numeric(ts$estimate)))
})

sub.stat <- t(sub.stat)
sub.stat <- as.data.frame(sub.stat)
names(sub.stat) <- c("pvalue","cor")
sub.stat$celltype <- rownames(sub.stat)
sub.stat$celltype <- factor(sub.stat$celltype,levels=cluster.order)

options(repr.plot.width=21,repr.plot.height=3)

px3 <- ggplot(mat.comb3,aes(tmCG,BSmCG)) + 
       geom_point_rast(size = .1,alpha=.1,raster.dpi = 300) + 
       geom_smooth(method = "lm", se = FALSE) +
       geom_text(data = sub.stat,
                 aes(x = .1,y=0.45,label = formatC(pvalue,digits = 2)),color = "red",hjust = 0) +
       geom_text(data = sub.stat,
                 aes(x = .1,y=0.35,label = round(cor,2)),color = "red",hjust = 0) +
       facet_grid(~celltype) + theme_bw() #+ theme_cowplot()
px3

options(repr.plot.width=18,repr.plot.height=9)

plot_grid(px1,px2,px3,nrow=3)

#old analysis 
options(repr.plot.width=18,repr.plot.height=9)

plot_grid(px1,px2,px3,nrow=3)

ggsave("20211029_scatterplot_with_stat.pdf",width=18,height=9)

dat1 <- mat_mCG_celltype_6g_raw_scatter
#dat2 <- log1p(mat_rna_celltype_6g_raw_scatter)
dat2 <- log2(mat_rna_celltype_6g_raw_scatter)



gene.common <- intersect(rownames(dat1),rownames(dat2))

dat1.2 <- dat1[gene.common,]
dat2.2 <- dat2[gene.common,]

mat.comb <- cbind(dat1.2,dat2.2)
mat.comb2 <- reshape2::melt(mat.comb,by="row.names")
names(mat.comb2) <- c("gene","type","value")
mat.comb2$celltype <- apply(mat.comb2,1,function(x){
    unlist(strsplit(x[2],"_"))[1]
})
mat.comb2$tech <- apply(mat.comb2,1,function(x){
    unlist(strsplit(x[2],"_"))[2]
})
mat.comb3 <- reshape2::dcast(mat.comb2,gene+celltype~tech,value.var = "value")

cluster.order <- c("Ex","Inh","Striatum","Astro","Oligo","MG")

mat.comb3$celltype <- factor(mat.comb3$celltype,levels=cluster.order)

sub.stat <- sapply(cluster.order,function(i){
    tmp <- mat.comb3 %>% dplyr::filter(celltype == i) %>% droplevels
    ts <- cor.test(tmp$BSmCG, tmp$RNA)
    return(c(as.numeric(ts$p.value),as.numeric(ts$estimate)))
})

sub.stat <- t(sub.stat)
sub.stat <- as.data.frame(sub.stat)
names(sub.stat) <- c("pvalue","cor")
sub.stat$celltype <- rownames(sub.stat)
sub.stat$celltype <- factor(sub.stat$celltype,levels=cluster.order)

options(repr.plot.width=21,repr.plot.height=3)

px4 <- ggplot(mat.comb3,aes(RNA,BSmCG)) + 
       geom_point_rast(size = .5,alpha=.1,raster.dpi = 300) + 
       geom_smooth(method = "lm", se = FALSE) +
       geom_text(data = sub.stat,
                 aes(x = .1,y=0.45,label = formatC(pvalue,digits = 2)),color = "red",hjust = 0) +
       geom_text(data = sub.stat,
                 aes(x = .1,y=0.35,label = round(cor,2)),color = "red",hjust = 0) +
       facet_grid(~celltype) + theme_bw() #+ theme_cowplot()
px4

dat1 <- mat_tmCG_celltype_6g_raw_scatter
#dat2 <- log1p(mat_rna_celltype_6g_raw_scatter)
dat2 <- log2(mat_rna_celltype_6g_raw_scatter)

gene.common <- intersect(rownames(dat1),rownames(dat2))

dat1.2 <- dat1[gene.common,]
dat2.2 <- dat2[gene.common,]

mat.comb <- cbind(dat1.2,dat2.2)
mat.comb2 <- reshape2::melt(mat.comb,by="row.names")
names(mat.comb2) <- c("gene","type","value")
mat.comb2$celltype <- apply(mat.comb2,1,function(x){
    unlist(strsplit(x[2],"_"))[1]
})
mat.comb2$tech <- apply(mat.comb2,1,function(x){
    unlist(strsplit(x[2],"_"))[2]
})
mat.comb3 <- reshape2::dcast(mat.comb2,gene+celltype~tech,value.var = "value")

cluster.order <- c("Ex","Inh","Striatum","Astro","Oligo","MG")

mat.comb3$celltype <- factor(mat.comb3$celltype,levels=cluster.order)

sub.stat <- sapply(cluster.order,function(i){
    tmp <- mat.comb3 %>% dplyr::filter(celltype == i) %>% droplevels
    ts <- cor.test(tmp$RNA, tmp$tmCG)
    return(c(as.numeric(ts$p.value),as.numeric(ts$estimate)))
})

sub.stat <- t(sub.stat)
sub.stat <- as.data.frame(sub.stat)
names(sub.stat) <- c("pvalue","cor")
sub.stat$celltype <- rownames(sub.stat)
sub.stat$celltype <- factor(sub.stat$celltype,levels=cluster.order)

options(repr.plot.width=21,repr.plot.height=3)

px5 <- ggplot(mat.comb3,aes(RNA,tmCG)) + 
       geom_point_rast(size = .5,alpha=.1,raster.dpi = 300) + 
       geom_smooth(method = "lm", se = FALSE) +
       geom_text(data = sub.stat,
                 aes(x = .1,y=0.45,label = formatC(pvalue,digits = 2)),color = "red",hjust = 0) +
       geom_text(data = sub.stat,
                 aes(x = .1,y=0.35,label = round(cor,2)),color = "red",hjust = 0) +
       facet_grid(~celltype) + theme_bw() #+ theme_cowplot()
px5

dat1 <- mat_hmCG_celltype_6g_raw_scatter
#dat2 <- log1p(mat_rna_celltype_6g_raw_scatter)
dat2 <- log2(mat_rna_celltype_6g_raw_scatter)



gene.common <- intersect(rownames(dat1),rownames(dat2))

dat1.2 <- dat1[gene.common,]
dat2.2 <- dat2[gene.common,]

mat.comb <- cbind(dat1.2,dat2.2)
mat.comb2 <- reshape2::melt(mat.comb,by="row.names")
names(mat.comb2) <- c("gene","type","value")
mat.comb2$celltype <- apply(mat.comb2,1,function(x){
    unlist(strsplit(x[2],"_"))[1]
})
mat.comb2$tech <- apply(mat.comb2,1,function(x){
    unlist(strsplit(x[2],"_"))[2]
})
mat.comb3 <- reshape2::dcast(mat.comb2,gene+celltype~tech,value.var = "value")

cluster.order <- c("Ex","Inh","Striatum","Astro","Oligo","MG")

mat.comb3$celltype <- factor(mat.comb3$celltype,levels=cluster.order)

sub.stat <- sapply(cluster.order,function(i){
    tmp <- mat.comb3 %>% dplyr::filter(celltype == i) %>% droplevels
    ts <- cor.test(tmp$RNA, tmp$hmCG)
    return(c(as.numeric(ts$p.value),as.numeric(ts$estimate)))
})

sub.stat <- t(sub.stat)
sub.stat <- as.data.frame(sub.stat)
names(sub.stat) <- c("pvalue","cor")
sub.stat$celltype <- rownames(sub.stat)
sub.stat$celltype <- factor(sub.stat$celltype,levels=cluster.order)

options(repr.plot.width=21,repr.plot.height=3)

px6 <- ggplot(mat.comb3,aes(RNA,hmCG)) + 
       geom_point_rast(size = .5,alpha=.1,raster.dpi = 300) + 
#       stat_density_2d(geom='polygon', aes(fill=..level..), base='identity', bins=100) +
       geom_smooth(method = "lm", se = T) +
       geom_text(data = sub.stat,
                 aes(x = .1,y=0.45,label = formatC(pvalue,digits = 2)),color = "red",hjust = 0) +
       geom_text(data = sub.stat,
                 aes(x = .1,y=0.35,label = round(cor,2)),color = "red",hjust = 0) +
       facet_grid(~celltype) + theme_bw() #+ theme_cowplot()
px6

dat1 <- mat_mCH_celltype_6g_raw_scatter
#dat2 <- log1p(mat_rna_celltype_6g_raw_scatter)
dat2 <- log2(mat_rna_celltype_6g_raw_scatter)


gene.common <- intersect(rownames(dat1),rownames(dat2))

dat1.2 <- dat1[gene.common,]
dat2.2 <- dat2[gene.common,]

mat.comb <- cbind(dat1.2,dat2.2)
mat.comb2 <- reshape2::melt(mat.comb,by="row.names")
names(mat.comb2) <- c("gene","type","value")
mat.comb2$celltype <- apply(mat.comb2,1,function(x){
    unlist(strsplit(x[2],"_"))[1]
})
mat.comb2$tech <- apply(mat.comb2,1,function(x){
    unlist(strsplit(x[2],"_"))[2]
})
mat.comb3 <- reshape2::dcast(mat.comb2,gene+celltype~tech,value.var = "value")

cluster.order <- c("Ex","Inh","Striatum","Astro","Oligo","MG")

mat.comb3$celltype <- factor(mat.comb3$celltype,levels=cluster.order)

sub.stat <- sapply(cluster.order,function(i){
    tmp <- mat.comb3 %>% dplyr::filter(celltype == i) %>% droplevels
    ts <- cor.test(tmp$RNA, tmp$BSmCH)
    return(c(as.numeric(ts$p.value),as.numeric(ts$estimate)))
})

sub.stat <- t(sub.stat)
sub.stat <- as.data.frame(sub.stat)
names(sub.stat) <- c("pvalue","cor")
sub.stat$celltype <- rownames(sub.stat)
sub.stat$celltype <- factor(sub.stat$celltype,levels=cluster.order)

options(repr.plot.width=21,repr.plot.height=3)

px7 <- ggplot(mat.comb3,aes(RNA,BSmCH)) + 
       geom_point_rast(size = .5,alpha=.1,raster.dpi = 300) + 
       geom_smooth(method = "lm", se = FALSE) +
       coord_cartesian(xlim = c(0,5),ylim = c(0,0.1)) +
       geom_text(data = sub.stat,
                 aes(x = .1,y=0.09,label = formatC(pvalue,digits = 2)),color = "red",hjust = 0) +
       geom_text(data = sub.stat,
                 aes(x = .1,y=0.07,label = round(cor,2)),color = "red",hjust = 0) +
       facet_grid(~celltype) + theme_bw() #+ theme_cowplot()
px7

dat1 <- mat_ratio_celltype_6g_raw_scatter
#dat2 <- log1p(mat_rna_celltype_6g_raw_scatter)
dat2 <- log2(mat_rna_celltype_6g_raw_scatter)


gene.common <- intersect(rownames(dat1),rownames(dat2))

dat1.2 <- dat1[gene.common,]
dat2.2 <- dat2[gene.common,]

mat.comb <- cbind(dat1.2,dat2.2)
mat.comb2 <- reshape2::melt(mat.comb,by="row.names")
names(mat.comb2) <- c("gene","type","value")
mat.comb2$celltype <- apply(mat.comb2,1,function(x){
    unlist(strsplit(x[2],"_"))[1]
})
mat.comb2$tech <- apply(mat.comb2,1,function(x){
    unlist(strsplit(x[2],"_"))[2]
})
mat.comb3 <- reshape2::dcast(mat.comb2,gene+celltype~tech,value.var = "value")

cluster.order <- c("Ex","Inh","Striatum","Astro","Oligo","MG")

mat.comb3$celltype <- factor(mat.comb3$celltype,levels=cluster.order)

sub.stat <- sapply(cluster.order,function(i){
    tmp <- mat.comb3 %>% dplyr::filter(celltype == i) %>% droplevels
    ts <- cor.test(tmp$RNA, tmp$ratio)
    return(c(as.numeric(ts$p.value),as.numeric(ts$estimate)))
})

sub.stat <- t(sub.stat)
sub.stat <- as.data.frame(sub.stat)
names(sub.stat) <- c("pvalue","cor")
sub.stat$celltype <- rownames(sub.stat)
sub.stat$celltype <- factor(sub.stat$celltype,levels=cluster.order)

options(repr.plot.width=21,repr.plot.height=3)

px8 <- ggplot(mat.comb3,aes(RNA,ratio)) + 
       geom_point_rast(size = .5,alpha=.1,raster.dpi = 300) + 
       geom_smooth(method = "lm", se = FALSE) +
       geom_text(data = sub.stat,
                 aes(x = .1,y=0.45,label = formatC(pvalue,digits = 2)),color = "red",hjust = 0) +
       geom_text(data = sub.stat,
                 aes(x = .1,y=0.35,label = round(cor,2)),color = "red",hjust = 0) +
       facet_grid(~celltype) + theme_bw() #+ theme_cowplot()

px8

options(repr.plot.width=18,repr.plot.height=24)
# plot_grid(px1,px2,px3,px4,px5,px6,px7,px8,nrow=8)
plot_grid(px1,px2,px4,px5,px6,px7,px8, nrow=7)

ggsave("HW-2022.11.06_scatterplot_log2RNA.pdf",width=18,height=24)

options(repr.plot.width=18,repr.plot.height=24)
# plot_grid(px1,px2,px3,px4,px5,px6,px7,px8,nrow=8)
plot_grid(px1,px2,px4,px5,px6,px7,px8, nrow=7)

ggsave("HW-2022.01.06_scatterplot_with_stat.pdf",width=18,height=24)


