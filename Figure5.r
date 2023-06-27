#!/usr/bin/Rscript
# or open it with jupyter notebook
# Title: "NATURE BIOTECHNOLOGY README"
# Author: "Peng"

library("Seurat")

#==========
# Fig5c
#==========

# load the CpG clsutering results
pbmc.snmc <- readRDS(file="20210921_with_final_annotation_pbmc.snmc.rds")

pbmc.snmc@meta.data$res.name3 <- plyr::mapvalues(pbmc.snmc@meta.data$res.name2,
                from = c("L2/3","L4","L4/5","L6","DL","Vip/Ndnf","Pv","Sst","Striatum","In"),
                to = c("Ex-up","Ex-mid","Ex-mid","Ex-deep","Ex-deep","Inh-Vip/Ndnf","Inh-Pv","Inh-Sst","Striatum","In"))

Idents(pbmc.snmc) <- "res.name3"

pbmc.snmc@meta.data$res.name.6g <- plyr::mapvalues(pbmc.snmc@meta.data$res.name2,
                from = c("L2/3","L4","L4/5","L6","DL","Vip/Ndnf","Pv","Sst","Striatum","In"),
                to = c("Ex","Ex","Ex","Ex","Ex","Inh","Inh","Inh","Striatum","In"))

Idents(pbmc.snmc) <- "res.name.6g"

dff <- readRDS(file="20210921_metacell_snhmCseq.rds")
dff2 <- dff[,c('cell_snhmCseq','res.name2','cell_snmCseq2','cell_in_snmCseq2_matrix','cell_in_snhmCseq_matrix')]
dff2$cellID <- paste("Paired",1:545,sep="")
dff$cellID <- dff2$cellID
table(dff2$res.name)
table(dff2$res.name2)

# get the common cell id here
pbmc.snmc2 <- subset(pbmc.snmc,cells = dff2$cell_snmCseq2)
pairedID <- dff$cellID
names(pairedID) <- dff$cell_snmCseq2
# add paired ID here
pbmc.snmc2@meta.data$cellID <- pairedID[rownames(pbmc.snmc2@meta.data)]

Idents(pbmc.snmc2) <- 'res.name2'
pbmc.snmc3 <- subset(pbmc.snmc2,idents = c("In"),invert=T)
pbmc.snmc3@meta.data <- droplevels(pbmc.snmc3@meta.data)

table(pbmc.snmc3@meta.data$res.name2)

mC_rate_rename <- readRDS(file="20210930_GeneBody_mC_rate_rename.rds")
mC_nsite_rename <- readRDS(file="20210930_GeneBody_mC_nsite_rename.rds")
hmC_rate_rename <- readRDS(file="20210930_GeneBody_hmC_rate_rename.rds")
hmC_nsite_rename <- readRDS(file="20210930_GeneBody_hmC_nsite_rename.rds")

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


tmC_rate_post_impute <- filter_by_nsite_neg_percent_impute(mC_nsite=mC_nsite_rename,
                                                           hmC_nsite=hmC_nsite_rename,
                                                           mC_rate=mC_rate_rename,
                                                           hmC_rate=hmC_rate_rename,
                                                           percent=.9,
                                                           nsite_cutoff=5)

tmC_rate_post_impute2 <- filter_by_nsite_neg_percent_impute(mC_nsite=mC_nsite_rename,
                                                           hmC_nsite=hmC_nsite_rename,
                                                           mC_rate=mC_rate_rename,
                                                           hmC_rate=hmC_rate_rename,
                                                           percent=.8,
                                                           nsite_cutoff=2)

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

hmC_rate_post_impute <- filter_impute3(hmC_nsite_rename,
                                       hmC_rate_rename,
                                       nSite = 2,
                                       rownames(tmC_rate_post_impute2))
dim(hmC_rate_post_impute) # 12583 545

mC_rate_post_impute <- filter_impute3(mC_nsite_rename,
                                       mC_rate_rename,
                                       nSite = 2,
                                       rownames(tmC_rate_post_impute2))
dim(mC_rate_post_impute) # 12583 545

## load gene body CpH Here
GB_mCH_rate <- readRDS("20210915_552nuclei_genebody_CH_comb_rate.rds")
GB_mCH_nsite <- readRDS("20210915_552nuclei_genebody_CH_comb_nsite.rds")

mCH_rate_post_impute <- filter_impute3(GB_mCH_nsite, 
                                       GB_mCH_rate, 
                                       nSite =2, 
                                       rownames(tmC_rate_post_impute2))

colnames(mCH_rate_post_impute) <- gsub("CpH","CpG",colnames(mCH_rate_post_impute))

hmC_matrix <- hmC_rate_post_impute
mC_matrix <- mC_rate_post_impute
tmC_matrix <- tmC_rate_post_impute2
mCH_matrix <- mCH_rate_post_impute
dim(hmC_matrix)
dim(mC_matrix)
dim(tmC_matrix)
dim(mCH_matrix)



saveRDS(hmC_matrix,file="20230622_hmC_matrix.rds")
saveRDS(mC_matrix,file="20230622_mC_matrix.rds")
saveRDS(tmC_matrix,file="20230622_tmC_matrix.rds")
saveRDS(mCH_matrix,file="20230622_mCH_matrix.rds")

replace_cell_name_vector <- dff2$cell_snmCseq2
names(replace_cell_name_vector) <- dff2$cellID
replace_cell_name_vector[1:3]

colnames(hmC_matrix) <- replace_cell_name_vector[colnames(hmC_matrix)]
colnames(mC_matrix) <- replace_cell_name_vector[colnames(mC_matrix)]
colnames(tmC_matrix) <- replace_cell_name_vector[colnames(tmC_matrix)]

# load the previous pbmc.snmc item
pbmc.snmc <- readRDS(file="20210917_pbmc.snmc.rds")

pbmc.snmc

DefaultAssay(pbmc.snmc) <- "gbch"

pbmc.snmc <- ScaleData(pbmc.snmc)

options(repr.plot.width = 5, repr.plot.height = 4)
DimPlot(pbmc.snmc,pt.size = .5,group.by = "res.name.new",label = TRUE)

Idents(pbmc.snmc) <- "res.name.new"

pbmc.snmc2 <- subset(pbmc.snmc,idents = "In",invert=T)

options(repr.plot.width = 5, repr.plot.height = 4)
DimPlot(pbmc.snmc2,pt.size = .5,group.by = "res.name.new",label = TRUE)

pbmc.snmc2@meta.data <- droplevels(pbmc.snmc2@meta.data)

table(pbmc.snmc2@meta.data$res.name.new)

colnames(hmC_matrix)[1:3]

colnames(hmC_matrix) <- gsub("CpG","CpH",colnames(hmC_matrix))
colnames(mC_matrix) <- gsub("CpG","CpH",colnames(mC_matrix))
colnames(tmC_matrix) <- gsub("CpG","CpH",colnames(tmC_matrix))
colnames(mCH_matrix) <- gsub("CpG","CpH",colnames(mCH_matrix))

head(pbmc.snmc2@meta.data,n=2)

pbmc.snmc3 <- subset(pbmc.snmc2,cells=intersect(colnames(hmC_matrix),rownames(pbmc.snmc2@meta.data)))

hmC_matrix <- hmC_matrix[,rownames(pbmc.snmc3@meta.data)]
mC_matrix <- mC_matrix[,rownames(pbmc.snmc3@meta.data)]
tmC_matrix <- tmC_matrix[,rownames(pbmc.snmc3@meta.data)]
mCH_matrix <- mCH_matrix[,rownames(pbmc.snmc3@meta.data)]

dim(hmC_matrix)
dim(mC_matrix)
dim(tmC_matrix)
dim(mCH_matrix)

# add the 4 assays
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

gene.plot <- c("Arpp21","Rorb","Tle4","Cux2","Slc6a1","Adarb2","Cacna2d2")

pbmc.snmc3

library(ggplot2)
library(RColorBrewer)

DefaultAssay(pbmc.snmc3) <- "BSmCG"
#pbmc.snmc3
options(repr.plot.width = 28,repr.plot.height = 4)
FeaturePlot(pbmc.snmc3,features = gene.plot,
#            min.cutoff = -2,
#            max.cutoff = 2,
            raster = TRUE,ncol = 7,pt.size = 4,slot = "scale.data")& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

DefaultAssay(pbmc.snmc3) <- "hmCG"
#pbmc.snmc3
options(repr.plot.width = 28,repr.plot.height = 4)
FeaturePlot(pbmc.snmc3,features = gene.plot,
#            min.cutoff = -2,
#            max.cutoff = 2,
            raster = TRUE,ncol = 7,pt.size = 4,slot = "scale.data")& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

DefaultAssay(pbmc.snmc3) <- "tmCG"
#pbmc.snmc3
options(repr.plot.width = 28,repr.plot.height = 4)
FeaturePlot(pbmc.snmc3,features = gene.plot,
#            min.cutoff = -2,
#            max.cutoff = 2,
            raster = TRUE,ncol = 7,pt.size = 4,slot = "scale.data")& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))




