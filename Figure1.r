#!/usr/bin/Rscript
# or open it with jupyter notebook
# Title: "NATURE BIOTECHNOLOGY README"
# Author: "Peng"

# load previous seurat object
# load clustering results
library("Seurat")
library(dplyr)
library(tidyr)

s.tmp <- readRDS(file="1228_5_CpGs_70cell_3000HVG_50pc_res1_s.tmp.with.normalization.clustering.rds")

s.tmp

chr.in <- sapply(rownames(s.tmp),function(x){
    as.character(unlist(strsplit(x,"-"))[1])
})

chr.in.uniq <- unique(chr.in)
chr.in.uniq

CG_final <- readRDS(file="0111_ourData_1mb_bin_5hmCG.rds")

# using 5 CpG here
mat_CpG <- CG_final$rate
mat_CpG[CG_final$nsite < 5] <- NA

# replcate CpH- to CpG to match the orignal data set
colnames(mat_CpG) <- gsub(".gz","",colnames(mat_CpG))
#colnames(mat_CpH) <- gsub("CpG_NeuN","NeuN",colnames(mat_CpG))

dim(mat_CpG)

na_num_per_bin <- apply(mat_CpG,1,function(x){
    sum(is.na(x))
})

keep.features <- names(na_num_per_bin)[na_num_per_bin < ncol(mat_CpG) * (1-0.6)] # less than 60%
length(keep.features) # 2610

# impute by rowMeans
mat_CpG2 <- mat_CpG[keep.features,]

m <- mat_CpG2
## from https://stackoverflow.com/questions/6918086/replace-na-values-by-row-means
k <- which(is.na(m), arr.ind=TRUE)
m[k] <- rowMeans(m, na.rm=TRUE)[k[,1]]
mat_CpG3 <- m

dim(mat_CpG3)

mat_CpG4 <- mat_CpG3[,rownames(s.tmp@meta.data)]

s.tmp[["CpG"]] <- CreateAssayObject(counts = mat_CpG3[,rownames(s.tmp@meta.data)])
DefaultAssay(s.tmp) <- "CpG"
s.tmp <- NormalizeData(s.tmp,assay = "CpG")
s.tmp <- FindVariableFeatures(s.tmp,assay = "CpG")
s.tmp <- ScaleData(s.tmp, assay = "CpG",features = rownames(s.tmp))

# add annotation here
s.tmp@meta.data$Cluster <- plyr::mapvalues(s.tmp@meta.data$`RNA_snn_res.1`,c(0,1,2),c("0_Ex","1_nonNeu","2_Inh"))

s.tmp@meta.data$Cluster2 <- plyr::mapvalues(s.tmp@meta.data$Cluster,
                                            c("0_Ex","1_nonNeu","2_Inh"),
                                            c("Neu","nonNeu","Neu"))

plot_list <- lapply(chr.in.uniq[1:19],function(i){
    j <- paste0(i,"-")
    features.in <- grep(j,rownames(s.tmp),value=T)
    # order features by first
    coord.start <- sapply(features.in,function(x){
        as.integer(unlist(strsplit(x,"-"))[2])
    })
    features.in2 <- names(sort(coord.start))
    
    mat <- as.matrix(s.tmp@assays$CpG@counts[features.in2,])
    
    if (nrow(mat) > 1){
    
    mat[mat > .4] <- .4
    
    NeuID <- rownames(s.tmp@meta.data)[s.tmp@meta.data$Cluster2 == "Neu"]
    nonNeuID <- rownames(s.tmp@meta.data)[s.tmp@meta.data$Cluster2 == "nonNeu"]
    
    # cellid as row, 1kb bin as column
    mat <- t(mat)
        
    mat.p1 <- mat[NeuID,]
    mat.p1 <- mat.p1[order(rowSums(mat.p1),decreasing = T),]
    mat.p2 <- mat[nonNeuID,]
    mat.p2 <- mat.p2[order(rowSums(mat.p2),decreasing = T),]    
        
    mat2 <- rbind(mat.p1,mat.p2)
    
    hf <- pheatmap::pheatmap(mat2,scale = "none",cluster_rows = F,cluster_cols = F,
                   show_colnames = F, show_rownames = F,
                       gaps_row = length(NeuID),
                   main = i)
    return(hf[[4]])
        
    }
})
