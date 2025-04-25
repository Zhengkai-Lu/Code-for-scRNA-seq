library(limma)
library(dplyr)
library(magrittr)
library(Seurat)
library(cowplot)
library(patchwork)
library(ggplot2)
library(metap)
library(ggpubr)

### setting working directory
setwd("F:\\E-MTAB-7376")  

### Set up Sham object
Sham <- Read10X("F:\\E-MTAB-7376\\Sham_day7_TIP")

Sham <- CreateSeuratObject(counts = Sham, project = "Sham", min.cells = 0, min.features = 0, names.delim = "_")
Sham$stim <- "Sham"

Sham[["percent.mt"]] <- PercentageFeatureSet(object = Sham, pattern = "^mt-")
VlnPlot(object = Sham, pt.size = 0.05, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

Sham <- subset(Sham, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA <20000 & percent.mt < 20) 
Sham <- FindVariableFeatures(Sham, selection.method = "vst", nfeatures = 3000)

### Set up MID3 object
MID3 <- Read10X("F:\\E-MTAB-7376\\MI_day3_TIP")

MID3 <- CreateSeuratObject(counts = MID3, project = "MI D3", min.cells = 0, min.features = 0, names.delim = "_")
MID3$stim <- "MI D3"

MID3[["percent.mt"]] <- PercentageFeatureSet(object = MID3, pattern = "^mt-")
VlnPlot(object = MID3, pt.size = 0.05, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

MID3 <- subset(MID3, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA <20000 & percent.mt < 20)   
MID3 <- NormalizeData(MID3)
MID3 <- FindVariableFeatures(MID3, selection.method = "vst", nfeatures = 3000)

### Set up MID7 object
MID7 <- Read10X("F:\\20230411;E-MTAB-7376\\MI_day7_TIP")

MID7 <- CreateSeuratObject(counts = MID7, project = "MI D7", min.cells = 0, min.features = 0, names.delim = "_")
MID7$stim <- "MI D7"

MID7[["percent.mt"]] <- PercentageFeatureSet(object = MID7, pattern = "^mt-")
VlnPlot(object = MID7, pt.size = 0.05, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

MID7 <- subset(MID7, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA <20000 & percent.mt < 20)    
MID7 <- NormalizeData(MID7)
MID7 <- FindVariableFeatures(MID7, selection.method = "vst", nfeatures = 3000)

### identify anchors and Perform integration
features <- SelectIntegrationFeatures(object.list = list(Sham, MID3, MID7), nfeatures = 3000)
anchors <- FindIntegrationAnchors(object.list = list(Sham, MID3, MID7), anchor.features = features, dims = 1:30)
combined <- IntegrateData(anchorset = anchors, dims = 1:30)

### Perform an integrated analysis
DefaultAssay(combined) <- "integrated"

### Run the standard workflow for visualization and clustering
all.genes <- rownames(combined)
combined <- ScaleData(combined, features = all.genes)   
combined <- RunPCA(combined, npcs = 50)

Idents(combined) ="stim"
DimPlot(object = combined, reduction = "pca")
DimHeatmap(object = combined, dims = 1:50, cells = 500, balanced = TRUE, nfeatures = 60, ncol=2) 
ElbowPlot(combined, ndims = 50, reduction = "pca")

combined <- FindNeighbors(combined, reduction = "pca", dims = 1:36)
combined <- FindClusters(combined, resolution = 0.3)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:36)
combined <- RunTSNE(combined, reduction = "pca", dims = 1:36)

color=DiscretePalette(25, palette = "alphabet2")
DimPlot(combined, pt.size = 0.5, reduction = "tsne", label = TRUE, repel = T, cols = color) + NoAxes()

combined$stim <- factor(combined$stim, levels = c("Sham","MI D3","MI D7"))

### Identify conserved cell type markers
DefaultAssay(combined) <- "RNA"

lapply(levels(combined), FUN = function(x) {
  M <- FindConservedMarkers(combined, ident.1 = x, grouping.var = "stim")
  path <- paste("cluster", x , sep="")
  write.table(M, file = paste(path,"xls", sep="."),sep="\t",row.names=T,quote=F)
})

combined <- RenameIdents(combined, `0` = "Endothelial cells", `1` = "Fibroblasts", `2` = "Macrophages",
                         `3` = "Granulocytes", `4` = "B cells", `5` = "Monocytes", 
                         `6` = "Endothelial cells", `7` = "T cells", `8` = "Endothelial cells", 
                         `9` = "Mural cells",`10` = "Dendritic-like cells", `11` = "Fibroblasts", `12` = "Endothelial cells",
                         `13` = "Fibroblasts", `14` = "Natural killer cells", `15` = "B cells", 
                         `16` = "Macrophages", `17` = "Cardiomyocytes",`18` = "Glial cells")

combined$celltype <- Idents(combined)
Idents(combined) <- "celltype"

color=DiscretePalette(25, palette = "alphabet2")
DimPlot(combined, pt.size = 0.5, reduction = "tsne", label = F, repel = T, cols = color)+ NoAxes()
DimPlot(combined, reduction = "tsne", pt.size = 0.5, repel = T, label = F, split.by = "stim", cols = color) & NoAxes()
DotPlot(object = combined, cols = c("lightgray", "red"),
        features = c("Emcn", "Cdh5", "Kdr","Pecam1", "Egfl7", 
                     "Col1a1", "Col5a1","Mmp2", "Fstl1", "Dcn",
                     "C1qa", "C1qb", "Cd68", "Lgals3", "Csf1r",
                     "S100a8", "S100a9", "Il1b","Clec4d", "Csf3r",
                     "Cd79a", "Cd79b", "Ly6d", "Ms4a1","Bank1",
                     "Ms4a6c", "Plac8","Ccr2","Clec4a3","Ifitm6",
                     "Cd3d", "Cd3g", "Cd3e", "Lef1","Thy1",
                     "Rgs5", "Pdgfrb", "Notch3", "Abcc9", "Cspg4",
                     "Plbd1", "H2-DMb1", "Ctnnd2", "Ciita", "Flt3", 
                     "Gzma", "Ccl5", "Nkg7", "Klrb1c", "Klrk1",
                     "Myl2", "Actc1", "Tnnt2", "Tnni3", "Tnnc1",
                     "Csmd1", "Prnp", "Plp1", "Gfra3", "Kcna1"
        ), col.min = 0, dot.scale = 10
)  + RotatedAxis()
dev.off() 

saveRDS(combined,file = "celltype.rds")


  
### cell number of Mki67+ ECs and mural cells
MC = subset(combined,idents = "Mural cells")
EC = subset(combined,idents = "Endothelial cells")

pro1=subset(MC,Mki67>0, slot="counts")
table(pro1$orig.ident) 
table(MC$orig.ident)

pro2=subset(EC,Mki67>0, slot="counts")
table(pro2$orig.ident) 
table(EC$orig.ident)

### expression levels of proliferation-associated genes in ECs and mural cells 
VlnPlot(combined, pt.size = 0.1, features = c("Cdk1","Top2a","Mki67","Birc5","Ube2c",
                                              "Prc1","Cdc20","Ccna2","Ccnb2","Ccnb1"),
        idents = c.factor("Endothelial cells","Mural cells"),  split.by = "stim", combine = T, ncol = 4) + RestoreLegend() 

### Expression of Cxcr4 and Cxcr7 in ECs
FeaturePlot(EC, reduction = "tsne", features = "Cxcr4", split.by = "stim",min.cutoff = 0, max.cutoff = 3,
            pt.size = 0.5, repel = T, label = F, order = T,
            cols = c("lightgrey","red")) + RestoreLegend() & NoAxes()
FeaturePlot(EC, reduction = "tsne", features = "Ackr3", split.by = "stim",min.cutoff = 0, max.cutoff = 3,
            pt.size = 0.5, repel = T, label = F, order = T,
            cols = c("lightgrey","red")) + RestoreLegend() & NoAxes()



### Expression of Wt1 and Krt19
FeaturePlot(combined, reduction = "tsne", features = c("Wt1","Krt19"), min.cutoff = 0, max.cutoff = 2,
            pt.size = 0.3, repel = T, label = F,ncol = 2,
            cols = c("lightgrey","red")) + RestoreLegend() & NoAxes()


### Expression of smooth muscle cell makers in mural cells
DotPlot(object = combined, cols = c("lightgray", "red"),
        features = c("Myl9", "Mylk", "Sncg","Acta2","Tagln", "Mustn1","Myh11"), col.min = 0, dot.scale = 12)+ 
  RotatedAxis()+ xlab(NULL)+ ylab(NULL)


### the expression of Pten in mural cells
VlnPlot(combined, pt.size = 0, features = "Pten", 
        idents = "Mural cells", ncol = 3,
        group.by = "stim", combine =  T ) & geom_boxplot(width=0.1,col="black",fill="white") &
  stat_compare_means(comparison= compared,label="p.signif", method = "wilcox.test") & ylim(-0.5,6)




sessionInfo()

# R version 4.2.2 (2022-10-31 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 7 x64 (build 7601) Service Pack 1

# Matrix products: default

# locale:
#   [1] LC_COLLATE=Chinese (Simplified)_People's Republic of China.936 
# [2] LC_CTYPE=Chinese (Simplified)_People's Republic of China.936   
# [3] LC_MONETARY=Chinese (Simplified)_People's Republic of China.936
# [4] LC_NUMERIC=C                                                   
# [5] LC_TIME=Chinese (Simplified)_People's Republic of China.936    

# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
#   [1] ggpubr_0.6.0       ggplot2_3.4.3      metap_1.8          patchwork_1.1.1    cowplot_1.1.1     
# [6] sp_1.5-0           SeuratObject_4.1.0 Seurat_4.1.1       magrittr_2.0.3     dplyr_1.1.3       
# [11] limma_3.50.3      

# loaded via a namespace (and not attached):
#   [1] TH.data_1.1-1         Rtsne_0.16            colorspace_2.0-3      ggsignif_0.6.3       
# [5] deldir_1.0-6          ellipsis_0.3.2        ggridges_0.5.4        rstudioapi_0.13      
# [9] spatstat.data_2.2-0   leiden_0.4.2          listenv_0.8.0         ggrepel_0.9.3        
# [13] fansi_1.0.3           mvtnorm_1.1-3         mathjaxr_1.6-0        codetools_0.2-18     
# [17] splines_4.2.2         mnormt_2.1.0          TFisher_0.2.0         polyclip_1.10-0      
# [21] jsonlite_1.8.0        broom_1.0.0           ica_1.0-3             cluster_2.1.4        
# [25] png_0.1-7             rgeos_0.5-9           uwot_0.1.11           shiny_1.7.1          
# [29] sctransform_0.3.3     spatstat.sparse_2.1-1 BiocManager_1.30.20   compiler_4.2.2       
# [33] httr_1.4.3            backports_1.4.1       Matrix_1.5-1          fastmap_1.1.0        
# [37] lazyeval_0.2.2        cli_3.6.1             later_1.3.0           htmltools_0.5.2      
# [41] tools_4.2.2           igraph_1.3.2          gtable_0.3.0          glue_1.6.2           
# [45] RANN_2.6.1            reshape2_1.4.4        Rcpp_1.0.10           carData_3.0-5        
# [49] Biobase_2.54.0        scattermore_0.8       vctrs_0.6.3           multtest_2.54.0      
# [53] nlme_3.1-160          progressr_0.10.1      lmtest_0.9-40         spatstat.random_2.2-0
# [57] stringr_1.5.0         globals_0.15.1        rbibutils_2.2.8       mime_0.12            
# [61] miniUI_0.1.1.1        lifecycle_1.0.3       irlba_2.3.5           rstatix_0.7.2        
# [65] goftest_1.2-3         future_1.26.1         qqconf_1.2.3          MASS_7.3-58.1        
# [69] zoo_1.8-10            scales_1.2.0          spatstat.core_2.4-4   promises_1.2.0.1     
# [73] spatstat.utils_2.3-1  sandwich_3.0-2        parallel_4.2.2        RColorBrewer_1.1-3   
# [77] reticulate_1.25       pbapply_1.5-0         gridExtra_2.3         rpart_4.1.19         
# [81] stringi_1.7.6         mutoss_0.1-12         plotrix_3.8-2         BiocGenerics_0.40.0  
# [85] Rdpack_2.3.1          rlang_1.1.1           pkgconfig_2.0.3       matrixStats_0.62.0   
# [89] lattice_0.20-45       ROCR_1.0-11           purrr_1.0.2           tensor_1.5           
# [93] htmlwidgets_1.5.4     tidyselect_1.2.0      parallelly_1.32.0     RcppAnnoy_0.0.20     
# [97] plyr_1.8.7            R6_2.5.1              generics_0.1.3        multcomp_1.4-19      
# [101] withr_2.5.0           pillar_1.9.0          mgcv_1.8-41           sn_2.0.2             
# [105] fitdistrplus_1.1-8    survival_3.4-0        abind_1.4-5           tibble_3.2.1         
# [109] future.apply_1.9.0    car_3.1-2             KernSmooth_2.23-20    utf8_1.2.2           
# [113] spatstat.geom_2.4-0   plotly_4.10.0         grid_4.2.2            data.table_1.14.2    
# [117] digest_0.6.29         xtable_1.8-4          numDeriv_2016.8-1.1   tidyr_1.3.0          
# [121] httpuv_1.6.5          stats4_4.2.2          munsell_0.5.0         viridisLite_0.4.2    

