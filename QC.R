

# chargement des librairies
library(ggplot2)
library(ggpubr)
library(Seurat)<
library(scater)
library(dplyr)

library(babelgene)
setwd("K:/BioInfo/SC_Analyses/Florian_Sophie")

lst_rds = list.files("matrix")
lst_rds

# Basic function to convert human to mouse gene names


m.s.genes <- orthologs(genes = cc.genes.updated.2019$s.genes, species = "mouse")
m.s.genes <- m.s.genes$symbol

m.g2m.genes <- orthologs(genes = cc.genes.updated.2019$g2m.genes, species = "mouse")
m.g2m.genes <- m.g2m.genes$symbol
m.g2m.genes


for (l in lst_rds[2:6]){
l
    rawData <- Read10X(paste("matrix/", l, "/sample_filtered_feature_bc_matrix/", sep=""))
    seuratObjet <- CreateSeuratObject(counts = rawData$`Gene Expression`, project = l, min.cells = 3,assay = "RNA")
    seuratObjet <- NormalizeData(seuratObjet, assay = "RNA")

   
    seuratObjet <- CellCycleScoring(seuratObjet, s.features = m.s.genes, g2m.features =  m.g2m.genes, set.ident = TRUE)
    seuratObjet <- FindVariableFeatures(seuratObjet, nfeatures = 2000)
    seuratObjet <- ScaleData(seuratObjet, vars.to.regress = c("S.Score", "G2M.Score"))
    seuratObjet <- RunPCA(seuratObjet)
    seuratObjet <- RunUMAP(seuratObjet, dims = 1:30)
    seuratObjet <- FindNeighbors(seuratObjet, dims = 1:30)
    seuratObjet <- FindClusters(seuratObjet, resolution = c(0,seq(0.1, 1, by=0.1)))

    OUtliernCount <-isOutlier(seuratObjet$nCount_RNA, nmads=3, type="both", log=FALSE)
    seuratObjet$QC_ncount <- seuratObjet$nCount_RNA < attributes(OUtliernCount)$thresholds[2]
    p1 <- (ggplot(as.data.frame(seuratObjet$nCount_RNA), aes(x=seuratObjet$nCount_RNA)) + 
        geom_density(color="darkblue", fill="lightblue") +   labs(x="nCount_RNA", y = "Density") + 
        geom_vline(xintercept = c(attributes(OUtliernCount)$thresholds[2],500))  + 
        ggtitle(gsub("E8-","",l)) )  + DimPlot(seuratObjet, group.by = "QC_ncount")


    OUtliernFeature <-isOutlier(seuratObjet$nFeature_RNA, nmads=3, type="both", log=FALSE)
    seuratObjet$QC_nfeature <- seuratObjet$nFeature_RNA < attributes(OUtliernFeature)$thresholds[2] & seuratObjet$nFeature_RNA > attributes(OUtliernFeature)$thresholds[1]
    p2 <- (ggplot(as.data.frame(seuratObjet$nFeature_RNA), aes(x=seuratObjet$nFeature_RNA)) + 
            geom_density(color="plum4", fill="thistle1") +   
            labs(x="nFeature_RNA", y = "Density") + 
            geom_vline(xintercept =c(attributes(OUtliernFeature)$thresholds[2],max(c(attributes(OUtliernFeature)$thresholds[1],200))))) + 
            DimPlot(seuratObjet, group.by = "QC_nfeature")

    seuratObjet$log10GenesPerUMI <- log10(seuratObjet$nFeature_RNA) / log10(seuratObjet$nCount_RNA)
    OUtlierlog10GenesPerUMI <-isOutlier(seuratObjet$log10GenesPerUMI, nmads=2, type="both", log=FALSE)
    seuratObjet$QC_log10GenesPerUMI <- seuratObjet$log10GenesPerUMI > attributes(OUtlierlog10GenesPerUMI)$thresholds[1]
    p3 <- (ggplot(as.data.frame(seuratObjet$log10GenesPerUMI), aes(x = seuratObjet$log10GenesPerUMI))+
        geom_density(color="darkgreen", fill="palegreen") +   
        labs(x="log10GenesPerUMI", y = "Density") + 
        geom_vline(xintercept = attributes(OUtlierlog10GenesPerUMI)$thresholds[1])) + 
        DimPlot(seuratObjet, group.by = "QC_log10GenesPerUMI")

    seuratObjet <- PercentageFeatureSet(seuratObjet, pattern = "^mt-", col.name = "percent.mt")
    OUtlierPerMT <- isOutlier(seuratObjet$percent.mt, nmads=3, type="both", log=FALSE)
    seuratObjet$QC_percent.mt <- seuratObjet$percent.mt < max(attributes(OUtlierPerMT)$thresholds[2],7)
    p4 <- (ggplot(as.data.frame(seuratObjet$percent.mt), aes(x=seuratObjet$percent.mt)) + 
        geom_density(color="grey27", fill="grey50") +   labs(x="percent.mt", y = "Density") + 
        geom_vline(xintercept = max(attributes(OUtlierPerMT)$thresholds[2],7))) + 
        DimPlot(seuratObjet, group.by = "QC_percent.mt") + 
        FeaturePlot(seuratObjet, features = "percent.mt", max.cutoff = 25)
    
    ### doublets detection
    # scDblFinder
    dir_n = paste("QC_", l, sep = "")
    if(!dir.exists(dir_n)){
        dir.create(dir_n)
    }


    sce <- as.SingleCellExperiment(DietSeurat(seuratObjet))
    library(scDblFinder)
    sce <- scDblFinder(sce)
    seuratObjet$scDblFinder_score <- sce$scDblFinder.score
    seuratObjet$scDblFinder_class <- sce$scDblFinder.class
    d <- scDblFinder(sce, verbose=FALSE, returnType="table")
    th <- doubletThresholding(d)

    ##### QC plots #####
    ### doublet detection
    # scDblFinder
    OutlierscDblFinder <- th
    seuratObjet$QC_scDblFinder <- seuratObjet$scDblFinder_score < th
    p5 <- (ggplot(as.data.frame(seuratObjet$scDblFinder_score), aes(x = seuratObjet$scDblFinder_score))+geom_density(color="#92A9BD", fill="#92A9BD") +   labs(x="scDblFinder score", y = "Density") + geom_vline(xintercept = th))+  DimPlot(seuratObjet, group.by = "QC_scDblFinder")


    #png(paste("/SCRATCH-BIRD/users/jrondineau/QC-",gsub("E8-","",args[i]),".png", sep=""))
    pdf(paste(dir_n,"/",l,".pdf", sep=""), width=12, height = 6)
    print(p1)
    print(p2)
    print(p3)
    print(p4)
    print(p5)

    Idents(seuratObjet) <- "RNA_snn_res.0.1"
    markers <- FindAllMarkers(seuratObjet, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
    print((DimPlot(seuratObjet, group.by = "RNA_snn_res.0.1") + labs(title="",subtitle="Clustering res.0.1") + NoLegend()) + DoHeatmap(seuratObjet, features = top10$gene))

    seuratObjet$to_discard <- seuratObjet$QC_ncount & seuratObjet$QC_nfeature & seuratObjet$QC_log10GenesPerUMI & seuratObjet$QC_percent.mt & seuratObjet$QC_scDblFinder

    print(DimPlot(seuratObjet, group.by = "to_discard") + ggtitle(paste("Percentage of cells failing QC : ",round(table(seuratObjet$to_discard)[1]/sum(table(seuratObjet$to_discard))*100,2),"%", sep="")))

    ##### filter cells
    seuratObjet <- subset(seuratObjet, subset = to_discard)
    ##### reprocess with good cells

    seuratObjet <- FindVariableFeatures(seuratObjet, nfeatures = 2000)
    seuratObjet <- ScaleData(seuratObjet, vars.to.regress = c("S.Score", "G2M.Score"))
    seuratObjet <- RunPCA(seuratObjet)
    seuratObjet <- RunUMAP(seuratObjet, dims = 1:30)
    seuratObjet <- FindNeighbors(seuratObjet, dims = 1:30)
    seuratObjet <- FindClusters(seuratObjet, resolution = c(0,seq(0.1, 1, by=0.1)))

    Idents(seuratObjet) <- "RNA_snn_res.0.1"
    markers <- FindAllMarkers(seuratObjet, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
    print((DimPlot(seuratObjet, group.by = "RNA_snn_res.0.1") + labs(title="",subtitle=paste("Clustering after QC - res.0.1 - n=",sum(table(seuratObjet$to_discard)),sep="")) + NoLegend()) + DoHeatmap(seuratObjet, features = top10$gene))
    print(FeaturePlot(seuratObjet, features = "Luc2") | VlnPlot(seuratObjet, features = 'Luc2') + geom_boxplot(width = 0.3, fill="white") )

    dev.off()
    saveRDS(seuratObjet, paste("RDS/",l,".rds", sep=""))
    
}
