# SpotSweeper<img src="https://user-images.githubusercontent.com/122006615/236639340-d832f708-5478-499a-9bfb-7bcb0dd7c89f.png" align="right" alt="" width="200" />
#### 
- The package provides basic visualization functions and convenient interactive tools for manually annotating spots
- It can select all nearby tumor spots and assign specific groups.
- It can estimate the score of immune cell infiltration in spots and compare the score difference of immune cell infiltration between specific tumor spots and adjacent spots.
- It can calculate the transcriptional heterogeneity of samples.
- It supports communication analysis between specific tumor spots and adjacent spots.

## Install
```
devtools::install_github('Biocxifu/SpotSweeper')
```  
## Usage
Most functions in this package support multiple samples. For convenience, we only use one sample for demonstration here.

data source：https://www.10xgenomics.com/resources/datasets/human-breast-cancer-ductal-carcinoma-in-situ-invasive-carcinoma-ffpe-1-standard-1-3-0

### Load example data
```
tissue <- Load10X_Spatial(data.dir = 'example/tissue1/',
                              filename = 'Visium_FFPE_Human_Breast_Cancer_filtered_feature_bc_matrix.h5')
tissue[["percent.mt"]] <- PercentageFeatureSet(tissue, pattern = "^[MT-]")
tissue <- SCTransform(tissue, assay = "Spatial", verbose = FALSE)

#filter
#VlnPlot(tissue, features = c("nCount_Spatial","nFeature_Spatial","percent.mt"), pt.size = 0.1) + NoLegend()
#tissue<- subset(tissue, subset = nFeature_Spatial > 200 & nFeature_Spatial <8000 &
#                nCount_Spatial > 1000 & nCount_Spatial < 40000 & percent.mt < 20)

##reduction
tissue <- RunPCA(tissue, assay = "SCT", verbose = FALSE)
tissue <- FindNeighbors(tissue, reduction = "pca", dims = 1:30)
tissue <- FindClusters(tissue, verbose = FALSE,resolution = 0.9)
tissue <- RunUMAP(tissue, reduction = "pca", dims = 1:30)
tissue <- RunTSNE(tissue, dims = 1:30)
```

### Visualization
#### The package provide visualizable function for Spatial Transcriptomics data
##### *For metadata*
```
cols <- c('#0083c4','#ffdd00','#6aa692','#c27874','#ffaabf','#ffa500','#97cf16',
          '#8ddcc3','#69008c','#d9d9d9','#ffbbff','#00bfff','#ff4040','#ffff00')
p1 <- HexSpatialPlot(object = tissue,group.by = 'seurat_clusters',color = col,
               legend = 'top',plot.image = T)
p2 <- HexSpatialPlot(object = tissue,group.by = 'seurat_clusters',color = col,
                     legend = 'top',plot.image =T,size = 1.4,alpha = 0.3)
p3 <- HexSpatialPlot(object = tissue,group.by = 'seurat_clusters',color = col,
                     legend = 'top',plot.image = F,size = 1.5)

p1+p2+p3
```
![image](https://github.com/Biocxifu/SpotSweeper/assets/122006615/ff52eba5-2505-4845-89de-0bb89764398e)

##### *Or gene expression*

```
p1 <- HexSpatialPlot(object = tissue,group.by = 'KRT18',
                     legend = 'top',plot.image = T,size = 1.2)
p2 <- HexSpatialPlot(object = tissue,group.by = 'KRT18',
                     legend = 'top',plot.image = F,size = 1.5)
p1+p2
```
![image](https://github.com/Biocxifu/SpotSweeper/assets/122006615/d7265313-4862-4535-a102-f34fbb01c9d0)

### Roughly annotate the spots
The grouping here is casual, just for demonstration purposes
```
celltype <- c("Tumor0"=c(5,6,9,11,12,13),
              "Fibroblasts0"=c(0,1,2,4,7),
              'Unknown0'=c(3,8,10)
) %>% as.data.frame() 
celltype %<>%mutate(celltype=str_split(rownames(celltype),'0',simplify = T)[,1]) %>% 
  rownames_to_column(var = 'dele') %>% select('celltype','.') %>% 
  column_to_rownames(var = '.') 
sort(unique(rownames(celltype)))
table(rownames(celltype))

tissue@meta.data$celltype='NA'
for(i in 1:nrow(tissue@meta.data)){
  index=as.character(tissue@meta.data$SCT_snn_res.0.9[i])
  tissue@meta.data$celltype[i]=celltype[rownames(celltype)==index,1]
}
```


### Modify spot annotation
We develop a Rshiny to modify spot annotation accurately. Rshiny can support the mapping of gene expression levels and meta information to assist in the manual annotation. Spot annotation is displayed in the first line and the information you select to show is in the second line.The label color is equivalent to color bar. 
```
colon <- SpotAnnotation(object = colon,celltype_var = 'celltype')
```
<img width="849" alt="shiny" src="https://github.com/Biocxifu/SpotSweeper/assets/122006615/0431b498-8172-43a7-96b6-ea042dcb2e81">

We recommend setting the transparency of points to 0 and the shape of points to 1 or 16, and then manually annotating them.
In addition to referencing the expression values of some marker genes, the malignancy score calculated using spaCET is also recommended as a reference for correction, which may even be important

<img width="863" alt="shiny2" src="https://github.com/Biocxifu/SpotSweeper/assets/122006615/47d7e858-0224-43c1-a5ec-261529253c2b">

### Detect spot
After annotating the spots, tumor grouping and detection of adjacent spots can be carried out

This process requires setting the tumor group first. If our research is based on gene expression for grouping, the simplest method is to directly group all spots without worrying about the group information of non tumor spots being calibrated, as the group information of non tumor spots will be ignored during function operation
```
tissue$group <- ifelse(tissue@assays$SCT@data['MKI67',]>0,'MKI67_pos','MKI67_neg')
tissue <- SpotSweeper(tissue,group_var = 'group',celltype_var = 'celltype',
                      tumor_name = c('Tumor'), nearby_name = c('Fibroblasts'))
```
Visualization for result
```
col2 <- c('#6d419c','#f8766d','#d7ebff','#cab3d6','#ffc0cb','#cccccc','grey69')
p1 <- HexSpatialPlot(object = tissue,group.by = 'Sweeper_type',color = col2,
                     legend = 'top',plot.image = F,size = 1.2)

p2 <- HexSpatialPlot(object = tissue,group.by = 'MKI67',
                     legend = 'top',plot.image = F,size = 1.2)
p1+p2
```
![image](https://github.com/Biocxifu/SpotSweeper/assets/122006615/83e8be1a-0f1a-4696-8e10-f85e36934500)


### Estimate the immune cell infiltration score
In fact, this applies not only to the gene set of immune cells, but also to other gene sets, depending on the research purpose
```
tissue <- CellScore(object = tissue,method = 'AddModuleScore',genelist = NULL)
tissue <- CellScore(object = tissue,method = 'ssGSEA',genelist = NULL,parallel.sz = 4)
```
### Visualization for immune cell infiltration score
```
colnames(tissue@assays[["ImmuneScore"]]$ssGSEA)
colnames(tissue@assays[["ImmuneScore"]]$AddModuleScore)
p1 <- HexScorePlot(object = tissue,ScoreType = 'ssGSEA',
                   type = 'Activated_CD8_T_cell',legend = 'top')
p2 <- HexScorePlot(object = tissue,ScoreType = 'AddModuleScore',
                   type = 'AddModuleScore_Activated_CD8_T_cell1',legend = 'top')
p1+p2
```
![image](https://github.com/Biocxifu/SpotSweeper/assets/122006615/e6973f54-dfd7-48a8-9a0e-f96aadaa0e61)


### Compare the scores of immune cell infiltration between groups and calculate enrichment fold
```
table(tissue$Sweeper_type)
ctr.group=list(c('MKI67_pos','Nearby_MKI67_pos'),c('MKI67_neg','Nearby_MKI67_neg'))
DeScore <- DeScore(object = tissue,ScoreType = 'ssGSEA',hide.ns = T,
               ctr.group = ctr.group)
DeScore$Fold
DeScore$plot
```
![image](https://github.com/Biocxifu/SpotSweeper/assets/122006615/59b0befa-9d27-4011-ab2d-20e485a27d2e)

### Evaluate transcriptome Heterogeneity of Samples

```
heterogeneity <- Heterogeneity(object = colon,TumorName = 'Tumor',celltype_var = 'celltype')
heterogeneity <-as.data.frame(heterogeneity)
```

### Cellular communication network
```
SpotChat(object = colon,subSample = 'colon1',cores = 8,
         json.path = '../data/DataSource/ST-colon1/outs/spatial/scalefactors_json.json') 
```




## Citations
- 	Gao S, Shi Q, Zhang Y, Liang G, Kang Z, Huang B, Ma D, Wang L, Jiao J, Fang X, et al. Identification of HSC/MPP expansion units in fetal liver by single-cell spatiotemporal transcriptomics. Cell Res (2022) 32:38–53. doi: 10.1038/s41422-021-00540-7
- Ru B, Huang J, Zhang Y, Aldape K, Jiang P. Estimation of cell lineages in tumors from spatial transcriptomics data. Nat Commun (2023) 14:568. doi: 10.1038/s41467-023-36062-6
- Wu R, Guo W, Qiu X, Wang S, Sui C, Lian Q, Wu J, Shan Y, Yang Z, Yang S, et al. Comprehensive analysis of spatial architecture in primary liver cancer. Sci Adv (2021) 7:eabg3750. doi: 10.1126/sciadv.abg3750
