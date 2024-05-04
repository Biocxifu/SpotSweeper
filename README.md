![Static Badge](https://img.shields.io/badge/build-Seurat%20v4-brightgreen?label=Base)

# SpotSweeper<img src="https://user-images.githubusercontent.com/122006615/236639340-d832f708-5478-499a-9bfb-7bcb0dd7c89f.png" align="right" alt="" width="200" />
#### 
- Utilize fundamental visualization functions and user-friendly interactive tools for manual annotation of spots.
- Identify spots proximate to tumors and categorize them into specific groups.
- Estimate the score of immune cell infiltration in spots and compare the score difference between specific tumor spots and adjacent spots.
- Enable the computation of transcriptional heterogeneity within samples.
- Facilitate communication analysis between specific tumor spots and adjacent spots.
- Execute spatially resolved niche classification and calculate the infiltration of characteristic cells in the micro-ecotype.

## Install
```
devtools::install_github('Biocxifu/SpotSweeper')
```  
## Usage
The majority of functions within this software package are designed to accommodate multiple samples. However, for the purpose of clarity and convenience, we demonstrate the functionality using a single sample in this instance.

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
The categorization presented herein is informal and is solely employed for illustrative purposes.
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
We have developed an R Shiny application to facilitate precise spot annotation. R Shiny enables the mapping of gene expression levels and meta-information to aid in the manual annotation process. Spot annotations are presented in the first line, and the selected information for display is showcased in the second line. The label color corresponds to the color bar.
```
tissue <- SpotAnnotation(object = tissue,celltype_var = 'celltype')
```
<img width="849" alt="shiny" src="https://github.com/Biocxifu/SpotSweeper/assets/122006615/0431b498-8172-43a7-96b6-ea042dcb2e81">

We suggest configuring the transparency of points to 0 and selecting either shape 1 or 16 for points, followed by manual annotation.
In addition to referencing the expression values of specific marker genes, we advocate considering the malignancy score computed by spaCET as an additional reference for correction.

<img width="863" alt="shiny2" src="https://github.com/Biocxifu/SpotSweeper/assets/122006615/47d7e858-0224-43c1-a5ec-261529253c2b">

### Detect spot
Following the annotation of the spots, subsequent steps involve tumor grouping and the identification of adjacent spots.

This procedural sequence necessitates the initial establishment of the tumor group. When our research relies on gene expression for grouping, the most straightforward approach is to directly group all spots, without the need to calibrate the group information of non-tumor spots, as such information will be disregarded during the functioning of the process.
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
This assertion holds true not only for gene sets related to immune cells but extends to other gene sets as well, contingent on the objectives of the research.
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
heterogeneity <- Heterogeneity(object = tissue,TumorName = 'Tumor',celltype_var = 'celltype')
heterogeneity <-as.data.frame(heterogeneity)
```

### Cellular communication network
```
SpotChat(object = tissue,subSample = 'tissue1',cores = 8,
         json.path = '../data/DataSource/ST-tissue1/outs/spatial/scalefactors_json.json') 
```
### Spatially resolved niche classification (Kmeans method)
##### *Determine the number of subtypes*
```
Cluster_p <- NicheKmeans(object = tissue,ncol = 2)
Cluster_p
```
![image](https://github.com/Biocxifu/SpotSweeper/assets/122006615/44707b72-0351-4bbc-b340-cc51bdee2078)



##### *Niche subtypes*
```
Ecotype_P <- NicheCluster(object = tissue,centers = 6,
                         color = c("#87ceeb" ,'#ffa8ad' ,
                                   "#ffc37f",'#fc8d62',"#33b24a",'#ad71b5'))
Ecotype_P[1]
Ecotype_P[2]
Ecotype_P[3]
```
![image](https://github.com/Biocxifu/SpotSweeper/assets/122006615/e822b2ef-f665-43a4-85d2-3b74e5f15f4e)

##### *Niche plot*
```
NichePlot <- NicheCalculate(Ecotype_P = Ecotype_P)
NichePlot[1]
NichePlot[2]
```
![image](https://github.com/Biocxifu/SpotSweeper/assets/122006615/aa685b25-71de-4a9c-92fc-4bade17af4b1)



## Citations
- 	Gao S, Shi Q, Zhang Y, Liang G, Kang Z, Huang B, Ma D, Wang L, Jiao J, Fang X, et al. Identification of HSC/MPP expansion units in fetal liver by single-cell spatiotemporal transcriptomics. Cell Res (2022) 32:38–53. doi: 10.1038/s41422-021-00540-7
- Ru B, Huang J, Zhang Y, Aldape K, Jiang P. Estimation of cell lineages in tumors from spatial transcriptomics data. Nat Commun (2023) 14:568. doi: 10.1038/s41467-023-36062-6
- Wu R, Guo W, Qiu X, Wang S, Sui C, Lian Q, Wu J, Shan Y, Yang Z, Yang S, et al. Comprehensive analysis of spatial architecture in primary liver cancer. Sci Adv (2021) 7:eabg3750. doi: 10.1126/sciadv.abg3750
