# SpotSweeper<img src="https://user-images.githubusercontent.com/122006615/236639340-d832f708-5478-499a-9bfb-7bcb0dd7c89f.png" align="right" alt="" width="200" />
#### 
- The package provides basic visualization functions and convenient interactive tools for manually annotating spots
- It can select all nearby tumor spots and assign specific groups.
- It can estimate the score of spot immune cell infiltration and compare the differences in immune cell infiltration between specific tumor spots and adjacent spots.
- It supports communication analysis between specific tumor spots and adjacent spots

## Install
```
devtools::install_github('Biocxifu/SpotSweeper')
```  
## Usage

### Load example data
```

```
### Modify spot annotation
We develop a Rshiny to modify spot annotation accurately. Rshiny can support the mapping of gene expression levels and meta information to assist in the manual annotation. Spot annotation is displayed in the first line and the information you select to show is in the second line.The label color is equivalent to color bar. 
```
colon <- SpotAnnotation(object = colon,celltype_var = 'celltype')
```
<img width="1000" alt="shniy_example" src="https://user-images.githubusercontent.com/122006615/236654377-cd4c73a1-04da-49da-973b-d6ae4d1862fb.png">

We recommend setting the transparency of points to 0 and the shape of points to 1 or 16, and then manually annotating them.
![shiny_example2](https://user-images.githubusercontent.com/122006615/236639369-adc3824f-8432-4b80-9436-8bcf306b2106.gif)


### Visualization
#### The package provide basic visualizable function for Spatial Transcriptomics data
##### *For metadata*
```
HexSpatialPlot(object = colon.sct,group.by ='celltype',align ='v',
               color = c('grey','blue','green',col,'skyblue','indianred1','khaki1')
```
><img width="852" alt="b6178ecf05e799832c7e10923153269" src="https://user-images.githubusercontent.com/122006615/233295037-a748a6d2-826c-4bf2-8a45-c2040a9783ac.png">

##### *Or gene expression*

```
HexSpatialPlot(colon.sct,group.by = 'MUC2',legend = 'top',align = 'hv',ncol = 2)
```
>![MUC2](https://user-images.githubusercontent.com/122006615/233325374-3d4ae00c-d3ac-4cdc-97f1-1f312dd5cf81.png)
### Detect spot
```

```
### Visualization for result
```

```
>![image](https://user-images.githubusercontent.com/122006615/235993210-3a841544-c772-4191-b2b1-55e456680756.png)


### Estimate the immune cell infiltration score
```
colon1 <- RunssGSEA(object = colon1,genelist = immune_list)

```

### DEboxplot
```
DEboxplot(colon)

```
>![image](https://user-images.githubusercontent.com/122006615/233399369-7799a532-706a-4568-8cd6-ad1148519c84.png)

### Visualization for ssGSEA score
```
HexssGSEAplot(colon,type = 'C1QC_TAM')

```
![image](https://user-images.githubusercontent.com/122006615/233401762-76ee47aa-95de-4b4f-babf-c0da475f9752.png)


### Calculate enrichment fold
```
```
### Evaluate transcriptome Heterogeneity of Samples

```
heterogeneity <- Heterogeneity(object = colon,TumorName = 'Tumor',celltype_var = 'celltype')
heterogeneity <-as.data.frame(heterogeneity)
```
![image](https://github.com/Biocxifu/SpotSweeper/assets/122006615/b25edad2-bfe7-4db0-b4ee-b5879faae6b5)


### Cellular communication network

```
SpotChat(object = colon,subSample = 'colon1',cores = 8,
         json.path = '../data/DataSource/ST-colon1/outs/spatial/scalefactors_json.json') 
```
![image](https://github.com/Biocxifu/SpotSweeper/assets/122006615/2c8618e2-35eb-4ac1-8867-9575faebc999)



## Citations
- Gao S, Shi Q, Zhang Y, et al. Identification of HSC/MPP expansion units in fetal liver by single-cell spatiotemporal transcriptomics. Cell Res. 2022;32(1):38-53. doi:10.1038/s41422-021-00540-7
- https://github.com/data2intelligence/SpaCET
