# SpotSweeper
#### Estimate the immune cell infiltration score of adjacent stromal spots in specific groups of tumor spots
![image](https://user-images.githubusercontent.com/122006615/233296103-12f4cde2-51f3-4b52-9826-357f8d6bd913.png)



## Install
```
devtools::install_github('')
```  
## Usage

### Load example data
```

```

### Visualization
#### The package provide visualizable function for Spatial Transcriptomics data
##### *For metadata*
```
HexSpatialPlot(object = colon.sct,group.by ='celltype',align ='v',
               color = c('grey','blue','green',col,'skyblue','indianred1','khaki1')
```
<img width="852" alt="b6178ecf05e799832c7e10923153269" src="https://user-images.githubusercontent.com/122006615/233295037-a748a6d2-826c-4bf2-8a45-c2040a9783ac.png">

##### *Or gene expression*

```
HexSpatialPlot(colon.sct,group.by = 'MUC2',legend = 'top',align = 'hv',ncol = 2)
```
![MUC2](https://user-images.githubusercontent.com/122006615/233325374-3d4ae00c-d3ac-4cdc-97f1-1f312dd5cf81.png)
### Detect spot
```

```
### Estimate the immune cell infiltration score
```
colon1 <- RunssGSEA(object = colon1,genelist = immune_list)

```

### DEboxplot
```
DEboxplot(colon)

```

### Visualization for ssGSEA score
```
HexssGSEAplot(colon,type = 'C1QC_TAM')

```
### Calculate enrichment fold
```
HexssGSEAplot(colon,type = 'C1QC_TAM')

```

## Citations
- Gao S, Shi Q, Zhang Y, et al. Identification of HSC/MPP expansion units in fetal liver by single-cell spatiotemporal transcriptomics. Cell Res. 2022;32(1):38-53. doi:10.1038/s41422-021-00540-7
