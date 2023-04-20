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
library(SeuratData)
library(Get.CibersortxTpm)
data("pbmc3k")
pbmc3k <- pbmc3k[,1:100]
```

### Visualization
##### *The package provide visualizable function for Spatial Transcriptomics data*
```
HexSpatialPlot(object = colon.sct,group.by ='celltype',align ='v',
               color = c('grey','blue','green',col,'skyblue','indianred1','khaki1')
```
<img width="852" alt="b6178ecf05e799832c7e10923153269" src="https://user-images.githubusercontent.com/122006615/233295037-a748a6d2-826c-4bf2-8a45-c2040a9783ac.png">

##### *Or you can select specific cells*

```

```
### Detect spot

### Estimate the immune cell infiltration score
