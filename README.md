# SpotSweeper
#### Estimate the immune cell infiltration score of adjacent stromal spots in specific groups of tumor spots
![image](https://user-images.githubusercontent.com/122006615/233296103-12f4cde2-51f3-4b52-9826-357f8d6bd913.png)



## Install
```
devtools::install_github('Biocxifu/SpotSweeper')
```  
## Usage

### Load example data
```

```
### Contact to SpaCET
SpaCET also provide the mailngant score and immune infiltration score.
This function can Run SpaCET deconvolution and import the result from SpaCET object into Seurat object.
```

```
### manual annotation
<img width="882" alt="shniy" src="https://user-images.githubusercontent.com/122006615/236437705-821f551d-3853-4e28-beed-6f1a2c445be0.png">

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
### Visualization for result
```

```
![image](https://user-images.githubusercontent.com/122006615/235993210-3a841544-c772-4191-b2b1-55e456680756.png)


### Estimate the immune cell infiltration score
```
colon1 <- RunssGSEA(object = colon1,genelist = immune_list)

```

### DEboxplot
```
DEboxplot(colon)

```
![image](https://user-images.githubusercontent.com/122006615/233399369-7799a532-706a-4568-8cd6-ad1148519c84.png)

### Visualization for ssGSEA score
```
HexssGSEAplot(colon,type = 'C1QC_TAM')

```
![image](https://user-images.githubusercontent.com/122006615/233401762-76ee47aa-95de-4b4f-babf-c0da475f9752.png)


### Calculate enrichment fold
```
```
### cytokine differential analysis

```

```



## Citations
- Gao S, Shi Q, Zhang Y, et al. Identification of HSC/MPP expansion units in fetal liver by single-cell spatiotemporal transcriptomics. Cell Res. 2022;32(1):38-53. doi:10.1038/s41422-021-00540-7
- https://github.com/data2intelligence/SpaCET
