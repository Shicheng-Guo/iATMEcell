# iATMEcell

A software R package for identification of abnormal tumor microenvironment cells

# Introduce

> Here, we propose a network-based calculated method, iATMEcell, to identify abnormal tumor microenvironment cells with gene expression data, and cell crosstalk network. Furthermore, the package can perform regression analysis to verify cells prognostic efficacy. There are also some functions used to visualize the results of survival analysis.

This network-based method consists of two major parts:

  - 1.Calculate the DEscore. We conducted a statistical comparison of gene expression values between case and control groups (e.g. disease and normal, dead and alive). In this case, we use Student’s t-test method to compute the gene differentially expressed level for each gene between dead and alive samples, and convert the t-test p-value of each gene to z-score. The z-score is defined as DEscore, and a larger DEscore indicates the gene regulated by survival status to a greater extent.

  - 2.Constructing network and randomization. In our method, we fist constructed a cell-GO bipartite network. We define an edge between a cell and a Go term, if they have a common gene, and give the weight of this edge that calculate by the Jaccard index and DEscore. Next, we made the cell-GO network convert to cell-cell network, similarly, we define an edge between two cells, if they have common biology function, and the edge weights will be larger for pairs of cells that relate more to GO and survival. Then, we use eigenvector centrality measure to calculate how central each cell is in this network. Finally, the significance of these centrality scores is assessed using a bootstrap-based randomization method.

   This package provides the **GetExampleSet** function to return example data set and environment variables, such as the gene expression profile and so on.

# How to install
> Installation method：
```
library(devtools); 
install_github("hanjunwei-lab/iATMEcell", build_vignettes = TRUE)

Use：
library(iATMEcell)
```
# Example 1 
Calculate the DEscore, constructing network and randomization
>  The function **iTMEcell** including the two major parts of this network-based method used to find significant abnormal cells. Users need to input a gene expression profile of specific disease and corresponding survival data of samples. And the last variable is nperm, that representative number of disturbances, usually nperm = 1000 or bigger.
```
# load depend package
library(igraph)
#Obtain input data
GEP<-GetExampleSet('GEP')
clinicaldata<-GetExampleSet('clinicaldata')
#Run the function
iTMEcellresult<-iTMEcell(ExpData=GEP,clinical=clinicaldata,nperm=1000)
```
# Example 2
Constructing the cox risk regression model with cell marker genes
> The function **RiskRegressModel** is used to construct the cox risk regression model. Users can specify a cell, and then the function will perform cox regression analysis on expression of marker genes of the cell and survival data. Statistical significant genes will be selected, and with these genes, a gene risk score model was constructed using a formula derived from the expression of the genes weighted by their cox proportional hazards regression coefficient.
```
library(survival)
#Run the function
R.result<-RiskRegressModel(cellname='M1 Macrophages',ExpData=GEP,clinical=clinicaldata,p.cutoff=0.01,method = 'univ')
```
# Visualize 1 
Draw a forest plot
> The function **plotforest** can visualize the result of cox regression analysis through forest plot. 
```
library(forestplot)
library(survival)
#Run the function
plotforest(Regress.list=R.result,p.cutoff=0.01)
```
# Visualize 2
Draw a Kaplan-Meier curve
> The function **plotKMcurve** is used to draw the Kaplan-Meier curve according to the riskscore of samples from function **RiskRegressModel**.
```
library(survminer)
library(survival)
#Run the function
plotKMcurve(Regress.list=R.result,ExpData=GEP)
```
# Visualize 3
Draw a heat map
> The function **plotHeatmap** is used to draw a heat map of marker genes.
```
library(pheatmap)
#Run the function
plotHeatmap(Regress.list=R.result,ExpData=GEP,p.cutoff=0.01)
```
# Visualize 4
Draw a split violin plot
> The function **plotSplitViolin** is used to draw a split violin plot of gene expression.
```
library(ggplot2)
library(reshape2)
library(plyr)
library(ggpubr)
#Run the function
plotSplitViolin(Regress.list=R.result,ExpData=GEP,gene.name="CD96")
```

