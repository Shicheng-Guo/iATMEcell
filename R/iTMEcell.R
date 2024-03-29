random_network<-function(cell_random){
  adj.final1<-as.matrix(cell_random)
  graph1 = graph.adjacency(adj.final1,mode=c("undirected"), weighted=TRUE,add.rownames=T)
  temp1 = page.rank(graph1, vids=V(graph1), directed=FALSE, damping=0.90, weights=NULL)
  rank1 = temp1$vector
  rank2 = as.matrix(rank1)
  return(rank2)
}

go_cell_score_row<-function(genecount,descore,cell_size){
  score_row<-rep(0,cell_size)
  for(j in 1:length(genecount)){
    if(genecount[j]!=""){
      gene<-unlist(strsplit(genecount[j], split = ","))
      location<-match(gene,descore[,1])
      de1<-descore[location,2]
      de_score1<-median(as.numeric(de1))
      if (!is.na(de_score1)) {
        score_row[j]<-de_score1
      }
    }
  }
  return(score_row)
}

CalDEscore <- function(ExpData, clinical){
  samplename<-intersect(clinical[,"sample"],colnames(ExpData))
  ExpData<-ExpData[,samplename]
  clinicaldata<-clinical[samplename,]

  Label<-clinicaldata[,2]
  expdata<-as.matrix(ExpData)
  test<-matrix(nrow=nrow(expdata),ncol=1)
  rownames(test)<-rownames(expdata)
  colnames(test)<-c("Zscore")
  Z.vector<-apply(expdata, 1, function(x){
    ind1 <- which(Label == 1)
    ind2 <- which(Label == 0)
    m <- length(ind1 <- which(Label == 1))
    n <- length(ind2 <- which(Label == 0))
    expdata1 <- x[ind1]
    expdata2 <- x[ind2]
    rmean1 <- mean(expdata1)
    rmean2 <- mean(expdata2)
    ss1 <- sum((expdata1 - rmean1)^2)
    ss2 <- sum((expdata2 - rmean2)^2)
    Tscore <- (m + n - 2)^0.5*(rmean2 - rmean1)/((1/m + 1/n)*(ss1 + ss2))^0.5
    Pvalue <- pt(abs(Tscore),lower.tail=F,df=m+n-2)
    Zvalue <- qnorm(Pvalue, mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)
    return(Zvalue)
  })
  test[,1]<-Z.vector
  test <- as.matrix(test[!is.infinite(test[,1]),])
  test <- as.matrix(test[!is.na(test[,1]),])
  test <- as.matrix(test[!is.nan(test[,1]),])
  #test<-as.matrix(test[,1]/(sum(test[,1]^2))^0.5)
  #test[,1] <- scale(test[,1])
  colnames(test)<-c("Zscore")
  return(test)
}

##' @title Identification of abnormal tumor microenvironment (TME) cells
##' @description The function "iTMEcell" is used to calculate the eigenvector centrality of TME cells and identify abnormal TME cells.
##' @param ExpData A gene expression profile of interest (rows are genes, columns are samples).
##' @param clinical A dataframe with three columns which are "sample" (sample id),"status" (survival status of samples, "0" represents live and "1" represents dead) and "time" (survival time of samples).
##' @param nperm Number of random permutations (default: 1000).
##' @return A dataframe with seven columns those are cell names, marker source, marker size, marker genes, centrality (eigenvector centrality), P-value and FDR.
##' @importFrom igraph graph.adjacency
##' @importFrom igraph V
##' @importFrom igraph page.rank
##' @importFrom stats median
##' @importFrom stats qnorm
##' @importFrom stats pt
##' @importFrom stats p.adjust
##' @usage iTMEcell(ExpData,clinical,nperm=1000)
##' @export
##' @examples
##' library(igraph)
##' #Obtain input data
##' GEP<-GetExampleSet('GEP')
##' clinicaldata<-GetExampleSet('clinicaldata')
##' #Run the function
##' iTMEcellresult<-iTMEcell(ExpData=GEP,clinical=clinicaldata,nperm=1000)

iTMEcell <- function(ExpData,clinical,nperm=1000){


  TMEcell<-GetExampleSet("TMEcellinfo")
  Go<-GetExampleSet("TME_related_Goterm")
  Jaccard<-GetExampleSet("Jaccardscore")
  Go_TMEcell_gene<-GetExampleSet("GoCellconGene")

  DEscore<-CalDEscore(ExpData,clinical)
  cell<-cbind(as.character(TMEcell[,"TMEcells"]),as.character(TMEcell[,'Genes']))
  cell_size<-length(cell[,1])
  go<-as.matrix(Go)
  go_size<-length(go[,"Go_BP"])
  jaccard<-as.matrix(Jaccard)
  go_path_gene<-as.matrix(Go_TMEcell_gene)
  DEscore<-cbind(rownames(DEscore),abs(as.numeric(DEscore[,1])))
  median_score<-matrix(0,nrow=go_size,ncol=cell_size)
  for(k in 1:go_size){
    con_gene<-go_path_gene[k,]
    row<-go_cell_score_row(con_gene,DEscore,cell_size)
    median_score[k,]<-row
  }
  go_cell<-median_score*jaccard
  colnames(go_cell)<-cell[,1]
  rownames(go_cell)<-go[,1]
  #######构建cell-cell网
  edge<-as.matrix(go_cell)
  edget<-t(edge)
  cell_cell<-edget%*%edge
  #######计算centra,
  adj.final<-as.matrix(cell_cell)
  graph = graph.adjacency(adj.final,mode=c("undirected"), weighted=TRUE,add.rownames=T)
  temp = page.rank(graph, vids=V(graph), directed=FALSE, damping=0.90, weights=NULL)
  rank = temp$vector
  rank1 = as.matrix(rank)
  ########扰动网络
  iter<-nperm
  Centrality_Scores<-matrix(nrow=cell_size,ncol=iter+1)
  real.centra<-rank1
  Centrality_Scores[,1]<-real.centra
  real.subname<-rownames(real.centra)
  cells<-colnames(cell_cell)
  for(i in 1:iter){
    per_TMEcell<-sample(cells,replace = F)
    per_cell_cell<-cell_cell
    colnames(per_cell_cell)<-per_TMEcell
    rownames(per_cell_cell)<-per_TMEcell
    per_centra<-random_network(per_cell_cell)
    location<-match(real.subname,per_TMEcell)
    per_centra1<-per_centra[location,]
    Centrality_Scores[,i+1]<-per_centra1
  }
  adj = as.matrix(Centrality_Scores)
  perm_rank = adj[,2:(iter+1)]
  perm_rank<-as.matrix(perm_rank)
  orig_rank = adj[,1]
  pval = matrix(data=NA, nrow=cell_size, ncol=1)
  for ( j in 1:cell_size ) {
    pval[j] = sum(perm_rank[j,] > orig_rank[j])/iter
  }

  p_padjust<-p.adjust(pval,method = "fdr")
  pa<-as.numeric(p_padjust)
  fdr<-round(pa,3)
  allresult<-cbind(TMEcell,rank1)
  allresult<-cbind(allresult,pval,fdr)
  colnames(allresult)<-c("TMEcells","Source","MarkerSize","gene","Centrality","P-value","FDR")
  allresult<-as.data.frame(allresult)
  allresult[,3]<-as.numeric(as.character(allresult[,3]))
  allresult[,5]<-as.numeric(as.character(allresult[,5]))
  allresult[,6]<-as.numeric(as.character(allresult[,6]))
  allresult[,7]<-as.numeric(as.character(allresult[,7]))
  allresult[,5]<-round(allresult[,5],5)

  ###rank
  rankresult<-allresult[order(allresult$`P-value`), ]
  rownames(rankresult)<-c(1:dim(cell)[1])
  return(rankresult)
}


