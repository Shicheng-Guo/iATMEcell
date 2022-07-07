
##' @title Constructing the cox risk regression model with cell marker genes
##' @description Users can specify a cell, and then the function will perform cox regression analysis on expression of marker genes of the cell and survival data. Statistical significant genes will be selected, and with these genes, a gene risk score model was constructed using a formula derived from the expression of the genes weighted by their cox proportional hazards regression coefficient.
##' @param cellname A cell whose marker genes will be used to perform regression analysis. The format of the entered cell name should refer to the cell information we provide.
##' @param ExpData A gene expression profile of interest (rows are genes, columns are samples).
##' @param clinical A dataframe with three columns which are "sample" (sample id),"status" (survival status of samples, "0" represents live and "1" represents dead) and "time" (survival time of samples).
##' @param marker A character vector composed with marker genes. If you does not want to use the marker genes provided by us, you can specify the marker genes you need with this parameter.
##' @param p.cutoff Statistical significance threshold for regression analysis (default: 0.05).
##' @param method Users can specify the univariate cox regression analysis (default, method = "univ") or multivariate cox regression analysis (method = "muti") in this function.
##' @return A list with two dataframes which are riskscores of samples and result of cox regression analysis respectively.
##' @importFrom survival coxph
##' @importFrom survival Surv
##' @importFrom stats as.formula
##' @usage RiskRegressModel(cellname,ExpData,clinical,marker=NULL,p.cutoff=0.01,method='univ')
##' @export
##' @examples
##' library(survival)
##' #Obtain input data
##' GEP<-GetExampleSet('GEP')
##' clinicaldata<-GetExampleSet('clinicaldata')
##' #Run the function
##' R.result<-RiskRegressModel(cellname='NK cells',ExpData=GEP,
##'   clinical=clinicaldata,p.cutoff=0.01,method = 'univ')


RiskRegressModel<-function(cellname,ExpData,clinical,marker=NULL,p.cutoff=0.01,method='univ'){
  if(is.null(marker)){
    TMEcellinfo<-GetExampleSet('TMEcellinfo')
    markergenes<-as.data.frame(TMEcellinfo)[which(TMEcellinfo[,1]== cellname),4]
    markergenes<-as.character(markergenes)
    markergenes<-unlist(strsplit(markergenes,','))
  }else{
    markergenes<-marker
  }

  congenes<-intersect(markergenes,rownames(ExpData))
  markergep<-ExpData[congenes,]
  markergep<-t(markergep)

  samplename<-intersect(clinical[,"sample"],rownames(markergep))
  markergep<-markergep[samplename,]
  clinicaldata<-clinical[samplename,]

  markergep<-markergep[match(clinicaldata[,"sample"],rownames(markergep)),]
  survframe<-cbind(clinicaldata,markergep)

  if(method=='univ'){
    covariates <- colnames(survframe[,4:length(survframe[1,])])

    univ_formulas <- sapply(covariates,function(x) as.formula(paste('Surv(time, status) ~',"`",x,"`",sep = "")))

    univ_models <- lapply( univ_formulas, function(x){coxph(x, survframe)})

    univ_results <- lapply(univ_models,
                           function(x){
                             x <- summary(x)
                             p.value<-signif(x$wald["pvalue"], digits=2)
                             #wald.test<-signif(x$wald["test"], digits=2)
                             beta<-signif(x$coef[1], digits=2)
                             HR <-signif(x$coef[2], digits=2)
                             HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                             HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                             #HR <- paste0(HR, " (",
                             #            HR.confint.lower, "-", HR.confint.upper, ")")
                             res<-c(beta, HR,HR.confint.lower,HR.confint.upper, p.value)
                             names(res)<-c("beta", "HR", "HR.95L", "HR.95H","p.value")
                             return(res)
                             #return(exp(cbind(coef(x),confint(x))))
                           })

    result_cox <- data.frame(t(as.data.frame(univ_results, check.names = FALSE)))
    result_cox[,'beta']<-as.numeric(as.character(result_cox[,'beta']))
    result_cox[,'p.value']<-as.numeric(as.character(result_cox[,'p.value']))
    result_cox<-result_cox[which(result_cox[,'p.value']<p.cutoff),]
    sig.genes<-rownames(result_cox)
    sig.survframe<-cbind(survframe[,1:3],survframe[,which(colnames(survframe)%in%sig.genes)])

    riskscore<-c()
    for (i in 1:dim(sig.survframe)[1]) {
      rs<-sum(sig.survframe[i,4:dim(sig.survframe)[2]]*result_cox[match(colnames(sig.survframe)[4:dim(sig.survframe)[2]],rownames(result_cox)),'beta'])
      riskscore<-c(riskscore,rs)
    }

    risk.result<-sig.survframe[,1:3]
    risk.result$riskscore<-riskscore

    result.list<-list()
    result.list[[1]]<-risk.result
    result.list[[2]]<-result_cox

  }

  if(method=='muti'){
    muti_cox<-summary(coxph(Surv(time,status)~.,data =survframe[,-1]))
    muti_cox_frame<-cbind(as.data.frame(muti_cox[[7]])[,c(1,2)],as.data.frame(muti_cox[[8]])[,3:4],as.data.frame(muti_cox[[7]])[,5])
    muti_cox_frame<-round(muti_cox_frame,2)
    colnames(muti_cox_frame)<-c("beta", "HR", "HR.95L", "HR.95H","p.value")

    muti_cox_frame<-muti_cox_frame[which(muti_cox_frame[,'p.value']<p.cutoff),]
    sig.genes<-rownames(muti_cox_frame)
    sig.survframe<-cbind(survframe[,1:3],survframe[,which(colnames(survframe)%in%sig.genes)])
    riskscore<-c()
    for (i in 1:dim(sig.survframe)[1]) {
      rs<-sum(sig.survframe[i,4:dim(sig.survframe)[2]]*muti_cox_frame[match(colnames(sig.survframe)[4:dim(sig.survframe)[2]],rownames(muti_cox_frame)),'beta'])
      riskscore<-c(riskscore,rs)
    }

    risk.result<-sig.survframe[,1:3]
    risk.result$riskscore<-riskscore

    result.list<-list()
    result.list[[1]]<-risk.result
    result.list[[2]]<-muti_cox_frame

  }

  names(result.list)<-c('riskscore.frame','cox.result')
  return(result.list)
}
