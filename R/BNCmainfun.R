#'@import mclust
#' @title {Bayesian Node Classification with Missing Values on The Biological Network}
#' @aliases BNC
#' @name BNC
#' @usage
#' BNC(net,test.stat,pvalue.stat=FALSE,candidate.z.set=c(-1,0,1),
#' seed.main=1024,na.action=c("NN","Bayes","na.remove"),niter.densupd=5,niter=10,
#' paras=list(tau=c(2,10,2),alpha=NULL,gamma=NULL,xi=NULL, beta=rep(10,3),
#' rho=c(1.003,0.479,0.988,0.000),pivec=c(0.15,0.7,0.15),densAcc=0.001,
#' null.quantile=c(0.25, 0.75),null.method="biGaussianModeSplit",
#' transitionMatrix.Z.11=0.6,miss.stat=2,min.node=5),
#' para.DPM=NULL,para.HODC=NULL,para.DMH=NULL)
#' @description Main function. Two steps: Given density specification, update selection indicator z by Swendsen- Wang algorithm; and given selection indicator z, update density specification by DPM fitting algorithm.
#' @param net The binary adjacent matrix representing the gene network configuration. The element at (i,j) location of value 0 indicating "gene i and gene j are connected" or "gene i and gene j are not directly connected"
#' @param test.stat The observed test statistics. Missing values are represented as NA. If the test statistics are p-values, the boolean label indicating the statistics are p-values: pvalue.stat should be TRUE otherwise FALSE
#' @param pvalue.stat Logical. Indicating whether the test statistics are generated as pvalues or not. The default is FALSE
#' @param candidate.z.set A set of all the possible values of the regulation type, the values of indicators. The default is (-1,0,1), 1=down-regulated genes, 2=not differentially expressed genes, 3=up-regulated genes
#' @param seed.main The random seed used for generating reproducible results. The default is 1024
#' @param na.action The method used for handling genes with missing values. It can be "NN" for nearest neighbor imputation method, "Bayes" for fully Bayesian inference for imputation, or "na.remove" for removing all the missing nodes and their edges in the network
#' @param niter.densupd The total number of iterations for updating density. The default is 5
#' @param niter The total number of iterations for main algorithm. The default is 10
#' @param paras A list contains hyper-parameters and other parameters used before main algorithm runs
#' \itemize{
#' \item niter.densupd The iteration is from 1 to the maximum steps when we update density specification by DPM. The default is 20
#' \item tau A vector of length 3, the default is c(2,10,2)
#' \item alpha A vector of length 3. The default is NULL
#' \item gamma A vector of length 3. The default is NULL
#' \item xi A vector of length 3. The default is NULL
#' \item beta A vector of length 3. The default is c(10,10,10)
#' \item rho A vector of length 4. The default is c(1.003,0.479,0.988,0.000), indicating local smoothness for Potts prior. Note: the default value is calculated based on data(net) strucutre by DMH
#' \item pivec A vector of length 3. The default is c(0.15,0.7,0.15). Contains prior knowledge globally about selection indicator z
#' \item densityAcc A number, need to specify precision for K-L integration when to use the numerical approximation. The default is 0.001
#' \item null.quantile A vector of length 2, representing the lower quantile and the upper quantile for calculating prior null density if lack of prior biology knowledge. The default is c(0.25, 0.75)
#' \item null.method A string. The name of method we used to estimate the null density. It can be  "biGaussianMean0"-- the null density is formed by two half normals. They share the same mean of 0. "biGaussianModeSplit"-- split data from their median value, then flip each part to the other side to estimate normal distribution
#' \item transitionMatrix.Z.11 The [1,1] element in transition matrix for z. The default is 0.6
#' \item miss.stat To impute NAs in the test statistics when apply Double Metropolis-Hastings sampler (DMH) to find hyperparameters rho and pi
#' \item min.node The minimum number of nodes in class of regulations
#' }
#' @param para.DPM A list of parameters used for fitting Dirichlet process mixture (DPM) model. The default is NULL when all members will use their built-in default values
#' \itemize{
#' \item niter The number of iterations. The default is 10
#' \item nsample The number of samples. The default is 10
#' \item KLrange A vector of length 2, the first element is integration lower bound and the second element is integration upper bound in calculating K-L distance of the proposed null density vs. prior null density. The default is c(-6,6), usually we suggest wider range than the boundaries of observed test statistics
#' \item KLprecision A number. The numerical precision for integration calculating the K-L distance. The default is 0.001
#' \item KLNullmethod The name of the method to calculate the proposal density. The default is "biGaussianMean0"
#' \item mcmc A list of the MCMC parameters for DPM. nburn is the number of burn-in iterations, nskip is the thinning interval, nsave is the total number of iterations to be saved, and ndisplay is the number of saved iterations to be displayed on screen. The default value is nburn=10000, nsave=100, nskip=0, ndisplay=10000
#' \item prior A list of the prior information for DPM, alpha is the value of the precision, m1 is the mean of the normal part of the baseline distribution, nu1 and psiinv1 is the hyperparameters of the inverted Wishart part of the baseline distribution, tau1 and tau2 is the hyperparameters for the gamma prior distribution. The default is alpha=3, m1=rep(0,1), psiinv1=diag(0.5,1), nu1=4, tau1=1, tau2=100
#'}
#' @param para.HODC A list of the parameters used for KL-HODC algorithm to find the proper initial values for fast convergence. The default is NULL when all members will use their built-in default values
#' \itemize{
#' \item nsample The number of samples. The default is 10
#' \item KLrange A vector of length 2, the first element is integration lower bound and the second element is integration upper bound in calculating K-L distance of the proposed null density vs. prior null density. The default is c(-6,6), usually we suggest wider range than the boundaries of observed test statistics
#' \item KLprecision A number. The numerical precision for integration calculating the K-L distance. The default is 0.001
#' \item KLNullmethod The name of the method to calculate the proposal density. The default is "biGaussianMean0"
#' \item mcmc A list of the MCMC parameters for DPM. nburn is the number of burn-in iterations, nskip is the thinning interval, nsave is the total number of iterations to be saved, and ndisplay is the number of saved iterations to be displayed on screen. The default value is nburn=1000, nsave=100, nskip=0, ndisplay=1000
#' \item prior A list of the prior information for DPM, alpha is the value of the precision, m2 and s2 is the mean and the covariance of the normal prior for the mean, nu2 and psiinv2 are the hyperparameters of the inverted Wishart prior distribution for the scale matrix, which is calculated from data along the way, nu1 is the hyperparameter of the inverted Wishart part of the baseline distribution, tau1 and tau2 is the hyperparameters for the gamma prior distribution. The default is alpha=3, m2=rep(0,1), s2=diag(100000,1), nu1=4, nu2=4, tau1=1, tau2=100
#'}
#' @param para.DMH If inside paras, the rho and pivec is not given, Double Metropolis-Hastings algorithm (DMH) is used for pre-calculating rho and pivec. The default parameter list is used which contains
#' \itemize{
#' \item niter The number of iterations in total. The default is 1000
#' \item pistat The initial values of pivec. The default is c(0.25,0.5,0.25)
#' \item pisd The vector of standard deviations when sampling pivec. The default is c(0.03,0.03,0.03)
#' \item rhostat The initial values of rho. The default is c(1,0.5,1,0)
#' \item rhosd The vector of standard deviations when sampling rho. The default is c(0.03,0.03,0.03,0.03)
#' \item rhoLowB The vector of lower bounds in sampling rho. The default is c(0,0,0,0)
#' \item rhoUpB The vector of upper bounds in sampling rho. The default is c(1.5,1.5,1.5,1.5)
#' \item piLowB The vector of lower bounds in sampling pi. The default is c(0,0,0)
#' \item piUpB The vector of upper bounds in sampling pi. The default is c(1,1,1)
#' \item niter.z The number of iterations within each iteration updating rho and pivec to update z. The default is 1
#' \item replaceInf The value used to replace Inf. The default is -99999
#' \item DMHplot Logical. Whether to plot. The default is FALSE
#' }
#'
#' @return A list contains
#' \item{initialValue}{The initial parameter list used in main algorithm}
#' \item{zTrack}{The trace for z}
#' \item{FinalValue}{The final parameter list after convergence}
#' \item{iters}{The total number of iterations}
#' \item{rmisTrack}{The trace for imputed values of the missing values in test statistics. Only applicable when there are missing values in the test statistics and imputation methods are used}
#' @details The fully Bayesian updating algorithm is executed in the following procedure
#' \itemize{
#' \item Input data r and gene network represented by graph G=<V,E>
#' \item Update z|theta via Swendsen-Wang
#' \item Update theta|z via DPM Fitting
#' }
#' @examples
#' \dontrun{
#' ## The simulation settings based on real gene network (takes time)
#' data("geneNetwork")
#' data("testStatistics")
#' ## three different regulation types- three class classification example
#' res=BNC(net=geneNetwork,test.stat=testStatistics,niter=300,na.action="NN")
#' ## three different regulation types- two class classification example
#' res=BNC(net=geneNetwork,test.stat=pnorm(testStatistics),pvalue.stat=TRUE,
#' candidate.z.set=c(0,1),na.action="NN",niter=300,
#' paras=list(tau=c(2,10),alpha=NULL,gamma=NULL,xi=NULL, beta=rep(10,2),rho=c(1,0.5,0),
#' pivec=c(0.2,0.8),densAcc=0.001,null.quantile=c(0.25, 1),
#' null.method="biGaussianModeSplit",transitionMatrix.Z.11=0.6,miss.stat=2,min.node=5))
#'
#' ## A toy example
#' simdata=simulatedDatasetGenerator(nnode=100,missing=TRUE,missrate=0.1,dist="norm",
#' plot=FALSE,nbin=c(20,20,10),rng=1024)
#' res=BNC(net=simdata$net,test.stat=simdata$testcov,niter=100,na.action="NN",
#' paras=list(tau=c(2,10,8),alpha=NULL,gamma=NULL,xi=NULL,beta=rep(10,3),
#' rho=c(1.003,0.479,0.988,0.000),pivec=c(0.25,0.40,0.35),densAcc=0.001,
#' null.quantile=c(0.25, 0.75),null.method="biGaussianModeSplit",transitionMatrix.Z.11=0.5,
#' miss.stat=2,min.node=5))
#' classLabelEst=summaryClassLabels(simdata$net,simdata$testcov,res$zTrack,
#' method="MajorVote",nburn=50)
#' print(table(classLabelEst))
#' }
#' @export
BNC=function(net,test.stat,pvalue.stat=FALSE,candidate.z.set=c(-1,0,1),seed.main=1024,na.action=c("NN","Bayes","na.remove"),niter.densupd=5,niter=10,paras=list(tau=c(2,10,2),alpha=NULL,gamma=NULL,xi=NULL, beta=rep(10,3),rho=c(1.003,0.479,0.988,0.000),pivec=c(0.15,0.70,0.15),densAcc=0.001,null.quantile=c(0.25, 0.75),null.method="biGaussianModeSplit",transitionMatrix.Z.11=0.6,miss.stat=2,min.node=5),para.DPM=NULL,para.HODC=NULL,para.DMH=NULL){

  if(!requireNamespace("mclust", quietly = TRUE)) {
    stop("Please load package: mclust as required ")
  }
  if(!requireNamespace("dirichletprocess", quietly = TRUE)) {
    stop("Please load package: dirichletprocess as required ")
  }

  set.seed(seed.main)
  # library(DPpackage,verbose = F,quietly = T)
  # library(mclust,verbose = F,quietly = T)
  if(pvalue.stat){
    test.stat=-stats::qnorm(test.stat)
  }

  #### if NA.action=NA.remove:
  if(sum(is.na(test.stat))>0){
    method.na.impute=match.arg(na.action)
    if(method.na.impute=="na.remove"){
      na.loc=which(is.na(test.stat))
      net=net[-na.loc,-na.loc]
      test.stat=test.stat[-na.loc]
    }
  }else{
    method.na.impute="NoImpute"
  }
  #### create a parameter list:
  numz=length(candidate.z.set)
  paraList=list(zSet=seq(1,numz),numZ=numz,pivec=paras$pivec,rho=paras$rho, rStat=test.stat, misLoc=which(is.na(test.stat)),numNodes=length(test.stat), zLab=NULL, tau=paras$tau, beta=paras$beta, alpha=rep(3,numz), gamma=NULL, xi=NULL, net=net, densityAccu=paras$densAcc,null.quantile=paras$null.quantile,null.method=paras$null.method,TransZ11=paras$transitionMatrix.Z.11,densIter=list(to=niter.densupd,niter=niter),min.node=paras$min.node)
  if(length(paraList$misLoc)>0){
    Imputeflag=TRUE
  }else{
    Imputeflag=FALSE
  }
  paraList$priorNull=biGaussianNull(paraList$rStat,paraList$null.quantile,method=paraList$null.method)
  # temp=mixtools::normalmixEM(paraList$rStat[!is.na(paraList$rStat)],k=numz)
  # paraList$mulist=as.list(temp$mu)
  # paraList$sdlist=as.list(temp$sigma)
  temp=mclust::Mclust(paraList$rStat[!is.na(paraList$rStat)],G=length(candidate.z.set),modelName="E")
  paraList$mulist=as.list(temp$parameter$mean)
  temp.sdlist=temp$parameter$variance$sigmasq
  if(length(temp.sdlist)!=length(candidate.z.set)){
    temp.sdlist=rep(sqrt(temp.sdlist),length(candidate.z.set))
  }
  paraList$sdlist=as.list(temp.sdlist)
  # paraList$zLab=temp$classification
  paraList$proplist=as.list(rep(1,length(candidate.z.set)))

  if(is.null(paras$alpha)){
    paraList$alpha=paraList$beta/unlist(paraList$sdlist)+1
  }else{
    paraList$alpha=paras$alpha
  }
  if(is.null(paras$gamma)){
    paraList$gamma=unlist(paraList$mulist)
  }else{
    paraList$gamma=paras$gamma
  }
  if(is.null(paras$xi)){
    paraList$xi=unlist(paraList$sdlist)*sqrt(2)
  }else{
    paraList$xi=paras$xi
  }
  paraList$TransZ=TransZ(paraList$TransZ11)

  # mcmc.default=list(nburn=1000,nsave=100,nskip=0,ndisplay=1000)
  mcmc.default=list(nburn=10000,nsave=1000,nskip=0,ndisplay=100000)
  # prior.default=list(alpha=3,m2=rep(0,1),s2=diag(100000,1),psiinv2=diag(100000,1),nu1=4,nu2=4,tau1=1,tau2=100)
  prior.default=list(alpha=3,m1=0,psiinv1=0.5,nu1=4,tau1=1,tau2=100)
  para.DPM.default=list(niter=10,nsample=10,KLrange=c(floor(min(test.stat,na.rm=TRUE)),ceiling(max(test.stat,na.rm=TRUE))),KLprecision=0.001,KLNullmethod="biGaussianMean0")
  vars=setdiff(names(mcmc.default),names(para.DPM$mcmc))
  para.DPM$mcmc[vars]=lapply(vars,function(x){mcmc.default[[x]]})
  vars=setdiff(names(prior.default),names(para.DPM$prior))
  para.DPM$prior[vars]=lapply(vars,function(x){prior.default[[x]]})
  vars=setdiff(names(para.DPM.default),names(para.DPM))
  para.DPM[vars]=lapply(vars,function(x){para.DPM.default[[x]]})
  paraList$DPM=para.DPM


  paraList=BNC::DPdensCluster(paraList,twogrps=FALSE)
  utils::flush.console()
  paraList$zLab[!is.na(paraList$rStat)]=temp$classification
  paraList$zLab[is.na(paraList$rStat)]=rep(2,length(paraList$misLoc))
  paraList$mk=table(paraList$zLab)[paraList$zSet]

  mcmc.default=list(nburn=1000,nsave=100,nskip=0,ndisplay=1000)
  vars=setdiff(names(mcmc.default),names(para.HODC$mcmc))
  para.HODC$mcmc[vars]=lapply(vars,function(x){mcmc.default[[x]]})
  prior.default=list(alpha=3,m2=rep(0,1),s2=diag(100000,1),psiinv2=diag(100000,1),nu1=4,nu2=4,tau1=1,tau2=100)
  para.HODC$prior=lapply(1:numz,function(x){
    vars=setdiff(names(prior.default),names(para.HODC$prior[[x]]))
    para.HODC$prior[[x]][vars]=lapply(vars,function(x){prior.default[[x]]})
    return(para.HODC$prior[[x]])})
  vars=setdiff(names(para.DPM.default),names(para.HODC))
  para.HODC[vars]=lapply(vars,function(x){para.DPM.default[[x]]})
  paraList$HODCPara=para.HODC

  if(is.null(paraList$rho) | is.null(paraList$pivec)){
    para.DMH.default=list(niter=1000,pistat=c(0.25,0.5,0.25),pisd=rep(0.03,3),rhostat=c(1,0.5,1,0),rhosd=rep(0.03,4),rhoLowB=c(0,0,0,0),rhoUpB=c(1.5,1.5,1.5,1.5),piLowB=c(0,0,0),piUpB=c(1,1,1),niter.z=1,replaceInf=-99999,DMHplot=FALSE,mean.z=temp$parameters$mean,sd.z=temp.sdlist)

    vars=setdiff(names(para.DMH.default),names(para.DMH))
    para.DMH[vars]=lapply(vars,function(x){para.DMH.default[[x]]})
    mypara=paraList
    mypara$hyperPara=para.DMH
    mypara$rStat[is.na(paraList$rStat)]=rep(mean(mypara$rStat,na.rm=TRUE),sum(is.na(mypara$rStat)))
    hyperRes=hyperParaDMH(mypara,plot=mypara$hyperPara$DMHplot)
    paraList$rho=hyperRes$rho
    paraList$pivec=hyperRes$pi
  }

  paraList$rho=c(paraList$rho[2],paraList$rho)
  paraList$pivec=c(paraList$pivec,paraList$pivec[1])
  print("Complete creating the parameter list")

  #### Main Loop for Fast S-W Updates
  simRes=list()
  simRes$initialValue=paraList
  zTrack=paraList$zLab
  iter=0
  if(Imputeflag){
    paraList$rStat[paraList$misLoc]=0
    rmisTrack=paraList$rStat[paraList$misLoc]
    if(method.na.impute=="NN"){
      paraList$nbrMis=lapply(paraList$misLoc,function(k){
        which(paraList$net[k,]>0)
      })
      names(paraList$nbrMis)=paraList$misLoc
    }
  }

  #### Main Loop + density updates:
  kk=1
  while(iter<paraList$densIter$to){
    ### updating z
    paraList=GCut.z(paraList)
    currtSWCut=SWCut(paraList)
    paraList=z.update.fastSW(paraList,currtSWCut)
    ### updating mu/sd/etc
    paraList=BNC::DensityDPOne(paraList)
    if(Imputeflag & method.na.impute=="NN"){
      paraList=r.knnImpute(paraList)
      rmisTrack=cbind(rmisTrack,paraList$rStat[paraList$misLoc])
    }else if(Imputeflag & method.na.impute=="Bayes"){
      paraList=r.BayesImpute(paraList)
      rmisTrack=cbind(rmisTrack,paraList$rStat[paraList$misLoc])
    }
    ## keep track of z changes:
    zTrack=cbind(zTrack,paraList$zLab)
    kk=kk+1
    iter=iter+1
  }
  print("Stop updating density, only update the class indicators")
  #### Main Loop - density updates:
  while(iter<paraList$densIter$niter){
    ### updating z
    paraList=GCut.z(paraList)
    currtSWCut=SWCut(paraList)
    paraList=z.update.fastSW(paraList,currtSWCut)
    if(Imputeflag & method.na.impute=="NN"){
      paraList=r.knnImpute(paraList)
      rmisTrack=cbind(rmisTrack,paraList$rStat[paraList$misLoc])
    }else if(Imputeflag & method.na.impute=="Bayes"){
      paraList=r.BayesImpute(paraList)
      rmisTrack=cbind(rmisTrack,paraList$rStat[paraList$misLoc])
    }
    ## keep track of z changes:
    zTrack=cbind(zTrack,paraList$zLab)
    iter=iter+1
    # print(iter)
  }
  print("Algorithm completed")
  #### Results
  colnames(zTrack)=NULL
  simRes$zTrack=zTrack
  simRes$FinalValue=paraList
  simRes$iters=iter
  if(Imputeflag){
    simRes$rmisTrack=rmisTrack
  }
  return(simRes)
}

#' @title {biGaussianNull}
#' @description Compute a prior null distribution without given any prior biological knowledge. Methods used to estimate this prior include: "biGaussian","biGaussianWith0Cutoff","biGaussianMean0","biGaussianModeSplit" etc.
#' @param rstat The observed test statistics.
#' @param null.quantile Data quantiles, pre-fixed so that we could reduce the impact of extreme values for estimating prior null distribution. Default value is c(0.25, 0.75).
#' @param method A char. The methods we used to do estimation: either "biGaussian"-- EM algorithm for mixtures of two univariate normals, or "biGaussianWith0Cutoff"-- assume all negative test statistics forms one normal and all positive test statistics forms the other one normal. And proportion parameter is proportional to a number of observations each class, or "biGaussianMean0"-- null is formed by two half normals. "biGaussianModeSplit"-- split data from median value, then flip each part to the other side to estimate normal distribution.
#' @return A list with element
#'  \item{alpha}{cutoff value (if given)}
#'  \item{mu}{mean positon for two normals}
#'  \item{sd}{sd value for two normals}
#'  \item{prop}{proportion for each normal}
#' @keywords internal
#' @export
biGaussianNull=function(rstat,null.quantile=c(0.25, 0.75),method=c("biGaussianMean0","biGaussianModeSplit")){
  if(sum(is.na(rstat))>0){
    rstat=rstat[!is.na(rstat)]
  }
  rr=stats::quantile(rstat,prob=null.quantile)
  rr.center=rstat[which(rstat>rr[1] & rstat<rr[2])]
  methods=match.arg(method)
  if(methods=="biGaussianMean0"){
    rr.pos=rr.center[rr.center>0]
    rr.neg=rr.center[rr.center<=0]
    sd=sqrt(c(stats::var(rr.pos),stats::var(rr.neg)))
    res=list(alpha=0,mu=c(0,0), sd=sd, prop=sd/sum(sd))
  }else if(methods=="biGaussianModeSplit"){
    cutoff=stats::quantile(rr,0.5)
    rr.pos=rr.center[rr.center>cutoff]
    rr.pos=c((2*cutoff-rr.pos),rr.pos)
    rr.neg=rr.center[rr.center<=cutoff]
    rr.neg=c(rr.neg,(2*cutoff-rr.neg))
    res=list(alpha=cutoff,mu=c(cutoff,cutoff), sd=sqrt(c(stats::var(rr.pos),stats::var(rr.neg))), prop=c(length(rr.neg),length(rr.pos))/(length(rr.neg)+length(rr.pos)))
  }
  return(res)
}

#' @title {TransZ}
#' @description Compute a special transition matrix of Z: symmetric, down-regulated genes are not able to jump directly to up-regulated genes; a probability of a down-regulated gene jumping to null gene equals the probability of an up-regulated gene jumping to a null gene.
#' @param p The probability of a down-regulated gene jumping to null gene
#' @return A 3*3 transition matrix for Z
#' @keywords internal
#' @export
TransZ=function(p){
  m=matrix(c(p,(1-p)/2,0,(1-p),p,(1-p),0,(1-p)/2,p),nrow=3)
  return(m)
}

#' @import dirichletprocess
#' @title {DPdensCluster}
#' @description First, apply DP-package to split data into several small normal groups based on their test statistics, then merge until we attain total three groups based on K-L distance between proposed null density vs. prior null density from biological background knowledge. Finally, return each group's mixture normal parameters as a list.
#' @param para The parameter list
#' @return The updated parameter list
#' @keywords internal
#' @export
DPdensCluster=function(para,twogrps=FALSE){
  #### initialized L-1 L0 L1 as sum(L_k)=fit$state$ncluster
  # set.seed(para$DPM$seed)
  #fit=DPpackage::DPdensity(y=para$rStat[!is.na(para$rStat)],prior=para$DPM$prior,mcmc=para$DPM$mcmc,state=NULL,status=TRUE)
  dp <- dirichletprocess::DirichletProcessGaussian(y = para$rStat[!is.na(para$rStat)],
                                                   goPriors = c(para$DPM$prior$m1,
                                                   para$DPM$prior$k0,
                                                   para$DPM$prior$nu1,
                                                   1/para$DPM$prior$psiinv1))
  dp <- dirichletprocess::Fit(dp,para$DPM$mcmc)
  utils::flush.console()
  #### sort sum(L_k)
  #nclust=fit$state$ncluster
  nclust = dp$numberClusters
  #prop=table(fit$state$ss)/sum(table(fit$state$ss))
  prop = table(dp$clusterLabels)/length(dp$clusterLabels)
  #mu=fit$state$muclus[1:nclust]
  mu = as.vector(dp$clusterParameters[[1]])
  #sigma=fit$state$sigmaclus[1:nclust]
  sigma = as.vector(dp$clusterParameters[[2]])
  index=order(mu,decreasing=FALSE)
  mu=mu[index]
  sigma=sigma[index]
  prop=prop[index]

  KLdist=c()
  for(id in c(1:nclust)){
    # print(id)
    KLdist=c(KLdist,BNC::myKLdivergenece(para$priorNull,list(mu=mu[id],sd=sigma[id],prop=prop[id]),integral=para$DPM$KLrange,precision=para$DPM$KLprecision,para$DPM$KLNullmethod))
  }
  nullclassID=which(KLdist==min(KLdist,na.rm=TRUE))[1]

  if(max(nullclassID)==1){nullclassID=2}
  if(min(nullclassID)==nclust){nullclassID=nclust-1}
  KLD=KLdist[nullclassID][1]
  KLdiff=KLD

  while(nclust>3 & KLdiff>0){
    KLdist=c()
    # print(nclust)
    for(id in c(min(nullclassID)-1,min(max(nullclassID)+1,length(mu)))){
      # print(id)
      KLdist=c(KLdist,BNC::myKLdivergenece(para$priorNull,list(mu=mu[c(nullclassID,id)],sd=sigma[c(nullclassID,id)],prop=prop[c(nullclassID,id)]),integral=para$DPM$KLrange,precision=para$DPM$KLprecision,para$DPM$KLNullmethod))
    }

    ll = which(KLdist==min(KLdist,na.rm=TRUE))[1]

    KLdiff=KLD-KLdist[ll]
    if(KLdiff>0){
      nullclassID=c(nullclassID,c(min(nullclassID)-1,max(nullclassID)+1)[ll])
      nullclassID=unique(nullclassID)
      KLD=KLdist[ll]
      nclust=nclust-1
    }
  }

  if(twogrps){
    idlist=list(c(1:(min(nullclassID)-1)),sort(c(nullclassID,(max(nullclassID)+1):length(mu))))
  }else{
    idlist=list(c(1:(min(nullclassID)-1)),sort(nullclassID),c((max(nullclassID)+1):length(mu)))
  }

  if(max(table(unlist(idlist)))==1){
    ### empty para$mulist:
    para$mulist=lapply(1:para$numZ,function(k){list()})
    para$sdlist=lapply(1:para$numZ,function(k){list()})
    para$proplist=lapply(1:para$numZ,function(k){list()})

    para$zLab=rep(2,length(para$rStat))
    for(k in 1:length(idlist)){
      para$mulist[[k]]=mu[idlist[[k]]]
      para$sdlist[[k]]=sigma[idlist[[k]]]
      para$proplist[[k]]=prop[idlist[[k]]]/sum(prop[idlist[[k]]])
      ll=index[idlist[[k]]]
      para$zLab[which(is.element(dp$clusterLabels,ll))]=k
    }
  }
  return(para)
}

#' @title {myKLdivergenece}
#' @description Compute Kullback Leibler divergence between biologic null distribution vs. proposed null distribution
#' @param nulldens A list, mu: mean location for each of normal density; sd: sd value for each of normal density; proportion: proportion for each normal density.
#' @param dens A list, mu: mean location for each of normal density; sd: sd value for each of normal density; proportion: proportion for each normal density.
#' @param integral A two element vector, the first element is integration lower bound; the second element is integration upper bound.
#' @param precision A number. The numerical precision for integration calculating. Default=0.001
#' @param method The way proposal density is calculated from wether it is "MixNorm" or "biGaussianMean0".
#' @return A number, KL distance between nulldens and dens
#' @keywords internal
#' @export
myKLdivergenece=function(nulldens,dens,integral=c(-6,6),precision=0.001,method=c("MixNorm","biGaussianMean0")){
  x=seq(integral[1],integral[2],by=precision)
  p=dens$prop/sum(dens$prop)
  y=sapply(1:length(p),function(pp){
    p[pp]*stats::dnorm(x,mean=dens$mu[pp],sd=dens$sd[pp])
  })
  y=rowSums(y,na.rm=TRUE)

  if(method=="biGaussianMean0"){
    ynull=x
    ll=which(x>=nulldens$alpha)
    ynull[ll]=stats::dnorm(x[ll],mean=nulldens$mu[2],sd=nulldens$sd[2])
    ll=which(x<nulldens$alpha)
    ynull[ll]=stats::dnorm(x[ll],mean=nulldens$mu[1],sd=nulldens$sd[1])
  }else{
    ynull=sapply(1:length(nulldens$prop),function(pp){
      nulldens$prop[pp]*stats::dnorm(x,mean=nulldens$mu[pp],sd=nulldens$sd[pp])
    })
    ynull=rowSums(ynull)
  }

  KL=flexmix::KLdiv(cbind(ynull,y),eps =1e-9)
  return(KL[1,2])
}

#' @title {hyperParaDMH}
#' @description Apply Double Metropolis-Hastings sampler (DMH) to find hyperparameters: rho & pi
#' @param para The parameter list
#' @param plot Logical. T: plot trace plot; F: not plot output.
#' @return a list:
#' \item{rhoTrack}{saved rho value each iteration}
#' \item{piTrack}{saved pi value each iteration}
#' \item{zTrack}{saved z label value each iteration}
#' \item{rho}{posterior mean calculated based on second half of the iterations}
#' \item{zstat}{posterior z value calculated based on major vote method.}
#' @keywords internal
#' @export
hyperParaDMH=function(para,plot=TRUE){
  rho=para$hyperPara$rhostat
  pi=para$hyperPara$pistat
  z=para$zLab
  rhoTrack=rho
  piTrack=pi
  zTrack=z
  rho.L=length(rho)
  pi.L=length(pi)
  for(iter in 1:para$hyperPara$niter){
    #### revise initial values for rho & pi
    repeat{
      newrho=truncnorm::rtruncnorm(rho.L,mean=rho,sd=para$hyperPara$rhosd,a=para$hyperPara$rhoLowB,b=para$hyperPara$rhoUpB)
      newrho[rho.L]=0
      newpi=truncnorm::rtruncnorm(pi.L,mean=pi,sd=para$hyperPara$pisd,a=para$hyperPara$piLowB,b=para$hyperPara$piUpB)
      newpi[2]=1-sum(newpi[-2])
      if(min(newrho[-c(2,rho.L)])>newrho[2] & newpi[2]>0.5) break;
    }
    ### updating z
    newz=BNC::z.update.Mclustprior(para,z,newrho,newpi)
    AccRate=BNC::llikDMH(para,newz,rho,pi)+BNC::llikDMH(para,z,newrho,newpi)-BNC::llikDMH(para,z,rho,pi)-BNC::llikDMH(para,newz,newrho,newpi)
    if(is.nan(AccRate)){
      accrate=1
    }else{
      accrate=min(AccRate,0)
    }
    if(log(stats::runif(1))<accrate){
      ### accept this new state.
      rho=newrho
      pi=newpi
      z=newz
    }
    rhoTrack=cbind(rhoTrack,rho)
    piTrack=cbind(piTrack,pi)
    zTrack=cbind(zTrack,z)
  }
  if(plot){
    graphics::par(mfrow=c(2,2))
    for(rr in 1:nrow(rhoTrack)){
      plot(seq(1:ncol(rhoTrack)),rhoTrack[rr,],type="l",main=paste0("DMH for rho",rr),xlab="iterations",ylab=paste0("rho",rr))
    }
    #     dev.off()
    graphics::par(mfrow=c(1,3))
    for(rr in 1:nrow(piTrack)){
      plot(seq(1:ncol(piTrack)),piTrack[rr,],type="l",main=paste0("DMH for pi",rr),xlab="iterations",ylab=paste0("pi",rr))
    }
    #     dev.off()
  }
  pi=vector(length=pi.L)
  pi[-2]=rowMeans(matrix(piTrack[-2,round(0.5*para$hyperPara$niter):para$hyperPara$niter],nrow=(pi.L-1)))
  pi[2]=1-sum(pi[-2])
  rho=rowMeans(rhoTrack[,round(0.5*para$hyperPara$niter):para$hyperPara$niter])
  zsw=apply(zTrack,1,function(x){
    l=length(which(x==1))
    ll=length(which(x==2))
    lll=length(which(x==3))
    mm=which(ll==max(ll))
    return(mm)
  })
  res=list(rhoTrack=rhoTrack,piTrack=piTrack,zTrack=zTrack,pi=pi,rho=rho,zstat=zsw)
  return(res)
}

#' @title {z.update.Mclustprior}
#' @description z.update.Mclustprior
#' @param para The parameter list
#' @param oldz NULL
#' @param rho NULL
#' @param pi NULL
#' @return voted z
#' @keywords internal
#' @export
z.update.Mclustprior=function(para,oldz,rho,pi){
  zTrack=c()
  for(iter in 1:para$hyperPara$niter.z){
    newz=vector(length=length(oldz))
    for(node in 1:para$numNodes){
      nbr=which(para$net[node,]>0)
      loglike=sapply(1:length(para$zSet),function(z){
        log(pi[z])+rho[z]*sum(oldz[nbr]==z)+log(stats::dnorm(para$rStat[node],mean=para$hyperPara$mean.z[z],sd=para$hyperPara$sd.z[z]))
      })
      if(sum(is.infinite(loglike))>0){
        loglike[which(is.infinite(loglike))]=para$hyperPara$replaceInf
      }
      sampleprob=sapply(1:length(loglike),function(z){
        vv=loglike-loglike[z]
        return(1/sum(exp(vv)))
      })
      newz[node]=sample(para$zSet,1,prob=sampleprob)
      # print(node)
    }
    zTrack=cbind(zTrack,newz)
    oldz=newz
  }

  zvote=apply(zTrack,1,function(x){
    a=length(which(x==1))
    aa=length(which(x==2))
    aaa=length(which(x==3))
    mm=max(a,aa,aaa)
    return(ifelse(mm==a, 1, ifelse(mm==aa,2,3)))
  })
  return(zvote)
}

#' @title {llikDMH}
#' @description llikDMH
#' @param para NULL
#' @param z NULL
#' @param rho NULL
#' @param pi NULL
#' @return llik Value
#' @keywords internal
#' @export
llikDMH=function(para,z,rho,pi){
  res=sapply(1:para$numNodes,function(node){
    nbr=which(para$net[node,]>0)
    return(log(pi[z[node]])+rho[z[node]]*sum(z[nbr]==z[node])+stats::dnorm(para$rStat[node],mean=para$hyperPara$mean.z[z[node]],sd=para$hyperPara$sd.z[z[node]]))
  })
  return(sum(res))
}

#' @title {GCut.z}
#' @description Split network into three where nodes located within each group has the same class label (z value).
#' @param para The parameter list
#' @return The updated parameter list
#' @keywords internal
#' @export
GCut.z=function(para){
  GCut=list()
  for(z in para$zSet){
    ll=which(para$zLab==z)
    G=igraph::graph.adjacency(para$net[ll,ll],mode="undirected")
    G=igraph::set.vertex.attribute(G,"node",value=ll)
    GCut[[z]]=G
  }
  para$GCut=GCut
  return(para)
}


#' @title {SWCut}
#' @description Swendsen-Wang graph cut techniques. The edges located within each subgraph (same z) has a positive probability to stay "on" or "off". The edges with "off" labels are closed so that it is further cut to smaller subgraphs.
#' @param para parameter list
#' @param rhoScale a 4-element vector, used when rho values need to be scaled for a propor graph cut. Default=c(1,1,1,0)
#' @return  A list with three elements, representing z=1, z=2, or z=3. Each element is also a list, contains all subgraph with same z value.
#' @keywords internal
#' @export
SWCut=function(para){
  ProbEOn=1-exp(-para$rho[-length(para$rho)])
  gList=list()
  gl=1

  locvec=c()

  for(k in 1:length(para$GCut)){
    subg=para$GCut[[k]]
    #     plot(subg)
    E=igraph::E(subg)
    V=igraph::V(subg)
    vnode=igraph::get.vertex.attribute(subg,"node")
    if(length(vnode)==0){
      cat("Empty subgraph!","\n")
    }
    myz=para$zLab[vnode]
    if(length(E)>0){
      Eid=stats::rbinom(length(E),1,ProbEOn[myz[1]])
      newsubg=igraph::delete.edges(subg,(E[which(Eid==0)]))
      newsubg=igraph::clusters(newsubg)
      gList[gl:(gl+newsubg$no-1)]=lapply(1:newsubg$no,function(kk){
        igraph::induced.subgraph(graph=subg,vids=V[which(newsubg$membership==kk)])
      })

      gl=gl+newsubg$no
    }else{
      nnode=length(V)
      gList[gl:(gl+nnode-1)]=lapply(1:nnode,function(kk){
        mygraph=igraph::graph.empty(0, directed=FALSE)
        mygraph=igraph::add.vertices(mygraph,nv=1)
        mygraph=igraph::set.vertex.attribute(mygraph,"node",value=vnode[kk])
        return(mygraph)
      })
      gl=gl+nnode
    }
    locvec=c(locvec,vnode)
  }
  #   e=sapply(gList,function(x){length(E(x))})
  return(gList)
}

#' @title {z.update.fastSW}
#' @description Given density specification. Update regulation type (z value) within each subgraph cut by Swendsen-Wang graph cut techniques.
#' @param para The parameter list
#' @param swcutList a list, the direct output from SWCut.
#' @return An updated parameter list.
#' @keywords internal
#' @export
z.update.fastSW=function(para,swcutList){
  for(gph in 1:length(swcutList)){
    if((min(table(para$zLab))<para$min.node) | (length(unique(para$zLab))<para$numZ)){
      next;
    }
    # cat(gph,"\n")
    currt.gph=swcutList[[gph]]
    node.gph=igraph::get.vertex.attribute(currt.gph,"node")
    # plot(currt.gph)
    nnode=length(node.gph)
    oldz=para$zLab[node.gph][1]
    nbr=lapply(1:length(node.gph),function(kk){which(para$net[node.gph[kk],]>0)})
    prob.z=sapply(1:length(para$zSet),function(z){
      mk=sapply(nbr,function(nbrs){sum(para$zLab[nbrs]==z)})
      return(log(para$pivec[z])*length(node.gph)+(para$rho[z]*sum(mk))+BNC::logdensMixNorm(para$rStat[node.gph],mu=para$mulist[[z]],sd=para$sdlist[[z]],prop=para$proplist[[z]]))
    })
    if(sum(is.infinite(prob.z))>0){
      loc=which(is.infinite(prob.z))
      prob.z[loc]=1
    }
    sampleprob=sapply(1:length(prob.z),function(z){
      vv=prob.z-prob.z[z]
      return(1/sum(exp(vv)))
    })
    newz=sample(size=1,x=seq(1,length(prob.z)),prob=sampleprob)

    AccRate=exp(prob.z[newz]-prob.z[oldz])*para$TransZ[newz,oldz]/para$TransZ[oldz,newz]
    if(is.nan(AccRate)){
      accrate=0
    }else{
      accrate=min(AccRate,1)
    }

    if(stats::runif(1)<accrate){
      #       acclist[gph]=TRUE
      ### accept this new state.
      if((newz!=oldz)){
        if((para$mk[oldz]-nnode)>para$min.node){
          para$zLab[node.gph]=newz
          para$mk[c(oldz,newz)]=para$mk[c(oldz,newz)]+c(-1,+1)*nnode
        }
      }
    }
    #   return(list(para=para,ztrack=zz,problist=problist,acclist=acclist))
  }
  return(para)
}

#' @title {logdensMixNorm}
#' @description logdensMixNorm
#' @param dat NULL
#' @param mu NULL
#' @param sd NULL
#' @param prop NULL
#' @return log density value
#' @keywords internal
#' @export
logdensMixNorm=function(dat,mu,sd,prop){
  res=sapply(1:length(mu),function(l){
    prop[l]*stats::dnorm(dat,mean=mu[l],sd=sd[l])
  })
  res=as.matrix(res)
  if(ncol(res)!=length(mu)){
    res=t(res)
  }
  x=rowSums(res)
  loc=which(x==0)
  if(length(loc)>0){
    #     cat("inf log density:,",loc,"\n")
    #     x=x[-loc]
    return(-999999)
  }else{
    return(sum(log(x)))
  }
}

#' @title {DensityDPOne}
#' @description Given regulation type (z value). Update density specification based on DPM fitting.
#' @param para parameter list.
#' @return An updated parameter list.
#' @keywords internal
#' @export
DensityDPOne=function(para){
  for(z in para$zSet){
    ll=which(para$zLab==z)
    if(length(ll)>para$min.node){
      myprior=para$HODCPara$prior[[z]]
      # myprior$m1=para$mulist[[z]]
      # nn=length(para$sdlist[[z]])
      # myprior$psiinv1=diag(para$sdlist[[z]],nrow=nn,ncol=nn)
      res=BNC::DPdensitySubset(ll,subdat=para$rStat,subprior=myprior,submcmc=para$HODCPara$mcmc,substatus=TRUE)
      utils::flush.console()
      para$mulist[[z]]=res$mu
      para$sdlist[[z]]=sqrt(res$sigma)
      para$proplist[[z]]=res$prop
      para$gvec[ll]=res$gvec
      para$zgloc[[z]]=res$gloc
    }else{
      para$mulist[[z]]=mean(para$rStat[ll])
      para$sdlist[[z]]=sqrt(stats::var(para$rStat[ll]))
      para$proplist[[z]]=1
      para$gvec[ll]=rep(1,length(ll))
      para$zgloc[[z]]=list(ll)
    }
  }
  return(para)
}

#' @title {DPdensitySubset}
#' @description DPdensitySubset
#' @param ll NULL
#' @param subdat NULL
#' @param subprior NULL
#' @param submcmc NULL
#' @param substatus NULL
#' @return para
#' @keywords internal
#' @export
DPdensitySubset=function(ll,subdat,subprior,submcmc,substatus){
  #fit=DPpackage::DPdensity(subdat[ll],prior=subprior,mcmc=submcmc,status=substatus)
  dp <- dirichletprocess::DirichletProcessGaussian(subdat[ll],
                                                   g0Priors = c(subprior$m1,
                                                                subprior$k0,
                                                                subprior$nu1,
                                                                1/subprior$psiinv1))
  utils::flush.console()
  #nclust=fit$state$ncluster
  nclust = dp$numberClusters
  #prop=table(fit$state$ss)/sum(table(fit$state$ss))
  prop = table(dp$clusterLabels)/length(dp$clusterLabels)
  #mu=fit$state$muclus[1:nclust]
  mu = as.vector(dp$clusterParameters[[1]])
  #sigma=fit$state$sigmaclus[1:nclust]
  sigma = as.vector(dp$clusterParameters[[2]])
  index=order(mu,decreasing=FALSE)
  oldg=dp$clusterLabels
  newg=oldg
  k=1
  gloc=list()
  for(g in index){
    lll=which(oldg==g)
    newg[lll]=k
    gloc[[k]]=ll[lll]
    k=k+1
  }
  return(list(mu=mu[index],sigma=sqrt(sigma)[index],prop=prop[index],gvec=newg,gloc=gloc))
}

#' @title {r.knnImpute}
#' @description Impute missing node using their nearest neighbor test statistics.
#' @param para parameter list.
#' @return  An updated parameter list.
#' @keywords internal
#' @export
r.knnImpute=function(para){
  res=sapply(1:length(para$misLoc),function(k){
    nbr=para$nbrMis[[k]]
    rr=para$rStat[nbr]
    return(mean(rr))
  })
  para$rStat[para$misLoc]=res
  return(para)
}


#' @title {r.BayesImpute}
#' @description Impute missing node by sampling from posterior density.
#' @param para parameter list.
#' @return An updated parameter list.
#' @keywords internal
#' @export
r.BayesImpute=function(para){
  zloc=lapply(para$zgloc,function(x){unlist(x)})
  miszloc=sapply(para$misLoc,function(k){
    a=sapply(1:length(zloc),function(z){
      is.element(k,zloc[[z]])
    })
    which(a)
  })
  nn=length(para$misLoc)
  miszgloc=sapply(1:nn,function(k){
    tt=sapply(para$zgloc[[miszloc[k]]],function(x){is.element(para$misLoc[k],x)})
    return(which(tt))
  })
  miszgloc=cbind(miszloc,miszgloc)
  #   print(mizgloc[1,])
  mm=sapply(1:nn,function(k){c(para$mulist[[miszgloc[k,1]]][[miszgloc[k,2]]],para$sdlist[[miszgloc[k,1]]][[miszgloc[k,2]]])})
  rr=stats::rnorm(nn,mean=mm[1,],sd=mm[2,])
  #   print(rr[1])
  para$rStat[para$misLoc]=rr
  #   print(para$rStat[773])
  return(para)
}



