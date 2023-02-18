# Bayesian Node Classification (BNC)

## **Overview**

An R package conducts the nonparametric Bayesian node classification (BNC) with missing values on biological networks. The posterior computation algorithm is designed utilizing Swendsen-Wang algorithm for efficiently updating the gene class indicators as well as an updated Dirichlet process mixture model as a prior model. The algorithm can naturally handle missing values in the Bayesian modeling framework and provide two algorithms for imputation. Functions used in simulation studies are also provided for illustration purpose.

## **Installation**

* **Please make sure the following R Packages installed** 
```{r}
for(Rpackage in c("devtools","flexmix",
"GGally","igraph","mclust","truncnorm","dirichletprocess")){

  if(!is.element(Rpackage,installed.packages())){
   
    install.packages(Rpackage)
    
  }
  
}


```

* **Install and load the package**
```{r}
devtools::install_github("kangjian2016/BNC")
library(BNC)
```

## **Example**
```{r}
simdata=simulatedDatasetGenerator(nnode=100,missing=TRUE,missrate=0.1,dist="norm",plot=TRUE,nbin=c(20,20,10),rng=1024)
res=BNC(net=simdata$net,test.stat=simdata$testcov,niter=100,na.action="NN",
paras=list(tau=c(2,10,8),alpha=NULL,gamma=NULL,xi=NULL,beta=rep(10,3),
rho=c(1.003,0.479,0.988,0.000),pivec=c(0.25,0.40,0.35),densAcc=0.001,
null.quantile=c(0.25, 0.75),null.method="biGaussianModeSplit",transitionMatrix.Z.11=0.5,
miss.stat=2,min.node=5))
classLabelEst=summaryClassLabels(simdata$net,simdata$testcov,res$zTrack,
method="MajorVote",nburn=50)
print(table(classLabelEst))
```

