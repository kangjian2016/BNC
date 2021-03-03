# Bayesian Node Classification (BNC)

**An R package conduct the nonparametric Bayesian node classification (BNC) with missing values on biological networks. The posterior computation algorithm is designed utilizing Swendsen-Wang algorithm for efficiently updating the gene class indicators as well as an updated Dirichlet process mixture model as a prior model. The algorithm can naturally handle missing values in the Bayesian modeling framework and provide two algorithms for imputation. Functions used in simulation studies are also provided for illustration purpose.**

* **Please make sure the following R Packages installed** 
```{r}
for(Rpackage in c("devtools","flexmix","GGally","igraph","mclust","truncnorm")){

  if(!is.element(Rpackage,installed.packages())){
   
    install.packages(Rpackage)
    
  }
  
}

if(!is.element("DPpackage",installed.packages()){    

   devtools::install_github("konkam/DPpackage")
   
}
```

* **Install and load the package**
```{r}
devtools::install_github("kangjian2016/BNC")
library(BNC)
```

