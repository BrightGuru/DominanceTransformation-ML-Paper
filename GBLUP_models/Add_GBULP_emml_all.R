library(tidyverse)
library(data.table)
library(foreach)
library(doParallel)
library(doFuture)
library(doMC)
library(plyr)
library(dplyr)
library(AGHmatrix)
library(EMMREML)


####Genomic prediction with Addictive model ####

### load in blues of the phenotype data (already calculated on line mean basis)

trait <- commandArgs(trailingOnly = T)[1]

print(trait)

file <- paste0("/usr/users/osatohanmwen/GP_Dom_G2Fdata/emmreml/",trait ,"_BLUES_All_final",".csv")

Pheno_Blues_2018_2021 <- fread(file)

Pheno_Blues_2018_2021 <- Pheno_Blues_2018_2021 %>%
  mutate(Hybrid=gsub("Hybrid", "",Hybrid))

print(summary(unique(Pheno_Blues_2018_2021$Hybrid)))

##load in filtered and recoded genotype data and sort with the hybrid from phenotype data

SNP_data <- data.frame(fread(file = "/usr/users/osatohanmwen/GP_Dom_G2Fdata/emmreml/Genotype_True1.csv"))%>%
  filter(Hybrid %in% Pheno_Blues_2018_2021$Hybrid)

Pheno_Blues_2018_2021 <- Pheno_Blues_2018_2021 %>%
  filter(Hybrid %in% SNP_data$Hybrid)


SNP_data_G <- as.data.frame(SNP_data[,-1])%>%
  mutate_if(is.character, as.numeric)

## Van Raden(2008) Genomic relationship matrix (additive)##
Genomic_add_matrix <- Gmatrix(SNPmatrix=as.matrix(SNP_data_G),
                              maf=0.05, method="VanRaden")


#5 folds cross validation 20 times ####

#1. MEAN ABSOLUTE PERCENTAGE ERROR (MAPE)
MAPE = function(y_actual,y_predict){
  mean(abs((y_actual-y_predict)/y_actual))*100
}

#2. R SQUARED error metric -- Coefficient of Determination
RSQUARE = function(y_actual,y_predict){
  cor(y_actual,y_predict)^2
}

#3 ROOT MEAN SQUARE ERROR
RMSE = function(y_actual,y_predict){
  sqrt(mean((y_actual - y_predict)^2))
}


#number of codes
n.cpus <- Sys.getenv("SLURM_CPUS_PER_TASK")

n.cpus

class(n.cpus)

# we need this to be numeric below so:

n.cpus <- as.numeric(n.cpus)

n.cpus

class(n.cpus)

# register a parallel backend specifying the number of CPUs as the number we imported via Sys.getenv()

registerDoMC(cores = n.cpus) 


###Genomic prediction with Addictive effects only.

All_TEN <- list()

for (time in 1:20){
  seed  <- 1234
  nfolds <- 5
  set.seed(seed + time)
  
  y <- as.matrix(Pheno_Blues_2018_2021[[trait]])
  
  n=length(y) #Computing the length of the vector containing the phenotypic data
  
  Z=diag(n)
  
  #x=matrix(1,n,1) ##Computing the vector of fixed effects
  
  folds <- sample(rep(1:nfolds, each = n %/% nfolds))
  
  cross_predict_add <- foreach (
    i=1:nfolds, .multicombine=T,.combine="rbind"
  ) %dopar% {
    
    Train <- which(folds != i)
    Test  <- which(folds == i)
    
    funout<-emmreml(y=y[Train], X=matrix(rep(1, n)[Train], ncol=1), Z=Z[Train,], K=Genomic_add_matrix)
    
    va   <- funout$Vu
    vd   <- 0
    beta <- funout$betahat
    ve   <- funout$Ve
    
    
    yobs <- y[Test]
    x=matrix(1,length(yobs),1)

    yhat <- as.matrix(funout$uhat[Test]) + (x%*%beta) 
    
    list(Fold = i,
         Repeat = time,
         Beta =beta,
         Va = va,
         Ve =ve,
         Vd =vd,
         COR=cor(yobs,yhat),
         Rsquare=RSQUARE(yobs,yhat),
         Rmse=RMSE(yobs,yhat),
         Mape=MAPE(yobs,yhat)
         )
  }
  
  RESULT_DF  <- data.frame(cross_predict_add)
  RESULT_DF <- apply(RESULT_DF,2,as.character)
  
  All_TEN[[time]] <- cbind(RESULT_DF)
  
}

RESULT_add <- ldply(All_TEN,rbind)

write.table(RESULT_add,file = paste0(trait,"_Add_Result", ".csv"),col.names = T,row.names = F,sep = '\t',quote = F)