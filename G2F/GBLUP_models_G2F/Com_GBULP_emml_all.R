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
##Genomic prediction with Addictive and dominance 

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

#5 folds cross validation 10 times ####

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


####Genomic prediction with the combined model

All_TEN <- list()

for (time in 1:10){
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
    
    Train   <- which(folds != i)
    Test    <- which(folds == i)
    y_train <- y[Train]
    y_test  <- y[Test]
    
    ## Genotype transformation
    
    #Weights on heterozygous genotypes 
    SNP_data_G <- as.data.frame(SNP_data[,-1])
    
    genotypeAD_train <- SNP_data_G[Train, ]
    
    
    ## calculating the corrected snp mean
    snpMean=array(0,dim=c(ncol(genotypeAD_train),3))
    
    for (snp in 1:ncol(genotypeAD_train)) {
      
      snpMean[snp,1]=mean(y_train[genotypeAD_train[,snp]==0], is.na=T)
      snpMean[snp,2]=mean(y_train[genotypeAD_train[,snp]==1], is.na=T)
      snpMean[snp,3]=mean(y_train[genotypeAD_train[,snp]==2], is.na=T)
      
    }
    
    snpMean[is.na(snpMean)] <- 0
    
    ##if mean(A1A1)>mean(A2A2), A1A1 and A2A2 are recoded as 2 and 0
    
    snpMeanOri=snpMean
    
    mychoose=snpMeanOri[,1]>snpMeanOri[,3]
    
    snpMean[mychoose,1]=snpMeanOri[mychoose,3]
    
    snpMean[mychoose,3]=snpMeanOri[mychoose,1]
    
    genotypeAD_train[,mychoose]=abs(genotypeAD_train[,mychoose]-2)
    
    
    ####Calculate d, which is the degree of dominance
    
    d=(snpMean[,2]-snpMean[,1])/(snpMean[,3]-snpMean[,1])*2
    
    mychooseNa=is.na(d)
    d[mychooseNa]=1
    
    ##The boundaries of 0 and 2 for d, the definition of the weighted heterogeneous genotype.
    mychoose1=d>=2
    d[mychoose1]=2
    
    mychoose2=d<=0
    d[mychoose2]= 0
    
    d
    
    ##Weighted heterozygous genotypes 
    for (snp in 1:ncol(SNP_data_G)) {
      mychoose=SNP_data_G[,snp]==1
      SNP_data_G[mychoose,snp]= d[snp]
    }
    
    ### Van Raden(2008) Genomic relationship matrix (additive and dominance combined)
    
    alelleFreq <- function(x, y) {
      (2 * length(which(x == y)) + length(which(x == 1)))/(2 * 
                                                             length(which(!is.na(x))))
    }
    
    
    Frequency <- cbind(apply(SNP_data_G, 2, function(x) alelleFreq(x,0)), 
                       apply(SNP_data_G, 2, function(x) alelleFreq(x, 2)),
                       colMeans(SNP_data_G)/2)
    
    
    FreqP <- matrix(rep(Frequency[, 2], each = nrow(SNP_data_G)), 
                    ncol = ncol(SNP_data_G))
    
    MyTwoPQ <- sum((2 *Frequency[, 1]) * Frequency[, 2])
    
    SNP_data_G <- SNP_data_G - 2 *FreqP
    
    
    SNP_data_G[is.na(SNP_data_G)] <- 0
    
    SNP_data_G <- as.matrix(SNP_data_G)
    
    Genomic_com_matrix <- (tcrossprod(SNP_data_G, SNP_data_G)) /as.numeric(MyTwoPQ)
    
    
    ###Prediction usiong emmreml
    
    funout<-emmreml(y=y_train, X=matrix(rep(1, n)[Train], ncol=1), Z=Z[Train,], K=Genomic_com_matrix)
    
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

RESULT_com <- ldply(All_TEN,rbind)

write.table(RESULT_com,file = paste0(trait,"_Com_Result", ".csv"),col.names = T,row.names = F,sep = '\t',quote = F)
