library(tidyverse)
library(AlphaSimR)
library(gtools)
# Remove all objects from the environment
rm(list = ls())

#Simulate founder genome
set.seed(100)
FounderGenomeHybrid =runMacs(nInd = 1000,
                             nChr=10,
                             inbred = TRUE,
                             segSites =  3000,
                             species = "MAIZE")

#create the object holding simulation parameters
SP =SimParam$new(FounderGenomeHybrid)

SP$setSexes("yes_rand")

# Add Trait 
# set 1
set.seed(100)
names_vector <- paste0("Oligo", 1:30)

for (name in names_vector) {
  
  SP$altAddTraitAD(
    nQtlPerChr = 3,
    mean = 20,
    varA = 0.6,
    varD = 1,
    inbrDepr = 2,
    limMeanDD = c(0, 1.5),
    limVarDD = c(0, 0.5),
    silent = FALSE,
    force = FALSE,
    name = name
  )
}
#set 2
set.seed(100)
names_vector <- paste0("Oligo",31:60)

for (name in names_vector) {
  SP$altAddTraitAD(
    nQtlPerChr = 3,
    mean = 20,
    varA = 0.6,
    varD = 0.4,
    inbrDepr = 2,
    limMeanDD = c(0.5, 1),
    limVarDD = c(0.5, 1),
    silent = FALSE,
    force = FALSE,
    name = name
  )
}

#set 3
set.seed(100)
names_vector <- paste0("Oligo", 61:90)

for (name in names_vector) {
  SP$altAddTraitAD(
    nQtlPerChr = 3,
    mean = 20,
    varA = 0.6,
    varD = 0.4,
    inbrDepr = 2,
    limMeanDD = c(0, 1),
    limVarDD = c(0, 1),
    silent = FALSE,
    force = FALSE,
    name = name
  )
}

#set 4
set.seed(100)
names_vector <- paste0("Oligo", 91:96)

for (name in names_vector) {
  SP$altAddTraitAD(
    nQtlPerChr = 3,
    mean = 20,
    varA = 0.6,
    varD = 0,
    inbrDepr = 0,
    limMeanDD = c(0, 1),
    limVarDD = c(0, 1),
    silent = FALSE,
    force = FALSE,
    name = name
  )
}
#Number of traits

SP$nTraits

#Create a population of individuals
set.seed(100)
basePop_female=newPop(FounderGenomeHybrid)

basePop_male=newPop(FounderGenomeHybrid)

#Make crosses for full diallele

female <- basePop_female[basePop_female@id[basePop_female@sex=="F" ][1:300]]
male   <- basePop_male[basePop_male@id[basePop_male@sex=="M" ][1:8]]

#?hybridCross
set.seed(100)
pop2 <- hybridCross(
  female,
  male,
  crossPlan = "testcross",
  returnHybridPop = FALSE,
  simParam = SP
)

# Error variance
Pro_dom = data.frame(Trait=rownames(varG(pop2)), Pro_dom=diag(varD(pop2,simParam=SP)/varG(pop2)))

Pro_dom <- Pro_dom[order(Pro_dom$Pro_dom, decreasing = F), ]
Pro_dom$GroupLabel <- rep(1:3, length.out = nrow(Pro_dom))

VarianceE_0.8 <- ((diag(varG(pop2))[c(Pro_dom[Pro_dom$GroupLabel == 1, ]$Trait)])*(1-0.8)) / 0.8
VarianceE_0.6 <- ((diag(varG(pop2))[c(Pro_dom[Pro_dom$GroupLabel == 2, ]$Trait)])*(1-0.6)) / 0.6
VarianceE_0.3 <- ((diag(varG(pop2))[c(Pro_dom[Pro_dom$GroupLabel == 3, ]$Trait)])*(1-0.3)) / 0.3

VarianceE <- c(VarianceE_0.8,VarianceE_0.6,VarianceE_0.3)
traits <- names(VarianceE)

traits <- traits[order(as.numeric(sub("Oligo", "", traits)))]

VarianceE <- VarianceE[traits]

##Phenotype Values
basePop_Pheno=setPheno(pop2,
                       varE=VarianceE)


##Extract the genome of all individual
Geno_sim = data.frame(Genotype=basePop_Pheno@id, pullSegSiteGeno(basePop_Pheno))
Geno_sim$Genotype = paste("Geno_",Geno_sim$Genotype, sep="")


##Extract the phenome of all individual
Peno_sim = data.frame(Genotype=basePop_Pheno@id,pheno(basePop_Pheno))
Peno_sim$Genotype = paste("Geno_",Peno_sim$Genotype, sep="")

##dominance variance of all trait
varianceA = varA(basePop_Pheno,simParam=SP)
varianceD = varD(basePop_Pheno,simParam=SP)
varianceP = varP(basePop_Pheno)
varianceG = varG(basePop_Pheno)

Her_dominance = diag(varD(basePop_Pheno,simParam=SP)/varP(basePop_Pheno))
Her_narow = diag(varA(basePop_Pheno,simParam=SP)/varP(basePop_Pheno))
Her_broad = diag(varG(basePop_Pheno)/varP(basePop_Pheno))

Pro_dominance = diag(varD(basePop_Pheno,simParam=SP)/varG(basePop_Pheno))
Pro_additive = diag(varA(basePop_Pheno,simParam=SP)/varG(basePop_Pheno))

Her <- data.frame(Trait=rownames(varianceA), varianceA =diag(varianceA), varianceD =diag(varianceD),varianceP =diag(varianceP),
                  varianceG =diag(varianceG), Her_dominance=Her_dominance,Her_narow= Her_narow,Her_broad=Her_broad,
                  Pro_dominance = Pro_dominance,Pro_additive = Pro_additive)%>%
  merge(.,Pro_dom[,-2])


Geno_pheno<- merge(Geno_sim,Peno_sim)

write.table(Geno_pheno,file = "Geno_pheno_hybrid.csv",col.names = T,row.names = F,sep = '\t',quote = F)


write.table(Geno_sim,file = "Geno_sim_hybrid.csv",col.names = T,row.names = F,sep = '\t',quote = F)
write.table(Peno_sim,file = "Peno_sim_hybrid.csv",col.names = T,row.names = F,sep = '\t',quote = F)
write.table(Her,file = "Her_hybrid.csv",col.names = T,row.names = F,sep = '\t',quote = F)
