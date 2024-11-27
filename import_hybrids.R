### Packages 
library(tidyverse)
library(data.table)
library(ggplot2)

# Pheno_data_2014_21 <- fread(file = "/home/uni03/UANZ/osatohanmwen/G_2_F_competition/G_2_F-Final/1_Training_Trait_Data_2014_2021.csv")
# 
# Pheno_data_2014_21 <- as_tibble(Pheno_data_2014_21)
# 
# 
# Pheno_data_2018_2021 <- filter(Pheno_data_2014_21, Year %in% c("2018", "2019", "2020", "2021")) 
# 
# Hybrids_2018_2021 <- c(unique(Pheno_data_2018_2021$Hybrid))
# 
# write.table(Hybrids_2018_2021,file = "import_hybrid.txt",col.names = T,row.names = F,sep = '\t',quote = F)
# 
# ##used the hybrid names to subset in the SNP filtering script
# 
# #Examining statistics (before filtering) from the SNP filtering script
# 
# ####Variant based statistics####
# 
# ## 1- Variant quality
# 
# var_qual <- read_delim("vcfresults/Geno_18_21_fil.lqual", delim = "\t", col_names = c("chr", "pos", "qual"), skip = 1)
# summary(var_qual$qual)
# 
# a <- ggplot(var_qual, aes(qual)) +
#   geom_density(fill="dodgerblue1")+
#   theme_light()
# 
# #### 2- Variant mean depth
# 
# var_depth <- read_delim("vcfresults/Geno_18_21_fil.ldepth.mean", delim = "\t", 
#                         col_names = c("chr", "pos", "mean_depth","var_dept"), skip = 1)
# 
# a <- ggplot(var_depth, aes(mean_depth)) +
#   geom_density(fill="dodgerblue1")+
#   theme_light()+
#   xlim(0,100)
# 
# #### 3- Variant missingnes
# var_miss<- read_delim("vcfresults/Geno_18_21_fil.lmiss", delim = "\t", 
#                       col_names = c("chr", "pos", "nchr","nfiltered", "nmiss", "fmiss"), skip = 1)
# 
# a <- ggplot(var_miss, aes(fmiss)) +
#   geom_density(fill="dodgerblue1")+
#   xlim(0.0, 0.7)+
#   theme_light(base_size = 30)
# 
# summary(var_miss$fmiss)
# 
# #### 4- Minor allele frequency
# 
# var_freq <- read_delim("vcfresults/Geno_18_21_fil.frq", delim = "\t", 
#                        col_names = c("chr", "pos", "nalleles","nchr", "a1", "a2"), skip = 1)
# ###to get the minor allele frequency
# 
# var_freq$maf <- var_freq %>%
#   dplyr::select(a1,a2)%>%
#   apply(1, function(z) min(z))
# 
# a <- ggplot(var_freq, aes(maf)) +
#   geom_density(fill="dodgerblue1")+
#   theme_light()
# 
# summary(var_freq$maf)
# 
# 
# ####Individual based statistics####
# 
# #### 1- Variant mean depth
# 
# ind_depth <- read_delim("vcfresults/Geno_18_21_fil.idepth", delim = "\t", 
#                         col_names = c("ind", "depth"), skip = 1)
# 
# a <- ggplot(ind_depth, aes(depth)) +
#   geom_histogram(fill="dodgerblue1")+
#   theme_light()+
#   xlim(0,100)
# 
# ### 2 - Missingness per individual
# 
# ind_miss<- read_delim("vcfresults/Geno_18_21_fil.imiss", delim = "\t", 
#                       col_names = c("ind", "ndata","nfiltered", "nmiss", "fmiss"), skip = 1)
# 
# a <- ggplot(ind_miss, aes(fmiss)) +
#   geom_histogram(fill="dodgerblue1")+
#   theme_light()+
#   xlim (0,0.3)
# 
# summary(ind_miss$fmiss)
# 
# 
# ### 3 - Heterozygosity and inbreeding coefficient
# 
# ind_het <- read_delim("vcfresults/Geno_18_21_fil.het", delim = "\t", 
#                       col_names = c("ind", "ho","he", "nsites", "f"), skip = 1)
# 
# a <- ggplot(ind_het, aes(f)) +
#   geom_histogram(fill="dodgerblue1")+
#   theme_light()+
#   xlim (0,0.3)
# 
# summary(ind_het$f)
# 
# ##imputation of missing data



##load in filtered genotype data and recode 

SNP_data_subset <- tibble(fread(file = "Geno_18_21_ld.raw"))

SNP_data_subset_select <- SNP_data_subset %>%
  dplyr::select(colnames(SNP_data_subset),-c("FID","PAT","MAT","PHENOTYPE","SEX"))

#recode
# for (i in 2:ncol(SNP_data_subset_select)){
#   SNP_data_subset_select[,i]<- ifelse(SNP_data_subset_select[,i]== "0/0",0,
#                                       ifelse(SNP_data_subset_select[,i]== "0/1", 1,
#                                              ifelse(SNP_data_subset_select[,i]=="1/0",1,
#                                                     ifelse(SNP_data_subset_select[,i]=="1/1", 2, "NA"))))
# }
# 
# rownames(SNP_data_subset_select) <- SNP_data_subset_select$ID
# 
# SNP_data_subset_select<-as.data.frame(t((SNP_data_subset_select)))[-1,]
# 
# SNP_data_subset_select$Genotype <- rownames(SNP_data_subset_select)
# SNP_data_subset_select <- data_frame(SNP_data_subset_select) %>%
#   dplyr::select(Genotype, everything())

##Missing marker imputation, replacing missing marker with the mean values
SNP_data_subset_select_1 <- tibble(SNP_data_subset_select[,-1]) %>%
  mutate_if(is.character, as.numeric)

mean_val <- round(colMeans(SNP_data_subset_select_1, na.rm=TRUE))

for (i in colnames(SNP_data_subset_select_1)){
  SNP_data_subset_select_1[,i][is.na(SNP_data_subset_select_1[,i])] <- mean_val[i]
}

SNP_data_subset_select_1 <- data.frame(Hybrid =SNP_data_subset_select$IID, SNP_data_subset_select_1)

write.table(SNP_data_subset_select_1,file = 'Genotype_True2.csv',col.names = T,row.names = F,sep = '\t',quote = F)


