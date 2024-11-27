# Predicting Maize Hybrid Performance with Machine Learning and a Locus-Specific Degree of Dominance Transformation

The genetic architecture of a trait plays a vital role in the predictive ability of genomic models. While classical models like genomic best linear unbiased prediction (GBLUP) dominate plant breeding, machine learning (ML) methods are gaining traction for their superior handling of non-linear effects.

This study assessed two models integrating ML and classical statistical methods, incorporating a novel locus-specific weighted dominance effect transformation matrix for genomic prediction in hybrid maize.

## A total of five models were compared: 
1. XGBoost combined with locus-specific weighted dominance effects, 
2. XGBoost only, 
3. GBLUP model with additive effects only, 
4. GBLUP model, including additive and dominance effects, 
5. GBLUP combined with locus-specific weighted dominance effects.

The models were evaluated using two simulated maize populations (one polygenic scenario and the other, an oligogenic scenario; in both scenarios, the dominance varied from 0% to 40%) and a real maize hybrid population (G2F data (2018 - 2021)). Below is a step-by-step process for all the analyses done. **Please follow the steps below to reproduce the results from our work**.

We start with the Genotype analysis.
#steps to filter the VCF file for use for my degree of dominance analysis

1. download the genomic file from the g2f site

```shell
wget https://de.cyverse.org/anon-files//iplant/home/shared/commons_repo/curated/GenomesToFields_GenotypeByEnvironment_PredictionCompetition_2023/Training_data/5_Genotype_Data_All_Years.vcf.zip
```
2. Unzip the VCF file

unzip 5_Genotype_Data_All_Years.vcf.zip 

3. filter out the hybrid from year 2018-2021

bcftools view -S import_hybrid.txt 5_Genotype_Data_All_Years.vcf > Genotype_data_18_21.vcf --force-samples

4. compress and index the vcf to make it easier to access

bgzip Genotype_data_18_21.vcf

bcftools index Genotype_data_18_21.vcf.gz

5. Generating important statistics to make it easier to filter.

#Calculate allele frequency
vcftools --gzvcf Genotype_data_18_21.vcf.gz --freq2 --out Geno_18_21_fil --max-alleles 2
#Calculate mean depth per individual
vcftools --gzvcf Genotype_data_18_21.vcf.gz --depth --out Geno_18_21_fil
#Calculate mean depth per sites
vcftools --gzvcf Genotype_data_18_21.vcf.gz --site-mean-depth --out Geno_18_21_fil
#Calculate site quality
vcftools --gzvcf Genotype_data_18_21.vcf.gz --site-quality  --out Geno_18_21_fil
#Calculate proportion of missing data per individual
vcftools --gzvcf Genotype_data_18_21.vcf.gz --missing-indv --out Geno_18_21_fil
#Calculate proportion of missing data per site
vcftools --gzvcf Genotype_data_18_21.vcf.gz --missing-site --out Geno_18_21_fil
#Calculate heterozygosity and inbreeding coefficient per individual
vcftools --gzvcf Genotype_data_18_21.vcf.gz --het --out Geno_18_21_fil

#go to "import_hybrid script" to examine statistics

6. remove indels and max alleles of more than 2

vcftools --gzvcf Genotype_data_18_21.vcf.gz --remove-indels --max-alleles 2 --recode --stdout | gzip -c > Geno_18_21_filtered.vcf.gz

7. After examining the Stats, filter using ;

MAF =0.05
Variant missingness = 0.9

vcftools --gzvcf Geno_18_21_filtered.vcf.gz --maf 0.05 --max-missing 0.8 --remove-indels --max-alleles 2 --recode --stdout | gzip -c > Geno_18_21_ffinal.vcf.gz

8. filter with LD 

#convert the vcf to a plink file
vcftools --gzvcf Geno_18_21_ffinal.vcf.gz --plink --out ./plink/myplink
#chech the snps R^2 values
plink --file ./plink/myplink --r2 --out ./plink/ld_output
#prune out snps with R^2 values greater than 0.99
plink --noweb --file ./plink/myplink --indep-pairwise 100 5 0.99 --out ./plink/plink
#subset the original SNP matrix with the prune snp
plink --file ./plink/myplink --extract ./plink/plink.prune.in --make-bed --out ./plink/pruneddata
#convert back to VCF format.
plink --bfile ./plink/pruneddata --recode vcf --out Geno_18_21_ld
#convert to csv
plink --bfile ./plink/pruneddata --recode A --out Geno_18_21_ld

Then, the Phenotype analysis.
GBLUP combined with locus-specific weighted dominance effects.
XGBoost combined with locus-specific weighted dominance effects,

The assessment of the real maize population focused on their performance across several traits, including test weight (kg/m³), ear height (cm), plant height (cm), days after pollination for pollen (Pollen DAP), days after pollination for silk (Silk DAP), and grain moisture. This evaluation aimed to encompass a range of diverse genetic architectures. The results show that the transformation with dominance effects did not improve ML methods. While ML models performed better than the GBLUP model with additive effects only and the GBLUP model including additive and dominance effects, the GBLUP Transformed model outperformed all others, including ML models across the three levels of heritability in the polygenic scenarios, demonstrating the superiority of genomic prediction with transformed GBLUP for polygenic traits. This cannot be said for oligogenic scenarios, as ML models outperformed all other models in oligogenic scenarios. Our findings contribute to optimizing breeding strategies and advancing genomic prediction methodologies.



