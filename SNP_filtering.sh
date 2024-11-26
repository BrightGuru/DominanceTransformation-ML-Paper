#steps to filter the VCF file for use for my degree of dominance analysis

1. download the genomic file from the g2f site

wget https://de.cyverse.org/anon-files//iplant/home/shared/commons_repo/curated/GenomesToFields_GenotypeByEnvironment_PredictionCompetition_2023/Training_data/5_Genotype_Data_All_Years.vcf.zip

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