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

**Genotype analysis**
1. download the genomic file from the g2f site
``` shell
wget https://de.cyverse.org/anon-files//iplant/home/shared/commons_repo/curated/GenomesToFields_GenotypeByEnvironment_PredictionCompetition_2023/Training_data/5_Genotype_Data_All_Years.vcf.zip
```
2. Unzip the VCF file
3. filter out the hybrid from year 2018-2021
4. compress and index the vcf to make it easier to access
5. Generating important statistics to make it easier to filter.
6. remove indels and max alleles of more than 2
7. After examining the Stats, filter using; MAF =0.05, Variant missingness = 0.9
8. filter with LD 
**Phenotype analysis**



