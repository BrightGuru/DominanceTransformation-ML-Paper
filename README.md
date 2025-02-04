# Predicting Maize Hybrid Performance with Machine Learning and a Locus-Specific Degree of Dominance Transformation

The genetic architecture of a trait plays a vital role in the predictive ability of genomic models. While classical models like genomic best linear unbiased prediction (GBLUP) dominate plant breeding, machine learning (ML) methods are gaining traction for their superior handling of non-linear effects.

This study assessed two models integrating ML and classical statistical methods, incorporating a novel locus-specific weighted dominance effect transformation matrix for genomic prediction in hybrid maize.

## A total of five models were compared: 
1. XGBoost combined with locus-specific weighted dominance effects, 
2. XGBoost only, 
3. GBLUP model with additive effects only, 
4. GBLUP model, including additive and dominance effects, 
5. GBLUP combined with locus-specific weighted dominance effects.

The models were evaluated using two simulated maize populations (one polygenic scenario and the other, an oligogenic scenario; in both scenarios, the dominance varied from 0% to 40%) and a real maize hybrid population (G2F data (2018 - 2021)). Below is a step-by-step process for all the analyses done. 

## Please follow the steps below to reproduce the results from our work with the G2F data

**Genotypic analysis**
1. download the genomic file from the g2f site
``` shell
wget https://de.cyverse.org/anon-files//iplant/home/shared/commons_repo/curated/GenomesToFields_GenotypeByEnvironment_PredictionCompetition_2023/Training_data/5_Genotype_Data_All_Years.vcf.zip
```
2. Unzip the VCF file
3. filter out the hybrid from year 2018-2021
4. compress and index the VCF to make it easier to access
5. Generating important statistics to make it easier to filter.
6. remove indels and max alleles of more than 2
7. After examining the Stats, filter using **MAF =0.05**, **Variant missingness = 0.9**.
8. filter with LD
Note. For bash codes, see *SNP_filtering.sh*

**Phenotypic analysis**

The phenotypic analysis includes outlier analysis, best linear unbiased estimation (BLUES) with replicate as fixed effects, BLUES with Environment(Field location and Year) and hybrid as random effects, and Heritability and Variance. The codes for these analyses can be found in *phenotypic_analyses_G2F_Phenotype.R*.

**GBLUP models**

To run GBLUP models, use the scripts in the **GBLUP_models** directory.  
- Use *Add_GBULP_emml_all.R* for GBLUP with additive effects only use.  
- Use *Add_dom_GBULP_emml_all.R* for the GBLUP model, including additive and dominance effects. 
- Use *Com_GBULP_emml_all.R* for GBLUP combined with locus-specific weighted dominance effects.

**Machine learning Models**

To run ML models, use the scripts in the **ML_models** directory.
- Use *XGBoost_model_script.py* for XGBoost only.
- Use *XGBoostdom_model_script.py* for XGBoost combined with locus-specific weighted dominance effects.

## To reproduce the results from our work with the simulated data, Please follow the steps below

**Simulating the Polygenic and oligogenic scenarios**
To simulate the Polygenic scenarios use *HybridSim_Polyscript.R* and *HybridSim_Oligoscript.R* for the Oligogenic scenarios.

**GBLUP models**

To run GBLUP models, use the scripts in the **Simulated_GBLUP_Model** directory.  
- Use *Sim_add_GBLUP.R* for GBLUP with additive effects only use.  
- Use *Sim_adddom_GBLUP.R* for the GBLUP model, including additive and dominance effects. 
- Use *Sim_com_GBLUP.R* for GBLUP combined with locus-specific weighted dominance effects.

**Machine learning Models**

To run ML models, use the scripts in the **Simulated_ML_Model** directory.
- Use *sim_worker_XGBoost.py* for XGBoost only.
- Use *sim_worker_XGBoostDom.py* for XGBoost combined with locus-specific weighted dominance effects.

**Note:** To get similar ML results for simulated and real data, ensure a similar conda environment like mine is replicated using the ML environment file *ML_environment.yml*.
