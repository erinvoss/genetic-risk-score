# Polygenic Risk Score (PRS) Calculation

- **Date Written:** 22.04.2019
- **Date Last Updated:** 14.04.2026

## Description
Polygenic risk scores (AKA: genetic risk scores, polygenic  scores, or genome-wide scores) is a summary measure of a set of risk-associated genetic variants and can be easily calculated using PLINK.

## Requirements 
1. **PLINK** (v2.0 is used here)
- The most common and easiest way to do this in PLINK [v2.0](https://www.cog-genomics.org/plink/2.0/score)
- GP2 has developed an end-to-end PRS calculation pipeline using Nextflow, which is described [here](https://github.com/hirotaka-i/gp2-gwas-variants-cumulative-burden). 
- Other tools, including [PRSice](https://choishingwan.github.io/PRSice/) and [LDpred2](https://privefl.github.io/bigsnpr/articles/LDpred2.html), are available as well.

2. **PLINK Binary Files** (`.bed`, `.bim`, `.fam`) 
3. **Score File** 
- This is a file with a variant identifier, allele, and an associated score value
- These scores come from GWA Studies of your disease of interest.
  	- Parkinson's Disease score file: [insert here] 
  	- Alzheimer's Disease score file: [insert here] 
  	- Additional diseases: [insert here] 
- Important: your plink files must use the same  reference genome build used to generate the score file (GRCh38). 
- The format of the file should match this example: 
```
1:154898185	C	0.2812
1:155135036	A	0.6068
1:155205634	T	-0.7467
1:161469054	C	0.065
```



## PLINK Commands

```bash 
module load plink2 #if on Biowulf
plink2 --bfile $yourfile --score $scorefile --out $outputfilename
```
Where: 
- `$yourfile` = standard binary file prefix (will point to `.bed`, `.bim`, and `.fam` files)
- `$outputfilename` = whatever you want it to be, the output will have the extension `.sscore` (or `.profile` if using plink v1.9)
- `$scorefile` = file with variant-name, allele and score-value

## Optional Follow-Up Analyses (Using Case-Control Data)

The code below provides R and Python steps to compare PRS scores between case and control individuals in your PLINK dataset. 
First, individual PRS scores are converted to Z-scores before generating a boxplot to visually compare PRS between groups. 
Lastly, we apply a GLM to the data to evaluate whether PRS differs significantly between cases and controls after accounting for covariates includge age, sex, and principal components. 

Option 1: R
```bash, R
module load R #if on Biowulf

# Load in the necessary packages 
library(dplyr)

# Load in data from the step above
data1 = read.table("$yourfile.sscore",header=T)

# Also load in covariates (this can be done in PLINK or using flashPCA)
data2 = read.table("$yourfile_covariates.txt",header=T)

# Drop the IID column to prevent double columns
data2$IID <- NULL

# Merge by FID column, can also make shorter in this case -> by='FID'
MM = merge(data2,data1,by.x='FID',by.y='FID')
temp <- anti_join(data1, data2, by='FID')
data <- rbind(temp, MM)

# Convert score values to Z-score
meanPRS <- mean(data$SCORE1_SUM) 
sdPRS <- sd(data$SCORE1_SUM)
data$SCOREZ <- (data$SCORE1_SUM - meanPRS)/sdPRS

# Generate boxplot case-control data 
pdf("$FILENAME.pdf",width=4)
boxplot(data$SCOREZ~data$PHENO.x,col=c('powderblue', 'pink1'),xlab="0 control, 1 PD-case",ylab="PRS Z-score")
grid()
dev.off()

# Calculate if the PRS is statistically significant between cases and controls using AGE, SEX and first 5 principal components (PC)
thisFormula1 <- formula(paste("PHENO ~ SCOREZ + SEX_COV + AGE + PC1 + PC2 + PC3 + PC4 + PC5"))
model1 <- glm(thisFormula1, data = data, ,family=binomial)
print(summary(model1))

# note sometimes the P-value shows this < 2e-16
# then use this : summary(model1)$coefficients[,4] 
# this will output the real P-value
```

Option 2: Python 
```bash, python3
module load python

# Load in the necessary packages 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm

# Load in data from PRS calculation step above (file ending in .sscore
# PLINK2 PRS output
data1 = pd.read_csv("yourfile.sscore", sep=r"\s+")

# If needed, reformat plink output header
# data1.rename(columns={"#FID": "FID"}, inplace=True)

# Covariates
data2 = pd.read_csv("yourfile_covariates.txt", sep=r"\s+")

# Drop the IID column to prevent columns duplication
if "IID" in data2.columns:
    data2 = data2.drop(columns=["IID"])

# Merge PRS and covariate dataframes by FID column
data = pd.merge(data1, data2, on="FID", how="left")

# Convert PRS score values to Z-score
data["SCOREZ"] = (
    data["SCORE1_SUM"] - data["SCORE1_SUM"].mean()
) / data["SCORE1_SUM"].std(ddof=1)

# Generate boxplot case-control data 
plt.figure(figsize=(4,4))

plt.boxplot(
    [
        data.loc[data["PHENO"] == 0, "SCOREZ"], # Verify that Plink file encodes control = 0
        data.loc[data["PHENO"] == 1, "SCOREZ"]  # Verify that Plink file encodes case = 1
    ],
    labels=["0 control", "1 PD-case"]
)

plt.xlabel("0 control, 1 PD-case")
plt.ylabel("PRS Z-score")
plt.grid(True)

plt.savefig("FILENAME.pdf")
plt.close()

# Logistic regression (glm family=binomial)
# Calculate if the PRS is statistically significant between cases and controls using AGE, SEX and first 5 principal components (PC)

predictors = [
    "SCOREZ", "SEX_COV", "AGE",
    "PC1", "PC2", "PC3", "PC4", "PC5"
]

model_data = data.dropna(subset=["PHENO"] + predictors)

X = model_data[predictors]
X = sm.add_constant(X)
y = model_data["PHENO"]

model1 = sm.Logit(y, X).fit()

print(model1.summary())
print("\nCoefficient p-values:")
print(model1.pvalues)
```
