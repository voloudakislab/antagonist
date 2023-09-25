# antagonist
Gene target prioritization and computational drug repurposing


# Expected input file format
The input df usually is a TWAS/GFI output file. If another file is used thene some column name changes are needed to work as expected. The required columns are as follows (with default names):
Column      | description
---         | ---
`feature`   | expected input is ENSEMBL ID and currently only works with genes
`zscore`    | z-score statistic
`gwas`      | trait name
`model_ID`  | disease signature source, e.g. microglia, DLPFC, meta-analysis, etc.

# Resources disclaimer



# STEP 1: Run antagonism
For one trait-tissue combination, it takes about 23,800 thread-minutes on an Intel 10th gen core.
## Interactive version
```
library(antagonist)

# Step 1
perform_antagonism()

# Step 2
```
## LSF version
For now only works on minerva. Potentially adding ability to run in other platforms in the future.
```
library=/sc/arion/projects/roussp01a/sanan/Rlibs/230919_R_4.2.0_MultiWAS_Antagonist
recipe=/sc/arion/projects/va-biobank/PROJECTS/2023_09_microglia_DGE_gtp_cdr/project.recipe.csv

# Step 1 mothership
bsub -J Ant_S01_mothership -P acc_va-biobank -q premium -n 20 -R span[hosts=1] \
-R rusage[mem=3000] -W 1440 --oo logs/S01mothership.out -oe logs/S01mothership.err \
-L /bin/bash $library/exec/antagonist_S01_wrapper.R --recipe $recipe
# $library/exec/antagonist_S01_wrapper.R --recipe $recipe --prototyping 2 # if you want to run with just two signature files for troubleshooting

# Step 2 mothership
# Waits for the first step to be completed


```



# Parameters for external resource file location
## This should be idea
parent.signature.dir = "/sc/arion/projects/va-biobank/resources/CMap/cmap_l1000_2021-11-20/"
signature.dir        = "/sc/arion/projects/va-biobank/resources/CMap/cmap_l1000_2021-11-20/eachDrug/"
gene.anno.file       = "/sc/arion/projects/va-biobank/resources/CMap/cmap_l1000_2021-11-20/geneinfo_beta.txt"
sig.annotation.file  = "/sc/arion/projects/va-biobank/resources/CMap/cmap_l1000_2021-11-20/siginfo_beta.txt"
