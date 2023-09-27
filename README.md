# antagonist
Gene target prioritization and computational drug repurposing


# Expected inputs 

## Disease file format (csv, csv.gz or RDS)
The input df usually is a TWAS/GFI output file. If another file is used thene some column name changes are needed to work as expected. The required columns are as follows (with default names):
Column      | description
---         | ---
`feature`   | expected input is ENSEMBL ID and currently only works with genes
`zscore`    | z-score statistic
`gwas`      | trait name
`model_ID`  | disease signature source, e.g. microglia, DLPFC, meta-analysis, etc.

## Recipe file format (csv)

# Installation

## Local

## Minerva (interactive), run the code line by line.
```
ml pigz/2.3.1
ml git

R

envName    <- "/sc/arion/projects/roussp01a/sanan/Rlibs/230919_R_4.2.0_MultiWAS_Antagonist"
user.name  <- readline("What is your github username?\n")
user.email <- readline("What is your github email?\n")
user.PAT   <- readline("What is your github token?\n")

libs <- .libPaths()
libs[3] <- envName # replaces user path /hpc/users/[user]/.Rlib
libs2 <- libs[c(3,1)] # excludes bioconductor (2) and shared R libraries 
# .libPaths(libs) # doesn't work
.libPaths(libs2) # works

usethis::use_git_config(user.name = user.name, user.email = user.email)
credentials::set_github_pat(user.PAT,force_new =T) # Will need to reenter PAT in command line. Minerva often defaults to environmental PAT which can cause conflicts

options(timeout=9999999)

#Repeat following until libraries are installed. In some cases certain libraries are known to cause issues and may need to be independently installed with another command
.libPaths(libs)
remotes::install_github("voloudakislab/MultiWAS", build = FALSE)
remotes::install_github("voloudakislab/antagonist", build = FALSE)

.libPaths(libs2)
remotes::install_github("voloudakislab/MultiWAS", build = FALSE)
remotes::install_github("voloudakislab/antagonist", build = FALSE)
```

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
Todo:
- make logs folder

```
library=/sc/arion/projects/roussp01a/sanan/Rlibs/230919_R_4.2.0_MultiWAS_Antagonist
# library=/home/georgios/R/x86_64-pc-linux-gnu-library/4.2
recipe=/sc/arion/projects/va-biobank/PROJECTS/2023_09_microglia_DGE_gtp_cdr/project.recipe.csv

# Step 1 mothership
bsub -J Ant_S01_mothership -P acc_va-biobank -q premium -n 20 -R span[hosts=1] \
-R rusage[mem=3000] -W 1440 --oo logs/S01mothership.out -oe logs/S01mothership.err \
-L /bin/bash Rscript $library/antagonist/exec/antagonist_S01_wrapper.R --recipe $recipe \
--prototyping 2 --dryrun TRUE
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
