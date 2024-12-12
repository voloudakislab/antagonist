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
```
MultiWAS::gv_install_packages(
cran.packages = c("optparse"),
bioc.packages = c("cmapR")
)
```
## Minerva (interactive), run the code line by line.
```
ml pigz/2.3.1
ml git
ml R
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

# Updating the package

## Minerva
Interactive bash
```
ml pigz/2.3.1
ml git
ml R
R
```
Interactive R
```
# Run this line by line
user.name  <- readline("What is your github username?\n")
user.email <- readline("What is your github email?\n")
user.PAT   <- readline("What is your github token?\n")

# 
envName    <- "/sc/arion/projects/roussp01a/sanan/Rlibs/230919_R_4.2.0_MultiWAS_Antagonist"
libs <- .libPaths()
libs[3] <- envName # replaces user path /hpc/users/[user]/.Rlib
libs2 <- libs[c(3,1)] # excludes bioconductor (2) and shared R libraries 
.libPaths(libs2) # works
usethis::use_git_config(user.name = user.name, user.email = user.email)

# Run this and provide the PAT
credentials::set_github_pat(user.PAT,force_new =T) # Will need to reenter PAT in command line. Minerva often defaults to environmental PAT which can cause conflicts

# Run this
options(timeout=9999999)
.libPaths(libs2)
remotes::install_github("voloudakislab/antagonist", build = FALSE, upgrade = "never")
```
To reload the package in an interactive R session
```
MultiWAS::reload_package("antagonist")
```

# Resources disclaimer




# STEP 1: Run antagonism (outdated)
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
- sh run ml R in bsub

```
bsub -P acc_va-biobank -q premium -n 30 -W 24:00 -Ip /bin/bash
ml R
library=/sc/arion/projects/roussp01a/sanan/Rlibs/230919_R_4.2.0_MultiWAS_Antagonist
# library=/home/georgios/R/x86_64-pc-linux-gnu-library/4.2
recipe=/sc/arion/projects/va-biobank/PROJECTS/2023_09_microglia_DGE_gtp_cdr/project.recipe.csv

# Step 1 mothership

bsub -J Ant_S01_mothership -P acc_va-biobank -q premium -n 20 -R span[hosts=1] \
-R rusage[mem=3000] -W 1440 --oo logs/S01mothership.out -oe logs/S01mothership.err \
-L /bin/bash Rscript --verbose $library/antagonist/exec/antagonist_S01_wrapper.sh --recipe $recipe \
--prototyping 2

bsub -J Ant_S01_mothership -P acc_va-biobank -q premium -n 20 -R span[hosts=1] \
-R rusage[mem=3000] -W 1440 --oo logs/S01mothership.out -oe logs/S01mothership.err \
-L /bin/bash Rscript --verbose $library/antagonist/exec/antagonist_S01_wrapper.R --recipe $recipe \
--prototyping 2



# $library/exec/antagonist_S01_wrapper.R --recipe $recipe --prototyping 2 # if you want to run with just two signature files for troubleshooting

# Step 2 mothership
# Waits for the first step to be completed


```

# Run Pipeline (only works in Minerva for now)
## Setup
- Create a directory to store all your recipes e.g. /path/to/recipes/ 
- Store in the above directory file my_recipe.csv
- Edit the antagonist_wrapper_V2.R in the following way
```
path_to_recipes = '/path/to/recipes/'

recipe.file = paste0(path_to_recipes, 'my_recipe.csv')

```
## Execution
Execute the Rscript in 'nohup' mode because the execution itself may take quite some time:
```
nohup Rscript antagonist_wrapper_V2.R > output.log 2>&1 &
```

# Pipeline Explained
## Summary / TLDR
In general you can run the pipeline both step-wise (each part executed on its own) and all at once.
I discuss this in detail further below, but briefly you can do this by switching run.me to TRUE/FALSE.

In my efforts to keep the code as modulated/readable as possible I organized the job submission operations
in different objects. In this README; objects that are used to submit jobs will be referenced as "Job-submission objects"

Some parameters of the Job-submission objects are determined inside the R script, while some others are passed to the object through the recipe.

Each Job-submission object submits a job that will execute a script (which can be either R or python), which I will from now on reference as "Job-script".

I recommend that you keep open the Recipe example as you continue further in this guide.

## Recipe Columns explained
Most of the parameters of the submitted jobs are passed to the Job-submission object through the Recipe.
Recipe is pretty self-explanatory but here goes nothing; 
elements of the 'variable' column will be "same name variables" in the wrapper script,
while elements of the 'value' column will be the values to be assigned to the variables of the row they are found in.

## Classes explained
As of 12/12/2024 there are 2 classes:
	class Antag (antagonizes signatures)
	class postAntag (processes results of class Antag or of class postAntag; objects of this class will usually be dependent on their "parent Job" 
	meaning they will wait for the completion of the parent Job before proceeding in executing their own job (through bsub -w))

This way you can have "chains of Job objects" executing one after the other.

## Explain variables for each class

For each object of class Antag, recipe contains 3 rows, e.g. for fiveRankJob:
	fiveRankJob
	parameters.fiveRankJob
	core.fiveRankJob


For each object of class postAntag, recipe contains 4 rows in , e.g. for avgRankJob:
	avgRankJob
	parameters.avgRankJob
	core.avgRankJob
	parent.avgRankJob

## Values for each of the Job-submission objects explained (todo: insert a screenshot img)

Example for avgRankJob (each of the below responds to an element of the recipe's 'variable' column with its corresponding value from the 'value' column):

avgRankJob: 
	- run.me = TRUE (regulates whether you want to submit that job)
	- walltime = 00:10 (sets walltime of the job)
	- n.threads = 1 (sets the number of cores)
	- mem = 400 (sets memory per core)

parameters.avgRankJob:
	- these are parameters we want to pass to the recipe 

core.avgRankJob:
	- "/path/to/script.R" (the Job-submission object submits a job that executes this script)

parent.avgRankJob:
	- fiveRankJob (determines the which "Parent Job" this job will be referring to, in this example avgRankJob will wait for fiveRankJob to complete, names must be     consistent with the actual name of the parent Job)

## Adding a new job (todo: give coding examples here)

All the above being said, to add a job you have to:
- Use another job as template; for example copy-paste the rows avgRankJob, parameters.avgRankJob , core.avgRankJob, parent.avgRankJob.
- Then edit each of the above to xyzJob, parameters.xyzJob, core.xyzJob, parent.xyzJob
- Assign the appropriate values to the value column.
- Create the object in the wrapper script (I recommend the object to be created after the object to which it will depend on, e.g. avgRankJob after fiveRankJob)
- Add the object to the pipeline.






# Parameters for external resource file location
## This should be idea
parent.signature.dir = "/sc/arion/projects/va-biobank/resources/CMap/cmap_l1000_2021-11-20/"
signature.dir        = "/sc/arion/projects/va-biobank/resources/CMap/cmap_l1000_2021-11-20/eachDrug/"
gene.anno.file       = "/sc/arion/projects/va-biobank/resources/CMap/cmap_l1000_2021-11-20/geneinfo_beta.txt"
sig.annotation.file  = "/sc/arion/projects/va-biobank/resources/CMap/cmap_l1000_2021-11-20/siginfo_beta.txt"
