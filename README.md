# antagonist

A multithreaded R package and wrapper for gene target prioritization and computational drug repurposing.

*If you use this package for gene target prioritization (GTP), cite::*

Voloudakis G, Vicari JM, Venkatesh S, Hoffman GE, Dobrindt K, Zhang W, Beckmann ND, Higgins CA, Argyriou S, Jiang S, Hoagland D, Gao L, Corvelo A, Cho K, Lee KM, Bian J, Lee JS, Iyengar SK, Luoh SW, Akbarian S, Striker R, Assimes TL, Schadt EE, Lynch JA, Merad M, tenOever BR, Charney AW; Mount Sinai COVID-19 Biobank; VA Million Veteran Program COVID-19 Science Initiative; Brennand KJ, Fullard JF, Roussos P. A translational genomics approach identifies IL10RB as the top candidate gene target for COVID-19 susceptibility. NPJ Genom Med. 2022 Sep 5;7(1):52. doi: [10.1038/s41525-022-00324-x](https://doi.org/10.1038/s41525-022-00324-x). PMID: [36064543](https://pubmed.ncbi.nlm.nih.gov/36064543/); PMCID: PMC9441828.

*If you use this package for computational drug repurposing (CDR), cite:*

Voloudakis G, Lee KM, Vicari JM, Zhang W, Hoagland D, Venkatesh S, Bian J, Anyfantakis M, Wu Z, Rahman S, Gao L, Cho K, Lee JS, Iyengar SK, Luoh S-W21,22, Themistocles L. Assimes16,17, Gabriel E. Hoffman1-4, Benjamin R. tenOever9-11, John F. Fullard1-4‡, Julie A. Lynch8,12‡, Panos Roussos1-4,6,7‡†.

*If you use the 5-method rank (default and only option at this moment), please cite the relevant paper:*


# Computational environment requirements
1. A linux computer (package has been developed and tested in linux; may work in other operating systems but it hasn't been tested)
2. R>=4.0
3. Dependencies as per `DESCRIPTION` file
4. Lots of RAM if you want it to run fast over a higher number of threads

# Installation
```
devtools::install_github("DiseaseNeuroGenomics/antagonist") # link for the center's repository
# devtools::install_github("voloudakislab/antagonist") # Voloudakis lab development link
```

# Overview of the inputs
1. [The perturbagen signature library](#-perturbagen-signature-library): a data.frame with known transcriptional signatures for compounds/shRNAs, etc, in this case LINCS
2. A disease signature: a data.frame with genes and their respective changees (can be logFC, z-score, effect sizes, etc.)
3. A recipe file: this is only required with job schedulers such as IBM's LSF; the wiki will be updated in the future for such applications.

# Perturbagen signature library
We are currently using the Expanded CMap LINCS Resource 2020 signature files from clue.io. For installation of the perturbagen library, a total of 25GB are required (11.5 GB after deleting intermediate files). We are using level5 signatures (see picture with different levels below).
![LINCS signature level overview](/data-raw/readme.images/L1000_Lvl5.png)

## Download the Expanded CMap LINCS Resource 2020 signature files from clue.io:
We are currently using the version last updated on 11/23/201 (created on 11/20/2020) which can be downloaded [here](https://clue.io/data/CMap2020#LINCS2020).

The following files are required:

|File name                              | Data matrix  | File size | MD5 Check sum
|---                                    | ---          | ---       | ---
|geneinfo_beta.txt                      | N/A          |   1.09 MB | 45c725d17ce6c377f1e7de07b821a5f0
|siginfo_beta.txt                       | N/A          | 443.69 MB | ab609fc04fab21180b07833119d1c7b6
|level5_beta_ctl_n58022x12328.gctx      | 58022x12328  |   2.66 GB | 4c70a6939637670d185775f2a2d98d67
|level5_beta_trt_cp_n720216x12328.gctx  | 720216x12328 |  33.08 GB | 9a82806e2aba6ec2a866cba77bd57fda
|level5_beta_trt_misc_n8283x12328.gctx  | 8283x12328   | 389.61 MB | b58bcaa628f9f2afeadbb815b49a7684
|level5_beta_trt_oe_n34171x12328.gctx   | 34171x12328  |   1.57 GB | a068e259e2d71676e8fa46a2d7a2de86
|level5_beta_trt_sh_n238351x12328.gctx  | 238351x12328 |  10.95 GB | 16952edbdc39756370a075b25f874029
|level5_beta_trt_xpr_n142901x12328.gctx | 142901x12328 |   6.07 GB | c852ca26affaa144f1b042463036702b

Download and save these in a folder that will be used for storage of these resource, for the purposes of this tutorial we will use `ExpandedCMapLINCS2020/`. 

Out of the 1,201,944 signatures the vast majority are not considered reproducible or distinct (`is_gold`; n=237,922; distil_cc_q75 >= 0.2 and pct_self_rank_q25 <= 0.05), or may be further subsetted to remove replicates as much as possible (`is_exemplar`; n=423,422) as described in the [GEO CMap LINCS User Guide v2.1](https://docs.google.com/document/d/1q2gciWRhVCAAnlvF2iRLuJ7whrGP6QjpsCMq1yWz7dU). The breakdown is as follows:

|                  | is_exemplar (FALSE) | is_exemplar (TRUE)
|---               | ---                 | ---
|is_gold (FALSE)   | 668,957             | 295,065
|is_gold (TRUE)    | 109,565             | 128,357

For our projects we include all `is_gold` signatures and do no filtering based on `is_exemplar` status. However, filtering for both will half the computational costs.


## Chunk the signature files in `.RDS` objects
This is done for easier batch processing and improved IO performance when running multiple versions.
```
signature.dir <- "ExpandedCMapLINCS2020/" # where the .gctx files were downloaded.
antagonist::split_gctx(parent.signature.dir = signature.dir)
```

## Disease file format (csv, csv.gz or RDS)
The input data frame usually is a TWAS/GFI/DGE output file. If another file is used then some column name changes are needed to work as expected. The required columns are as follows (with default names):

|Column      | description
|---         | ---
|`feature`   | expected input is ENSEMBL ID and currently only works with genes
|`zscore`    | z-score statistic
|`gwas`      | trait name
|`model_ID`  | disease signature source, e.g. microglia, DLPFC, meta-analysis, etc.

Of note, human transcripts are expected. If you need to use data from another species, genes have to be mapped to their othologs first. 

For the purposes of the tutorial we will use 

!!! Add file for testing !!!

# Run the analyses
Setting up the variables and loading the package

```
library(antagonism)
signature.dir <- "ExpandedCMapLINCS2020/"
```


# STEP 1: Run antagonism (outdated)
For one trait-tissue combination, it takes about 23,800 thread-minutes on an Intel 10th gen core.
<<<<<<< HEAD
=======
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
>>>>>>> 63a8b6435046fff5029c992bce6fc3e8c0863e58

```
perform_antagonism(
signature.dir  = paste0(signature.dir, "eachDrug/",
gene.anno.file = paste0(signature.dir, "geneinfo_beta.txt"
)
```

<<<<<<< HEAD
# Step 2 Aggregate and prioritize()

```
aggregate_and_prioritize(
signature.dir  = paste0(signature.dir, "eachDrug/",
gene.anno.file = paste0(signature.dir, "geneinfo_beta.txt"
)
```

> Please note that the default way of meta-analyzing the results is pulling all the tissues or cell types (whatever is in model_ID) together.

# Inspect the output

# Additional figures
=======
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

## Variables of each class (in the "Variables" column) (both classes are explained for clarity)

For each object of **class Antag**, recipe contains 3 rows, e.g. for fiveRankJob:
- fiveRankJob
- parameters.fiveRankJob
- core.fiveRankJob


For each object of **class postAntag**, recipe contains 4 rows in , e.g. for avgRankJob:
- avgRankJob
- parameters.avgRankJob
- core.avgRankJob
- parent.avgRankJob

## Values for each variable (see right above) of class postAntag (containing both classes would be redundant) (todo: insert a screenshot img)

Example for avgRankJob (each of the below responds to an element of the recipe's 'variable' column with its corresponding value from the 'value' column):

**avgRankJob**: 
- run.me = TRUE (regulates whether you want to submit that job)
- walltime = 00:10 (sets walltime of the job)
- n.threads = 1 (sets the number of cores)
- mem = 400 (sets memory per core)

**parameters.avgRankJob**:
- these are parameters we want to pass to the recipe 

**core.avgRankJob**:
- "/path/to/script.R" (the Job-submission object submits a job that executes this script)

**parent.avgRankJob**:
- fiveRankJob (determines the which "Parent Job" this job will be referring to, in this example avgRankJob will wait for fiveRankJob to complete, names must be     consistent with the actual name of the parent Job)

## Adding a new job (todo: give coding examples here)

All the above being said, to add a job you have to:
- Use another job as template; for example copy-paste the rows avgRankJob, parameters.avgRankJob , core.avgRankJob, parent.avgRankJob.
- Then edit each of the above to xyzJob, parameters.xyzJob, core.xyzJob, parent.xyzJob
- Assign the appropriate values to the value column.
- Create the object in the wrapper script (I recommend the object to be created after the object to which it will depend on, e.g. avgRankJob after fiveRankJob)
- Add the object to the pipeline.



>>>>>>> 63a8b6435046fff5029c992bce6fc3e8c0863e58



