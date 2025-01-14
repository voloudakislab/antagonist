<!-- TOC start (generated with https://github.com/derlin/bitdowntoc) -->

- [How to cite this manuscript](#how-to-cite-this-manuscript)
- [Computational environment requirements](#computational-environment-requirements)
- [Installation](#installation)
- [Overview of the inputs](#overview-of-the-inputs)
- [Perturbagen signature library](#perturbagen-signature-library)
   * [Download the Expanded CMap LINCS Resource 2020 signature files from clue.io:](#download-the-expanded-cmap-lincs-resource-2020-signature-files-from-clueio)
   * [Chunk the signature files in `.RDS` objects](#chunk-the-signature-files-in-rds-objects)
   * [Disease file format (csv, csv.gz or RDS)](#disease-file-format-csv-csvgz-or-rds)
- [Run the analyses](#run-the-analyses)
- [STEP 1: Run antagonism](#step-1-run-antagonism)
- [STEP 2: Aggregate and prioritize](#step-2-aggregate-and-prioritize)
- [STEP 3: The output](#step-3-the-output)
- [STEP 4: Additional figures](#step-4-additional-figures)
   * [Showcasing a signature](#showcasing-a-signature)
   * [Generating a gene-target prioritization plot](#generating-a-gene-target-prioritization-plot)

<!-- TOC end -->

antagonist
==========

A multithreaded R package and wrapper for gene target prioritization and computational drug repurposing.

<!-- TOC --><a name="how-to-cite-this-manuscript"></a>
# How to cite this manuscript

*If you use this package for gene target prioritization (GTP), cite::*

Voloudakis G, Vicari JM, Venkatesh S, Hoffman GE, Dobrindt K, Zhang W, Beckmann ND, Higgins CA, Argyriou S, Jiang S, Hoagland D, Gao L, Corvelo A, Cho K, Lee KM, Bian J, Lee JS, Iyengar SK, Luoh SW, Akbarian S, Striker R, Assimes TL, Schadt EE, Lynch JA, Merad M, tenOever BR, Charney AW; Mount Sinai COVID-19 Biobank; VA Million Veteran Program COVID-19 Science Initiative; Brennand KJ, Fullard JF, Roussos P. A translational genomics approach identifies IL10RB as the top candidate gene target for COVID-19 susceptibility. NPJ Genom Med. 2022 Sep 5;7(1):52. doi: [10.1038/s41525-022-00324-x](https://doi.org/10.1038/s41525-022-00324-x). PMID: [36064543](https://pubmed.ncbi.nlm.nih.gov/36064543/); PMCID: PMC9441828.

*If you use this package for computational drug repurposing (CDR), cite:*

Voloudakis G, Lee KM, Vicari JM, Zhang W, Hoagland D, Venkatesh S, Bian J, Anyfantakis M, Wu Z, Rahman S, Gao L, Cho K, Lee JS, Iyengar SK, Luoh S-W21,22, Themistocles L. Assimes16,17, Gabriel E. Hoffman1-4, Benjamin R. tenOever9-11, John F. Fullard1-4‡, Julie A. Lynch8,12‡, Panos Roussos1-4,6,7‡†.

*If you use the 5-method rank (default and only option at this moment), please cite the relevant paper:*

So H-C, Chau CK-L, Chiu W-T, Ho K-S, Lo C-P, Yim SH-Y, et al. Analysis of genome-wide association data highlights candidates for drug repositioning in psychiatry. Nat Neurosci. 2017;20:1342–9. PMID:[28805813](https://pubmed.ncbi.nlm.nih.gov/28805813/)

<!-- TOC --><a name="computational-environment-requirements"></a>
# Computational environment requirements
1. A linux computer (package has been developed and tested in linux; may work in other operating systems but it hasn't been tested)
2. R>=4.0
3. Dependencies as per `DESCRIPTION` file
4. Lots of RAM if you want it to run fast over a higher number of threads

If you are using a standard Ubuntu installation, you will need to install the following packages:
```
sudo apt-get install libmpfr-dev 
```


<!-- TOC --><a name="installation"></a>
# Installation
```
devtools::install_github("DiseaseNeuroGenomics/antagonist") # link for the center's repository
# devtools::install_github("voloudakislab/antagonist") # Voloudakis lab development link
```
If there is an error with the installation of cytolib please try the [following](https://stackoverflow.com/questions/40721182/error-in-installing-flowcore-package-r):
```
install.packages('BiocManager')
BiocManager::install('flowCore', Ncpus = 8) # will fail but still necessary
install.packages('remotes')
remotes::install_github('RGLab/cytolib')
```

<!-- TOC --><a name="overview-of-the-inputs"></a>
# Overview of the inputs
1. [The perturbagen signature library](#-perturbagen-signature-library): a data.frame with known transcriptional signatures for compounds/shRNAs, etc, in this case LINCS
2. A disease signature: a data.frame with genes and their respective changees (can be logFC, z-score, effect sizes, etc.)
3. A recipe file: this is only required with job schedulers such as IBM's LSF; the wiki will be updated in the future for such applications.

<!-- TOC --><a name="perturbagen-signature-library"></a>
# Perturbagen signature library
We are currently using the Expanded CMap LINCS Resource 2020 signature files from clue.io. For installation of the perturbagen library, a total of ~67GB are required (11.5 GB after deleting intermediate files). We are using level5 signatures (see picture with different levels below).
![LINCS signature level overview](/data-raw/readme.images/L1000_Lvl5.png)

<!-- TOC --><a name="download-the-expanded-cmap-lincs-resource-2020-signature-files-from-clueio"></a>
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


<!-- TOC --><a name="chunk-the-signature-files-in-rds-objects"></a>
## Chunk the signature files in `.RDS` objects
This is done for easier batch processing and improved IO performance when running multiple versions.
```
library(antagonist)
signature.dir <- "ExpandedCMapLINCS2020/" # where the .gctx files were downloaded.
antagonist::split_gctx(parent.signature.dir = signature.dir)
```
The `level5_*.gctx` files can now be safely deleted.

<!-- TOC --><a name="disease-file-format-csv-csvgz-or-rds"></a>
## Disease file format (csv, csv.gz or RDS)
The input data frame usually is a TWAS/GFI/DGE output file. If another file is used then some column name changes are needed to work as expected. The required columns are as follows (with default names):

|Column         | description
|---            | ---
|`feature`      | expected input is ENSEMBL ID and currently only works with genes
|`feature_name` | name of the feature, e.g. gene_symbol *PSEN1*
|`zscore`       | z-score statistic
|`pvalue`       | p value
|`gwas`         | trait name
|`model_ID`     | disease signature source, e.g. microglia, DLPFC, meta-analysis, etc.

Of note, human transcripts are expected. If you need to use data from another species, genes have to be mapped to their othologs first. 

For the purposes of the tutorial we will use 

!!! Add file for testing !!!

<!-- TOC --><a name="run-the-analyses"></a>
# Run the analyses
Setting up the variables and loading the package
> For the tutorial, we will use the genetically regulated gene expression for Rheumatoid Arthritis 

```
library(antagonist)
signature.dir    <- "ExpandedCMapLINCS2020/"
disease.sig.file <- system.file("extdata", "sample.datasets/RA.epixcan.csv.gz", package="antagonist")
```

<!-- TOC --><a name="step-1-run-antagonism"></a>
# STEP 1: Run antagonism
For one trait-tissue combination, it takes about 23,800 thread-minutes on an Intel 10th gen core.
> For testing if the pipeline is running, setting `prototyping = 10`, for example, which means that only 10/300 signatures will be used from each signature file would allow to see if there are any errors.
> In the case that one is not interested in whether individual signatures significantly reverse the disease signature, but want to use just the aggregation and prioritiziation pipeline, they can set `noperm=3` (default is 100) as in this tutorial. This decreases runtime and RAM usage significantly
> Default `n.threads` setting is to use all available logical CPUs minus 2 which assumes server-grade RAM availability. In consumer-grade hardware, this may need to be decreased.

```
perform_antagonism(
df             = disease.sig.file,
signature.dir  = paste0(signature.dir, "eachDrug/"),
gene.anno.file = paste0(signature.dir, "geneinfo_beta.txt"),
noperm         = 3 # this is just for the tutorial to reduce run times
)
```

<!-- TOC --><a name="step-2-aggregate-and-prioritize"></a>
# STEP 2: Aggregate and prioritize

```
aggregate_and_prioritize()
```

> Please note that the default way of meta-analyzing the results is pulling all the tissues or cell types (whatever is in model_ID) together.

<!-- TOC --><a name="step-3-the-output"></a>
# STEP 3: The output
This is the out put folder structure, if there are more than one tissues parsed, then 
```
.
└── results
    └── GTP_CDR
        ├── intermediate.files
        │   ├── all.signatures.csv.gz
        │   ├── signature.inventory.csv
        │   └── signature.location.csv
        ├── RA
        │   ├── ALL 
        │   │   ├── RA_cdr_all_compound_level.csv
        │   │   ├── RA_cdr_AvgRank_distribution_landscape.pdf
        │   │   ├── RA_cdr_AvgRank_distribution_square.pdf
        │   │   ├── RA_cdr_disease_area_level.csv
        │   │   ├── RA_cdr_indication_level.csv
        │   │   ├── RA_cdr_launched_compound_level.csv
        │   │   ├── RA_cdr_moa_level.csv
        │   │   ├── RA_cdr_P3_and_launched_compound_level.csv
        │   │   ├── RA_cdr_signature_level.csv.gz
        │   │   ├── RA_cdr_target_level.csv
        │   │   ├── RA_gtp_AvgRank_distribution_landscape.pdf
        │   │   ├── RA_gtp_AvgRank_distribution_square.pdf
        │   │   ├── RA_gtp_compound_level.csv
        │   │   └── RA_gtp_signature_level.csv.gz
        │   └── STARNET_BLD
        │       ├── RA_cdr_all_compound_level.csv
        │       ├── RA_cdr_AvgRank_distribution_landscape.pdf
        │       ├── RA_cdr_AvgRank_distribution_square.pdf
        │       ├── RA_cdr_disease_area_level.csv
        │       ├── RA_cdr_indication_level.csv
        │       ├── RA_cdr_launched_compound_level.csv
        │       ├── RA_cdr_moa_level.csv
        │       ├── RA_cdr_P3_and_launched_compound_level.csv
        │       ├── RA_cdr_signature_level.csv.gz
        │       ├── RA_cdr_target_level.csv
        │       ├── RA_gtp_AvgRank_distribution_landscape.pdf
        │       ├── RA_gtp_AvgRank_distribution_square.pdf
        │       ├── RA_gtp_compound_level.csv
        │       └── RA_gtp_signature_level.csv.gz
        ├── trt_cp_all.signatures.AvgRank.csv.gz
        ├── trt_cp_all.signatures.AvgRank.readme
        ├── trt_oe_all.signatures.AvgRank.csv.gz
        ├── trt_oe_all.signatures.AvgRank.readme
        ├── trt_sh_all.signatures.AvgRank.csv.gz
        ├── trt_sh_all.signatures.AvgRank.readme
        ├── trt_xpr_all.signatures.AvgRank.csv.gz
        └── trt_xpr_all.signatures.AvgRank.readme

```
Most of the files are intuitively named. The main abbreviations to understand the output are:
Treatment type abbreviations: `trt_cp`:compounds; `trt_sh`:shRNA; `trt_oe`:over-expression; `trt_xpr`:CRISPR
Type of analysis: `cdr`:computational drug repurposing using `trt_cp`; `gtp`: gene target prioritization using `trt_sh`, `trt_oe` and `trt_xpr`.
> Please note that the signatures are inverted for `trt_oe` for the summarization and this should be taken into account when interpreting `trt_oe` files.

Primary outputs are the `*_compound_level.csv` files that provide rankings and statistics summarized at the level of the perturbagen, compound and gene for the cdr and gtp analyses, respectively. The `*_signature_level.csv.gz` files allow for exploration of the performance of specific signatures for downstream applications (e.g. to determine optimal concentration).

By leveraging compound annotations for the cdr analysis from PMID:[28388612](https://pubmed.ncbi.nlm.nih.gov/28388612/), we are able to perform additional ranking/enrichment analyses based on compounds' mechanism of action (moa), indications, and applied disease area. Further prioritization can be done based on preclinical/clinical testing stage, e.g. `P3_and_launched_compound` refers to compounds that are either launched (may not be FDA-approved but used in other countries or safe for veterinarian use) or in phase III clinical trials.

A few useful column definitions across files are the following:

|Column|Description
|---|---|
|Compound|Compound name|
|Rank|Compound/moa rank|
|\*.MW.p|Compound/moa Mann-Whitney U test p value; this considers the full distribution prior any filtering, e.g. for launched compounds|
|\*.MW.FDR|FDR-adjusted \*.MW.P values|
|\*.pseudo.zscore|Compound/MOA pseudo zscore, a positive score means that it is predicted to antagonize the trait: -(Hedges-Lehmann estimator for compound / SD average ranks of all compounds) while considering compounds from all clinical phases (not just the “launched: ones)|
|target|Target gene(s)|
|disease_area|Disease area for compound indication|
|indication|Compound indication|
|N_experiments|Number of experiments/signatures prior to any filtering|

Average rank (`AvgRank`) distribution plots are also generated for diagnostic purposes. Please note that when using cluster comparison signatures, the z-scores of the disease signatures and the `AvgRank` distributions may not be normally distributed, but the pipeline is non-parametric and can accomodate for that.

![CDR AvgRank Distribution](/data-raw/readme.images/RA_cdr_AvgRank_distribution_landscape.png)

![GTP AvgRank Distribution](/data-raw/readme.images/RA_gtp_AvgRank_distribution_landscape.png)

<!-- TOC --><a name="step-4-additional-figures"></a>
# STEP 4: Additional figures
Additional figures can be prepared

<!-- TOC --><a name="showcasing-a-signature"></a>
## Showcasing a signature
For example, how does actinomycin D transcriptional signature antagonize the RA disease signature?

```
showcase_method_cdr_gtp(
  "actinomycin-d",
  twas.model = "STARNET_BLD",
  twas.trait = "RA",
  twas = disease.sig.file,
  gene.anno.file = paste0(signature.dir, "geneinfo_beta.txt")
)

```

![actinomycin D in RA](/data-raw/readme.images/RA.STARNET_BLD.actinomycin-d.CRCGN004_PC3_6H.BRD-A42383464-001-04-8.10.png)

<!-- TOC --><a name="generating-a-gene-target-prioritization-plot"></a>
## Generating a gene-target prioritization plot
```
gtp_pvalue_qqplot(    ### Parameters
  thistrait       = "RA",
  thistissue      = "STARNET_BLD",
  stouffer.ma     = TRUE,
  z.score.scaling = "before",
  twas.df         = disease.sig.file,
  figure.name     = "gtp_actions", # main body of the name of the figure
)
```
In this there is no actionable gene that demonstrates a FDR-significant compbined score while also having FDR-significant dysregulation in the disease signatures and antagonism via shRNA/overexpression/CRISPR.

![actinomycin D in RA](/data-raw/readme.images/prioritization.stouffer.scaled.before.png)

