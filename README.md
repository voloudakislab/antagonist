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

# Parameters for external resource file location
## This should be idea
parent.signature.dir = "/sc/arion/projects/va-biobank/resources/CMap/cmap_l1000_2021-11-20/"
signature.dir        = "/sc/arion/projects/va-biobank/resources/CMap/cmap_l1000_2021-11-20/eachDrug/"
gene.anno.file       = "/sc/arion/projects/va-biobank/resources/CMap/cmap_l1000_2021-11-20/geneinfo_beta.txt"
sig.annotation.file  = "/sc/arion/projects/va-biobank/resources/CMap/cmap_l1000_2021-11-20/siginfo_beta.txt"
