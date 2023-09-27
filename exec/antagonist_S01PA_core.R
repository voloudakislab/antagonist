######################
# Perform antagonism #
######################
# consider using an alias for shortcut.

# STEP 1: run the 5-rank method and permutations
#     Part A: run 5-rank method for each signature batch
#     Part B: join results and split by disease-source combination
#     Part C: Summarize and run permutation analysis for target

######################
# HARDCODED PARAMETERS
## Look here if the script fails
R.lib.dir <- "/sc/arion/projects/roussp01a/sanan/Rlibs/230919_R_4.2.0_MultiWAS_Antagonist"


###################################################
# R library settings (minerva and package specific)
env  <- R.lib.dir
libs <- .libPaths()
libs[3] <- env
.libPaths(libs)

###################
# Handle parameters
suppressMessages(library(optparse))
suppressMessages(library(antagonist))
option_list = list(
  make_option(c("-r", "--recipe"),
              type="character",
              default=NULL,
              help = "Point to recipe file for this project",
              metavar = "character"),
  make_option(c("-c", "--cmapfile"),
              type="character",
              default=NULL,
              help = "Use a number to only run a subset of the signatures, e.g. 2 works well",
              metavar = "character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
if (is.null(opt$recipe)){
  print_help(opt_parser)
  stop("A recipe file must be provided", call.=FALSE)
}
if (is.null(opt$cmapfile)){
  print_help(opt_parser)
  stop("A cdr file must be provided", call.=FALSE)
}

##############################################
# Parse the recipe file and optparse arguments
recipe <- parse_recipe(opt$recipe)
#  need to get directory information
setwd(recipe$working.directory)
# Perform first step of the function
i <- opt$cmapfile
file.prefix <- sub("\\.RDS$", "", basename(i))


##########################
# Load the necessary files
results.dir <- recipe$results.dir
n.threads   <- recipe$n.threads
df <- fread(paste0(results.dir, "intermediate.files/df.shaped.csv.gz"))


##########################
# Run the five-rank method

## Preparing the signature library
## the splitted drug-induced gene expressions from CMap database: need to be splitted due to large number of perturbations
## it is a matrix of expression t-statistics/z-statistics with entrez gene id as row names and drugs as column names (we assume a drug expression was compared with controls)
cmap <- readRDS(i)
if (grepl("trt_oe", i)) cmap <- cmap * -1 # So that the signatures can be integrated with the knockdown experiments.
cmap <- as.data.table(cmap, keep.rownames = "entrezgene")
cmap[, entrezgene := as.integer(entrezgene)]
# cmap=data.frame(cmap) # change the type back
# Load gene annotation
geneAnn <- fread(recipe$gene.anno.file)
geneAnn <- unique(geneAnn[, c("gene_id", "ensembl_id")])
geneAnn <- geneAnn[!is.na(ensembl_id) & ensembl_id != ""]
setnames(geneAnn, "gene_id", "entrezgene")
cmap <- data.table::as.data.table(suppressMessages(dplyr::left_join(geneAnn, cmap)))
cmap <- cmap[, entrezgene:=NULL]
cmap <- cmap[!duplicated(ensembl_id)]
setnames(cmap, "ensembl_id", "feature")
sig_id <- names(cmap)[2:(dim(cmap)[2])]
# table(duplicated(df$gene_name)) ensembl are the best intersection of IDs between cmap and our TWAS approach.
payload <- as.data.table(suppressMessages(dplyr::inner_join(df, cmap)))
# Order by absolute z-score.
payload <- payload[order(abs(zscore), decreasing = T)]
# rm(list = c("cmap", "df", "geneAnn"))
to.process <- as.data.table(tidyr::expand_grid(
  unique(payload[, c("gwas","model_ID") ]),
  sig_id)) #, "thres.N.Vector" = c(NA,thres.N.vector)))

# Run loop for each signature in the file
signature <- mclapply(
  seq(nrow(to.process)),
  # seq(5),
  FUN = function(j){

    # limit the payload
    colofinterest <- c("feature","zscore",to.process[j]$sig_id)
    thispayload   <- payload[
      gwas == to.process[j]$gwas &
        model_ID == to.process[j]$model_ID,
      ..colofinterest, with = F]


    # Prime output
    output.left <- to.process[j]
    output.right <- lapply(
      stats::setNames(c(NA,thres.N.vector), c("ALL",thres.N.vector)),
      five_rank_method, x = thispayload, scramble = F)
    output.perm <- lapply(
      stats::setNames(
        rep(c(NA,thres.N.vector), noperm),
        paste0("Perm_", rep(1:noperm, each = length(c(NA,thres.N.vector))),"_", c("ALL",thres.N.vector))),
      five_rank_method, x = thispayload, scramble = T)
    output.right <- do.call(cbind, output.right)
    output.right[, ALL.ks.signed:=NULL]
    output.perm <- do.call(cbind, output.perm)
    drop.cols   <- colnames(output.perm)[grep("ALL.ks.signed", colnames(output.perm))]
    output.perm[, (drop.cols) := NULL]
    return(cbind(cbind(output.left, output.right), output.perm))
  },
  mc.cores = ifelse(n.threads<1,1,n.threads)
) # mclapply loop for all signatures

fwrite(
  signature,
  paste0(results.dir, "/intermediate.files/5rank/S01A/",
         file.prefix,".csv.gz")
  )
