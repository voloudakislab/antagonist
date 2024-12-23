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
R.lib.dir <- "/sc/arion/projects/va-biobank/software/Georgios_dev/240702_R_4.2.0_MultiWAS_Antagonist/"

libs <- .libPaths()
libs[3] <- "/sc/arion/projects/roussp01a/sanan/Rlibs/231221_R_4.2.0_MultiWAS"
.libPaths(libs)
# Load MultiWAS
library(MultiWAS)

###################################################
# R library settings (minerva and package specific)
env  <- R.lib.dir
libs <- .libPaths()
libs[6] <- env
.libPaths(libs)


#######################
#' Parse LSF project recipe file
#'
#' @param recipe.file recipe file
#' @return A list file with all included parameters
#' @keywords antagonism step1 LSF
#' @export

five_rank_method <- function(
    thres.N     , # threshold NA is all values otherwise consider top values
    x           , # the payload of interest
    scramble    = F # scramble is used for permutations
) {
  
  # Sample if scramble for permutations
  if(scramble) x$zscore <- sample(x$zscore)
  
  # Do these once:
  x$cp.z.rank <- rank(-x[[3]])
  
  # Determine disease and signature vectors
  if (is.na(thres.N)) { # considers all genes
    genesoi          <- x$feature
    disease.zscore   <- as.numeric(x$zscore)
    signature.zscore <- as.numeric(x[[3]])
  } else { # considers top genes
    
    ###########
    # KS method
    # For the ones that are up
    geneset2 <- sort(x[zscore>=0][seq(thres.N)]$cp.z.rank)
    a.up = max(   (1:thres.N)/thres.N - geneset2/nrow(x)      )
    b.up = max(   geneset2/nrow(x) -  (1:thres.N-1)/thres.N )
    ks.test.obj.up = ifelse(a.up>b.up, a.up , -b.up)
    # For the ones that are down
    geneset2 <- sort(x[zscore<0][seq(thres.N)]$cp.z.rank)
    a.down = max(   (1:thres.N)/thres.N - geneset2/nrow(x)      )
    b.down = max(   geneset2/nrow(x) -  (1:thres.N-1)/thres.N )
    ks.test.obj.down = ifelse(a.down>b.down, a.down , -b.down)
    
    ks.signed <- ifelse(
      ks.test.obj.up != ks.test.obj.down,
      ks.test.obj.up - ks.test.obj.down, 0                 )
    # rm(list = c("geneset2", "a.up", "b.up", "a.down", "b.down",
    #             "ks.test.obj.up", "ks.test.obj.down"))
    
    disease.zscore   <- rbind(
      x[zscore>=0][seq(thres.N)],
      x[zscore<0][seq(thres.N)])
    genesoi          <- disease.zscore$feature
    signature.zscore <- as.numeric(disease.zscore[[3]])
    disease.zscore   <- as.numeric(disease.zscore$zscore)
  }
  
  # make the calculations
  spearman.obj = cor.test(
    signature.zscore, disease.zscore,
    method="spearman", use="na.or.complete", exact = FALSE)
  pearson.obj  = cor.test(
    signature.zscore, disease.zscore,
    method = "pearson",use="na.or.complete")
  
  
  # populate the table
  if (is.na(thres.N)) {
    preoutput <- data.table(
      cor.spearman = spearman.obj$estimate,
      cor.pearson  = pearson.obj$estimate,
      p.spearman   = spearman.obj$p.value,
      p.pearson    = pearson.obj$p.value,
      ks.signed    = ifelse(is.na(thres.N), NA, ks.signed)
    )
  } else {
    preoutput <- data.table(
      cor.spearman = spearman.obj$estimate,
      cor.pearson  = pearson.obj$estimate,
      ks.signed    = ifelse(is.na(thres.N), NA, ks.signed)
    )
  }
  
  # names(preoutput) <- paste0(names(preoutput), "_", ifelse(is.na(thres.N), "ALL", thres.N))
  return(preoutput)
}



#suppressMessages(library(antagonist))

message('Proceeding in Execution')
option_list = list(
  make_option(c("-r", "--recipe.file"),
              type="character",
              default=NULL,
              help = "Point to recipe file for this project",
              metavar = "character"),
  make_option(c("-g", "--thisgwas"),
              type="character",
              default=NULL,
              help = "Point to recipe file for this project",
              metavar = "character"),
  make_option(c("-m", "--thismodelID"),
              type="character",
              default=NULL,
              help = "Use a number to only run a subset of the signatures, e.g. 2 works well",
              metavar = "character"),
  make_option(c("-d", "--results.dir"),
              type="character",
              default=NULL,
              help = "Use a number to only run a subset of the signatures, e.g. 2 works well",
              metavar = "character"),
  make_option(c("-s", "--results.subdir"),
              type="character",
              default=NULL,
              help = "Use a number to only run a subset of the signatures, e.g. 2 works well",
              metavar = "character"),
  make_option(c("-p", "--parent.subdir"),
              type="character",
              default=NULL,
              help = "Use a number to only run a subset of the signatures, e.g. 2 works well",
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
message(' I loaded the bash environment arguments succesfully')

############# MARIOS MODIFICATION
##################################
###### FOR DEBUG START HERE ######
#opt = list(
#  thisgwas = 'AD',
#  thismodelID = 'Microglia_Genes_.._Microglia_FANS_.._Genes_.._PrediXcan_.._EUR',
#  cmapfile = '/sc/arion/projects/va-biobank/PROJECTS/marios_temp/parallel_antagonism/materials/compounds/test_run/LINCS_test/trt_cp_1.300.RDS',
#  recipe = '/sc/arion/projects/va-biobank/PROJECTS/marios_temp/parallel_antagonism/recipes/test_run_recipes/V1/S4_test_run_V2.csv',
#  results.subdir = '/intermediate.files/five.rank'
#)

message('I will now parse recipe')

#################################

##############################################
# Parse the recipe file and optparse arguments
recipe <- parse_recipe(opt$recipe)

#
#i <- readRDS(file("stdin"))

#  need to get directory information
setwd(recipe$working.directory)
results.subdir <- opt$results.subdir
thisgwas <- opt$thisgwas
thismodelID <- opt$thismodelID
i <- opt$cmapfile
file.prefix <- sub("\\.RDS$", "", basename(i))


##########################
# Load the necessary files
results.dir <- recipe$results.dir
df          <- fread(paste0(results.dir, "/intermediate.files/df.shapes/", thisgwas, '/', thismodelID, '.shaped.csv.gz'))

## Other script parameters
path.to.gene.anno.file = recipe$gene.anno.file
gene.anno.file = ma_paste0(file.path(system.file('extdata', package = 'MultiWAS'), path.to.gene.anno.file))


noperm.pair = unlist(strsplit(recipe$parameters.fiveRankJob, split = ', '))[1]
noperm.str = unlist(strsplit(noperm.pair, split = ' = '))[2]
noperm = eval(parse(text = noperm.str))

thres.N.vector.pair = unlist(strsplit(recipe$parameters.fiveRankJob, split = ', '))[2]
thres.N.vector.str = unlist(strsplit(thres.N.vector.pair, split = ' = '))[2]
thres.N.vector = eval(parse(text = thres.N.vector.str))

n.threads.pair = unlist(strsplit(recipe$parameters.fiveRankJob, split = ', '))[3]
n.threads.str = unlist(strsplit(n.threads.pair, split = ' = '))[2]
n.threads = eval(parse(text = n.threads.pair))


##########################
# Run the five-rank method

## Preparing the signature library
## the splitted drug-induced gene expressions from CMap database: need to be splitted due to large number of perturbations
## it is a matrix of expression t-statistics/z-statistics with entrez gene id as row names and drugs as column names (we assume a drug expression was compared with controls)
cmap <- readRDS(i)

message('All external variables have been successfully improted')


if (grepl("trt_oe", i)) cmap <- cmap * -1 # So that the signatures can be integrated with the knockdown experiments.
cmap <- as.data.table(cmap, keep.rownames = "entrezgene")
cmap[, entrezgene := as.integer(entrezgene)]
# cmap=data.frame(cmap) # change the type back
# Load gene annotation

message('I modified the cmap file and now will modify geneAnn')
geneAnn <- fread(recipe$gene.anno.file)
geneAnn <- unique(geneAnn[, c("gene_id", "ensembl_id")])
geneAnn <- geneAnn[!is.na(ensembl_id) & ensembl_id != ""]
setnames(geneAnn, "gene_id", "entrezgene")

message('Now will perform more modifications to cmap variable')
cmap <- data.table::as.data.table(suppressMessages(dplyr::left_join(geneAnn, cmap)))
cmap <- cmap[, entrezgene:=NULL]
cmap <- cmap[!duplicated(ensembl_id)]
setnames(cmap, "ensembl_id", "feature")

message('Now will extract the sig_id')
sig_id <- names(cmap)[2:(dim(cmap)[2])]
# table(duplicated(df$gene_name)) ensembl are the best intersection of IDs between cmap and our TWAS approach.
payload <- as.data.table(suppressMessages(dplyr::inner_join(df, cmap)))
# Order by absolute z-score.
payload <- payload[order(abs(zscore), decreasing = T)]
# rm(list = c("cmap", "df", "geneAnn"))

message("Will create the to.process file")
to.process <- as.data.table(tidyr::expand_grid(
  unique(payload[, c("gwas","model_ID") ]),
  sig_id)) #, "thres.N.Vector" = c(NA,thres.N.vector)))
#print(paste0('the old to.proceess has dimension: ', dim(to.process)))

### prototyping:
#to.process <- to.process[, .SD[sample(.N, min(20, .N))], by = .(gwas, model_ID)]
#to.process1 <- as.data.frame(to.process1)
#print(paste0('the new to.proceess has dimension: ', dim(to.process)))

# add tryCatch in pbmclapply that will raise an error

message('I will do: signature <- pbmclapply using ', parallel::detectCores() - 2, ' cores.')

# Run loop for each signature in the file
signature <- pbmclapply(
  seq(nrow(to.process)),
  #seq(5),
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
  mc.cores = parallel::detectCores() - 2
) # mclapply loop for all signatures

signature <-do.call(rbind,signature)

# NO NEED FOR LAPPLY HERE
# Easy way to go recursively
MultiWAS::gv_dir.create(ma_paste0(file.path(results.dir, results.subdir, "5rank/S01A/",
                                            MultiWAS::make_java_safe(thisgwas),
                                            MultiWAS::make_java_safe(thismodelID))))


fwrite(signature, ma_paste0(file.path(results.dir, results.subdir, "5rank/S01A/",
                                      MultiWAS::make_java_safe(thisgwas),
                                      MultiWAS::make_java_safe(thismodelID),
                                      paste0(file.prefix,".csv.gz"))))

#THIS IS WILL BE OVERWRITTEN BY EACH SET OF COMPOUNDS..
#Save master file (needed for checks)
#fwrite(
#  signature,
#  paste0(results.dir, "/intermediate.files/5rank/S01A/",
#         file.prefix,".csv.gz")
#)
