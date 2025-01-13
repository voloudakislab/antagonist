libs <- .libPaths()
libs[3] <- "/sc/arion/projects/roussp01a/sanan/Rlibs/231221_R_4.2.0_MultiWAS"
libs[4] <- "/sc/arion/projects/va-biobank/software/Georgios_dev/240702_R_4.2.0_MultiWAS_Antagonist/"
.libPaths(libs)
# Load MultiWAS
library(MultiWAS)
library(antagonist)

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
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#################################
#####         DEBUG         #####
#opt = list(
#  recipe = '/sc/arion/projects/va-biobank/PROJECTS/marios_temp/parallel_antagonism/recipes/test_run_recipes/V1/S4_test_run_V2.csv',
#  thisgwas = 'AD',
#  thismodelID = 'Microglia_Genes_.._Microglia_FANS_.._Genes_.._PrediXcan_.._EUR',
#  results.subdir = '/intermediate.files/avgRank',
#  parent.subdir = 'intermediate.files/five.rank'
#)


# From BASH
thisgwas = opt$thisgwas
thismodelID = opt$thismodelID
parent.subdir = opt$parent.subdir
results.subdir = opt$results.subdir

# FROM RECIPE
recipe = parse_recipe(opt$recipe)
results.dir = recipe$results.dir
grep.sig.pattern = recipe$grep.sig.pattern
sig.annotation = tryCatch(eval(parse(text = recipe$sig.annotation)),
                          error = function(e) {return(MultiWAS::return_df(recipe$sig.annotation))})

noperm.pair = unlist(strsplit(recipe$parameters.fiveRankJob, split = ', '))[1]
noperm.str = unlist(strsplit(noperm.pair, split = ' = '))[2]
noperm = eval(parse(text = noperm.str))


setwd(recipe$working.directory)


mylist <- readRDS(ma_paste0(file.path(results.dir, "intermediate.files/signature.inventory.list.RDS")))

# also created in 5_rank_part
types.of.data <- data.table(
  text.pattern = c("trt_cp", 	"trt_sh", "trt_oe", "trt_xpr", "trt_misc", "ctl"),
  data.name    = c("compounds", "shRNA", "over.expression", "CRISPR",
                   "other.treatments", "negative.controls"),
  perturbagen.group = c("drug", "gene", "gene", "gene", NA, NA))

signature.inventory <- fread(ma_paste0(file.path(results.dir, "intermediate.files/signature.inventory.csv")))

message("Collecting the individual signature files")

# results.dir <- '/sc/arion/projects/va-biobank/PROJECTS/marios_temp/parallel_antagonism/working.directories/test_runs/test_AD_fixed_names_V1/results/GTP_CDR/CCL_SCZ_controls/'
path_to_gwas <- list.files(path = ma_paste0(file.path(results.dir, parent.subdir, "/5rank/S01A")), full.names = TRUE)

# this is not really needed, in the original version it is only used to create gwas.model.combos
df <- fread(ma_paste0(file.path(results.dir, "intermediate.files/df.shapes/", thisgwas, '/', paste0(thismodelID, ".shaped.csv.gz"))))# mm

signatures <- list.files(ma_paste0(file.path(results.dir, parent.subdir, "/5rank/S01A/", thisgwas, '/', thismodelID)), full.names = TRUE)
thesesignatures <- do.call(
  rbind,
  pbapply::pblapply(
    stats::setNames(signatures, basename(signatures)),
    FUN = function(i) {
      fread(i)
    }
  )
)

signatures.dir = ma_paste0(file.path(results.dir, parent.subdir, "/signatures/S01A/", thisgwas))
MultiWAS::gv_dir.create(signatures.dir)
fwrite(thesesignatures, ma_paste0(file.path(signatures.dir, '/', paste0(thismodelID, ".signatures.csv.gz"))))

gc()

#   ##############################################################################
#   
#   # Script for averaging the ranking
#
used.types <- unlist(strsplit(grep.sig.pattern, "\\|"))[[1]]  ### only run for trt_cp
used.types <- types.of.data[text.pattern %in% used.types]


for (i in seq(nrow(used.types))) {
  
  # writing code to be able to be run both for groups and individual data.sets in the future.
  message("Loading the master signature file...")
  
  message(paste0("Now working on ", used.types$data.name[i], " (", used.types$text.pattern[i],  ")"))
  
  thesesignatures <- thesesignatures[sig_id %in% unique(signature.inventory[Signature.type %in% used.types$text.pattern[i]]$sig_id) ]
  
  fixed.columns   <- c("gwas", "model_ID", "sig_id")
  
  ### rank extractor ###
  
  rank_extractor <- function(
prefix,  # "" for actual otherwise it is the permutation of interest
z # the signature dataset as above (e.g. thesesignatures) but preferrably split by unique gwas model_ID combinations.
  ) {
    ranks.actual <- data.table(
      "gwas"                      = z$gwas,
      "model_ID"                  = z$model_ID,
      "sig_id"                    = z$sig_ID,
      "rank.est.pearson"          = rank(z[[paste0(prefix, "ALL.cor.pearson")]]),
      "rank.est.spearman"         = rank(z[[paste0(prefix, "ALL.cor.spearman")]]) )
    #       # ks
    cols <- names(z)[grep(paste0(ifelse(prefix == "", "^", prefix), "[[:digit:]]+\\.ks.signed"), names(z))]
    ranks.actual$mean.rank.ks <- rowMeans(z[,..cols][ , (cols) := lapply(.SD, rank), .SDcols = cols])
    #       # Extreme pearson
    cols <- names(z)[grep(paste0(ifelse(prefix == "", "^", prefix), "[[:digit:]]+\\.cor.pearson"), names(z))]
    ranks.actual$mean.rank.extreme.pearson <- rowMeans(z[,..cols][ , (cols) := lapply(.SD, rank), .SDcols = cols])
    #       # Extreme spearman
    cols <- names(z)[grep(paste0(ifelse(prefix == "", "^", prefix), "[[:digit:]]+\\.cor.spearman"), names(z))]
    ranks.actual$mean.rank.extreme.spearman <- rowMeans(z[,..cols][ , (cols) := lapply(.SD, rank), .SDcols = cols])
    #       # Avg rank
    cols <- c("rank.est.pearson", "rank.est.spearman", "mean.rank.ks", "mean.rank.extreme.pearson", "mean.rank.extreme.spearman")
    ranks.actual[[paste0(prefix, "AvgRank")]] <- rowMeans(ranks.actual[,..cols][ , (cols) := lapply(.SD, rank), .SDcols = cols])
    return(rowMeans(ranks.actual[,..cols][ , (cols) := lapply(.SD, rank), .SDcols = cols]))
  }
  
  message(paste0(thisgwas, '_', thismodelID))
  five.rank.all <- do.call(
    cbind,
    lapply(
      stats::setNames(
        unlist(c("", paste0("Perm_", seq(noperm), "_"))),
        unlist(c("Actual", paste0("Perm_", seq(noperm))))),
      rank_extractor,
      z = thesesignatures ) )
  
  #           ### Getting permutation p value
  message("Getting permutation p values...")
  perm.p <- unlist(lapply(
    seq(nrow(five.rank.all)),
    FUN = function(i) {
      sum(
        as.numeric(five.rank.all[i,2:ncol(five.rank.all)]) <=
          as.numeric(five.rank.all[i,1]) ) /
        (ncol(five.rank.all)-1)
    } ))
  #
  #           ### Compiling the per gwas, model_ID and sig_od results.
  final.rank <- data.table(
    "gwas"         = thesesignatures$gwas,
    "model_ID"     = thesesignatures$model_ID,
    "sig_id"       = thesesignatures$sig_id,
    "AvgRank"      = as.numeric(five.rank.all[,1]),
    "perm.p"       = perm.p
  ) # do.call end.
  
  output <- as.data.table(dplyr::left_join(final.rank, sig.annotation))
  
  MultiWAS::gv_dir.create(ma_paste0(file.path(results.dir, results.subdir, '/results', thisgwas)))
  
  message('avgRank calculation successful')
  
  fwrite(
    output,
    ma_paste0(file.path(results.dir, results.subdir, '/results', thisgwas, '/',
           paste0(used.types$text.pattern[i], '_', MultiWAS::make_java_safe(thismodelID), ".signatures.AvgRank.csv.gz"))))
}
      # sink caused 'sink stack is full error'
      #sink( paste0(results.dir, "/AvgRank/", thisgwas, '/readme/', used.types$text.pattern[i], '_', this_modelID, "all.signatures.AvgRank.readme"))
      #print(paste0("Please note that the permutation test is run for each gwas-model combination")) # this is irrelevant
      #











