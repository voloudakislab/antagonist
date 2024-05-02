######################
# Perform antagonism #
######################
# consider using an alias for shortcut.
# submit as coordinator with max walltime
# ~2GB per thread

# STEP 1: run the 5-rank method and permutations
#     Part A: run 5-rank method for each signature batch
#     Part B: join results and split by disease-source combination
#     Part C: Summarize and run permutation analysis for target


###################################################
# R library settings (minerva and package specific)
env  <- "/sc/arion/projects/roussp01a/sanan/Rlibs/230919_R_4.2.0_MultiWAS_Antagonist"
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
  make_option(c("-p", "--prototyping"),
              type="character",
              default="NA",
              help = "Use a number to only run a subset of the signatures, e.g. 2 works well",
              metavar = "character"),
  make_option(c("-d", "--dryrun"),
              type="character",
              default="FALSE",
              help = "If dryrun is TRUE then the scripts are generated but not submitted",
              metavar = "character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
if (is.null(opt$recipe)){
  print_help(opt_parser)
  stop("A recipe file must be provided", call.=FALSE)
}

# Debugging example
# opt <- list()
# opt$recipe <- "/sc/arion/projects/va-biobank/PROJECTS/2023_09_microglia_DGE_gtp_cdr/project.recipe.csv"
# opt$dryrun <- "TRUE"
# opt$prototyping <- "TRUE"



#######################
# Parse the recipe file
recipe <- parse_recipe(opt$recipe)
#  need to get directory information
setwd(recipe$working.directory)
# Perform first step of the function

########################
# Run the wrapper script
## I can add cluster files instead of specifying here
## Alternatively here use a compacted list so that changes in the main package can be propagated more easily.
perform_antagonism_lsf_S01_wrapper(
  # wrapper input
  working.directory      = recipe$working.directory,
  recipe.file            = opt$recipe,
  dryrun                 = eval(parse(text = opt$dryrun)),
  # Input
  df                     = recipe$df,
  column.feature         = recipe$column.feature,
  column.statistic       = recipe$column.statistic,
  column.trait           = recipe$column.trait,
  column.source          = recipe$column.source,
  # Output
  results.dir            = recipe$results.dir,
  # Parameters
  n.threads              = recipe$n.threads,
  signature.dir          = recipe$signature.dir,
  gene.anno.file         = recipe$gene.anno.file,
  grep.sig.pattern       = recipe$grep.sig.pattern,
  noperm                 = recipe$noperm,
  thres.N.vector         = recipe$thres.N.vector,
  sig.annotation         = recipe$sig.annotation,
  overwrite.intermediate = recipe$overwrite.intermediate,
  model.banlist.grep     = recipe$model.banlist.grep,
  prototyping            = eval(parse(text = opt$prototyping))
)


