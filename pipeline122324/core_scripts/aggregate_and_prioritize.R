############ LOAD PACKAGES

libs <- .libPaths()
libs[3] <- "/sc/arion/projects/roussp01a/sanan/Rlibs/231221_R_4.2.0_MultiWAS"
libs[4] <- "/sc/arion/projects/va-biobank/software/Georgios_dev/240702_R_4.2.0_MultiWAS_Antagonist/"
.libPaths(libs)
# Load MultiWAS
library(MultiWAS)
library(antagonist)

######################
# LOAD VARIABLES from shell

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

############################################################ 




# DEBUG variables
#opt <- list()
#opt$recipe.file <- '/sc/arion/projects/va-biobank/PROJECTS/marios_temp/parallel_antagonism/recipes/test_run_recipes/bash_test/S4_test_run_V3.csv'
#opt$thisgwas = 'AD'
#opt$thismodelID = 'Microglia_Genes_.._Microglia_FANS_.._Genes_.._PrediXcan_.._EUR'
#opt$results.dir = 'results/GTP_CDR'
#opt$results.subdir = '/wilcoxRank'
#opt$parent.subdir = '/intermediate.files/avgRank'

# end of DEBUG
#############################

# Default variables
# these can be handled from recipe if needed

# Unload recipe
recipe <- parse_recipe(opt$recipe.file)

#Unload opt list
thisgwas <- opt$thisgwas
thismodelID <- opt$thismodelID
results.dir <- opt$results.dir
results.subdir <- opt$results.subdir
parent.subdir <- opt$parent.subdir
###################

working.directory <- recipe$working.directory 
setwd(working.directory)

dfs.to.process = ma_paste0(file.path(results.dir, parent.subdir, '/results', thisgwas, thismodelID))
limit.dfs.to = recipe$grep.sig.pattern

#dfs.to.process = "results/GTP_CDR/.+_.+_all\\.signatures\\.AvgRank\\.csv\\.gz$"
#limit.dfs.to = "trt_cp|trt_sh|trt_oe|trt_xpr"
ref.drug = TRT_CP.INFO.20200324
ref.cell = CELL.LINE.INFO
min.experiments.n = 2
n.cores = parallel::detectCores() - 2

output.dir <- ma_paste0(file.path(getwd(), results.dir, results.subdir, '/results', thisgwas, thismodelID))


gc()
"%!in%" <- function(x, y) !(x %in% y)
ref.drug <- MultiWAS::return_df(ref.drug)
ref.cell <- MultiWAS::return_df(ref.cell)
message("Loading the data...")

# splits df to x.cdr and x.gtp based on their name (they are stored as /gwas/modelID/file.prefix.csv.gz)
dfs.to.process <- list.files(path = dirname(dfs.to.process), 
                             pattern = basename(dfs.to.process), full.names = T)
dfs.to.process <- dfs.to.process[grep(limit.dfs.to, dfs.to.process)]
#gwass <- unique(MultiWAS::return_df(dfs.to.process[1])$gwas)

# x.cdr and x.gtp are dataframes of ranked drug signatures, ex: output of either avgRank or VAE
x.cdr <- dfs.to.process[grep("trt_cp", dfs.to.process)]
x.gtp <- dfs.to.process[grep("trt_sh|trt_oe|trt_xpr", dfs.to.process)]

if (length(x.cdr) != 0) {
  x.cdr <- MultiWAS::return_df(x.cdr)
  x.cdr <- as.data.table(dplyr::left_join(x.cdr, ref.drug))
  x.cdr$readily.repurposable <- !is.na(x.cdr$clinical_phase)
  x.cdr[, `:=`(N_experiments, length(AvgRank)), by = pert_iname]
  x.cdr <- x.cdr[N_experiments >= min.experiments.n]
}
if (length(x.gtp) != 0) {
  x.gtp <- do.call(rbind, lapply(x.gtp, FUN = function(z) {
    MultiWAS::return_df(z)
  }))
  x.gtp[, `:=`(N_experiments, length(AvgRank)), by = pert_iname]
  x.gtp <- x.gtp[N_experiments >= min.experiments.n]
}

# mylist may contain x.cdr and/or x.gtp if each one of them are NOT NULL
mylist = lapply(list(x.cdr, x.gtp), function(x) if(!is.null(nrow(x))) return(x))
mylist = mylist[sapply(seq_along(mylist), function(i) ifelse(is.null(mylist[[i]]), FALSE, TRUE))]

if(length(mylist) == 0) stop('After trying to load the dataframes from: ',
                             ma_paste0(file.path(results.dir, results.subdir, '/results', thisgwas)),
                             ' for model_ID: "', thismodelID, '" using the recipe-provided grep.sig.pattern = "',
                             limit.dfs.to, '" , I grouped them in a list but noticed that it has length 0. 
                             \nConsider checking:\n1. dataframes exist\n
                             2. their number of rows is > 0\n
                             3. "grep.sig.pattern" in recipe makes sense\n
                             Have a good one!')

lapply(mylist, FUN = function(x) {
      
      #if (!is.na(limit.models.to[1])) x <- x[model_ID %in% limit.models.to]
      #modelprefix <- as.character(MultiWAS::make_java_safe(paste(limit.models.to, 
      #                                                           collapse = "_X_")))
      #if (modelprefix == "NA") modelprefix <- "ALL"
      
      MultiWAS::gv_dir.create(output.dir)
      
      d.plot <- ggpubr::ggdensity(data.frame(Type = paste(x$pert_type, 
                                                          collapse = "|"), AvgRank = as.numeric(x$AvgRank)), 
                                  x = "AvgRank", add = "mean", xlab = "Signature AvgRank")
      cowplot::ggsave2(paste0(output.dir, "/", ifelse(unique(x$pert_type)[1] == 
                                                             "trt_cp", "cdr", "gtp"), "_AvgRank_distribution_square.pdf"), 
                       d.plot, height = 4, width = 4)
      cowplot::ggsave2(paste0(output.dir, "/", ifelse(unique(x$pert_type)[1] == "trt_cp", "cdr", "gtp"), "_AvgRank_distribution_landscape.pdf"), 
                       d.plot, height = 3, width = 8)
      
      ########################
      # Calculating statistics
      # hist(x$AvgRank) # Seems to be following the gaussian distribution
      # avg.rank.mean <- mean(x$AvgRank)
      start_time = Sys.time()
      avg.rank.sd   <- sd(x$AvgRank)
      message("Calculating the MW statistic (WRST)...")
      
      # Calculate MW statistic (WRST)
      x <- as.data.table(dplyr::left_join(
        x,
        do.call(
          rbind,
          pbmclapply(
            unique(x$pert_iname),
            FUN = function(i){
              mw.res <- wilcox.test(
                x[pert_iname == i]$AvgRank,
                x[pert_iname != i]$AvgRank,
                paired = FALSE,
                conf.int = T) # conf.int = T is required to get the estimate
              data.frame(
                pert_iname             = i,
                Compound.MW.p          = mw.res$p.value,
                Compound.pseudo.zscore = -1*(mw.res$estimate/avg.rank.sd), # the estimate is a deviation from the median/mean AvgRank and we are dividing that by the sd of AvgRank
                Compound.HL.estimate   = mw.res$estimate,
                Compound.avg.rank.sd   = avg.rank.sd
              )
            }, mc.cores = n.cores ))
      ))
      fwrite(x, paste0(output.dir, "/", ifelse(unique(x$pert_type)[1] == "trt_cp", 
                                                    "cdr", "gtp"), "_signature_level.csv.gz"))
      
      end_time = Sys.time()
      elapsed_time = end_time - start_time
      message(paste0('wrote _signature_level.csv.gz ; Elapsed time: ', elapsed_time))
      
      # Aggregate at the level of compound
      x$Rank <- as.numeric(NA) # reset
      x[, AvgRank := mean(AvgRank), by = pert_iname]
      x[, perm.p.all:= paste(sort(perm.p), collapse=";"), by = pert_iname ] # save permutation values
      
      
      if (unique(x$pert_type)[1] == "trt_cp") {
        ## CDR specific code
        compound.level = unique(x[, c(
          "gwas",
          "pert_iname", "Rank", "AvgRank", "Compound.MW.p", "Compound.pseudo.zscore",
          "Compound.HL.estimate", "Compound.avg.rank.sd",
          "clinical_phase",   "moa", "target", "disease_area", "indication",
          "N_experiments", "perm.p.all")])
      } else {
        compound.level = unique(x[, c(
          "gwas",
          "pert_iname", "Rank", "AvgRank", "Compound.MW.p",
          "Compound.pseudo.zscore", "Compound.HL.estimate", "Compound.avg.rank.sd",
          "N_experiments", "perm.p.all")])
        ## GTP specific code
      }
      
      compound.level$Rank <- rank(compound.level$AvgRank) # replace Rank
      compound.level      <- compound.level[order(Rank)]
      
      if (unique(x$pert_type)[1] == "trt_cp") {
        ## CDR specific code for moa (Mechanism of action)
        
        #### MOA ####
        # Aggregate at the level of moa
        
        start_time = Sys.time()
        
        moa.repurp      <- compound.level[, c("moa", "AvgRank")]
        # https://stackoverflow.com/questions/61684587/r-data-table-split-a-row-into-multiple-rows-based-on-string-values
        # Exclude those that don't have moa
        moa.repurp <- moa.repurp[!is.na(moa)]
        moa.repurp      = moa.repurp[moa != "unknown"]
        # Break multiple moas
        # moa.repurp <- moa.repurp[, .(moa = unlist(tstrsplit(moa, "\\|", type.convert = TRUE))), by = "AvgRank"] doesn't handle NULL well
        #moa.repurp <- moa.repurp[, .(moa = {
        #  allmoa <- unlist(tstrsplit(moa, "\\|", type.convert = TRUE))
        #  class(allmoa) <- 'character' # to handle NULL and unable to deduct class issue
        #  allmoa
        #}), by = "AvgRank"] # to handle NULLs
        moa.repurp <- moa.repurp[, .(
          moa = unlist(strsplit(moa, "\\|")),
          AvgRank = rep(AvgRank, sapply(strsplit(moa, "\\|"), length))
        )]
        moa.repurp[, N_compounds_moa := length(AvgRank), by = moa]
        moa.repurp      = moa.repurp[N_compounds_moa > 1] # only if there are more than 2
        moa.avg.rank.sd = sd(moa.repurp$AvgRank)
        # 96 unique MOA
        
        message("Running mechanism of action enrichment...")
        moa.repurp <- as.data.table(dplyr::left_join(
          moa.repurp,
          do.call(rbind,pbmclapply(
            unique(moa.repurp$moa),
            FUN = function(i){
              mw.res <- wilcox.test(
                moa.repurp[moa == i]$AvgRank,
                moa.repurp[moa != i]$AvgRank,
                paired = FALSE,
                conf.int = T) # conf.int = T is required to get the estimate
              data.frame(
                moa               = i,
                MOA.MW.p          = mw.res$p.value,
                MOA.pseudo.zscore = -1*(mw.res$estimate/moa.avg.rank.sd), # the estimate is a deviation from the median/mean AvgRank and we are dividing that by the sd of AvgRank
                MOA.HL.estimate   = mw.res$estimate,
                MOA.avg.rank.sd   = avg.rank.sd
              )
            }, mc.cores = detectCores()-2 ))
        ))
        moa.repurp[, AvgRank_moa := mean(AvgRank), by = moa]
        moa.repurp = unique(moa.repurp[
          , c("moa", "AvgRank_moa", "MOA.pseudo.zscore",
              "MOA.HL.estimate", "MOA.avg.rank.sd",
              "N_compounds_moa", "MOA.MW.p")])
        moa.repurp$Rank_moa <- rank(moa.repurp$AvgRank_moa)
        # moa.repurp$Rank_moa_percentile <- my_percent(moa.repurp$Rank_moa_percentile/length(unique(moa.repurp$moa)))
        # moa.repurp$Rank_moa_percentile <- moa.repurp$Rank_moa_percentile/length(unique(moa.repurp$moa))
        moa.repurp <- moa.repurp[order(Rank_moa)]
        moa.repurp$MOA.MW.FDR <- p.adjust(moa.repurp$MOA.MW.p, method = "fdr")
        fwrite(moa.repurp, paste0(output.dir, 
                                  "/", thisgwas, "_", ifelse(unique(x$pert_type)[1] == 
                                                               "trt_cp", "cdr", "gtp"), "_moa_level.csv"))
        
        end_time = Sys.time()
        elapsed_time = end_time - start_time
        message(paste0('Wrote _moa_level.csv ; Time elapsed from start to finish: ', elapsed_time))
        
        start_time = Sys.time()
        
        #### TARGET ####
        # Aggregate at the level of target
        target.repurp      <- compound.level[, c("target", "AvgRank")]
        # https://stackoverflow.com/questions/61684587/r-data-table-split-a-row-into-multiple-rows-based-on-string-values
        # Exclude those that don't have target
        target.repurp <- target.repurp[!is.na(target)]
        target.repurp      = target.repurp[target != "unknown"]
        # target.repurp <- target.repurp[, .(target = unlist(tstrsplit(target, "\\|", type.convert = TRUE))), by = "AvgRank"]
        #target.repurp <- target.repurp[, .(target = {
        #  alltargets <- unlist(tstrsplit(target, "\\|", type.convert = TRUE))
        #  class(alltargets) <- 'character' # to handle NULL and unable to deduct class issue
        #  alltargets
        #}), by = "AvgRank"] # to handle NULLs
        target.repurp <- target.repurp[, .(
          target = unlist(strsplit(target, "\\|")),
          AvgRank = rep(AvgRank, sapply(strsplit(target, "\\|"), length))
        )]
        target.repurp[, N_compounds_target := length(AvgRank), by = target]
        target.repurp      = target.repurp[N_compounds_target > 1] # only if there are more than 2
        target.avg.rank.sd = sd(target.repurp$AvgRank)
        
        message("Running target enrichment...")
        target.repurp <- as.data.table(dplyr::left_join(
          target.repurp,
          do.call(rbind,pbmclapply(
            unique(target.repurp$target),
            FUN = function(i){
              mw.res <- wilcox.test(
                target.repurp[target == i]$AvgRank,
                target.repurp[target != i]$AvgRank,
                paired = FALSE,
                conf.int = T) # conf.int = T is required to get the estimate
              data.frame(
                target               = i,
                target.MW.p          = mw.res$p.value,
                target.pseudo.zscore = -1*(mw.res$estimate/target.avg.rank.sd), # the estimate is a deviation from the median/mean AvgRank and we are dividing that by the sd of AvgRank
                target.HL.estimate   = mw.res$estimate,
                target.avg.rank.sd   = avg.rank.sd
              )
            }, mc.cores = detectCores()-2 ))
        ))
        target.repurp[, AvgRank_target := mean(AvgRank), by = target]
        target.repurp = unique(target.repurp[
          , c("target", "AvgRank_target", "target.pseudo.zscore",
              "target.HL.estimate", "target.avg.rank.sd",
              "N_compounds_target", "target.MW.p")])
        target.repurp$Rank_target <- rank(target.repurp$AvgRank_target)
        target.repurp <- target.repurp[order(Rank_target)]
        target.repurp$target.MW.FDR <- p.adjust(target.repurp$target.MW.p, method = "fdr")
        fwrite(target.repurp, paste0(output.dir, 
                                     "/", thisgwas, "_", ifelse(unique(x$pert_type)[1] == 
                                                                  "trt_cp", "cdr", "gtp"), "_target_level.csv"))
        end_time = Sys.time()
        elapsed_time = end_time - start_time
        message(paste0('Wrote _target_level.csv ; Time elapsed from start to finish: ', elapsed_time))
        
        start_time = Sys.time()
        
        #### disease_area ####
        # Aggregate at the level of disease_area
        disease_area.repurp      <- compound.level[, c("disease_area", "AvgRank")]
        # https://stackoverflow.com/questions/61684587/r-data-table-split-a-row-into-multiple-rows-based-on-string-values
        # Exclude those that don't have disease_area
        disease_area.repurp <- disease_area.repurp[!is.na(disease_area)]
        disease_area.repurp      = disease_area.repurp[disease_area != "unknown"]
        # disease_area.repurp <- disease_area.repurp[, .(disease_area = unlist(tstrsplit(disease_area, "\\|", type.convert = TRUE))), by = "AvgRank"]
        #disease_area.repurp <- disease_area.repurp[, .(disease_area = {
        #  disease_areas <- unlist(tstrsplit(disease_area, "\\|", type.convert = TRUE))
        #  class(disease_areas) <- 'character' # to handle NULL and unable to deduct class issue
        #  disease_areas
        #}), by = "AvgRank"] # to handle NULLs
        ### Zhenyi fixed:
        disease_area.repurp <- disease_area.repurp[, .(
          disease_area = unlist(strsplit(disease_area, "\\|")),
          AvgRank = rep(AvgRank, sapply(strsplit(disease_area, "\\|"), length))
        )]
        
        disease_area.repurp[, N_compounds_disease_area := length(AvgRank), by = disease_area]
        disease_area.repurp      = disease_area.repurp[N_compounds_disease_area > 1] # only if there are more than 2
        disease_area.avg.rank.sd = sd(disease_area.repurp$AvgRank)
        
        message("Running disease area enrichment...")
        disease_area.repurp <- as.data.table(dplyr::left_join(
          disease_area.repurp,
          do.call(rbind,pbmclapply(
            unique(disease_area.repurp$disease_area),
            FUN = function(i){
              mw.res <- wilcox.test(
                disease_area.repurp[disease_area == i]$AvgRank,
                disease_area.repurp[disease_area != i]$AvgRank,
                paired = FALSE,
                conf.int = T) # conf.int = T is required to get the estimate
              data.frame(
                disease_area               = i,
                disease_area.MW.p          = mw.res$p.value,
                disease_area.pseudo.zscore = -1*(mw.res$estimate/disease_area.avg.rank.sd), # the estimate is a deviation from the median/mean AvgRank and we are dividing that by the sd of AvgRank
                disease_area.HL.estimate   = mw.res$estimate,
                disease_area.avg.rank.sd   = avg.rank.sd
              )
            }, mc.cores = detectCores()-2 ))
        ))
        disease_area.repurp[, AvgRank_disease_area := mean(AvgRank), by = disease_area]
        disease_area.repurp = unique(disease_area.repurp[
          , c("disease_area", "AvgRank_disease_area", "disease_area.pseudo.zscore",
              "disease_area.HL.estimate", "disease_area.avg.rank.sd",
              "N_compounds_disease_area", "disease_area.MW.p")])
        disease_area.repurp$Rank_disease_area <- rank(disease_area.repurp$AvgRank_disease_area)
        disease_area.repurp <- disease_area.repurp[order(Rank_disease_area)]
        disease_area.repurp$disease_area.MW.FDR <- p.adjust(disease_area.repurp$disease_area.MW.p, method = "fdr")
        
        fwrite(disease_area.repurp, paste0(output.dir, 
                                           "/", thisgwas, "_", ifelse(unique(x$pert_type)[1] == 
                                                                        "trt_cp", "cdr", "gtp"), "_disease_area_level.csv"))
        end_time = Sys.time()
        elapsed_time = end_time - start_time
        message(paste0('Wrote _disease_area_level.csv ; Time elapsed from start to finish: ', elapsed_time))
        
        start_time = Sys.time()
        
        #### indication ####
        # Aggregate at the level of indication
        indication.repurp      <- compound.level[, c("indication", "AvgRank")]
        # https://stackoverflow.com/questions/61684587/r-data-table-split-a-row-into-multiple-rows-based-on-string-values
        # Exclude those that don't have indication
        indication.repurp <- indication.repurp[!is.na(indication)]
        indication.repurp      = indication.repurp[indication != "unknown"]
        # indication.repurp <- indication.repurp[, .(indication = unlist(tstrsplit(indication, "\\|", type.convert = TRUE))), by = "AvgRank"]
        indication.repurp <- indication.repurp[, .(
          indication = unlist(strsplit(indication, "\\|")),
          AvgRank = rep(AvgRank, sapply(strsplit(indication, "\\|"), length))
        )] # to handle NULLs
        indication.repurp[, N_compounds_indication := length(AvgRank), by = indication]
        indication.repurp      = indication.repurp[N_compounds_indication > 1] # only if there are more than 2
        indication.avg.rank.sd = sd(indication.repurp$AvgRank)
        
        message("Running indication enrichment...")
        indication.repurp <- as.data.table(dplyr::left_join(
          indication.repurp,
          do.call(rbind,pbmclapply(
            unique(indication.repurp$indication),
            FUN = function(i){
              mw.res <- wilcox.test(
                indication.repurp[indication == i]$AvgRank,
                indication.repurp[indication != i]$AvgRank,
                paired = FALSE,
                conf.int = T) # conf.int = T is required to get the estimate
              data.frame(
                indication               = i,
                indication.MW.p          = mw.res$p.value,
                indication.pseudo.zscore = -1*(mw.res$estimate/indication.avg.rank.sd), # the estimate is a deviation from the median/mean AvgRank and we are dividing that by the sd of AvgRank
                indication.HL.estimate   = mw.res$estimate,
                indication.avg.rank.sd   = avg.rank.sd
              )
            }, mc.cores = detectCores()-2 ))
        ))
        indication.repurp[, AvgRank_indication := mean(AvgRank), by = indication]
        indication.repurp = unique(indication.repurp[
          , c("indication", "AvgRank_indication", "indication.pseudo.zscore",
              "indication.HL.estimate", "indication.avg.rank.sd",
              "N_compounds_indication", "indication.MW.p")])
        indication.repurp$Rank_indication <- rank(indication.repurp$AvgRank_indication)
        indication.repurp <- indication.repurp[order(Rank_indication)]
        indication.repurp$indication.MW.FDR <- p.adjust(indication.repurp$indication.MW.p, method = "fdr")
        fwrite(indication.repurp, paste0(output.dir, "/", thisgwas, "_",
                                         ifelse(unique(x$pert_type)[1] == "trt_cp", "cdr","gtp"),
                                         "_indication_level.csv"))
        fwrite(indication.repurp, paste0(output.dir, 
                                         "/", thisgwas, "_", ifelse(unique(x$pert_type)[1] == 
                                                                      "trt_cp", "cdr", "gtp"), "_indication_level.csv"))
        
        end_time = Sys.time()
        elapsed_time = end_time - start_time
        message(paste0('Wrote _indication_level.csv ; Time elapsed from start to finish: ', elapsed_time))
        
        # Compile and save
        # This is imperfect for compounds that have more than one moa
        compound.level <- as.data.table(dplyr::left_join(
          compound.level, moa.repurp[, c(
            "moa", "Rank_moa", "MOA.pseudo.zscore",
            "MOA.HL.estimate", "MOA.avg.rank.sd",
            "MOA.MW.p", "MOA.MW.FDR")]))
        fwrite(compound.level, paste0(output.dir, 
                                      "/", thisgwas, "_", ifelse(unique(x$pert_type)[1] == 
                                                                   "trt_cp", "cdr", "gtp"), "_all_compound_level.csv"))
        message('wrote _all_compound_level')
        
        # Filter to only Phase 3 and lauched
        compound.level$Compound.MW.FDR <- as.numeric(NA)
        compound.level[clinical_phase %in% c("Phase 3", "Launched")]$Compound.MW.FDR <-
          p.adjust(compound.level[clinical_phase %in% c("Phase 3", "Launched")]$Compound.MW.p, method = "fdr")
        # compound.level$Compound.MW.FDR <- p.adjust(compound.level$Compound.MW.p, method = "fdr")
        compound.level <- compound.level[
          , c("pert_iname", "clinical_phase", "Rank", "AvgRank", "Compound.MW.p",
              "Compound.MW.FDR", "Compound.pseudo.zscore",
              "Compound.HL.estimate", "Compound.avg.rank.sd",
              "moa", "Rank_moa", "MOA.pseudo.zscore",
              "MOA.HL.estimate", "MOA.avg.rank.sd",
              "MOA.MW.p", "MOA.MW.FDR", "target", "disease_area",
              "indication", "N_experiments", "perm.p.all")]
        
        fwrite(compound.level, paste0(output.dir, 
                                      "/", thisgwas, "_", ifelse(unique(x$pert_type)[1] == 
                                                                   "trt_cp", "cdr", "gtp"), "_P3_and_launched_compound_level.csv"))
        
        
        message(paste0('Wrote _P3_and_launched_compound_level.csv'))
        
        compound.level <- compound.level[clinical_phase=="Launched"] # 201 compounds
        compound.level$Rank <- rank(compound.level$Rank) # rerank
        column.order <- names(compound.level)
        moa.ranks <- unique(compound.level[,c("moa", "Rank_moa")])
        moa.ranks <- moa.ranks[!is.na(Rank_moa)]
        moa.ranks$MOA.Rank <- rank(moa.ranks$Rank_moa) # rerank
        moa.ranks[, Rank_moa := NULL]
        compound.level[, Rank_moa := NULL] # remove
        compound.level <- as.data.table(dplyr::left_join(compound.level, moa.ranks))
        fwrite(compound.level, paste0(output.dir, "/", thisgwas, "_",
                                      ifelse(unique(x$pert_type)[1] == "trt_cp", "cdr","gtp"),
                                      "_launched_compound_level.csv"))
        fwrite(compound.level, paste0(output.dir, 
                                      "/", thisgwas, "_", ifelse(unique(x$pert_type)[1] == 
                                                                   "trt_cp", "cdr", "gtp"), "_launched_compound_level.csv"))
        
        message(paste0('Wrote _launched_compound_level.csv'))
      }
      else{
        compound.level$Rank <- rank(compound.level$AvgRank)
        compound.level$Compound.MW.FDR <- p.adjust(compound.level$Compound.MW.p, 
                                                   method = "fdr")
        compound.level <- compound.level[order(Rank)][, 
                                                      c("pert_iname", "Rank", "AvgRank", "Compound.MW.p", 
                                                        "Compound.pseudo.zscore", "Compound.MW.FDR", 
                                                        "N_experiments", "perm.p.all")]
        fwrite(compound.level, paste0(output.dir, "/", ifelse(unique(x$pert_type)[1] == 
                                                                     "trt_cp", "cdr", "gtp"), "_compound_level.csv"))
      }
    
  })

