#' Antagonism core script
#'
#' Called by:
#' - exec/antagonist_S01_wrapper.R
#' Calls:
#' - exec/antagonist_S01PA_core.R
#' This is the first part of step 1 of the pipeline. The file with the disease
#' signature is provided and may have multiple traits (e.g. TWASs) and sources
#' (e.g. tissues). The parallelization at this step is at the level of the
#' signature.
#'
#'General information below:
#' This is intended to run in an LSF custer. There is a version that is optimized for running locally.
#' Runs the 5 methods for each GWAS / model_ID combination.
#' The following 5 variables/methods need to be saved in splitting the drug matrices:
#' 1. cor.pearson : Pearson (all)
#' 2. cor.spearman : Spearman (all)
#' 3. extreme.cor.spearman : Spearman (most differentially expressed genes)
#' 4. extreme.cor.pearson : Pearson (most differentially expressed genes)
#' 5. ks.signed : KS method (most differentially expressed genes)
#' This is a multithreaded data.table based re-write of So HC et al. Analysis of
#'  genome-wide association data highlights candidates for drug repositioning in
#'  psychiatry. Nat Neurosci. 2017 Oct;20(10):1342-1349. PMID: 28805813.
#' It has been adapted to run for signatures (individual parameters) vs.
#' summarized compound scores.
#'
#' @param working.directory provide working directory
#' @param recipe.file provide recipe.file location
#' @param dryrun if true the jobs are not submitted
#' @param df twas data frame
#' @param column.feature expected input is ENSEMBL ID and currently only works with genes (default: feature)
#' @param column.statistic z-score statistic (default: zscore)
#' @param column.trait trait name (default: gwas)
#' @param column.source disease signature source, e.g. microglia, DLPFC, meta-analysis, etc. (default: model_ID)
#' @param results.dir Results main dir. Subdirectories will be created
#' @param n.threads number f threads for multicore strategy, default is detectCores()-2
#' @param signature.dir Location of signature files. Consider saving locally to speed things up (e.g. /scratch/cmap_l1000_2021_01_28/). Look into the split_gctx on how to generate these files.
#' @param gene.anno.file Gene annotation file (provided from L1000)
#' @param grep.sig.pattern Regular exrpression pattern to grep from signature name files. Default is to process sets for computational drug repurposing and gene target prioritization.
#' @param noperm Number of permutations (the final number of permutation data points is noperm × number of drugs)
#' @param thres.N.vector defining the threshold(s) K such that only top K items are included for comparison of expressions: use the same thresholds as So et al.
#' @param sig.annotation signature annotation. Default is SIG.INFO.20211120
#' @param overwrite.intermediate default is not to overwrite the all.signatures.csv.gz (FALSE). This doesn't change from analysis to analysis.
#' @param model.banlist.grep regular expression to exclude non-gene-level models. Default is: ":: Transcripts ::|:: H3K4me3 ::|:: H3K27ac ::|:: CA ::"
#' @param prototyping number of signature files (each one has around 300) to test again for pipeline testing (default is NA which is all)
#' @return N/A. Saves required files
#' @keywords antagonism step1 LSF
#' @export
#'
perform_antagonism_lsf_S01_wrapper <- function(
  working.directory      ,
  recipe.file            ,
  dryrun                = FALSE,
  # Input
  df                     ,
  column.feature         = "feature",
  column.statistic       = "zscore",
  column.trait           = "gwas",
  column.source          = "model_ID",
  # Output
  results.dir            = "results/GTP_CDR/",
  # Parameters
  n.threads              = parallel::detectCores()-2,
  signature.dir          = "/sc/arion/projects/va-biobank/software/Georgios_dev/forJamie/antagonist_files/drugs", # /scratch/cmap_l1000_2021_01_28/
  gene.anno.file         = "/sc/arion/projects/va-biobank/resources/CMap/cmap_l1000_2021-11-20/geneinfo_beta.txt",
  grep.sig.pattern       = "trt_cp|trt_sh|trt_oe|trt_xpr",
  noperm                 = 100, #
  thres.N.vector         = c(50,100,250,500),
  sig.annotation         = SIG.INFO.20211120,
  overwrite.intermediate = TRUE,
  model.banlist.grep     = ":: Transcripts ::|:: H3K4me3 ::|:: H3K27ac ::|:: CA ::",
  prototyping            = NA
) {


  ######################
  # HARDCODED PARAMETERS
  ## Look here if the script fails
  setwd(working.directory)
  R.lib.dir <- "/sc/arion/projects/va-biobank/software/Georgios_dev/240702_R_4.2.0_MultiWAS_Antagonist/"
  types.of.data <- data.table(
    text.pattern = c("trt_cp", 	"trt_sh", "trt_oe", "trt_xpr", "trt_misc", "ctl"),
    data.name    = c("compounds", "shRNA", "over.expression", "CRISPR",
                        "other.treatments", "negative.controls"),
    perturbagen.group = c("drug", "gene", "gene", "gene", NA, NA))

  ##############################################################################
  # RUN STATISTICS FOR ALL SIGNATURES ##########################################
  ##############################################################################

  #######################
  # Preparing Directories
  gv_dir.create(results.dir)
  gv_dir.create(paste0(results.dir, "/intermediate.files/"))
  gv_dir.create(paste0(results.dir, "/intermediate.files/5rank/S01A"))
  gv_dir.create(paste0(results.dir, "/intermediate.files/scripts/S01A"))
  gv_dir.create(paste0(results.dir, "/intermediate.files/logs/S01A"))
  gv_dir.create(paste0(results.dir, "/intermediate.files/5rank/S01B"))
  gv_dir.create(paste0(results.dir, "/intermediate.files/scripts/S01B"))
  gv_dir.create(paste0(results.dir, "/intermediate.files/logs/S01B"))
  # Preparing TWAS
  df               <- return_df(df)
  # Preparing columns
  df$feature   <- df[[column.feature]]
  df$zscore    <- df[[column.statistic]]
  df$gwas      <- df[[column.trait]]
  df$model_ID  <- df[[column.source]]
  df           <- df[zscore!=-Inf & zscore!=Inf]
  # gwas.model.combs <- unique(df[, c("gwas","model_ID") ])
  df <- df[, c("feature", "zscore", "gwas", "model_ID")]
  # Keeping only genes:
  if (length(grep(model.banlist.grep, df$model_ID))>0) {
    message(paste0("Currently expression banlist is: ", model.banlist.grep))
    message("Intention is to only keep genes")
    df <- df[-grep(model.banlist.grep, model_ID)] }
  fwrite(df, paste0(results.dir, "intermediate.files/df.shaped.csv.gz"))

  ############################################
  # Build project-specific signature inventory
  # Identifying signatures
  cmap.list <- list.files(signature.dir, full.names = T)
  # cmap=readRDS("/sc/hydra/projects/roussp01b/Wen/LINCS_Level5/eachDrugPhase2/1-300_cp.RDS")
  # to.process <- as.data.table(tidyr::crossing(to.process, data.table("cmap.file" = cmap.list)))
  if (!is.na(grep.sig.pattern[1])) mylist <- cmap.list[grep(grep.sig.pattern, cmap.list)] else mylist <- cmap.list
  # Create a signature inventory
  if (!file.exists(paste0(results.dir, "intermediate.files/signature.inventory.csv")) | overwrite.intermediate) {
  message("Building signature inventory by perturbagen type")
  signature.inventory <- pbmclapply(
    stats::setNames(mylist, sub("_[[:digit:]]+\\.[[:digit:]]+\\.RDS", "", basename(mylist)) ),
    FUN = function(x) {colnames(readRDS(x))},
    mc.cores = n.threads
  )
  message('First signature.inventory works')
  signature.inventory <- do.call(rbind, lapply(
    unique(unlist(names(signature.inventory))),
    FUN = function(x) {
      data.table(
        "Signature.type" = x,
        "sig_id" = unique(unlist(signature.inventory[grep(x, names(signature.inventory))]))
      ) } ) )
  fwrite(signature.inventory, paste0(results.dir, "intermediate.files/signature.inventory.csv"))
  } else {
    message("Loading signature inventory by perturbagen type")
    signature.inventory <- fread(paste0(results.dir, "intermediate.files/signature.inventory.csv"))
  }
  if (!is.na(prototyping)) {
    mylist <- unlist(
      lapply(
        unique(signature.inventory$Signature.type),
        FUN = function(x) {
          mylist[grep(x, mylist)][seq(prototyping)]
    }))
  }
  saveRDS(mylist, paste0(results.dir, "intermediate.files/signature.inventory.list.RDS"))

################################################################################


  message("Submitting jobs for running five methods across each signature")
  if (!file.exists(paste0(paste0(results.dir, "intermediate.files/all.signatures.csv.gz"))) | overwrite.intermediate) {
    # FIXME: this script only outputs one model ID..
    # Multicore strategy split by RDS file
    pbapply::pblapply(
      stats::setNames(mylist, basename(mylist)),
      FUN = function(i) {
        # Pass recipe file and signature file.
        # add ml R
        job.name   <- paste0(sub("\\.RDS$", "", basename(i)))
        ## Populate the files.info with this new information
        ### always remove the old outputs and generate new ones for rbinding the signatures.
        # Usually takes 30' for about 50 gwas - model combinations
        b.sub <- paste0('bsub -P acc_va-biobank -q premium -n ', n.threads,
                        ' -W 2:00 -J ', job.name, ' -R span[hosts=1] -R rusage[mem=3000] -oo '
                        ,results.dir, '/intermediate.files/logs/S01A/', job.name, '.out ',
                        '-eo ',results.dir, '/intermediate.files/logs/S01A/', job.name , '.err ',
                        '-L /bin/bash', ' < ')
        # prepare script
        thiscommand <- paste0(
          "cd ", working.directory, "\n",
          "ml R", "\n",
          "Rscript --verbose  ", R.lib.dir, "/antagonist/exec/antagonist_S01PA_core.R", " ",
          "--recipe ", recipe.file, " ",
          "--cmapfile ", i
          )
        writeLines(
          paste0(
            "#!/bin/bash", "\n",
            thiscommand),
          paste0(results.dir, "/intermediate.files/scripts/S01A/",
                 job.name, ".sh")
          )
        submission.command <- paste0(
          b.sub, results.dir, "/intermediate.files/scripts/S01A/", job.name, ".sh")

        if (!dryrun) { # if dry run then don't execute the scripts
          system(submission.command)
          } else {
            writeLines(
              submission.command,
              paste0(results.dir, "/intermediate.files/scripts/S01A/",
                     job.name, ".dryrun.txt") )
            }
        #while (!file.exists(signature.file)) {
        #  Sys.sleep(10)  # Check every 10 seconds
        #}
        #message(paste0('Signature file for job', job.name, ' is finished'))
        #signature <- fread(signature.file)
        #return(signature)
      #sink()
        } # signature loop ends
      ) # pblapply ends
  } # this only runs if the output doesn't exist

  ##############################################################################
  # Collect all the signatures from 5-rank method

  ###
  message("Waiting for the rank to be completed...")
  percent.complete <- 0
  # Initializes the progress bar
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = 100, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  while (percent.complete < 100) {
    # Sys.sleep(5)
    Sys.sleep(300) # wait for 5' before checking the files
    percent.complete <- 100 * sum(as.numeric(
      sub("\\.RDS","",basename(mylist)) %in%
        sub("\\.csv\\.gz", "", basename(list.files(paste0(results.dir, "/intermediate.files/5rank/S01A/"))))
    )) / length(mylist) # calculate percent complete
    setTxtProgressBar(pb, percent.complete)
  }
  close(pb) # Close the connection

  ###
  message("Collecting the individual signature files")

  # list.dirs(paste0(results.dir, "/intermediate.files/5rank/S01A/"))
  message('All jobs finished and writing all signatures to all.signatures.csv.gz file')
  signatures <- list.files(path = paste0(results.dir, "intermediate.files/5rank/S01A/"), pattern = "\\.csv\\.gz$", full.names = TRUE)
  all.signatures <- do.call(
    rbind,
    pbapply::pblapply(
      # pbmcapply::pbmclapply(
      stats::setNames(signatures, basename(signatures)),
      FUN = function(i) {
        fread(i)
      }
    )
  )
  fwrite( all.signatures, paste0(paste0(results.dir, "intermediate.files/all.signatures.csv.gz")) )
  #rm(all.signatures)
  gc()



#   ##############################################################################
#   # Now collecting across all the results
#
#   ###
#   # Script for averaging the ranking
#
#
  used.types <- unlist(strsplit(grep.sig.pattern, "\\|"))[[1]]  ### only run for trt_cp
  used.types <- types.of.data[text.pattern %in% used.types]
#
#
#
# #### !!!!!!!!!!!!!!!!!!!!!!!!!! split by GWAS and model?
#
#
  for (i in seq(nrow(used.types))) {
     # writing code to be able to be run both for groups and individual data.sets in the future.
    message("Loading the master signature file...")
    thesesignatures <- all.signatures
    #thesesignatures <- fread(paste0(paste0(results.dir, "intermediate.files/all.signatures.csv.gz")))
    message(paste0("Now working on ", used.types$data.name[i], " (", used.types$text.pattern[i],  ")"))
    thesesignatures <- thesesignatures[sig_id %in% unique(signature.inventory[Signature.type %in% used.types$text.pattern[i]]$sig_id) ]
    fixed.columns   <- c("gwas", "model_ID", "sig_id")
#
#     ### rank extractor ###
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

    gwas.model.combos <- unique(df[, c("gwas","model_ID") ])
    output <- do.call(
      rbind,
      pbmcapply::pbmclapply(
        seq(nrow(gwas.model.combos)),
        FUN = function (j) {
          thisz <- thesesignatures[gwas == gwas.model.combos$gwas[j] & model_ID == gwas.model.combos$model_ID[j]]
#
#           ### Getting Avg Rank
          message("Summarizing the five rank method...")
          five.rank.all <- do.call(
            cbind,
            lapply(
              stats::setNames(
                unlist(c("", paste0("Perm_", seq(noperm), "_"))),
                unlist(c("Actual", paste0("Perm_", seq(noperm))))),
              rank_extractor,
              z = thisz ) )
#
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
            "gwas"         = gwas.model.combos$gwas[j],
            "model_ID"     = gwas.model.combos$model_ID[j],
            "sig_id"       = thisz$sig_id,
            "AvgRank"      = as.numeric(five.rank.all[,1]),
            "perm.p"       = perm.p
            )
          return(final.rank)
          }, mc.cores = n.threads# more conservative.
) # pbmclapply gwas.model.combos end
) # do.call end.
#
#     # Annotate
    sig.annotation <- return_df(sig.annotation)
    output <- as.data.table(dplyr::left_join(output, sig.annotation))
#
#     # Save
    fwrite(
      output,
      paste0(results.dir, "/", used.types$text.pattern[i],
             "_", "all.signatures.AvgRank.csv.gz"))
    sink( paste0(results.dir, "/", used.types$text.pattern[i],
                 "_", "all.signatures.AvgRank.readme") )
    print(paste0("Please note that the permutation test is run for each gwas-model combination"))
#
    } # save for all data.types.

} # function ends
