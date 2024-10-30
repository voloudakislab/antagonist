#' Antagonism core script
#'
#' This is intended to run in your personal computer. There is a version that is optimized for LSF based clusters.
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
#' Analysis is run for each "gwas" (trait) - "model_ID" (imputation model or specific analysis) pair
#' It has been adapted to run for signatures (individual parameters) vs.
#' summarized compound scores.
#'
#'
#' @param df twas data frame
#' @param column.feature expected input is ENSEMBL ID and currently only works with genes (default: feature)
#' @param column.statistic z-score statistic (default: zscore)
#' @param column.trait trait name (default: gwas)
#' @param column.source disease signature source, e.g. microglia, DLPFC, meta-analysis, etc. (default: model_ID). For non-twas inputs make sure there is a column describing the input.
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
#' @return N/A. Saves the figures
#' @keywords visualization regional miami plot
#' @export
#'
perform_antagonism <- function(
  # Input
  df                     = paste0("output/2.METAXCAN/", basename(getwd()),"_df.all.annotated.onlytranscripts.csv.gz"),
  column.feature         = "feature",
  column.statistic       = "zscore",
  column.trait           = "gwas",
  column.source          = "model_ID",
  # Output
  results.dir            = "results/GTP_CDR/",
  # Parameters
  n.threads              = parallel::detectCores()-2,
  signature.dir          = "/sc/arion/projects/va-biobank/resources/CMap/cmap_l1000_2021-11-20/eachDrug/", # /scratch/cmap_l1000_2021_01_28/
  gene.anno.file         = "/sc/arion/projects/va-biobank/resources/CMap/cmap_l1000_2021-11-20/geneinfo_beta.txt",
  grep.sig.pattern       = "trt_cp|trt_sh|trt_oe|trt_xpr",
  # grep.cdr.pattern       = "Launched|Phase 3|Phase 2", # NA for no filtering
  noperm                 = 100, #
  thres.N.vector         = c(50,100,250,500),
  sig.annotation         = SIG.INFO.20211120,
  # trt_cp.annotation      = TRT_CP.INFO.20200324,
  overwrite.intermediate = TRUE,
  model.banlist.grep     = ":: Transcripts ::|:: H3K4me3 ::|:: H3K27ac ::|:: CA ::",
  prototyping            = NA
) {

  types.of.data <- data.table(
    text.pattern = c("trt_cp", 	"trt_sh", "trt_oe", "trt_xpr", "trt_misc", "ctl"),
    data.name    = c("compounds", "shRNA", "over.expression", "CRISPR",
                        "other.treatments", "negative.controls"),
    perturbagen.group = c("drug", "gene", "gene", "gene", NA, NA))

  ##############################################################################
  # RUN STATISTICS FOR ALL SIGNATURES ##########################################
  ##############################################################################
  # TODO: make a for loop that goes from trt_cp etc... Combine trt_sh and trt_xpr

  # Prototyping
  # signature.dir <- "/scratch/cmap_l1000_2021_01_28/"

  # suppressMessages(library(mygene))
  # suppressMessages(library(dplyr))
  # suppressMessages(library(ggplot2))
  # suppressMessages(library(Hmisc))

  # library(data.table)

  # Preparing TWAS
  df               <- MultiWAS::return_df(df)

  # Preparing columns
  df$feature   <- df[[column.feature]]
  df$zscore    <- df[[column.statistic]]
  df$gwas      <- df[[column.trait]]
  df$model_ID  <- df[[column.source]]
  df           <- df[zscore!=-Inf & zscore!=Inf]
  # gwas.model.combs <- unique(df[, c("gwas","model_ID") ])
  df <- df[, c("feature", "zscore", "gwas", "model_ID")]

  # Keeping only genes - this is based on banlinst:
  if (length(grep(model.banlist.grep, df$model_ID))>0) {
    message(paste0("Currently expression banlist is: ", model.banlist.grep))
    message("Intention is to only keep genes")
    df <- df[-grep(model.banlist.grep, model_ID)] }


  # Preparing Directories
  MultiWAS::gv_dir.create(results.dir)
  MultiWAS::gv_dir.create(paste0(results.dir, "intermediate.files/"))

  # Identifying signatures

  cmap.list <- list.files(signature.dir, full.names = T)
  # cmap=readRDS("/sc/hydra/projects/roussp01b/Wen/LINCS_Level5/eachDrugPhase2/1-300_cp.RDS")
  # to.process <- as.data.table(tidyr::crossing(to.process, data.table("cmap.file" = cmap.list)))
  if (!is.na(grep.sig.pattern[1])) mylist <- cmap.list[grep(grep.sig.pattern, cmap.list)] else mylist <- cmap.list


  # Create a signature inventory
  if (!file.exists(paste0(results.dir, "intermediate.files/signature.inventory.csv")) | overwrite.intermediate) {
  # This part saves a list of signature file name table.
    message("Building signature location inventory")
  fwrite(do.call(rbind, pbmclapply(
    stats::setNames(mylist, mylist),
    FUN = function(x) {
      data.table(
        sig_id   = colnames(readRDS(x)),
        filename = x )
      }, mc.cores = n.threads
  )), paste0(results.dir, "intermediate.files/signature.location.csv"))


  # This part makes a signature inventory
    message("Building signature inventory by perturbagen type")
  signature.inventory <- pbmclapply(
    stats::setNames(mylist, sub("_[[:digit:]]+\\.[[:digit:]]+\\.RDS", "", basename(mylist)) ),
    FUN = function(x) {colnames(readRDS(x))},
    mc.cores = n.threads
  )
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

  # TODO: Split by GWAS and model_ID

  message("Running five methods across each signature")
  if (!file.exists(paste0(paste0(results.dir, "intermediate.files/all.signatures.csv.gz"))) | overwrite.intermediate) {
    # FIXME: this script only outputs one model ID..
    # Multicore strategy split by RDS file
    all.signatures <- do.call(
      rbind,
      pbapply::pblapply(
        # pbmcapply::pbmclapply(
        stats::setNames(mylist, basename(mylist)),
        FUN = function(i) {
          # message(paste0("Now processing: ", basename(i)))
          # i = 100
          # Preparing the signature library
          ## the splitted drug-induced gene expressions from CMap database: need to be splitted due to large number of perturbations
          ## it is a matrix of expression t-statistics/z-statistics with entrez gene id as row names and drugs as column names (we assume a drug expression was compared with controls)
          cmap <- readRDS(i)
          if (grepl("trt_oe", i)) cmap <- cmap * -1 # So that the signatures can be integrated with the knockdown experiments.
          cmap <- as.data.table(cmap, keep.rownames = "entrezgene")
          cmap[, entrezgene := as.integer(entrezgene)]
          # cmap=data.frame(cmap) # change the type back
          # Load gene annotation
          geneAnn <- fread(gene.anno.file)
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
          ) # mclapply loop for all signattures
          return(do.call(rbind,signature))
        }#,
        # mc.cores = ifelse(parallel::detectCores()-2<1,1,parallel::detectCores()-2)
      ) # lapply loop ends
    ) # do.call rbind ends

    # TODO: add the fwrite below.
    fwrite( all.signatures, paste0(paste0(results.dir, "intermediate.files/all.signatures.csv.gz")) )
    rm(all.signatures)
    gc()
    # fwrite(all.signatures, "inst/extdata/all.signatures.csv.gz")
  } # this only runs if the output doesn't exist



  ##############################################################################
  # Now collecting across all the results

  ###
  # Script for averaging the ranking


  used.types <- unlist(strsplit(grep.sig.pattern, "\\|"))
  used.types <- types.of.data[text.pattern %in% used.types]



#### !!!!!!!!!!!!!!!!!!!!!!!!!! split by GWAS and model?


  for (i in seq(nrow(used.types))) {
    # writing code to be able to be run both for groups and individual data.sets in the future.
    message("Loading the master signature file...")
    thesesignatures <- fread(paste0(paste0(results.dir, "intermediate.files/all.signatures.csv.gz")))
    message(paste0("Now working on ", used.types$data.name[i], " (", used.types$text.pattern[i],  ")"))
    thesesignatures <- thesesignatures[sig_id %in% unique(signature.inventory[Signature.type %in% used.types$text.pattern[i]]$sig_id) ]
    fixed.columns   <- c("gwas", "model_ID", "sig_id")

    ### rank extractor ###
    rank_extractor <- function(
    prefix,  # "" for actual otherwise it is the permutation of interest
    z  # the signature dataset as above (e.g. thesesignatures) but preferrably split by unique gwas model_ID combinations.
    ) {
      ranks.actual <- data.table(
        "gwas"                      = z$gwas,
        "model_ID"                  = z$model_ID,
        "sig_id"                    = z$sig_ID,
        "rank.est.pearson"          = rank(z[[paste0(prefix, "ALL.cor.pearson")]]),
        "rank.est.spearman"         = rank(z[[paste0(prefix, "ALL.cor.spearman")]]) )
      # ks (extreme)
      cols <- names(z)[grep(paste0(ifelse(prefix == "", "^", prefix), "[[:digit:]]+\\.ks.signed"), names(z))]
      ranks.actual$mean.rank.ks <- rowMeans(z[,..cols][ , (cols) := lapply(.SD, rank), .SDcols = cols])
      # Extreme pearson
      cols <- names(z)[grep(paste0(ifelse(prefix == "", "^", prefix), "[[:digit:]]+\\.cor.pearson"), names(z))]
      ranks.actual$mean.rank.extreme.pearson <- rowMeans(z[,..cols][ , (cols) := lapply(.SD, rank), .SDcols = cols])
      # Extreme spearman
      cols <- names(z)[grep(paste0(ifelse(prefix == "", "^", prefix), "[[:digit:]]+\\.cor.spearman"), names(z))]
      ranks.actual$mean.rank.extreme.spearman <- rowMeans(z[,..cols][ , (cols) := lapply(.SD, rank), .SDcols = cols])
      # Avg rank
      cols <- c("rank.est.pearson", "rank.est.spearman", "mean.rank.ks", "mean.rank.extreme.pearson", "mean.rank.extreme.spearman")
      # ranks.actual[[paste0(prefix, "AvgRank")]] <- rowMeans(ranks.actual[,..cols][ , (cols) := lapply(.SD, rank), .SDcols = cols])
      return(rowMeans(ranks.actual[,..cols][ , (cols) := lapply(.SD, rank), .SDcols = cols]))
    }


    gwas.model.combos <- unique(df[, c("gwas","model_ID") ])
    output <- do.call(
      rbind,
      pbmcapply::pbmclapply(
        seq(nrow(gwas.model.combos)),
        FUN = function (j) {

          thisz <- thesesignatures[gwas == gwas.model.combos$gwas[j] & model_ID == gwas.model.combos$model_ID[j]]

          ### Getting Avg Rank
          message("Summarizing the five rank method...")
          five.rank.all <- do.call(
            cbind,
            lapply(
            stats::setNames(
              unlist(c("", paste0("Perm_", seq(noperm), "_"))),
              unlist(c("Actual", paste0("Perm_", seq(noperm))))),
            rank_extractor,
            z = thisz ) )

          ### Getting permutation p value
          message("Getting permutation p values...")
          perm.p <- unlist(lapply(
            seq(nrow(five.rank.all)),
            FUN = function(i) {
              sum(
                as.numeric(five.rank.all[i,2:ncol(five.rank.all)]) <=
                  as.numeric(five.rank.all[i,1]) ) /
                (ncol(five.rank.all)-1)
            } ))

          ### Compiling the per gwas, model_ID and sig_id results.
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

    # Annotate
    sig.annotation <- MultiWAS::return_df(sig.annotation)
    # hist(as.data.table(table(as.data.table(dplyr::left_join(output, sig.annotation))$pert_iname, useNA = "always"))$N)
    output <- as.data.table(dplyr::left_join(output, sig.annotation))

    # Save
    fwrite(
      output,
      paste0(results.dir, "/", used.types$text.pattern[i],
      "_", "all.signatures.AvgRank.csv.gz"))
    sink( paste0(results.dir, "/", used.types$text.pattern[i],
                 "_", "all.signatures.AvgRank.readme") )
    print(paste0("Please note that the permutation test is run for each gwas-model combination"))
    sink()
    # plyr::join_all(list(x,y,z), by='Flag', type='left')

  } # save for all data.types.

} # function ends
