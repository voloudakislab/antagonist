#' Aggregate and prioritize compounds
#'
#' Aggregates across models. Run separately if you wish to separate models. Perturbagen Library used. We use the LINKS Phase II L1000 dataset (GSE70138) perturbagen reference library1. All inferred genes (AIG; n=12,328) are considered. Only “gold” signatures are considered. Imputed transcriptomes used. Analysis is limited to 17 imputed transcriptomes: (1) the B2 phenotype, and (2) EpiXcan tissue models that have at least one FDR-significant finding (FDR significance takes into account all COVID phenotypes and all 42 tissue models considered): Adipose: subcutaneous (GTEx), Adipose: subcutaneous (STARNET), Adipose: visceral (GTEx), Adipose: visceral (STARNET), Artery: Aorta (GTEx), Artery: Aorta (STARNET), Artery: Mammary (STARNET), Blood (STARNET), GI: esophagus, GE junction (GTEx), GI: esophagus, mucosa (GTEx), GI: muscularis (GTEx), GI: pancreas (GTEx), Muscle: skeletal (GTEx), Muscle: skeletal (STARNET), Reproductive: mammary tissue (GTEx), Respiratory: lung (GTEx), Skin: sun exposed lower leg (GTEx). Summarization and prioritization approach: From the original CDR pipeline2 (using 100 permutations), applied to the imputed transcriptome of each tissue, we obtain the average rank of the compound in antagonizing the GReX and the permutation p values. After pulling together all the results, we perform a Mann-Whithney U test for each candidate compound/ shRNA against all other compounds to see if the candidate’s rankings significantly deviate from the median rank. For each candidate we also estimate a GReX antagonism pseudo z-score, which is defined as the negative Hodges-Lehmann estimator (of the median difference between that specific candidate vs. the other candidates) divided by the standard deviation of the ranks of the compounds (-Hodges-Lehmann estimatorperturbagenSD average ranks of all perturbagens); a positive pseudo z-score is interpreted as a potential therapeutic candidate whereas a negative pseudo z-score would suggest that the shRNA is not antagonizing the imputed transcriptome and is thus likely to exacerbate the phenotype. Of note is that at this stage each candidate is compared against the other candidates but we can confirm that the candidate is effectively antagonizing the GReX by looking at the original permutation p values. FDR is estimated with the Benjamini–Hochberg procedure3. Additional information for chemical compound analyses: Analysis is limited to compounds eligible for drug repurposing (n=495). Drug information for the compounds under consideration (e.g. clinical phase, mechanism of action and molecular targets) was obtained from http://www.broadinstitute.org/repurposing (file date: 3/24/2020). For comparison with other studies; the compounds under question were compared with all the other compounds. For the mechanism of action comparison all compounds with a known mechanism of action represented with two or more candidates are considered. Final recommendations are for launched medications and FDR correction is applied only to launched compounds. Additional information for shRNA analyses: All shRNAs were considered.
#' To run models separately a loop will have to pass arguments to limit.models.to (and keep NA if needed). Requires up to 10GB per thread.
#' Consider "clinical_phase" column options are "Launched", "Phase 3", "Phase 2/Phase 3", "Phase 2", "Phase 1/Phase 2", "Phase 1", "Preclinical" and "Withdrawn". A sensible list for CDR is "Launched|Phase 3|Phase 2"
#'
#' @param dfs.to.process dfs to be loaded based on regular expression.
#' @param limit.dfs.to regular expression for inclusion rlist of datasets.
#' @param limit.models.to Limit model_IDs
#' @param limit.gwass.to Limit gwas
#' @param output.dir Default is "results/CDR_GTP/"
#' @param compound.type Options are CDR for computational drug repurposing/repositioning and GTP for gene target prioritization which considers trt_sh (shRNA), trt_oe ().
#' @param drug.file trt_cp specirific parametereCompound file; default is repurposing_drugs_20200324.txt from clue.io
#' @param cell.file Compound file; default is cell_lkines.csv from clue.io
#' @param ref.drug File with information for compounds (trt_cp) such as mechanism of action, indication etc. Look at data-raw if you want to generate a similar file
#' @param ref.cell File with informaiton about the cell lines used in the perturbagen library.
#' @param iterative.tissue when having more than one tissue/cell type aggregate for each tissue as well.
#' @param min.experiments.n Min experiment per perturbagen (default is 2)
#' @param n.cores Number of cores to use if multithreaded
#' @return Compound level summary
#' @export
#'

aggregate_and_prioritize = function(
  dfs.to.process                  = "results/GTP_CDR/.+_.+_all\\.signatures\\.AvgRank\\.csv\\.gz$",
  #dfs.to.process                  = '/sc/arion/projects/va-biobank/software/Georgios_dev/results/GTP_CDR_prototyping_2RDS_5com_Jamie_SCZ_070924/trt_cp_all.signatures.AvgRank.csv.gz',
  limit.dfs.to                    = "trt_cp|trt_sh|trt_oe|trt_xpr",
  limit.models.to                 = NA,
  limit.gwass.to                  = NA,
  output.dir                      = "results/GTP_CDR/",
  #output.dir                      = "/sc/arion/projects/va-biobank/software/Georgios_dev/results/GTP_CDR_prototyping_2RDS_5com_Jamie_SCZ_070924",
  # For CDR
  ref.drug                        = TRT_CP.INFO.20200324,
  # Cell experiments
  ref.cell                        = CELL.LINE.INFO,
  # Iterative
  iterative.tissue                = TRUE,
  # parameters
  min.experiments.n               = 2,
  n.cores                         = parallel::detectCores()-2
) {

  # FIXME: make sure that this runs for each GWAS separately.
  # FIXME: do both CDR and CTG if the load data are there.
  # TODO: multitissue vs. singletissue.

  gc()

  ### HELPER SCRIPTS ###
  '%!in%' <- function(x,y)!('%in%'(x,y))

  ### SETTING THE PARAMETERS ###
  ref.drug <- return_df(ref.drug)
  ref.cell <- return_df(ref.cell)

  ### LOAD THE DATA ###
  # TODO: Load the files from the save directory and then automatically do appropriate actions.
  message("Loading the data...")
  dfs.to.process <- list.files(
    path = dirname(dfs.to.process),
    pattern = basename(dfs.to.process),
    full.names = T )
  dfs.to.process <- dfs.to.process[grep(limit.dfs.to, dfs.to.process)]
  # Probe one file to get some info
  gwass  <- unique(return_df(dfs.to.process[1])$gwas)
  if (!is.na(limit.gwass.to[1])) {
    gwass <- gwass[gwass %in% limit.gwass.to] }

  x.cdr <- dfs.to.process[grep("trt_cp", dfs.to.process)]
  x.gtp <- dfs.to.process[grep("trt_sh|trt_oe|trt_xpr", dfs.to.process)]


  ### CDR ###
  if (length(x.cdr) != 0) {
    x.cdr <- return_df(x.cdr)
    x.cdr <- as.data.table(dplyr::left_join(x.cdr, ref.drug))
    x.cdr$readily.repurposable <- !is.na(x.cdr$clinical_phase)

    # Keep at least min.experiments
    x.cdr[, N_experiments := length(AvgRank), by = pert_iname]
    x.cdr <- x.cdr[N_experiments >= min.experiments.n]

  }

  ### GTP ###
  if (length(x.gtp) != 0) {
    x.gtp <- do.call(
      rbind,
      lapply(
        x.gtp,
        FUN = function(z) {
          return_df(z)
        }))

    # Keep at least min.experiments
    x.gtp[, N_experiments := length(AvgRank), by = pert_iname]
    x.gtp <- x.gtp[N_experiments >= min.experiments.n]

  }

  # #################################
  # # Density plot for IL10RB (paper)
  # x$shRNA <- paste0("other (n=",length(x[pert_iname != "IL10RB"]$shRNA),")")
  # x[pert_iname == "IL10RB"]$shRNA <- paste0("IL10RB (n=",length(x[pert_iname == "IL10RB"]$shRNA),")")
  # IL10RBplot <- ggpubr::ggdensity( # https://rpkgs.datanovia.com/ggpubr/reference/ggdensity.html
  #   x, x = "AvgRank", add = "mean", xlab = "Signature AvgRank",
  #   color = "shRNA", fill = "shRNA",
  #   palette = vector_to_colors(unique(x$shRNA)) # , rug = TRUE
  #   )
  # cowplot::ggsave2( paste0(output.dir, "/IL10RB_AvgRank_distribution_square.pdf") , IL10RBplot, height = 4, width = 4 )
  # cowplot::ggsave2( paste0(output.dir, "/IL10RB_AvgRank_distribution_landscape.pdf") , IL10RBplot, height = 3, width = 10 )


  # This is no longer used, since most signatures come only from specific cell lines.
  # # Remove those that were done in irrelevant cell lines - keeping lung and immune cell lines
  # if (discard.non.relevant.cell.lines) {
  #   x <- as.data.table(inner_join(x, cell_lines))
  #   x <- x[!is.na(clinical_phase)]
  # }

  # TODO: Use a combination approach like this to iterate
  # tempdf <- fread("output/2.METAXCAN/OUD_df.all.annotated.onlytranscripts.csv.gz")
  # ucomb  <- unique(tempdf[,c("gwas", "model_ID")]) # this is need for the ancestry specific approach
  # for (i in seq(nrow(ucomb)) ) {
  #   aggregate_and_prioritize(
  #     limit.gwass.to   = ucomb[i]$gwas,
  #     limit.models.to  = ucomb[i]$model_ID,
  #     n.cores          = parallel::detectCores()/2)
  # }


  lapply(
    gwass,
    FUN = function(thisgwas) {

      message(paste0("Now processing ", thisgwas))
      this.output.dir <- paste0(output.dir, "/", thisgwas, "/")
      gv_dir.create(this.output.dir)

      ### MAIN SCRIPT ###
      pbapply::pblapply(
        list(x.cdr, x.gtp),
        FUN = function(x) {
          if (length(x) != 0) {

            # Limit to GWAS of interest
            x <- x[gwas == thisgwas]

            # Limit to models of interest and implement iterative approach
            if (!is.na(limit.models.to[1])) x <- x[model_ID %in% limit.models.to]
            if (iterative.tissue) { # this is a bit complicated to accomodate for merging a list with a vector
              limit.models.to <- list(limit.models.to)
              existing.models <- unique(x$model_ID)
              for (z in seq(length(existing.models))) {
                limit.models.to[z+1] <- existing.models[z] }
              limit.models.to <- unique(limit.models.to)
              # If it is one model, don't run it twice
              if (length(limit.models.to)<=2) limit.models.to <- limit.models.to[-1]
              }



            # Now run all the model_IDs
            pbapply::pblapply(
              limit.models.to,
              FUN = function(limit.models.to) {
                message(paste0("Now working on: ", paste(limit.models.to, collapse = ", ")))
                if (!is.na(limit.models.to[1])) x <- x[model_ID %in% limit.models.to] # relimit to model under question
                modelprefix <- as.character(make_java_safe(paste(limit.models.to, collapse = "_X_")))
                if (modelprefix == "NA") modelprefix <- "ALL"
                this.output.dir <- paste0(
                  this.output.dir, "/",
                  modelprefix, "/")
                gv_dir.create(this.output.dir)

                ###########################
                # Density plot for Avg rank
                d.plot <- ggpubr::ggdensity( # https://rpkgs.datanovia.com/ggpubr/reference/ggdensity.html
                  data.frame(Type = paste(x$pert_type,collapse = "|"),
                             AvgRank = as.numeric(x$AvgRank)),
                  x = "AvgRank", add = "mean", xlab = "Signature AvgRank")
                cowplot::ggsave2( paste0(this.output.dir, "/", thisgwas,"_", ifelse(unique(x$pert_type)[1] == "trt_cp", "cdr", "gtp"), "_AvgRank_distribution_square.pdf") , d.plot, height = 4, width = 4 )
                cowplot::ggsave2( paste0(this.output.dir, "/", thisgwas,"_", ifelse(unique(x$pert_type)[1] == "trt_cp", "cdr", "gtp"), "_AvgRank_distribution_landscape.pdf") , d.plot, height = 3, width = 8 )


                ########################
                # Calculating statistics
                # hist(x$AvgRank) # Seems to be following the gaussian distribution
                # avg.rank.mean <- mean(x$AvgRank)
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

                fwrite(x, paste0(this.output.dir, "/", thisgwas,"_",
                                 ifelse(unique(x$pert_type)[1] == "trt_cp", "cdr","gtp"),
                                 "_signature_level.csv.gz"))

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
                  fwrite(moa.repurp, paste0(this.output.dir, "/", thisgwas, "_",
                                            ifelse(unique(x$pert_type)[1] == "trt_cp", "cdr","gtp"),
                                            "_moa_level.csv"))

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
                  fwrite(target.repurp, paste0(this.output.dir, "/", thisgwas, "_",
                                               ifelse(unique(x$pert_type)[1] == "trt_cp", "cdr","gtp"),
                                               "_target_level.csv"))


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
                  fwrite(disease_area.repurp, paste0(this.output.dir, "/", thisgwas, "_",
                                                     ifelse(unique(x$pert_type)[1] == "trt_cp", "cdr","gtp"),
                                                     "_disease_area_level.csv"))

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
                  fwrite(indication.repurp, paste0(this.output.dir, "/", thisgwas, "_",
                                                   ifelse(unique(x$pert_type)[1] == "trt_cp", "cdr","gtp"),
                                                   "_indication_level.csv"))


                  # Compile and save
                  # This is imperfect for compounds that have more than one moa
                  compound.level <- as.data.table(dplyr::left_join(
                    compound.level, moa.repurp[, c(
                      "moa", "Rank_moa", "MOA.pseudo.zscore",
                      "MOA.HL.estimate", "MOA.avg.rank.sd",
                      "MOA.MW.p", "MOA.MW.FDR")]))
                  fwrite(compound.level, paste0(this.output.dir, "/", thisgwas, "_",
                                                ifelse(unique(x$pert_type)[1] == "trt_cp", "cdr","gtp"),
                                                "_all_compound_level.csv"))

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

                  fwrite(compound.level, paste0(this.output.dir, "/", thisgwas, "_",
                                                ifelse(unique(x$pert_type)[1] == "trt_cp", "cdr","gtp"),
                                                "_P3_and_launched_compound_level.csv"))

                  compound.level <- compound.level[clinical_phase=="Launched"] # 201 compounds
                  compound.level$Rank <- rank(compound.level$Rank) # rerank
                  column.order <- names(compound.level)
                  moa.ranks <- unique(compound.level[,c("moa", "Rank_moa")])
                  moa.ranks <- moa.ranks[!is.na(Rank_moa)]
                  moa.ranks$MOA.Rank <- rank(moa.ranks$Rank_moa) # rerank
                  moa.ranks[, Rank_moa := NULL]
                  compound.level[, Rank_moa := NULL] # remove
                  compound.level <- as.data.table(dplyr::left_join(compound.level, moa.ranks))
                  fwrite(compound.level, paste0(this.output.dir, "/", thisgwas, "_",
                                                ifelse(unique(x$pert_type)[1] == "trt_cp", "cdr","gtp"),
                                                "_launched_compound_level.csv"))
                } else {

                  ### GTP SPECIFIC SCRIPTS ###
                  compound.level$Rank <- rank(compound.level$AvgRank) # replace Rank
                  compound.level$Compound.MW.FDR <- p.adjust(compound.level$Compound.MW.p, method = "fdr")
                  compound.level <- compound.level[order(Rank)][
                    , c("pert_iname", "Rank", "AvgRank", "Compound.MW.p",
                        "Compound.pseudo.zscore",
                        "Compound.HL.estimate", "Compound.avg.rank.sd",
                        "Compound.MW.FDR",
                        "N_experiments", "perm.p.all")]
                  fwrite(compound.level, paste0(this.output.dir, "/", thisgwas, "_",
                                                ifelse(unique(x$pert_type)[1] == "trt_cp", "cdr","gtp"),
                                                "_compound_level.csv"))
                }
              }) # Iterative model loop ends here.

          } # section of analysis specific scripts done.




        } # function for thisgwas loop ends here
      ) # lapply loop for CDR and GTP dataframes


    } # lapply function for all GWASs
  ) # lapply loop for all GWASs

  # TODO: Compile the results and save as xls add this description as you do in MultiWAS package
  # pert_iname: common name of compound
  # clinical_phase: clinical phase of compound as per 2018-09-07
  # Rank: Overall Rank of compound
  # AvgRank: average of all average ranks of compounds during permutation analysis
  # Compound.WRS.p: Wilcoxon rank sum test p value for placement out of all experiments and GWASs if applicable
  # moa: mechanism of action
  # Rank_moa_percentile: moa ranking percentile (only if 2 or more drugs per moa)
  # MOA.WRS.P: Wilcoxon rank sum test p value for placement out of all compounds (were already ranked at compound level)
  # target: molecular targets
  # disease_area:
  # indication:
  # N_experiments: number of experiments for each TWAS
  # perm.p.all: all perm.p values collapsed
  # return(compound.level)

} # function definition ends here.
