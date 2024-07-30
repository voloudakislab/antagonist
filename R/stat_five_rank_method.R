#' Distilled main function for five-rank method approach by So et al 2017
#'
#' Intention is to be used to run the 5 methods for each GWAS - model_ID combination.
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
#' Of note that in their paper, it is a bit ambiguous what is done, that's why
#' the code was migrated ("Permutation_test_compare_expression_corrected.R")
#'
#'
#' @param thres.N threshold NA is all values otherwise consider top values
#' @param x payload of interest
#' @param scramble scramble the z-scores (used for permutation analysis)
#' @return Returns five ranks
#' @keywords five rank method
#' @export
#'
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
