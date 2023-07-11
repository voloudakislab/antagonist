#' Antagonism core script
#'
#' The following 5 variables/methods need to be saved in splitting the drug matrices:
#' 1. cor.pearson : Pearson (all)
#' 2. cor.spearman : Spearman (all)
#' 3. extreme.cor.spearman : Spearman (most differentially expressed genes)
#' 4. extreme.cor.pearson : Pearson (most differentially expressed genes)
#' 5. ks.signed : KS method (most differentially expressed genes)
#' This is a re-write of So HC et al. Analysis of genome-wide
#' association data highlights candidates for drug repositioning in psychiatry.
#' Nat Neurosci. 2017 Oct;20(10):1342-1349. PMID: 28805813.
#' It has been adapted to run for signatures (individual parameters) vs.
#' summarized compound scores.
#'
#' @param df twas data frame
#' @param results.dir Results main dir. Subdirectories will be created
#' @param signature.dir Location of signature files. Consider saving locally to speed things up (e.g. /scratch/cmap_l1000_2021_01_28/)
#' @param gene.anno.file Gene annotation file (provided from L1000)
#' @param grep.sig.pattern Regular exrpression pattern to grep from signature name files. Default is to process sets for computational drug repurposing and gene target prioritization.
#' @param noperm Number of permutations (the final number of permutation data points is noperm × number of drugs)
#' @param thres.N.vector defining the threshold(s) K such that only top K items are included for comparison of expressions: use the same thresholds as So et al.
#'
#' @return N/A. Saves the figures
#' @keywords visualization regional miami plot
#' @export
#'
perform_antagonism <- function(
  # Input
  df              ,
  # Output
  results.dir      = "results/GTP_CDR/",
  # Parameters
  signature.dir    = "/sc/arion/projects/va-biobank/resources/CMap/cmap_l1000_2021-11-20/eachDrug/", # /scratch/cmap_l1000_2021_01_28/
  gene.anno.file   = "/sc/arion/projects/va-biobank/resources/CMap/cmap_l1000_2021-11-20/geneinfo_beta.txt",
  grep.sig.pattern = "trt_cp|trt_sh|trt_oe|trt_xpr",
  noperm           = 100, #
  thres.N.vector   = c(50,100,250,500)

) {

  types.of.data <- data.table(
    text.pattern = c("trt_cp", 	"trt_sh", "trt_oe", "trt_xpr", "trt_misc", "ctl"),
    data.name    = c("compounds", "shRNA", "over.expression", "CRISPR",
                        "other.treatments", "negative.controls"),
    perturbagen.group = c("drug", "gene", "gene", "gene", NA, NA))


  ##############################################################################
  # DISTILED MAIN FUNCTION FOR FIVE METHOD APPROACH BY SO ET AL 2017 ###########
  ##############################################################################
  five_rank_method <- function(
    thres.N     , # threshold NA is all values otherwise consider top values
    x           , # the payload of interest
    scramble    = F) {

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

  ##############################################################################
  # RUN STATISTICS FOR ALL SIGNATURES ##########################################
  ##############################################################################
  # TODO: make a for loop that goes from trt_cp etc... Combine trt_sh and trt_xpr

  # Prototyping
  # signature.dir <- "/scratch/cmap_l1000_2021_01_28/"

  suppressMessages(library(mygene))
  suppressMessages(library(dplyr))
  suppressMessages(library(ggplot2))
  suppressMessages(library(Hmisc))

  library(data.table)

  # Preparing TWAS
  df               <- MultiWAS::return_df(df)
  df               <- df[zscore!=-Inf & zscore!=Inf]
  # gwas.model.combs <- unique(df[, c("gwas","model_ID") ])

  # Preparing Directories
  MultiWAS::gv_dir.create(results.dir)
  MultiWAS::gv_dir.create(paste0(results.dir, "intermediate.files/"))

  # Identifying signatures

  cmap.list <- list.files(signature.dir, full.names = T)
  # cmap=readRDS("/sc/hydra/projects/roussp01b/Wen/LINCS_Level5/eachDrugPhase2/1-300_cp.RDS")
  # to.process <- as.data.table(tidyr::crossing(to.process, data.table("cmap.file" = cmap.list)))
  if (!is.na(grep.sig.pattern[1])) mylist <- cmap.list[grep(grep.sig.pattern, cmap.list)]

  # Create a signature inventory
  message("Building signature inventory by perturbagen type")
  signature.inventory <- pbmclapply(
    stats::setNames(mylist, sub("_[[:digit:]]+\\.[[:digit:]]+\\.RDS", "", basename(mylist)) ),
    FUN = function(x) {colnames(readRDS(x))},
    mc.cores = detectCores()-2
  )
  signature.inventory <- do.call(rbind, lapply(
    unique(unlist(names(signature.inventory))),
    FUN = function(x) {
      data.table(
        "Signature.type" = x,
        "sig_id" = unique(unlist(signature.inventory[grep(x, names(signature.inventory))]))
      ) } ) )



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
          mc.cores = ifelse(parallel::detectCores()-2<1,1,parallel::detectCores()-2)
        ) # mclapply loop for all signattures
        return(do.call(rbind,signature))
      }#,
      # mc.cores = ifelse(parallel::detectCores()-2<1,1,parallel::detectCores()-2)
    ) # lapply loop ends
    ) # do.call rbind ends

  # TODO: add the fwrite below.
  fwrite( all.signatures, paste0(paste0(results.dir, "intermediate.files/all.signatures.csv.gz")) )
  rm(all.signatures)
  # fwrite(all.signatures, "inst/extdata/all.signatures.csv.gz")

  ##############################################################################
  # Now collecting across all the results

  ###
  # Script for averaging the ranking


  used.types <- unlist(strsplit(grep.sig.pattern, "\\|"))
  used.types <- types.of.data[text.pattern %in% used.types]

  for (i in seq(nrow(used.types))) {
    # writing code to be able to be run both for groups and individual data.sets in the future.
    thesesignatures <- fread(paste0(paste0(results.dir, "intermediate.files/all.signatures.csv.gz")))
    message(paste0("Now working on ", used.types$data.name[i], " (", used.types$text.pattern[i],  ")"))
    thesesignatures <- unique(signature.inventory[Signature.type %in% used.types$text.pattern[i] ]$sig_id)
    thesesignatures <- all.signatures[sig_id %in% thesesignatures]
    fixed.columns   <- c("gwas", "model_ID", "sig_id")
    ranks.actual <- data.table(
      "gwas"                      = thesesignatures$gwas,
      "model_ID"                  = thesesignatures$model_ID,
      "sig_id"                    = thesesignatures$sig_ID,
      rank.est.pearson            = rank(thesesignatures$ALL.cor.pearson),
      rank.est.spearman           = rank(thesesignatures$ALL.cor.spearman),
      mean.rank.ks                = ,
      mean.rank.extreme.pearson   =,
      mean.rank.extreme.spearman  =
    )

    thesesignatures[ ,.(names(thesesignatures)[grep("^[[:digit:]]+\\.cor.pearson", names(thesesignatures))]), with = F]


  }




    ## The following 5 variables/methods need to be saved in splitting the drug matrices
    # cor.pearson : Pearson all
    # cor.spearman : Spearman all
    # extreme.cor.spearman : Spearman with most differentially expressed genes
    # extreme.cor.pearson : Pearson with most differentially expressed genes
    # ks.signed : KS method

    drugname=colnames(finalres)[14:ncol(finalres)]
    save(nodrugs,no.thres.N,drugname, cor.pearson,cor.spearman,ks.signed, extreme.cor.spearman,extreme.cor.pearson,file=paste("/sc/arion/projects/CommonMind/wen/covid19/ANA_A2_V2/results/drugRepurpose/EpiXcan/RData/", "Correlations_5methods_w_permutation_",
                                                                                                                              gsub(".csv", "_1-300_cp.RData", x=filenames[q]),sep="") )

    # above is the correlation calculations without permutations. After getting correlations of each of the splitted sub-drug sets
    # we combine them and calculate the Average ranks according to Five methods


    #*********************************************************************
    # permutation to determine the significance of the avg. rank obtained (repeats the above code)
    #*********************************************************************
    noperm=100
    cor.spearman=NULL
    cor.pearson=NULL

    ks.signed = matrix(ncol=no.thres.N,nrow=nodrugs)
    extreme.cor.spearman = matrix(ncol=no.thres.N,nrow=nodrugs)
    extreme.cor.pearson = matrix(ncol=no.thres.N,nrow=nodrugs)
    avgrank.perm <- matrix(nrow=nodrugs, ncol=noperm)

    for (r in 1:noperm) {
      disease.zscore = sample(finalres$zscore)
      disease.zscore=as.numeric(disease.zscore)
      for (i in 1:nodrugs) {
        spearman.obj = cor.test(as.numeric(finalres[,i+13]),disease.zscore, method="spearman",use="na.or.complete", exact = FALSE)
        cor.spearman[i] = spearman.obj$estimate
        pearson.obj  = cor.test(as.numeric(finalres[,i+13]),disease.zscore, method = "pearson",use="na.or.complete")
        cor.pearson[i]= pearson.obj$estimate

        nogenes = nrow(finalres)

        pos.zscore = disease.zscore[disease.zscore>=0]
        no.pos.genes = length(pos.zscore)
        neg.zscore = disease.zscore[disease.zscore<0]
        no.neg.genes = length(neg.zscore)


        for (j in 1:no.thres.N) {
          thres.N = thres.N.vector[j]
          ind.upreg = order(disease.zscore,decreasing=TRUE)[1:thres.N]

          #****************************************
          #   KS method
          #******************************************
          geneset2 = rank(-as.numeric(finalres[,i+13]))[ind.upreg]
          geneset2 = sort(geneset2)

          a.up = max(   (1:thres.N)/thres.N - geneset2/nogenes      )
          b.up = max(   geneset2/nogenes -  (1:thres.N-1)/thres.N )
          ks.test.obj.up = ifelse(a.up>b.up, a.up , -b.up)


          ind.downreg= order(disease.zscore,decreasing=FALSE)[1:thres.N]

          geneset2 = rank(-as.numeric(finalres[,i+13]))[ind.downreg]
          geneset2 = sort(geneset2)
          a.down = max(   (1:thres.N)/thres.N - geneset2/nogenes      )
          b.down = max(geneset2/nogenes - (1:thres.N-1)/thres.N )
          ks.test.obj.down = ifelse(a.down>b.down, a.down , -b.down)


          ks.signed[i,j]=0
          if (  sign(ks.test.obj.up)!=sign(ks.test.obj.down)  )
          {
            ks.signed[i,j] = ks.test.obj.up - ks.test.obj.down

          }

          #****************************************
          #   spearman and Pearson correlations (either use all observations or only compare the most extreme observations)
          #****************************************
          spearman.obj.extr = cor.test( as.numeric(finalres[,i+13][c(ind.upreg,ind.downreg)]), disease.zscore[c(ind.upreg,ind.downreg)]
                                        ,method="spearman",use="na.or.complete", exact = FALSE   )
          extreme.cor.spearman[i,j] =  spearman.obj.extr$estimate
          pearson.obj.extr  = cor.test(as.numeric(finalres[,i+13][c(ind.upreg,ind.downreg)]),disease.zscore[c(ind.upreg,ind.downreg)], method = "pearson",use="na.or.complete")
          extreme.cor.pearson[i,j]= pearson.obj.extr$estimate

        }

      } #end of looping over all drugs

      drugname=colnames(finalres)[14:ncol(finalres)]

      ## for each permutation, needs to store the five methods:
      save(nodrugs,no.thres.N,drugname, cor.pearson,cor.spearman,ks.signed, extreme.cor.spearman,extreme.cor.pearson,file=paste("/sc/arion/projects/CommonMind/wen/covid19/ANA_A2_V2/results/drugRepurpose/EpiXcan/RData/", "Correlations_5methods_permutation_",r,"_",
                                                                                                                                gsub(".csv", "_1-300_cp.RData", x=filenames[q]),sep="") )

      # then need to integrate all the above 100 permutation results


    }  # end of permutation loop

  } ## end of entire loop (looping over each brain region file)
  ## above code will generate the correlation files and permutations for all subsets

}
