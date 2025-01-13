#' Split gctx files in RDS files
#'
#' Prepare RDS files from LINCS gctx files for main pipeline. Need to run this for xpr
#' Requires 1) cmapR (either from bioconductor or github), 2) .gctx files, 3) siginfo_beta.txt
#'
#' Defaults to selecting "Golden" signatures, a heuristic for assessing whether
#' a signature is reproducible and distinct:
#' https://docs.google.com/document/d/1q2gciWRhVCAAnlvF2iRLuJ7whrGP6QjpsCMq1yWz7dU
#' Requires cmapR package either from bioconductor or github
#'
#' @param split.in.n How many signatures in each RDS file. Default is 300
#' @param parent.signature.dir Location of signature files. These are gw Consider saving locally to spead things up (e.g. /scratch/cmap_l1000_2021_01_28/)
#' @param distil_cc_q75_min 75th quantile of pairwise Spearman correlations in landmark space of replicate level 4 profiles. Default is 0.2 corresponding to gold signatures.
#' @param pct_self_rank_q25_max Self connectivity of replicates expressed as a percentage of total instances in a replicate set. Default is 0.05 corresponding to gold signatures.
#' @return N/A. Saves the figures
#' @keywords visualization regional miami plot
#' @export
#'

# TODO: run for xpr

split_gctx <- function(
    split.in.n    = 300,
    # File management
    parent.signature.dir = "/sc/arion/projects/va-biobank/resources/CMap/cmap_l1000_2021-11-20/",
    # Identifying Golden signatures
    # https://docs.google.com/document/d/1q2gciWRhVCAAnlvF2iRLuJ7whrGP6QjpsCMq1yWz7dU/edit#
    # is_gold: A heuristic for assessing whether a signature is reproducible and distinct. Requirements include: distil_cc_q75 >= 0.2 and pct_self_rank_q25 <= 0.05.
    distil_cc_q75_min     = 0.2,
    pct_self_rank_q25_max = 0.05,
    grep.sig.pattern      = "trt_cp|trt_sh|trt_oe|trt_xpr"
){
  ###############
  # Types of data

  types.of.data  <- c("trt_cp", 	"trt_sh", "trt_oe", "trt_xpr", "trt_misc", "ctl")
  names(types.of.data) <- c("compounds", "shRNA", "over.expression", "CRISPR",
                            "other.treatments", "negative.controls")

  ##################################
  # Identifying signature gctx files

  cmap.list <- list.files(
    parent.signature.dir,
    pattern = "gctx$",
    full.names = T)
  if (!is.na(grep.sig.pattern[1])) {
    gctx.files <- cmap.list[
      grep(
        paste(unlist(lapply(strsplit(grep.sig.pattern, "\\|"),
                            FUN = function(x) paste0("_", x, "_"))),
              collapse = "|"),
        cmap.list)]
  } else gctx.files <- cmap.list


  #############
  # Find golden
  metrics <- fread(paste0(parent.signature.dir, "/siginfo_beta.txt")) # 2021-11-20 version there are 1,201,944
  metrics$is_gold <- metrics$cc_q75 >= distil_cc_q75_min & metrics$pct_self_rank_q25 <= pct_self_rank_q25_max
  # only about 237,922 are gold and these are the ones we will save.
  gold.compounds <- metrics[(is_gold)]$sig_id


  #######################
  # Prepare the RDS files
  # name of genes

  for (thistype in types.of.data[grep(grep.sig.pattern, types.of.data)]) {
    message(paste0("Now processing: ", thistype))
    # thistype = as.character(types.of.data["other.treatments"])
    this.gctx.file = gctx.files[grep(thistype, gctx.files)]
    dsgene<-cmapR::read_gctx_ids(this.gctx.file,dim='row') # length(dsgene) #[1] 12328
    dsdrug<-cmapR::read_gctx_ids(this.gctx.file,dim='column') # name of drugs
    dsdrug <- dsdrug[dsdrug %in% gold.compounds]

    toprocess <- split(1:length(dsdrug), ceiling(seq_along(1:length(dsdrug))/split.in.n))
    names(toprocess) <- lapply(names(toprocess),FUN=function(x){
      start = ((as.numeric(x)-1)*split.in.n)+1
      stop  = start + length(toprocess[[x]]) - 1
      paste0(start, ".", stop)
    })
    # Save the matrices
    target.dir <- paste0(parent.signature.dir, "/eachDrug")
    if(!dir.exists(target.dir)) {
      dir.create(target.dir) }
    pbmclapply(names(toprocess), FUN = function(i){
      saveRDS(# get matrix and save it
        cmapR::parse_gctx(this.gctx.file, cid=dsdrug[toprocess[[i]]])@mat,
        file=paste0(target.dir, "/",thistype,"_",i,'.RDS')  )
    }, mc.cores = parallel::detectCores())
  }

}


