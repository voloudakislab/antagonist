#' Split gctx files in RDS files
#'
#' Prepare RDS files from LINCS gctx files for main pipeline. Need to run this for xpr
#'
#' @param split.in.n How many signatures in each RDS file. Default is 300
#' @param
#' @param noperm Number of permutations (the final number of permutation data points is noperm × number of drugs)
#' @param signature.dir Location of signature files. These are gw Consider saving locally to spead things up (e.g. /scratch/cmap_l1000_2021_01_28/)
#' @return N/A. Saves the figures
#' @keywords visualization regional miami plot
#' @export
#'

split_gctx <- function(
    split.in.n = 300,
    # File management

    # Identifying gold signatures


    )

#################################
# Split gctx files in RDS files #
#################################
#  Requires:
#    R >= 4
#    cmapR package either from bioconductor or github

############
# Parameters
split.in.n = 300 # split in 300 pertubagens
# Golden
# https://docs.google.com/document/d/1q2gciWRhVCAAnlvF2iRLuJ7whrGP6QjpsCMq1yWz7dU/edit#
# is_gold: A heuristic for assessing whether a signature is reproducible and distinct. Requirements include: distil_cc_q75 >= 0.2 and pct_self_rank_q25 <= 0.05.
distil_cc_q75_min     = 0.2
pct_self_rank_q25_max = 0.05
setwd("/sc/arion/projects/va-biobank/resources/CMap/cmap_l1000_2021_01_28/")
library(cmapR)
library(data.table)


#############
# Find golden
metrics <- fread("/sc/arion/projects/va-biobank/resources/CMap/cmap_l1000_2021_01_28/siginfo_beta.txt") # 1,202,656
metrics$is_gold <- metrics$cc_q75 >= distil_cc_q75_min & metrics$pct_self_rank_q25 <= pct_self_rank_q25_max
gold.compounds <- metrics[(is_gold)]$sig_id



###################
# Types of datasets
gctx.files     <- list.files(pattern = "gctx$")
types.of.data  <- c("trt_cp", 	"trt_sh", "trt_oe", "trt_xpr", "trt_misc", "ctl")
names(types.of.data) <- c("compounds", "shRNA", "over.expression", "CRISPR",
                          "other.treatments", "negative.controls")


##############
# Miscelaneous
# name of genes

for (thistype in types.of.data) {
  message(paste0("Now processing: ", thistype))
  # thistype = as.character(types.of.data["other.treatments"])
  this.gctx.file = gctx.files[grep(thistype, gctx.files)]
  dsgene<-read_gctx_ids(this.gctx.file,dim='row') # length(dsgene) #[1] 12328
  dsdrug<-read_gctx_ids(this.gctx.file,dim='column') # name of drugs
  dsdrug <- dsdrug[dsdrug %in% gold.compounds]

  toprocess <- split(1:length(dsdrug), ceiling(seq_along(1:length(dsdrug))/split.in.n))
  names(toprocess) <- lapply(names(toprocess),FUN=function(x){
    start = ((as.numeric(x)-1)*split.in.n)+1
    stop  = start + length(toprocess[[x]]) - 1
    paste0(start, ".", stop)
  })
  # Save the matrices
  if(!dir.exists("eachDrug")) dir.create("eachDrug")
  library(pbmcapply)
  pbmclapply(names(toprocess), FUN = function(i){
    saveRDS(# get matrix and save it
      parse_gctx(this.gctx.file, cid=dsdrug[toprocess[[i]]])@mat,
      file=paste0("eachDrug/",thistype,"_",i,'.RDS')  )
  }, mc.cores = 12)
}
