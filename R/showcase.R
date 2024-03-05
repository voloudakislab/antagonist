#' Showcase how the method works
#'
#' This script generates most of the figures needed for the Fig. S5 in [1]. As a
#' reminder the core method is the 5-rank method[2]. This function requires that
#' GTP or CDR has run already.
#' - Panel A: Scatter plot of disease and perturbagen effects (signature level)
#' - Panel B: Density plot of Signatures for perturbagen vs. baseline.
#' - Panel C:
#'    - Left: Scatter plot of disease and perturbagen effects (compound level)
#'    - Middle: histogram of combined z for all shRNAs (compound level)
#'    - Right: Combined z -score for prioritization.
#'
#' 1. Voloudakis, G., Vicari, J.M., Venkatesh, S. et al. A translational genomics
#' approach identifies IL10RB as the top candidate gene target for COVID-19
#' susceptibility. npj Genom. Med. 7, 52 (2022).
#' https://doi.org/10.1038/s41525-022-00324-x
#'
#' 2. So HC et al. Analysis of genome-wide association data highlights
#' candidates for drug repositioning in psychiatry. Nat Neurosci. 2017
#' Oct;20(10):1342-1349. PMID: 28805813
#'
#' TODO: Pending improvements:
#' - Add Panel C
#' - source the n from the column names for non stardard analyses
#' - use this script as guide: /lab/r.projects/COVID_drug_repurposing/4.drug_repurp_explanation_Fig.S5.R
#'
#' @param perturbagen pert_iname; e.g. "benactyzine"
#' @param twas.model GFI model name  (can also be DGE experiment type or cell type)
#' @param twas.trait Trait or GWAS
#' @param twas Disease signature data frame or file (default is the output of the MultiWAS pipeline)
#' @param twas.feature.column.name genomic feature column name (default: feature). Expects ENSEMBL IDs and filters out everything that doesn't start with ENSG.
#' @param twas.feature.name.column.name genomic feature column name (default: feature_name). Friendly name
#' @param twas.zscore.column.name z-score statistic column name (default: zscore)
#' @param twas.model.column.name GFI model (can also be DGE experiment type or cell type) column name (default: model_ID)
#' @param twas.trait.column.name Trait or GWAS column name (default: gwas)
#' @param sig.info File with perturbagen signature information (default: SIG.INFO.20211120)
#' @param gtp.cdr.dir Directory where GTP/CDR saves its results (default: "results/GTP_CDR/")
#' @param gene.anno.file Gene annotation file (provided from L1000)
#' @return N/A. Saves the figures and csv file
#' @keywords visualization cdr gtp showcase
#' @export

showcase_method_cdr_gtp <- function(
    # Perturbagen and context selection
  perturbagen                     , # e.g. "benactyzine"
  twas.model                    = NA,
  twas.trait                    = NA,
  perturbagen.sig_id            = NA,
  # Disease signature
  twas                          = "output/2.METAXCAN/CDR_Colllaborative_Merit_df.all.annotated.onlytranscripts.csv.gz", # conventionally TWAS
  twas.feature.column.name      = "feature", # Has to be ENSEMBL ids so that it can also keep only the genes
  twas.feature.name.column.name = "feature_name",
  twas.zscore.column.name       = "zscore",
  twas.model.column.name        = "model_ID",
  twas.trait.column.name        = "gwas",
  # Peturbagen signature
  # GTP/CDR
  sig.info                      = SIG.INFO.20211120,
  # signature.inventory           = "results/GTP_CDR/intermediate.files/signature.inventory.csv",
  gtp.cdr.dir                   = "results/GTP_CDR/",
  gene.anno.file                = "/sc/arion/projects/va-biobank/resources/CMap/cmap_l1000_2021-11-20/geneinfo_beta.txt"
){
  # finding what the signature is
  sig.info <- MultiWAS::return_df(sig.info)
  sig.info <- sig.info[pert_iname == perturbagen]
  perturbagen.type <- unique(sig.info$pert_type)

  # Process TWAS
  twas = MultiWAS::return_df(twas)
  twas$feature      <- twas[[twas.feature.column.name]]
  twas$feature_name <- twas[[twas.feature.name.column.name]]
  twas$zscore       <- twas[[twas.zscore.column.name]]
  twas$model_ID     <- twas[[twas.model.column.name]]
  twas$gwas         <- twas[[twas.trait.column.name]]
  twas              <- twas[grep("^ENSG", feature)] # limit analysis to genes
  ## Handle gwas
  if (length(unique(twas$gwas)) > 1) { stop(paste0(
    "There are more than one gwas/traits, please provide one of: \n",
    paste(unique(twas$gwas), collapse = "\n"), "\n",
    "Provide it as a twas.trait variable" )) } else {
      if (is.na(twas.trait)) twas.trait <- unique(twas$gwas)
      twas <- twas[gwas == twas.trait]
    }
  ## Handle model_ID
  if (length(unique(twas$model_ID)) > 1) { stop(paste0(
    "There are more than one models/conditions, please provide one of: \n",
    paste(unique(twas$model_ID), collapse = "\n"), "\n",
    "Provide it as a twas.model variable" )) } else {
      if (is.na(twas.model)) twas.model <- unique(twas$model_ID)
      twas <- twas[model_ID == twas.model]
    }

  # Handle sig_id
  if (length(perturbagen.sig_id)>1) stop("only one perturbagen.sig_id is supported") else {
    if (is.na(perturbagen.sig_id)) {
      message("No sig_id was provided, the best one will be used to showcase method.")
      # then choose the best signature
      sigs <- fread(paste0(gtp.cdr.dir, unique(twas$gwas),"/ALL/", unique(twas$gwas),
                           "_", ifelse(perturbagen.type == "trt_cp", "cdr", "gtp"),
                           "_signature_level.csv.gz" ))
      perturbagen.sig_id <- sigs[
        model_ID == twas.model & gwas == twas.trait & pert_iname == perturbagen][
          AvgRank == min(AvgRank)]$sig_id[1] # adding [1] in case of a tie
    } }

  # Get compound info
  signature.location <- fread(paste0(
    results.dir, "intermediate.files/signature.location.csv"))
  i <- signature.location[sig_id == perturbagen.sig_id]$filename
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
  cmap <- cmap[, c("feature", perturbagen.sig_id), with = F]
  sig.info[sig_id == perturbagen.sig_id]

  # Prepare the TWAS for comparison
  twas = twas[order(abs(zscore),decreasing = T)]
  twas$rank  = rev(rank(abs(twas$zscore)))
  twas$Group = ifelse(twas$rank <= 50, "top 050",
                      ifelse( twas$rank <= 100, "top 100",
                              ifelse(twas$rank <= 250, "top 250",
                                     ifelse(twas$rank <= 500, "top 500",
                                            "Top not"
                                     ) ) ) )
  twas.vs.compounds = as.data.table(dplyr::inner_join(twas, cmap))
  names(twas.vs.compounds) <- gsub(":|-", ".", names(twas.vs.compounds))

  library(ggpubr)
  p1 <- MultiWAS::my_ggscatter(twas.vs.compounds,
                               x              = "zscore",
                               y              = names(twas.vs.compounds)[ncol(twas.vs.compounds)] ,
                               color.group    = "Group",
                               custom.colors  = c("black", MultiWAS::vector_to_colors(unique(twas.vs.compounds$Group))),
                               # usethislabel   = "Text",
                               xlab           = unique(paste0("TWAS zscore for ", twas.vs.compounds$tissue)),
                               ylab           = paste0(perturbagen, " signature for ", sig.info[sig_id == perturbagen.sig_id]$cell_id))
  this.filename <- paste0(results.dir, "intermediate.files/",
                          MultiWAS::make_java_safe(twas.trait), ".",
                          MultiWAS::make_java_safe(twas.model), ".",
                          MultiWAS::make_java_safe(perturbagen), ".",
                          MultiWAS::make_java_safe(perturbagen.sig_id))
  ggsave((paste0(this.filename, ".pdf")), p1, width = 6, height = 6)
  ggsave((paste0(this.filename, ".png")), p1, width = 6, height = 6, units = "in")
  fwrite(twas.vs.compounds, paste0(this.filename, ".csv"))

}
