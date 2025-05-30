#' pvalue_qqplot
#'
#' Create a quantile-quantile plot with ggplot2. Script was adapted from: https://slowkow.com/notes/ggplot2-qqplot/. The null is the ordinates for probability plotting Generates the sequence of probability points (1:m - a)/(m + (1-a)-a) as generated by ppoints similarly to the one used in qqplot and qqnorm.
#'
#' Future improvements
#' - Stouffer's method or other zscore based meta-analysis
#' - Scaling prior to getting z combined
#' - Only label if both sig in TWAS and antagonism, otherwise grey
#'
#'
#' Assumptions:
#'   - Expected P values are uniformly distributed.
#'   - Confidence intervals assume independence between tests.
#'     We expect deviations past the confidence intervals if the tests are
#'     not independent.
#'     For example, in a genome-wide association study, the genotype at any
#'     position is correlated to nearby positions. Tests of nearby genotypes
#'     will result in similar test statistics.
#'
#' @param thistrait trait name (gwas column for TWAS)
#' @param thistissue tissue or cell type of choice (doesn't support multiple cell types at this point)
#' @param stouffer.ma perform Stouffer meta-analysis between TWAS and antagonism otherwise it is just z twas + z gtp which is just for visualization (not valid approach)
#' @param z.score.scaling scale the combined z scores. "before" MA, "after" MA or "FALSE" (after MA)
#' @param use.sample.data use sample data to play with plot parameters
#' @param twas.df data frame or file with the TWAS / disease signature data frame
#' @param twas.feature.name.column Column name for feature name in disease df. Default is "feature_name"
#' @param twas.trait.column Column name for trait/gwas in disease df. Default is "gwas"
#' @param twas.tissue.column Column name for tissue/cell type in disease df. Default is "model_ID"
#' @param twas.zscore.column Column name for zscore in disease df. Default is "zscore"
#' @param twas.pvalue.column Column name for pvalue in disease df. Default is "pvalue"
#' @param twas.FDR.column Column name for FDR in disease df. Default is "fdr_trait_model"
#' @param gtp.df data frame or file with the antagonist signature data frame
#' @param gtp.feature.name.column Column name for feature name in gtp df. Default is "pert_iname"
#' @param gtp.zscore.column Column name for zscore in gtp df. Default is "Compound.pseudo.zscore"
#' @param gtp.FDR.column Column name for FDR in gtp df. Default is "Compound.MW.FDR"
#' @param stat.prefix text for the method or title of the qqplot
#' @param ppoints.a set ppoints a parameter
#' @param figsDir Default is "figs/"
#' @param figure.name Default is "qqplot_TWAS_vs_FunWAS"
#' @param save_plot Save plots in various useful formats
#' @param plot.height plot heighth (inches)
#' @param plot.width plot width (inches)
#' @param add.lambda add lambda to the plot
#' @return returns ggplot graph
#' @export
#' @examples
#' library(ggplot2)
#' gg_qqplot(runif(1e2)) + theme_grey(base_size = 24)

gtp_pvalue_qqplot  <- function(
    ### Parameters
  thistrait       = NA,
  thistissue      = NA,
  stouffer.ma     = TRUE,
  z.score.scaling = "before",
  use.sample.data = FALSE,

  ### TWAS information
  twas.df            = paste0(
    "output/2.METAXCAN/", basename(getwd()),
    "_df.all.annotated.onlytranscripts.csv.gz"),
  twas.feature.name.column = "feature_name",
  twas.trait.column = "gwas",
  twas.tissue.column = "model_ID",
  twas.zscore.column = "zscore",
  twas.pvalue.column = "pvalue",
  twas.FDR.column    = "fdr_trait_model",

  ### Antagonist gtp output
  gtp.df                   = NA,
  gtp.feature.name.column  = "pert_iname",
  gtp.zscore.column        = "Compound.pseudo.zscore",
  gtp.FDR.column           = "Compound.MW.FDR",

  ### Graph settings
  ci                       = 0.95,
  stat.prefix              = "",
  ppoints.a                = NA,
  figsDir                  = "figs/",
  figure.name              = "gtp_actions", # main body of the name of the figure
  save_plot                = TRUE,
  plot.height              = 5,
  plot.width               = 5,
  add.lambda               = FALSE
) {
  library(ggplot2)
  # library(dplyr)
  # library(data.table)
  `%!in%` = Negate(`%in%`)

  # Prepare the TWAS df
  if (use.sample.data) { # for optimizing the plot
    twas.df <- data.table(
      "feature_name"       = random_gene_name(thisseed = 12345),
      "gwas"               = "sample.gwas.1",
      "model_ID"           = "sample.tissue.1",
      "zscore"             = c(rnorm(99, sd = 2), 10),
      "fdr_trait_model"    = c(ppoints(99), 1e-10)
    )
  } else {
    twas.df              <- return_df(twas.df)
  }
  setnames(twas.df, twas.feature.name.column,       "feature_name")
  setnames(twas.df, twas.trait.column,  "gwas")
  setnames(twas.df, twas.tissue.column, "model_ID")
  setnames(twas.df, twas.zscore.column, "zscore")
  # handle twas FDR column
  if (twas.FDR.column %!in% names(twas.df)) {
    comb.grid <- unique(twas.df[, c("gwas", "model_ID")])
    setnames(comb.grid, c("model_ID"), c("ID"))
    twas.df <- do.call(
      rbind,
      pbmclapply(
        1:nrow(comb.grid),
        FUN = function(i) {
          x   <- twas.df[gwas == comb.grid$gwas[i] & model_ID == comb.grid$ID[i]]
          x$fdr_trait_model <- p.adjust(x[[twas.pvalue.column]],method = "fdr")
          return(x)
        }
      ))
  } else {
    setnames(twas.df, twas.FDR.column,    "fdr_trait_model")
  }
  if (is.na(thistrait)) {
    if (length(unique(twas.df$gwas)) == 1) {
      thistrait <- unique(twas.df$gwas)
    } else {
    stop("The disease signature file has more than one tissues/cell types.\n
                You need to specify the one of interest...")
    } }
  if (is.na(thistissue)) {
    if (length(unique(twas.df$model_ID)) == 1) {
      thistissue <- unique(twas.df$model_ID)
    } else {
      stop("The disease signature file has more than one tissues/cell types.\n
                You need to specify the one of interest...")
    } }
  twas.df <- twas.df[gwas == thistrait & model_ID == thistissue]

  ### Prepare the gtp df ##
  if (use.sample.data) { # for optimizing the plot
    gtp.df = data.table(
      "pert_iname"              = random_gene_name(thisseed = 12345),
      "Compound.pseudo.zscore"  = c(rnorm(99, sd = 2), 10),
      "Compound.MW.FDR"         = c(ppoints(99), 1e-10)
    )
  } else {
    if (is.na(as.character(gtp.df [1])[1])) {
      # if it is NA retrieve it
      gtp.df <- return_df(paste0(
        "results/GTP_CDR/",
        thistrait, "/",
        make_java_safe(thistissue), "/",
        thistrait, "_gtp_compound_level.csv"
      ))
    } else {
      # if it is not NA just return the data.frame
      twas.df  <- return_df(twas.df)
    } }
  setnames(gtp.df, gtp.feature.name.column,  "feature_name") # for joining
  setnames(gtp.df, gtp.zscore.column,        "Compound.pseudo.zscore")
  setnames(gtp.df, gtp.FDR.column,           "Compound.MW.FDR")

  ### Prepare the merged df
  ps <- as.data.table(dplyr::inner_join(twas.df,gtp.df))

  if (z.score.scaling == "before") {
    ps[, zscore := scale(zscore)]
    ps[, Compound.pseudo.zscore:=scale(Compound.pseudo.zscore)]
  }
  if (stouffer.ma) {
    ps[, z.combined := (zscore + Compound.pseudo.zscore) / sqrt(2) ]
  } else {
    ps[, z.combined := zscore + Compound.pseudo.zscore]
  }
  if (z.score.scaling == "after") {
    ps[, z.combined := scale(z.combined)]
  }
  ps[, p.combined := z2p(z.combined)]
  ps[, FDR.combined := p.adjust(p.combined, method = "BH")]
  ps[, intervention := ifelse(z.combined*Compound.pseudo.zscore < 0, "unclear", ifelse(z.combined > 0, "downregulation", "upregulation"))]
  # 24: upregulation arrow up, 25 downregulation arrow down, 21 unclear circle
  ps <- ps[order(p.combined)]
  fwrite(ps, paste0(
    "results/GTP_CDR/",thistrait, "/",
    make_java_safe(thistissue), "/prioritization",
    ifelse(stouffer.ma, ".stouffer", ""),
    ifelse(isFALSE(z.score.scaling), "",
           ifelse(z.score.scaling == "before", ".scaled.before", "scaled.after")),
    ".csv"))

  ### Calculate inflation λ ###
  gv_lambda_pvalue <- function(pvalues.atomic.vector) {
    set.seed(12345)
    # chisq <- qchisq(1-pvalues.atomic.vector, 1) # change this, is not precise enough
    chisq <- qchisq(pvalues.atomic.vector, 1, lower.tail=FALSE) # use this
    lambda <- median(chisq) / qchisq(0.5, 1)
    return (lambda)   }

  ### Prepare the values ###
  stat.prefix <- paste0(stat.prefix, " ") # add some space

  # # P value sanity check
  # if (type[1] == "pvalue") {
  #   if(!all(ps >= 0 & ps <= 1))
  #     stop("You say you provided p values but the value range is not [0,1]")}

  # # Shape the df if we provided zscores
  # if (type[1] == "zscore"){
  #
  #   # z score to p value function
  #   zscore_to_pvalue <- function(z, min.p = 1e-300) {
  #     options(digits = 22)
  #     p <- 2*pnorm(-abs(z))
  #     p[p==0] = min.p # otherwise returns 0
  #     return(p) }
  #   ps             <- zscore_to_pvalue(ps)
  # }

  # Generate df for plotting

  n              <- nrow(ps)
  ps[, observed := -log10(p.combined)]
  ps[, expected := -log10(ppoints(
    n, a = ifelse(is.na(ppoints.a), ifelse(n <= 10, 3/8, 1/2), ppoints.a)))]
  ps[, clower := -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1))]
  ps[, cupper := -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))]
  df <- ps
  log10Pe <- bquote(paste("Expected ", .(stat.prefix), "-log"[10], plain(P), sep = ""))
  log10Po <- bquote(paste("Observed ", .(stat.prefix), "-log"[10], plain(P), sep = ""))

  ps[, FDR := colourvalues::colour_values(FDR.combined, palette = "magma") ]
  ps[, alpha.value := ifelse(fdr_trait_model < 0.05 & Compound.MW.FDR < 0.05, 1, 0.2)]
  # ps[, alpha.label := ifelse(fdr_trait_model < 0.05 & Compound.MW.FDR < 0.05, TRUE, FALSE)]
  ps[, label := ifelse(fdr_trait_model < 0.05 & Compound.MW.FDR < 0.05 & FDR.combined < 0.05, feature_name, NA)]

  df <- ps



  library(colourvalues)
  library(ggplot2)
  p1 <-  ggplot(df) +
    geom_ribbon(
      mapping = aes(x = expected, ymin = clower, ymax = cupper),
      alpha = 0.1) +
    geom_point(aes(expected, observed, fill = -log10(FDR.combined),
                   alpha = alpha.value, shape = intervention ), #shape = 16,
               size = 3) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    xlab(log10Pe) +
    ylab(log10Po) +
    theme_classic(base_size = 10) +
    scale_x_continuous(expand = c(0, .1)) + # expands only on the right side by 10%
    scale_y_continuous(expand = c(0, .15)) +
    scale_fill_viridis_c(
      bquote(-log[10]~"(GTP FDR)"),
      # "minus log10\nCombined FDR",
      option = "plasma", direction = -1) +
    scale_alpha_identity(name = "FDR-sig disease\nand anatagonism",
                         breaks = c(0.2,1),
                         labels = c("FALSE", "TRUE"),
                         guide = "legend") +
    scale_shape_manual(
      name   = "Recommended\nintervention",
      # breaks = c(24,25,21),
      # labels = c("upregulation", "downregulation", "unclear"),
      # guide  = "legend"
      values = c("upregulation" = 24, "downregulation" = 25, "unclear" = 21)#,
      # guide = guide_legend(order = 2)
      # 24: upregulation arrow up, 25 downregulation arrow down, 21 unclear circle
    )+
    ggrepel::geom_text_repel(
      aes(expected, observed, label = label), min.segment.length = 0, box.padding = 1,
      fontface = "italic", # for italic
      show.legend = FALSE, max.overlaps = Inf) + #this removes the 'a' from the legend
    theme(legend.position="right")
  if(add.lambda) p1 <- p1 +
    annotate(
      geom  = "text",
      x     = -Inf,
      y     = Inf,
      hjust = -0.15,
      vjust = 1 + 0.15 * 3,
      label = sprintf("λ = %.2f", gv_lambda_pvalue(df$p.combined)) #,
      # size = 8
    )
  cowplot::ggsave2(
    filename = paste0("results/GTP_CDR/",thistrait, "/",
                      make_java_safe(thistissue), "/prioritization",
                      ifelse(stouffer.ma, ".stouffer", ""),
                      ifelse(isFALSE(z.score.scaling), "",
                             ifelse(z.score.scaling == "before", ".scaled.before", "scaled.after")),
                      ".pdf"),
    plot   = p1,
    width  = plot.width,  height = plot.height )
  cowplot::ggsave2(
    filename = paste0("results/GTP_CDR/",thistrait, "/",
                      make_java_safe(thistissue), "/prioritization",
                      ifelse(stouffer.ma, ".stouffer", ""),
                      ifelse(isFALSE(z.score.scaling), "",
                             ifelse(z.score.scaling == "before", ".scaled.before", "scaled.after")),
                      ".png"),
    plot   = p1,
    width  = plot.width,  height = plot.height, dpi = 300 )
  return(p1)

}
