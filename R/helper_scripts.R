# This is a list from helper scripts from the MultiWAS package, they are added
# here to remove dependencies for public release. They are not exported to keep
# the MultiWAS version if there is a conflict

gv_dir.create <- function(x, ...) {
  invisible(lapply(x, FUN = function(y){
    if (!dir.exists(y)) dir.create(y, recursive = T, ...) }
  ))
}

return_df <- function(
    x
) {
  if (!is.character(x)) {
    x <- as.data.table(x, keep.rownames = T)
  } else {
    if (grepl(".RDS$", x, ignore.case = T)) {
      x <- as.data.table(readRDS(x))
    } else x <- fread(x)
  }
  return(x)
}

my_ggscatter <- function(
    df, # dataframe
    x,  # x-axis column name
    y,  # y-axis column name
    method         = c("pearson", "spearman", "kendall"),
    color.group    = "black",
    usethislabel   = NULL,
    custom.label.x = NA,
    custom.label.y = NA, # at some point it was 30
    custom.colors  = NULL,
    plot.repel     = TRUE,
    is.lm.and.passes.0.0 = FALSE,
    use.r2 = FALSE,
    ...
) {
  df <- return_df(df)
  if (!is.lm.and.passes.0.0) {
    p1 <- ggpubr::ggscatter(
      df,
      x          = x,
      y          = y,
      add        = "reg.line",                                 # Add regression line
      color      = color.group,
      palette    = custom.colors,
      conf.int   = TRUE,                                  # Add confidence interval
      add.params = list(color = "blue", fill = "lightgray"),
      label      = usethislabel,
      repel      = plot.repel,
      ...
    ) }
  if (is.lm.and.passes.0.0) {
    p1 <- ggpubr::ggscatter(
      df,
      x          = x,
      y          = y,
      add        = "none",                                 # Add regression line
      color      = color.group,
      palette    = custom.colors,
      conf.int   = TRUE,                                  # Add confidence interval
      add.params = list(color = "blue", fill = "lightgray"),
      label      = usethislabel,
      repel      = plot.repel,
      ...
    )
    p1 <- p1 + geom_smooth(method="lm", formula=y~0+x)
  }
  if (is.na(custom.label.x)) {
    custom.label.x = 1*min(na.omit(df[[x]])) }
  if (is.na(custom.label.y)) {
    custom.label.y = 1*max(na.omit(df[[y]])) }
  if (use.r2) {
    p1 <- p1 + ggpubr::stat_cor(
      method = method[1],
      aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
      label.x = custom.label.x,
      label.y = custom.label.y)  # Add correlation coefficient
  } else {
    p1 <- p1 + ggpubr::stat_cor(
      method = method[1],
      label.x = custom.label.x,
      label.y = custom.label.y)  # Add correlation coefficient
  }
  p1
}


vector_to_colors <- function(
    x, # vector for which we need colors
    additional.colors = "viridis") {
  x <- unique(x) # only need unique colors
  list.length <- length(x)
  # figPalette <- as.vector(character()) # initialize
  # cbPalette
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  figPalette <- cbPalette
  if (list.length <= length(figPalette) - 1) {
    figPalette <- figPalette[2:(list.length+1)]
  } else { # generate more colors if cbPalette is not enough
    if (additional.colors == "viridis" | (additional.colors == "random" & list.length > 74)) { # preferred option
      figPalette <- viridis::viridis_pal(option = "D")(list.length)
      # figPalette[2:(list.length+1)] <- # replace viridis
      # viridis_pal(option = "D")(list.length)
      # figPalette <- figPalette[2:(list.length+1)]
    }  # n = number of colors seeked
    if (additional.colors == "random") {
      qual_col_pals = RColorBrewer::brewer.pal.info[brewer.pal.info$category == 'qual',]
      col_vector = unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) # 74 colors
      more.colors <- sample(col_vector, (list.length - length(figPalette)+1))
      figPalette <- c(figPalette, more.colors) }  }#; return(cbPalette)
  # all.colors <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)] # sample(all.colors,1) to get a random color out of t >430colors
  scales::show_col(figPalette)
  return(figPalette)
}

make_java_safe <- function(x) {
  # function to make names compatible with paths
  x <- gsub("(\\(|\\))","",gsub(" ","_", gsub("/","_", gsub(",","", gsub("\\'", "", x)))))
  # function to remove : character
  x <- gsub("\\:", ".",x)
  x
}

z2p <- function(Z, method = c("pnorm")) {
  if (any(is.infinite(Z))) {
    warning(
      "The 'Z' vector contains infinite values. These will be turned into NAs, because no meaninful P value can be calculated from that."
    )
    is.na(Z) <- is.infinite(Z)
  }

  if (method == "pnorm") {
    P <- exp(pnorm(abs(Z), log.p = TRUE, lower.tail = FALSE)) * 2

    if (any(P == 0, na.rm = TRUE)) {
      message("Some P-values are equal to 0. Try using the option method = 'Rmpfr::pnorm'")
    }
  }

  if (method == "Rmpfr::pnorm") {
    ## from https://stackoverflow.com/questions/46416027/how-to-compute-p-values-from-z-scores-in-r-when-the-z-score-is-large-pvalue-muc
    message("using method Rmpfr::")

    ## remove all NAs
    P <- Rmpfr::mpfr(abs(Z), precBits = 100)

    P[!is.na(Z)] <- 2 * Rmpfr::pnorm(Rmpfr::mpfr(abs(Z[!is.na(Z)]),
                                                 precBits = 100
    ),
    lower.tail = FALSE, log.p = FALSE
    )
  }

  return(P)
}
