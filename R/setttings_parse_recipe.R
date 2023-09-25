#' Parse LSF project recipe file
#'
#' @param recipe.file recipe file
#' @return A list file with all included parameters
#' @keywords antagonism step1 LSF
#' @export

parse_recipe <- function(
    recipe.file
){
  x <- MultiWAS::return_df(recipe.file)
  recipe <- list()
  if(x[variable == "df"]$value == "NULL") {
    stop("A df value must be provided in the recipe file")
  } else { recipe$df <- x[variable == "df"]$value }
  if(x[variable == "working.directory"]$value == "NULL") {
    stop("A working.directory value must be provided in the recipe file")
  } else {  recipe$working.directory <- x[variable == "working.directory"]$value }
  recipe$column.feature         <- x[variable == "column.feature"]$value
  recipe$column.statistic       <- x[variable == "column.statistic"]$value
  recipe$column.trait           <- x[variable == "column.trait"]$value
  recipe$column.source          <- x[variable == "column.source"]$value
  recipe$results.dir            <- x[variable == "results.dir"]$value
  recipe$n.threads              <- eval(parse(text = x[variable == "n.threads"]$value))
  recipe$signature.dir          <- x[variable == "signature.dir"]$value
  recipe$gene.anno.file         <- x[variable == "gene.anno.file"]$value
  recipe$grep.sig.pattern       <- x[variable == "grep.sig.pattern"]$value
  recipe$noperm                 <- eval(parse(text = x[variable == "noperm"]$value))
  recipe$thres.N.vector         <- eval(parse(text = x[variable == "thres.N.vector"]$value))
  recipe$sig.annotation         <- eval(parse(text = x[variable == "sig.annotation"]$value))
  recipe$overwrite.intermediate <- eval(parse(text = x[variable == "overwrite.intermediate"]$value))
  recipe$model.banlist.grep     <- x[variable == "model.banlist.grep"]$value
  # recipe$prototyping            <- eval(parse(text = x[variable == "prototyping"]$value))
  return(recipe)
}
