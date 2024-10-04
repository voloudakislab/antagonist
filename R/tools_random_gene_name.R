#' Generate random gene names
#'
#' @param how.many.samples how many random genes to return
#' @param how.many.letters how many letters each gene
#' @param thisseed can use seed to avoid occurance of inappropriate words
#' @return returns vector
#' @export
#' @examples


random_gene_name <- function(
    how.many.samples = 100,
    how.many.letters = 3,
    thisseed = NULL)  {
  set.seed(thisseed)
  letter.vector <- lapply(
    1:how.many.samples,
    function(i) {
      toupper(sample(letters, size = how.many.letters))
    })
  unlist(lapply(letter.vector, FUN = function(i) paste0(i, collapse = "")))
}
