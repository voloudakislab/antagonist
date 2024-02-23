## code to prepare `CELL.LINE.INFO` dataset goes here

ref.cell.file <- system.file(
  "extdata", "resources.CDR", "cell_lines.csv",
  package="antagonist")

# Cell line information.
# TODO: is there a more updated list?
repurposing_cell_lines <- function(
    cell.file = ref.cell.file
){
  return(MultiWAS::return_df(cell.file)[, c("cell_id", "cell_desc")])
}

CELL.LINE.INFO <- repurposing_cell_lines()

usethis::use_data(CELL.LINE.INFO, overwrite = TRUE)
