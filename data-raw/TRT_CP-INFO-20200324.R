## code to prepare `TRT_CP.INFO.20200324` dataset goes here

TRT_CP.INFO.20200324                   = system.file(
  "extdata", "resources.CDR", "repurposing_drugs_20200324.txt",
  package="antagonist")

# CDR repurposing candidates
# TODO: is there a more updated list?
repurposing_candidates <- function(
    drug.file = ref.drug.file
){
  repurp_candidates <- fread(cmd=paste0("grep -v '^!' ", drug.file))
  repurp_candidates[moa == "", moa := "unknown"]
  repurp_candidates[duplicated(pert_iname)] # no duplicates
  final.repurp.candidates <- repurp_candidates
  return(final.repurp.candidates)
}
TRT_CP.INFO.20200324 <- repurposing_candidates(TRT_CP.INFO.20200324)

usethis::use_data(TRT_CP.INFO.20200324, overwrite = TRUE)
