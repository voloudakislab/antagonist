

## code to prepare `BBB.INFO.20240612` dataset goes here
BBB.regression                   = system.file(
  "extdata", "resources.CDR", "BBB", "B3DB_regression.tsv",
  package="antagonist")


BBB.classification                   = system.file(
  "extdata", "resources.CDR", "BBB", "B3DB_classification.tsv",
  package="antagonist")

sig.annotation         = SIG.INFO.20211120



### MAIN FUNCTIONS

compound_annotation <- function(
    BBB.classification = download.B3D3.classfication.file,
    BBB.regression = download.B3D3.regression.file,
    sig.annotation = SIG.INFO.20211120
){
  ##

  #BBB.classification <- BBB.classification[, c('compound_name', 'SMILES', 'CID', 'logBB', 'BBB+/BBB-')]
  #BBB.regression <- BBB.regression[, c('compound_name', 'SMILES', 'CID', 'logBB')]

  ### Add logBB ang BBB+/BBB- to the annotation
  sig.annotation.classification <- left_join(sig.annotation, BBB.classification[, c('compound_name', 'BBB+/BBB-', 'SMILES', 'CID')], by = c("pert_iname" = "compound_name"), suffix = c(".sig.annotation", ".BBB.classification"), relationship = 'many-to-many')
  sig.annotation.classification.regression <- left_join(sig.annotation.classification, BBB.regression[, c('compound_name', 'logBB', 'SMILES')], by = c('pert_iname' = 'compound_name', 'SMILES' = 'SMILES'), suffix = c(".sig.annotation.classification", ".BBB.regression"), relationship = 'many-to-many')
  #sig.annotation.classification.regression <- merge(sig.annotation.classification, BBB.regression[, c('compound_name', 'logBB', 'SMILES', 'CID')], by.x = 'pert_iname', by.y = 'compound_name', all.x = TRUE)


  ### BBB.classification also contains some logBB which are not included in BBB.regression
  different <- BBB.classification[!which(unique(BBB.classification$compound_name) %in% unique(BBB.regression$compound_name)), ]
  #sig.annotation.classification.regression <- merge(sig.annotation.classification.regression, different[, c('compound_name', 'logBB', 'SMILE', 'CID')], by.x = 'pert_iname', by.y = 'compound_name', all.x = TRUE)
  sig.annotation.classification.regression1 <- merge(sig.annotation.classification.regression, different[, c('compound_name', 'logBB', 'SMILES')], by.x = c('SMILES', 'pert_iname'), by.y = c('SMILES', 'compound_name'), all.x = TRUE)
  sig.annotation.classification.regression1$logBB.x[is.na(sig.annotation.classification.regression1$logBB.x)] <- sig.annotation.classification.regression1$logBB.y[is.na(sig.annotation.classification.regression1$logBB.x)]
  sig.annotation.classification.regression1$logBB.y <- NULL
  sig.annotation.classification.regression1 <- sig.annotation.classification.regression1[, !(names(sig.annotation.classification.regression1) %in% "logBB.y")]
  colnames(sig.annotation.classification.regression1)[colnames(sig.annotation.classification.regression1) == "logBB.x"] <- "logBB"

  #missing_logBB_SMILE <- sig.annotation.classification.regression[which(is.na(sig.annotation.classification.regression$logBB)), ]$SMILES
  #missing_logBB_index <- which(is.na(sig.annotation.classification.regression$logBB))
  #logBB_makeup <- different[which(different$SMILES %in% missing_logBB_SMILE), ]$logBB
  #for (i in missing_logBB_pert){
  #  if missing_logBB_pert[1]
  #}
  return(sig.annotation.classification.regression1)

}


BBB_ANNOTATION.INFO.20200423 <- compound_annotation(BBB.classification, BBB.regression, sig.annotation)

usethis::use_data(BBB_ANNOTATION.INFO.20200423, overwrite = TRUE)
