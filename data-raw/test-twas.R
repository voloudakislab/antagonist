## code to prepare `test.twas` dataset goes here

library(MultiWAS)
test.twas  <- return_df("/lab/r.projects/PolyXcan_v2/output/2.METAXCAN/PolyXcan_DLPFC_v1_df.all.annotated.csv.gz")
test.twas  <- test.twas[model_ID == "Brain: DLPFC (PEC) genes :: Homogenate :: Genes :: EpiXcan"]

usethis::use_data(test.twas, overwrite = TRUE)
