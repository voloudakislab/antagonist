## code to prepare `AD.MG.TWAS` dataset goes here

AD.MG.TWAS <- fread("data-raw/CDR_Colllaborative_Merit_df.all.annotated.onlytranscripts.csv.gz")
usethis::use_data(AD.MG.TWAS, overwrite = TRUE)
