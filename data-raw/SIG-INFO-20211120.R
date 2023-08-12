## code to prepare `SIG.INFO.20211120` dataset goes here

SIG.INFO.20211120 <- fread("/sc/arion/projects/va-biobank/resources/CMap/cmap_l1000_2021-11-20/siginfo_beta.txt")
SIG.INFO.20211120 <- SIG.INFO.20211120[,c("sig_id", "pert_id", "cmap_name", "cell_iname", "pert_dose", "pert_time", "pert_type")]
# For compatibility with previous versions
setnames(SIG.INFO.20211120,
         c("cmap_name", "cell_iname"),
         c("pert_iname", "cell_id"))


usethis::use_data(SIG.INFO.20211120, overwrite = TRUE)
