## code to prepare `SIG.INFO.20240729` dataset goes here

SIG.INFO.20240729 <- fread('/sc/arion/projects/va-biobank/PROJECTS/cdr.comparative.efficacy.marios/Resources/updated_readhead_siginfo.csv')
SIG.INFO.20240729 <- SIG.INFO.20240729[
  ,c("sig_id", "pert_id", "cmap_name", "cell_iname", "pert_dose", "pert_dose_unit",
     "pert_time", "pert_time_unit", "pert_type", "is_exemplar_sig", "is_ncs_sig",
     "is_null_sig")]
setnames(SIG.INFO.20240729,
         c("cmap_name", "cell_iname"),
         c("pert_iname", "cell_id"))

usethis::use_data(SIG.INFO.20240729, overwrite = TRUE)





