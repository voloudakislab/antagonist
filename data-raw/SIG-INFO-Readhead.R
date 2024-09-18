## code to prepare SIG.INFO fles for Readhead dataset goes here

## SIG.INFO.Readhead.CLL.control.SCZ


#
# SIG.INFO.Readhead.CLL.control.SCZ <- fread('/sc/arion/projects/va-biobank/PROJECTS/marios_temp/Resources/Readhead/no_shsy5y_siginfo/rh_siginfo_CLL_control_SCZ.csv')
# SIG.INFO.Readhead.CLL.control.SCZ <- SIG.INFO.Readhead.CLL.control.SCZ[
#   ,c("sig_id", "pert_id", "cmap_name", "cell_iname", "pert_dose", "pert_dose_unit",
#      "pert_time", "pert_time_unit", "pert_type", "is_exemplar_sig", "is_ncs_sig",
#      "is_null_sig")]
# setnames(SIG.INFO.Readhead.CLL.control.SCZ,
#          c("cmap_name", "cell_iname"),
#          c("pert_iname", "cell_id"))
#
# usethis::use_data(SIG.INFO.Readhead.CLL.control.SCZ, overwrite = TRUE)


## SIG.INFO.Readhead.CLL.control.SCZ



## SIG.INFO.Readhead.CLL.control.SCZ
