# code to prepare SIG.INFO fles for Readhead dataset goes here

# File locations
files <- list(
SIG.INFO.Readhead.CCL.NPC.no.SHSY5Y          = "/sc/arion/projects/va-biobank/PROJECTS/marios_temp/Resources/Readhead/091824_sig_info_files/SIG.INFO.Readhead.CCL.NPC.no.SHSY5Y.rda",
SIG.INFO.Readhead.CCL.NPC                    = "/sc/arion/projects/va-biobank/PROJECTS/marios_temp/Resources/Readhead/091824_sig_info_files/SIG.INFO.Readhead.CCL.NPC.rda",
SIG.INFO.Readhead.CCL.SCZ.controls.no.SHSY5Y = "/sc/arion/projects/va-biobank/PROJECTS/marios_temp/Resources/Readhead/091824_sig_info_files/SIG.INFO.Readhead.CCL.SCZ.controls.no.SHSY5Y.rda",
SIG.INFO.Readhead.CCL.SCZ.controls           = "/sc/arion/projects/va-biobank/PROJECTS/marios_temp/Resources/Readhead/091824_sig_info_files/SIG.INFO.Readhead.CCL.SCZ.controls.rda"
)

# Load and bring them to global environment
sapply(
  names(files),
  FUN = function(i) {
    load(as.character(files[[i]]))
    # message(i)
    x <- MultiWAS::return_df(get(i))
    # print(x)
    x <- x[
      ,c("sig_id", "pert_id", "pert_iname", "cell_id", "pert_dose", "pert_time",
         # "pert_dose_unit", "pert_time_unit",  "is_exemplar_sig", "is_ncs_sig", "is_null_sig", # these columns are not included
         "pert_type")]
    # get(i) <- x
    assign(i, x, envir=.GlobalEnv)
    # print(sub("^<environment: (.+)>$", "\\1", capture.output(print(environment()))))
    # usethis::use_data(i, overwrite = TRUE)
  }
)

# Save files
library(usethis)
for (i in names(files)) do.call("use_data", list(as.name(i)) )#, overwrite = TRUE)
