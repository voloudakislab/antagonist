# wrapper_main


# Script
# parse args..
opt = parse(recipe.file)

results.dir <- opt$results.dir
parent.subdir <- opt$parent.subdir

setwd(opt$working.directory)

logs_to_check = ma_paste0(file.path(results.dir, parent.subdir, '/logs/S01A/', thisgwas))

check.command = paste0(logs_to_check, ' -type f -name "*.out" -exec grep -iL "successfully completed" {} \;')
system(check.command)
