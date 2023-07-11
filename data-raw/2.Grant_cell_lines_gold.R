library(data.table)
x <- fread("Gold.signatures.csv")
x <- xp[pert_type == "trt_cp"]
nrow(x)
barplot(sort(table(x$cell_mfc_name), decreasing = T))
y <- as.data.table(sort(table(x$cell_mfc_name)))
y <- y[order(y$N, decreasing=T)]
nrow(y)
sum(y[1:25]$N)/nrow(x)
