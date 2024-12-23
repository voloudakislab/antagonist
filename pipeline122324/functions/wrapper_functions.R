
prepareTWAS.main = function(TWASList){
  extractColumnNames <- function(columnNamesStr){
    
    # clear quotas, quotas help in csv format. clear 'spaces' (better string manipulation)
    columnNamesStr = gsub('"', '', columnNamesStr)
    columnNamesStr = gsub(' ', '', columnNamesStr)
    pairedColumnNames = unlist(strsplit(columnNamesStr, split = ','))
    
    columnNamesList = list() 
    for(thispair in pairedColumnNames){
      sublistTitle = unlist(strsplit(thispair, split = '='))[1]
      columnNamesList[[sublistTitle]] = unlist(strsplit(thispair, split = '='))[2]
    }
    return(columnNamesList)
  }
  
  df = MultiWAS::return_df(TWASList$path)
  
  # input string: "feature = feature, trait = gwas, statistic = zscore, source = model_ID"
  # output list: columnNames where columnNames$feature = "feature", columnNames$trait = "gwas" , columnNames$statistic = "zscore"..
  columnNames = extractColumnNames(TWASList$columnNames)
  
  model.banlist.grep = TWASList$model.banlist
  
  # Check if the names in columnNames match the colnames(df)
  # ex: if columnNames$statistic is not a column name in df, then df[[columnNames$statistic]] will result in logical error.
  if(!all(c(columnNames$feature, columnNames$statistics, columnNames$trait, columnNames$source) %in% colnames(df))) stop('From prepareTWAS.main(): 
                                                                                                                         the "columnNames" variable by recipe.csv provided names that 
                                                                                                                         do not match the actual column names. Consider
                                                                                                                         checking BOTH the provided "df" (the twas/dge dataframe)
                                                                                                                         and the "columnNames" variable.')
  
  df$feature   <- df[[columnNames$feature]]
  df$zscore    <- df[[columnNames$statistic]]
  df$gwas      <- df[[columnNames$trait]]
  df$model_ID  <- df[[columnNames$source]]
  df           <- df[zscore!=-Inf & zscore!=Inf]
  

  df <- df[, c("feature", "zscore", "gwas", "model_ID")]
  
  if (length(grep(model.banlist.grep, df$model_ID))>0) {
    message(paste0("Currently expression banlist is: ", model.banlist.grep))
    message("Intention is to only keep genes")
    df <- df[-grep(model.banlist.grep, model_ID)]
  }
  return(df)
}

create_df.shape.main = function(recipe.obj, df, thisgwas, thismodelID){
  results.dir = recipe.obj@results.dir
  df = df[gwas == thisgwas & model_ID == thismodelID]
  
  thisgwas <- MultiWAS::make_java_safe(thisgwas)
  thismodelID <- MultiWAS::make_java_safe(thismodelID)
  
  df$gwas = thisgwas
  df$model_ID = thismodelID
  
  # this is init.job specific
  MultiWAS::gv_dir.create(ma_paste0(file.path(results.dir, "intermediate.files/df.shapes/", thisgwas)))
  
  # This is used as reference file that the core functions source.
  fwrite(df, paste0(results.dir, "/intermediate.files/df.shapes/", thisgwas, '/', thismodelID, '.shaped.csv.gz'))
}

ReturnExportSignatures.main <- function(siglist, recipe.object){
  
  results.dir <-      recipe.object@results.dir
  signature.dir <-    siglist$signature.dir
  grep.sig.pattern <- siglist$grep.sig.pattern
  overwrite.intermediate <- siglist$overwrite.intermediate
  prototyping <-      siglist$prototyping
  
  MultiWAS::gv_dir.create(ma_paste0(file.path(results.dir, "intermediate.files")))
  
  cmap.list <- list.files(signature.dir, full.names = T)
  # cmap=readRDS("/sc/hydra/projects/roussp01b/Wen/LINCS_Level5/eachDrugPhase2/1-300_cp.RDS")
  # to.process <- as.data.table(tidyr::crossing(to.process, data.table("cmap.file" = cmap.list)))
  if (!is.na(grep.sig.pattern[1])) mylist <- cmap.list[grep(grep.sig.pattern, cmap.list)] else mylist <- cmap.list
  # Create a signature inventory
  if (!file.exists(paste0(results.dir, "intermediate.files/signature.inventory.csv")) | overwrite.intermediate) {
    message("Building signature inventory by perturbagen type")
    signature.inventory <- pbmclapply(
      stats::setNames(mylist, sub("_[[:digit:]]+\\.[[:digit:]]+\\.RDS", "", basename(mylist)) ),
      FUN = function(x) {colnames(readRDS(x))},
      mc.cores = parallel::detectCores() - 2 
    )
    message('First signature.inventory works')
    
    signature.inventory <- do.call(rbind, lapply(
      unique(unlist(names(signature.inventory))),
      FUN = function(x) {
        data.table(
          "Signature.type" = x,
          "sig_id" = unique(unlist(signature.inventory[grep(x, names(signature.inventory))]))
        ) } ) )
    fwrite(signature.inventory, ma_paste0(file.path(results.dir, "intermediate.files/signature.inventory.csv")))
  } else {
    message("Loading signature inventory by perturbagen type")
    signature.inventory <- fread(ma_paste0(file.path(results.dir, "intermediate.files/signature.inventory.csv")))
  }
  if (!is.na(prototyping)) {
    mylist <- unlist(
      lapply(
        unique(signature.inventory$Signature.type),
        FUN = function(x) {
          mylist[grep(x, mylist)][seq(prototyping)]
        }))
  }
  saveRDS(mylist, ma_paste0(file.path(results.dir, "intermediate.files/signature.inventory.list.RDS")))
  
  return(mylist)
}
#
