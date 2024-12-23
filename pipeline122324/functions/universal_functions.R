setGeneric('submitJob.class', function(object, recipe.obj) standardGeneric('submitJob.class'))

setGeneric('create_bsub', function(object, ...) standardGeneric('create_bsub'))

setGeneric('create.core.command', function(object, ...) standardGeneric('create.core.command'))

setGeneric('collect_submitJob.class', function(object, ...) standardGeneric('collect_submitJob.class'))

#####################
# NON-SPECIFIC FUNCTIONS

# protects from 'double-dash' errors in bash
ma_paste0 = function(...) return(gsub('/+', '/', paste0(...)))

# c(1, 2, 3) %!in% c(2, 3, 4) returns TRUE FALSE FALSE
`%!in%` <- function(x,y) !(x %in% y)

parse_recipe <- function(recipe.file){
  df <- MultiWAS::return_df(recipe.file)
  recipe <- list()
  for(var in df$variable){
    if(var %in% c('overwrite.intermediate', 'dryrun', 'prototyping')){recipe[[var]] = eval(parse(text = df[variable == var]$value))
    } else recipe[[var]] <- df[variable == var]$value
  }
  return(recipe)
}

extractParameters <- function(JobSpecParameters){
  
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

injectCollectors.uni <- function(combo.table, type_filter = 'collector'){
  injectTotalModelsPerGwas <- function(object, models.per.gwas){
    object@total.models.per.gwas <- models.per.gwas
    return(object)
  }
  injectIDsCollected <- function(object, ids.collected){
    object@ids.collected <- ids.collected
    return(object)
  }
  
  injectModelsCollected <- function(object, models.collected){
    object@models.collected <- models.collected
    return(object)
  }
  
  # models.per.gwas and ids.collected are lists that contain the total number of model_IDs per gwas and the ids collected for each gwas
  models.per.gwas = list()
  ids.collected = list()
  models.collected = list()
  for(thisgwas in unique(combo.table$gwas)){
    models.per.gwas[[thisgwas]] = nrow(combo.table[gwas == thisgwas])
    ids.collected[[thisgwas]] = character()
    models.collected[[thisgwas]] = character()
  } 
  
  # all objects are sourced for global env, those of class 'collector' will be injected with models.per.gwas and ids.collected lists
  all_objects = ls(envir = .GlobalEnv)
  
  for(obj_name in all_objects){
    object = get(obj_name, envir = .GlobalEnv)
    if(class(object)[1] %in% type_filter){
      
      # models.per.gwas is a list of integers, where each corresponds to the number of model_IDs for a certain gwas
      object <- injectTotalModelsPerGwas(object, models.per.gwas)
      
      # ids.collected is initially a list of empty character vectors, each corresponding to a gwas
      object <- injectIDsCollected(object, ids.collected)
      
      # models.collected is also initially a list of empty character vectors, which will be populated with models
      object <- injectModelsCollected(object, models.collected)
      
      assign(obj_name, object, envir = .GlobalEnv)
    }
  }
}





# Injects Job settings (mem, n.threads, walltime) to the object
injectJobSettings.uni <- function(object, jobSettings){
  
  # ex: str1:"mem = 1500, n.threads = 20" becomes str2: "mem = 1500" and str3: "n.threads = 20"
  jobSettings = unlist(strsplit(jobSettings, split = ', '))
  
  # ex: iterates over str2 and str3 (see above comment)
  for(x in seq_along(jobSettings)){
    x = jobSettings[[x]]
    
    # returns 'run.me' from 'run.me = TRUE'
    variable = str_extract(x, '^[^" "]+')
    
    # return 'TRUE' from 'run.me = TRUE'
    value = str_extract(x, '[^" "]+$')
    
    # if 'variable' doesn't match the objects' slot-names you get error.
    if(variable %!in% slotNames(object)) stop(paste0('Variable \'', variable, '\' does not match, slot names of object ', object@name, '.'))
    
    # inject 'value' to the appropriate slot of object
    if(variable %in% c('run.me')){ slot(object, variable) <- eval(parse(text = value))
    } else slot(object, variable) <- value
  }
  return(object)
}

# At the beginning of the loop in the main section all the objects are injected with the gwas and modelID
injectGwasModelID.uni <- function(thisgwas, thismodelID, type_filter = list('antag', 'postAntag', 'recipeClass', 'collector')){
  # injects gwas and model ID in object
  inject.gwas.model  = function(object, thisgwas, thismodelID){
    slot(object, 'thisgwas')      <- thisgwas # this is 'i' in gv code
    slot(object, 'thismodelID')   <- thismodelID
    return(object)
  }
  
  
  all.objects = ls(envir = globalenv())
  
  # iterates over all objects, if class(object) %in% type_filter; gwas/modelID are injected.
  for(obj_name in all.objects){
    obj = get(obj_name, envir = globalenv())
    if(class(obj)[1] %in% type_filter){
      obj <-  inject.gwas.model(obj,thisgwas, thismodelID)
      assign(obj_name, obj, envir = globalenv())
    }
  }
}


create.dirs.uni = function(object, results.dir){
  
  results.subdir <- object@results.subdir
  thisgwas <- object@thisgwas
  thismodelID <- object@thismodelID
  # create directories
  if(class(object) == 'antag'){
    MultiWAS::gv_dir.create(ma_paste0(file.path(results.dir, results.subdir,'/intermediate.scripts/S01A/', thisgwas)))
    MultiWAS::gv_dir.create(ma_paste0(file.path(results.dir, results.subdir,'/intermediate.scripts/S01B/', thisgwas)))
    MultiWAS::gv_dir.create(ma_paste0(file.path(results.dir, results.subdir, "/scripts/S01A", '/', thisgwas, '/', thismodelID)))
    MultiWAS::gv_dir.create(ma_paste0(file.path(results.dir, results.subdir, "/scripts/S01B", '/', thisgwas, '/', thismodelID)))
    MultiWAS::gv_dir.create(ma_paste0(file.path(results.dir, results.subdir, "/logs/S01A", '/', thisgwas, '/', thismodelID)))
    MultiWAS::gv_dir.create(ma_paste0(file.path(results.dir, results.subdir, "/logs/S01B", '/', thisgwas, '/', thismodelID)))
    MultiWAS::gv_dir.create(ma_paste0(file.path(results.dir, results.subdir, "/jobs/S01A/", thisgwas, '/', thismodelID)))
    MultiWAS::gv_dir.create(ma_paste0(file.path(results.dir, results.subdir, "/jobs/S01B/", thisgwas, '/', thismodelID)))
  } else if(class(object) == 'postAntag'){
    MultiWAS::gv_dir.create(ma_paste0(file.path(results.dir, results.subdir, "/scripts/S01A", '/', thisgwas)))
    MultiWAS::gv_dir.create(ma_paste0(file.path(results.dir, results.subdir, "/scripts/S01B", '/', thisgwas)))
    MultiWAS::gv_dir.create(ma_paste0(file.path(results.dir, results.subdir, "/logs/S01A", '/', thisgwas)))
    MultiWAS::gv_dir.create(ma_paste0(file.path(results.dir, results.subdir, "/logs/S01B", '/', thisgwas)))
    MultiWAS::gv_dir.create(ma_paste0(file.path(results.dir, results.subdir, "/jobs/S01A/", thisgwas)))
    MultiWAS::gv_dir.create(ma_paste0(file.path(results.dir, results.subdir, "/jobs/S01B/", thisgwas)))
  } else if(class(object) == 'collector'){
    MultiWAS::gv_dir.create(ma_paste0(file.path(results.dir, results.subdir, "/scripts/S01A")))
    MultiWAS::gv_dir.create(ma_paste0(file.path(results.dir, results.subdir, "/scripts/S01B")))
    MultiWAS::gv_dir.create(ma_paste0(file.path(results.dir, results.subdir, "/logs/S01A")))
    MultiWAS::gv_dir.create(ma_paste0(file.path(results.dir, results.subdir, "/logs/S01B")))
    MultiWAS::gv_dir.create(ma_paste0(file.path(results.dir, results.subdir, "/jobs/S01A/")))
    MultiWAS::gv_dir.create(ma_paste0(file.path(results.dir, results.subdir, "/jobs/S01B/")))
  } else stop('Input class in create.dirs.uni() is neither antag nor postAntag nor collector')
}



