setClass('collector',
         slots = list(
           total.models.per.gwas = 'list',
           ids.collected = 'list',
           models.collected = 'list'
         ),
         prototype = list(
           total.models.per.gwas = NULL,
           ids.collected = NULL,
           models.collected = list()
         ),
         contains = 'postAntag')



########## the first if should be part of function injectCollectors

setMethod('initialize',
          'collector',
          function(.Object, name, results.subdir, recipe){
            .Object <- callNextMethod()
            .Object
          })

setMethod('collect_submitJob.class',
          'collector',
          function(object, recipe.obj){
            #if(!parent@run.me){
            #  check.file.existence
              
            #}
            
            parentJob <- objects.env[[object@parentJob]]
            parent.name = parentJob@name
            parent.job_id = parentJob@job_id
            parent.run.me = parentJob@run.me
            parent.subdir = parentJob@results.subdir
            
            
            thisgwas = object@thisgwas
            thismodelID = object@thismodelID
            
            # total.models.per.gwas is a list of integers, where each corresponds to the number of model_IDs for a certain gwas
            total.models.per.gwas <- object@total.models.per.gwas[[thisgwas]]
            
            ######### COLLECT MODEL IDS AND SET THEM AS THE CONDITION.. COLLECT IDS ONLY IF PARENT IS TRUE
            
            object <- injectModelsCollected(object, thismodelID)
            
            # 
            if(parent.run.me){
              parent.job_id <- gsub('\"', '', parent.job_id)
              
              theseIDs = object@ids.collected[[thisgwas]]
              
              if(length(theseIDs) == 0){ theseIDs <- paste0('\"', parent.job_id)
              }else theseIDs <-  paste0(theseIDs, ' && ', parent.job_id)
              
            }
            
            # when the last model is collected the job gets submitted :: total.models.per.gwas is an integer
            if(length(object@models.collected[[thisgwas]]) == object@total.models.per.gwas[[thisgwas]]){
              
              object <- injectTheseIDs(object, paste0(theseIDs, '\"'))
              
              # the "ended(...) && ended(...) && ..." is the new job.id which your submitted command will wait for..
              # FOR THIS TO WORK YOU NEED TO HAVE "DONE(...)" AVAILABLE FOR CLASS POSTANTAG TOO !!!!
              if(parent.run.me) updateParentJobID(object, parentJob)
              
              # class collector inherits from class postAntag
              submitJob.class(object, objRecipe)
              
            } else{
              object <- injectTheseIDs(object, theseIDs)
              assign(object@name, object, envir = .GlobalEnv)
            }
            
          })


injectTheseIDs <- function(object, theseIDs){
  thisgwas <- object@thisgwas
  object@ids.collected[[thisgwas]] <- theseIDs
  return(object)
}


injectModelsCollected <- function(object, thismodelID){
  thisgwas <- object@thisgwas
  object@models.collected[[thisgwas]] <- c(object@models.collected[[thisgwas]], thismodelID)
  return(object)
}



injectJobIDs <- function(object, thisgwas, ids.collected){
  object@ids.collected[[thisgwas]] <- ids.collected
  assign(object@name, object, envir = .GlobalEnv)
}


updateParentJobID <- function(object, parentJob){
  thisgwas <- object@thisgwas
  parentJob@job_id <- object@ids.collected[[thisgwas]]
  assign(parentJob@name, parentJob, envir = objects.env)
}




