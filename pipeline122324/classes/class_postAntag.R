setClass('postAntag',
         slots = list(
           parentJob =       'ANY'
         ),
         prototype = list(
           parentJob =       'ANY'
         ), validity = function(object){
           
           
           `%!in%` <- function(x,y) !(x%in%y)
           
           errors <- character()
           
           
           if(object@code %!in% c('R', 'py')) errors = c(errors, paste0('\nFrom postAntag validity(): The assigned core.function to ', name, ' has extension ',
                                                    paste0('.', object@code), ' :: Present version only recognizes extension ".R" and ".py"'))
           
           if(object@run.me && object@name == 'avgRankJob' && object@parentJob != 'fiveRankJob') errors = c(errors, paste0('\nFrom postAntag validity(): For object avgRankJob 
                                                                                       present version supports only fiveRankJob to be its parent.'))
           
           if(object@run.me && object@name == 'wilcoxRankJob' && object@parentJob %!in% c('avgRankJob', 'VAEJob')) errors = c(errors, paste0('\nFrom postAntag validity():
                                                                                                         For object wilcoxRankJob, present version
                                                                                                         supports only avgRankJob and VAEJob to be
                                                                                                         its parent.'))
           
           if(!file.exists(object@core.function)) errors = c(errors, paste0('\nFrom postAntag validity, object ', object@name, ' the provided path to its core script ',
                                                             object@core.function, ' does not exist. In present version of the package,
                                                the core functions are sourced from system.file("extdata", package = "MultiWAS").
                                                If your provided path looks like e.g. "/external.code/cdr/core.function.R" the initialize() function
                                                will look for it in system.file("extdata", package = "MultiWAS"), so make sure that it can actually be found there.'))
           
           tryCatch({
           parentJob = objects.env[[object@parentJob]]
           parent.subdir = parentJob@results.subdir
           parent.run.me = parentJob@run.me
           
           if(!file.exists(parent.subdir) && !parent.run.me) errors = c(errors, paste0('\nFrom postAntag validity, object ',
                                                                       object@name, ' has parent ', object@parentJob,
                                                                       ' . I notice that the parent.subdir "', parent.subdir,
                                                                       '" does NOT exist, while at the same time the parent.run.me
                                                                       is set to FALSE.'))
           },
           error = function(e){
             message(paste0('Warning from postAntag validity: I tried checking for the existence of the parent.subdir, however an unexpected error occured,
                     here is what the error message looked like:\n', e$message))
           })
           
           if(length(errors) == 0) TRUE else errors
           
         },
         contains = 'jobsParent'
)


setMethod(
  'initialize',
  'postAntag',
  function(.Object, name, results.subdir, recipe, parentJob = NULL, core.function = NULL, thisgwas = NULL, thismodelID = NULL){
    
    # for objects ensureSuccess, the input parameters are hardcoded and NOT from recipe
    if(grepl('ensureSuccess', name)){
      .Object@parentJob = parentJob
      .Object@name = name
      .Object@results.subdir = results.subdir
      .Object@core.function = core.function
      .Object@code = str_extract(.Object@core.function, '[^.]+$')
      
      # in objects that are not ensureSuccess these are injected in initialize through injectJobSettings.uni
      .Object@run.me = TRUE
      .Object@walltime = '11:00'
      .Object@mem = '500'
      .Object@n.threads = '1'
      
      # in objects that are not ensureSuccess these are injected at the beginning of the for loop
      .Object@thisgwas = thisgwas
      .Object@thismodelID = thismodelID
      objects.env[[name]] <- .Object
      return(.Object)
    }
    
    # returns a string with the name of the parent job from recipe ex: names(recipe) == 'parent.avgRankJob' then returns 'fiveRankJob' 
    find_parentJob = function() return(recipe[grepl(tolower(name), tolower(names(recipe))) & grepl('parent', tolower(names(recipe)))][[1]])
    
    # If user doesn't define parentJob then it will be 'found' based on names in recipe (parentJob is defined in the ensureSuccess context) 
    # spaces removed for convenience
    if(is.null(parentJob)) str_parentJob = gsub(' ', '', find_parentJob()) else str_parentJob = parentJob
    
    # if more than one parents are found, the function stops
    temp = unlist(strsplit(str_parentJob[[1]], split = ','))
    if(length(temp) > 1) stop('In object ', name, ' only a single parentJob is allowed.')
    
    .Object@parentJob = str_parentJob
    
    # string-name of the object is injected in the object
    .Object@name = name
    
    .Object@results.subdir = results.subdir
    
    if(is.null(core.function)){
      # Injects the /path/to/core.function.R in the object, sourcing it from the package
      core.path = recipe[[paste0('core.', name)]]
      if(!file.exists(core.path)){
        message('From initialize class.postAntag: the provided core.function path: ', core.path ,' is NOT an absolut path. I will check inside system.file("extdata", package = "MultiWAS"')
        core.path = ma_paste0(file.path(system.file('extdata', package = 'MultiWAS'), core.path))
        if(!file.exists(core.path)){
          stop('From initialize class.postAntag, after looking inside system.file("extdata", package = "MultiWAS")
                                       I could NOT find file ', core.path, ' . 
             Make sure core.function is either absolut path or relative to system.file("extdata", package = "MultiWAS")')
        } else message('core.function path was found in system.file("extdata", package = "MultiWAS")')
      }
    } else core.path = core.function
    
    
    
    .Object@core.function = core.path
    
    # Injects n.threads, mem and walltime in the object
    .Object <- injectJobSettings.uni(.Object, recipe[[name]])
    
    # the slot 'code' is determined. ex: my_core_script.R will return 'R'
    .Object@code = str_extract(.Object@core.function, '[^.]+$')
    
    # Attach object to the 'objects.env'. This environment is used as a dynamic link between objects
    objects.env[[name]] <- .Object
    .Object
  }
)

setMethod('submitJob.class',
          'postAntag',
          function(object, recipe.obj){
            # creates all necessary directories
            create.dirs.uni(object, recipe.obj@results.dir)
            tryCatch(
              {
                
                # Load the current job's "parent job" (which the current job waits before being submitted)
                parentJob <- objects.env[[object@parentJob]]
                parent.name = parentJob@name
                parent.job_id = parentJob@job_id
                parent.run.me = parentJob@run.me
                parent.subdir = parentJob@results.subdir
                
                
                if(!parent.run.me){
                  if(!file.exists(ma_paste0(file.path(recipe.obj@results.dir, parent.subdir)))) stop('custom error message: From submitJob.class for object ',
                                                                                                object@name,
                                                                                                ' the supposed "parent.subdir", ',
                                                                                                parent.subdir,
                                                                                                ' does NOT exist. Make sure you have set ', 
                                                                                                parent.name, ' to TRUE.')
                }
                
                # creates a string ex: 'bsub -w "job_id" -P acc_va-biobank -q premium ... -L /bin/bash '
                submission.command <- create_bsub(object, recipe.obj, parent.name, parent.run.me, parent.job_id)
                
                # creates a /path/to/scripts/gwas/model_id.sh script, which calls an R/python script e.g. avg_rank_core.R
                create.core.command(object, recipe.obj, parent.subdir)
                
                if(!recipe.obj@dryrun){
                  # ex: submits the job and returns 'Job <12345> submitted to queue premium'
                  bash.output <- system(submission.command, intern = TRUE)
                  
                  # extracts '12345' and attach it to object
                  job_id <- str_extract(bash.output, '[0-9]+')
                } else job_id <- '12345'
                
                
                # prepare job_id
                job_id <- paste0('"ended(', job_id, ')"')
                
                slot(object, 'job_id') <- job_id
                message('The job_id of ', object@name, ' is ', job_id)
                
                # we change the job priority, so that the chain is complete
                # Sys.sleep(10)
                #tryCatch({
                #temp = system(paste0('bmod -sp 100 ', job_ids[i]), intern = TRUE)
                # message('I tried changing the priority for job ID ', job_ids[i], ' and system returned: ',temp)
                #}, error = function(e){
                #  message('I caught this error while trying to change the permissions \n', e)
                #})
                
                # assign object to objects.env so the 'subsequent Job object' can extract its job ID
                objects.env[[object@name]] <- object
                
                # the condition will not be satisfied when the object itself is an ensureSuccess object, because then it would lead to a non-stop loop
                if(!grepl('ensureSuccess', object@name)){
                          # create a vector with all the jobs that are used in the context of ensureSuccess
                          String_jobObjects = c(object@name, paste0('ensureSuccess_V', 1:3))
                          
                          # multiple ensureSuccess job object are chained one behind the other so that each one guarantees the previous ones success 
                          for(i in seq_len(length(String_jobObjects) - 1)){
                            ensureSuccess <- new('postAntag',
                                                 name = String_jobObjects[i+1],
                                                 results.subdir = ma_paste0('/', object@results.subdir, '/ensureSuccess/', 'V', i),
                                                 recipe,
                                                 parentJob = String_jobObjects[i],
                                                 thisgwas = thisgwas,
                                                 thismodelID = thismodelID,
                                                 core.function = '/sc/arion/projects/va-biobank/PROJECTS/marios_temp/parallel_antagonism/local_antagonist_package/perform_antagonism_V3/setup/confirmSuccess_antag.R')
                            submitJob.class(ensureSuccess, recipe.obj)
                          }
                          
                          # for clarity. a bit complicated but this is accurate. first string_jobObject is 'fiveRankJob' so the second is gonna have V1
                          final_ind = i 
                          
                          # the ID of the last job object is assigned to be the job_ID that the 'subsequent job' will be waiting
                          success_job_id <- objects.env[[paste0('ensureSuccess_V', final_ind)]]@job_id
                          
                          slot(object, 'job_id') = success_job_id
                          
                          # assign object to objects.env so that the 'subsequent job' object can extract its job ID
                          objects.env[[object@name]] <- object
                }
                
              },
              error = function(e){
                message(paste0('In \'setMethod(\'submitJob\', \'postAntag\',...)\' for object: \'', object@name,'\' you got this error message: ',
                                e$message))
              })
          }
)

setMethod('create_bsub',
          'postAntag',
          function(object, recipe.obj, parent.name, wait, parent.job_id){
            tryCatch(
              {
                
                # from recipe
                account.name =        recipe.obj@account.name
                results.dir =         recipe.obj@results.dir
                queue =               recipe.obj@queue
                
                
                # from object
                results.subdir =      object@results.subdir
                thisgwas =            object@thisgwas
                thismodelID =         object@thismodelID
                n.threads =           object@n.threads
                walltime =            object@walltime
                mem.per.core =        object@mem
                
                # If the 'run.me' status of the 'parent Job' is TRUE, the 'current Job' will be on 'waiting status'. Else wise it is started immediately
                if(wait) message(object@name, ' is waiting for job with ID ', parent.job_id) else message('First job to be submitted is for object ', object@name)
                
                
                b.sub <- paste0(ifelse(wait, paste0('bsub -w ', parent.job_id), paste0('bsub -J ', thisgwas, '_', thismodelID)),
                                ' -P ', account.name,
                                ' -q ', queue,
                                ' -n ', n.threads,
                                ' -W ', walltime,
                                ' -R "span[hosts=1]"',
                                ' -R "affinity[core(1)]"',
                                ' -R \"rusage[mem=', mem.per.core, ']\"',
                                ' -oo ', ifelse(class(object) == 'collector',
                                                ma_paste0(file.path(results.dir, results.subdir,'/logs/S01A/', paste0(thisgwas, '.out'))),
                                                ma_paste0(file.path(results.dir, results.subdir,'/logs/S01A/', thisgwas, '/', paste0(thismodelID, '.out')))
                                                ),
                                ' -eo ', ifelse(class(object) == 'collector',
                                                ma_paste0(file.path(results.dir, results.subdir,'/logs/S01A/', paste0(thisgwas, '.err'))),
                                                ma_paste0(file.path(results.dir, results.subdir,'/logs/S01A/', thisgwas, '/', paste0(thismodelID, '.err')))
                                                ),
                                ' -L /bin/bash', ' < ')
                
                # This is the job that gets submitted. 'submission.command' will be submitted to bash. submission.command calls '< /path/to/core.function.sh' (created in create.core.command).
                submission.command <- paste0(
                  b.sub,  ifelse(class(object) == 'collector',
                                 ma_paste0(file.path(results.dir, results.subdir,'/scripts/S01A', paste0(thisgwas, ".sh"))),
                                 ma_paste0(file.path(results.dir, results.subdir,'/scripts/S01A/', thisgwas,'/', paste0(thismodelID, ".sh"))))# This file executes the core function.
                )
                
                # the bsub is also exported, because in case a job fails, you want to resubmit it.
                writeLines(
                  submission.command,
                  ifelse(class(object) == 'collector',
                         ma_paste0(file.path(results.dir, results.subdir,'/jobs/S01A', paste0(thisgwas, ".sh"))),
                         ma_paste0(file.path(results.dir, results.subdir,'/jobs/S01A', thisgwas,  paste0(thismodelID, ".sh")))
                         )
                )
                
                return(submission.command)
              },
              error = function(e){
                message(paste0('In \'setMethod(\'create_bsub\', \'postAntag\',...)\' ::: for object \'',object@name,'\' you got this error message:\n',e$message))
                
              }
            )
          }
)

setMethod('create.core.command',
          'postAntag',
          function(
            object,
            recipe.obj,
            parent.subdir
          ){
            tryCatch(
              {
                
                # from recipe
                working.directory =   recipe.obj@working.directory
                results.dir =         recipe.obj@results.dir
                recipe.file =         recipe.obj@recipe.file
                
                # from object
                thisgwas =            MultiWAS::make_java_safe(object@thisgwas)
                thismodelID =         MultiWAS::make_java_safe(object@thismodelID)
                results.subdir =      object@results.subdir
                core.function =       object@core.function
                code =                object@code
                
                # This gets executed by the 'current job' (submission.command).
                thiscommand <- paste0(
                  "cd ", working.directory, "\n",
                  
                  if(code == 'R'){ paste0("ml R", "\n")
                  }else if(code == 'py'){ paste0("ml python", "\n")
                  }else stop('From postAntag create.core.command: I cant handle the extension ', code),
                  
                  if(code == 'R'){ paste0("Rscript --verbose ", core.function, " ") 
                  }else if(code == 'py'){ paste0("python --verbose ", core.function, " ") 
                  }else stop('From postAntag create.core.command: I still cant handle extesion ', code, ' :: Please fix both "ml R" and "Rscript core.function.R"'),
                  
                  "--recipe.file ", recipe.file, " ",
                  "--thisgwas ", thisgwas, " ",
                  "--thismodelID ", thismodelID, " ",
                  "--results.dir ", results.dir, " ",
                  "--results.subdir ", results.subdir, " ",
                  "--parent.subdir ", parent.subdir, " "
                )
                
                # thiscommand is exported to its respective .sh file (to be executed by submission.command)
                writeLines(
                  paste0(
                    "#!/bin/bash", "\n",
                    thiscommand),
                  ifelse(class(object) == 'collector',
                         ma_paste0(file.path(results.dir, results.subdir, '/scripts/S01A/', paste0(thisgwas, ".sh"))),
                         ma_paste0(file.path(results.dir, results.subdir, '/scripts/S01A/', thisgwas, paste0(thismodelID, ".sh")))
                         )
                )
              },
              error = function(e){
                return(list(
                  message =paste0('In \'setMethod(\'create.core.command\', \'postAntag\',...)\' :\n',
                                  e$message)))
              }
            )
          }
)
