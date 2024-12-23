setClass(
  "antag",
  slots = list(
    compounds.chunk.path= 'character',  
    chunks.number =       'numeric'
  ),
  prototype = list(
    compounds.chunk.path =  '/path/to/compounds',
    chunks.number =         0
  ),
  validity = function(object){
    
    errors <- character()
    
    if(object@code %!in% c('R', 'py')) errors = c(errors, ' From antag validity: The assigned core.function to ', name, ' has extension ',
                                             paste0('.', object@code), ' :: Present version only recognizes extension ".R" and ".py"')
    
    if(!file.exists(object@core.function)) errors = c(errors, ' From antag validity, object ', object@name, ' the provided path to its core script ',
                                                object@core.function, ' does not exist. In present version of the package,
                                                the core functions are sourced from system.file("extdata", package = "MultiWAS").
                                                If your provided path looks like e.g. "/external.code/cdr/core.function.R" the initialize() function
                                                will look for it in system.file("extdata", package = "MultiWAS"), so make sure that it can actually be found there.')
    
    if(length(errors) == 0) TRUE else errors
  },
  contains = 'jobsParent'
)


setMethod(
  'initialize',
  'antag',
  function(.Object, name, results.subdir, recipe, core.function = NULL){
    

    if(name %!in% names(recipe)) stop('From initialize: For an object of class "antag" you provided the name ', name,
    ' which doesnt match any of the names in the "variable" column of the given recipe.')
    
    # string-name of the object is injected in the objected
    .Object@name = name
    
    # Injects n.threads, mem and walltime in the object
    .Object <- injectJobSettings.uni(.Object, recipe[[name]])
    
    
    
    if(is.null(core.function)){
      # Injects the /path/to/core.function.R in the object, sourcing it from the package
      core.path = recipe[[paste0('core.', name)]]
      
      if(!file.exists(core.path)){
        message('From initialize class.antag: the provided core.function path: ', core.path ,' is NOT an absolut path. I will check inside system.file("extdata", package = "MultiWAS"')
        core.path = ma_paste0(file.path(system.file('extdata', package = 'MultiWAS'), core.path))
        if(!file.exists(core.path)){
          stop('From initialize class.antag, after looking inside system.file("extdata", package = "MultiWAS")
                                       I could NOT find file ', core.path, ' . Make sure core.function is either absolut path or relative to system.file("extdata", package = "MultiWAS")')
        } else message('core.function path was found in system.file("extdata", package = "MultiWAS")')
      }
    } else object@core.function = core.function
    
    
    .Object@core.function = core.path
    
    # the slot 'code' is determined. ex: my_core_script.R will return 'R'
    .Object@code = str_extract(.Object@core.function, '[^.]+$')
    
    # all the dirs and files related to the object are stored in 'results.subdir'
    .Object@results.subdir = results.subdir
    
    # Attach object to the 'objects.env'. This environment is used as a dynamic link between objects
    objects.env[[name]] <- .Object
    
    # Return the object
    .Object
  }
)


setMethod('submitJob.class',
          'antag',
          function(object, recipe.obj){
            tryCatch(
              {
                # ex: object is fiveRankJob
                
                # mylist contains the full paths to each of compounds chunks (/path/to/trt_cp_1.300.RDS, /path/to/trt_cp_300.600.RDS , ...)
                mylist <- recipe.obj@compoundsList
                
                # creates all necessary directories (for results, logs, jobs (contain bsub commands), scripts (contain the 'core.commands'))
                create.dirs.uni(object, recipe.obj@results.dir)
                
                # chunks are the '.RDS batches of compounds' (trt_cp_1.300.RDS is a chunk). Their total number determines the size of the array.
                slot(object, 'chunks.number') <- length(mylist)
                
                for(ind in seq_along(mylist)){
                  
                  # inject compound information
                  # ex: compounds.chunk.path = /path/to/trt_cp.1.300.RDS
                  slot(object, 'compounds.chunk.path')  <- stats::setNames(mylist, basename(mylist))[ind] # this is 'i' in gv code
                  
                  # creates a /path/to/scripts/script_ind.sh script, which calls an R script that is specific for a .RDS compounds file(ex: trt_cp_1.300.RDS) : Each script will be run by run.script.sh (see below)
                  create.core.command(object, recipe.obj, ind)
                }
                
                # creates a bsub-string ex: 'bsub -J "myArray[1-10]" -P acc_va-biobank -q premium ... -L /bin/bash '
                b_sub <- create_bsub(object, recipe.obj)
                
                # creates a file which runs each of the scripts created by 'create.core.command()' using the $LSB_JOBINDEX
                # file is stored at the /intermediate.scripts/run.script.sh (actual name matches model_ID)
                path_to_run.scripts.sh <- create_run.scripts.sh(object, recipe.obj)
                
                # This is the job that gets submitted to the shell.
                # SUMMARY: 'b_sub' calls 'run.scripts.sh' which calls each of the /scripts/script_ind.sh (created by create.core.command).
                # Each script_ind.sh executes 'RScript core.function.R' for a specific trt_cp...RDS (e.g. 5-rank.core.R for trt_cp.1.300.RDS)
                submission.command <- paste0(b_sub, path_to_run.scripts.sh)

                # for debugging: job is exported
                writeLines(
                  submission.command,
                  ma_paste0(file.path(recipe.obj@results.dir, object@results.subdir,'/jobs/S01A/',
                                      object@thisgwas,'/',  paste0(object@thismodelID, '.txt'))))
                
                if(!recipe.obj@dryrun){
                  
                  # ex: submits the job and returns 'Job <12345> submitted to queue premium'
                  bash.output = system(submission.command, intern = TRUE) 
                  
                  # extract '12345' from above example
                  array_id <- str_extract(bash.output, '[0-9]+')
                  message('The array_id of ', object@name, ' is ', array_id)
                } else array_id <- '12345'
                
                # the pipeline accepts the final form (e.g. "done(12345)"). This is for a good reason (see class collector for more)
                array_id <- paste0('"ended(', array_id, ')"')
                
                # the job_id is passed to the object so that ensureSuccess/'subsequent job' will wait for it.
                slot(object, 'job_id') <- array_id
                
                # assign object to objects.env so the 'ensureSuccess' object can extract its job ID
                objects.env[[object@name]] <- object
                
                # This was supposed to be a way to 'save' failed jobs, but lsf probably handles this well with the '-retry 3 1' option.
                # I leave it here for now..
                if(TRUE){
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
                stop(paste0('For object ', object@name, ' in "setMethod(\'submitJob.class\', \'antag\',...)" you got the following message:\n',
                            e$message, '\n'))
              }
            )
          })



setMethod('create_bsub',
          'antag',
          function(object, recipe.obj){
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
                chunks.number =       object@chunks.number
                n.threads =           object@n.threads
                walltime =            object@walltime
                mem.per.core =        object@mem
                
                # ex: -J AD_excitatory_neurons[1-10] ; submits an 'lsf array' job with 10 cells in shell
                # in each cell of the 'lsf array' a script that loads the core.function loaded with a specific 'chunk of compounds' (e.g. trt_cp_1.300.RDS)
                name.of.array = paste0('"', thisgwas, '_', thismodelID, '[1-', chunks.number, ']"')
                name.of.array = gsub('\\.', '', name.of.array)
                
                b.sub <- paste0('bsub -J ', name.of.array,
                                ' -P ', account.name,
                                ' -q ', queue,
                                ' -n ', n.threads,
                                ' -W ', walltime,
                                ' -R "span[hosts=1]"',
                                ' -R "affinity[core(1)]"',
                                ' -R \"rusage[mem=', mem.per.core, ']\"',
                                ' -oo ', ma_paste0(file.path(results.dir, results.subdir,'/logs/S01A/', thisgwas, '/', thismodelID, '/', paste0('script_', '%I', '.out'))),
                                ' -eo ', ma_paste0(file.path(results.dir, results.subdir,'/logs/S01A/', thisgwas, '/', thismodelID, '/', paste0('script_', '%I', '.err'))),
                                ' -L /bin/bash ')
                return(b.sub)
              },
              error = function(e){
                message(paste0('In "setMethod(\'create_bsub\', \'antag\',...)" :\n',
                                  e$message, '\n'))
              }
            )
          })


setMethod('create.core.command',
          'antag',
          function(
    object,
    recipe.obj,
    ind
          ){
            tryCatch(
              {
                # from recipe
                working.directory =     recipe.obj@working.directory
                results.dir =           recipe.obj@results.dir
                recipe.file =           recipe.obj@recipe.file
                
                # from object
                thisgwas =              object@thisgwas
                thismodelID =           object@thismodelID
                results.subdir =        object@results.subdir
                core.function =         object@core.function
                code =                  object@code
                
                # /path/to/trt_cp_1.300.RDS
                compounds.chunk.path =  object@compounds.chunk.path
                
                # This is the command that executes the core.function
                thiscommand <- paste0(
                  "cd ", working.directory, "\n",
                  ifelse(code == 'R', paste0("ml R", "\n"), paste0("ml python", "\n")),
                  ifelse(code == 'R', paste0("Rscript --verbose ", core.function, " "), paste0("python --verbose ", core.function, " ")),
                  "--thisgwas ", thisgwas, " ",
                  "--thismodelID ", thismodelID, " ",
                  "--recipe ", recipe.file, " ",
                  "--cmapfile ", compounds.chunk.path, " ",
                  "--results.subdir ", results.subdir
                )
                
                # 'thiscommand' is exported to its respective .sh file (to be executed by run.scripts.sh)
                writeLines(
                  paste0(
                    "#!/bin/bash", "\n",
                    thiscommand),
                  ma_paste0(file.path(results.dir, results.subdir, '/scripts/S01A/', thisgwas, thismodelID, paste0('/script_', ind, '.sh'))) # /path/to/script_1.sh will be executed when LSB_JOBINDEX = 1
                )
              },
              error = function(e){
                return(list(
                  message = paste0('In \'setMethod(\'create.core.command\', \'antag\',...)\' : ',
                                   e$message)))
              }
            )
          })



create_run.scripts.sh <- function(object, recipe.obj){
  
  thismodelID <- object@thismodelID
  thisgwas <- object@thisgwas
  results.dir <- recipe.obj@results.dir
  results.subdir <- object@results.subdir 
  working.directory <- recipe.obj@working.directory
  
  # ex: /path/to/scripts/gwas/model_ID , contains the script_1.sh , script_2.sh etc (created by 'create.core.command()')
  # each of those scripts execute the 'core.function' for a specific .RDS batch of compounds
  path_to_core.scripts = ma_paste0(file.path(results.dir, results.subdir,'/scripts/S01A/', thisgwas,'/',  thismodelID))
  
  # 'run.scripts.sh' changes directory to the 'path_to_core.scripts' (see above line) and submits the containing scripts according to the $LSB_JOBINDEX
  run.scripts.sh <- paste0('#!/bin/bash\n',
                            'cd ', ma_paste0(file.path(working.directory, path_to_core.scripts)), '\n',
                            'scripts=($(ls))\n',
                            'thisscript=${scripts[$LSB_JOBINDEX-1]}\n',
                            'sh $thisscript')
  
  path_to_run.scripts.sh <- ma_paste0(file.path(getwd(), results.dir, results.subdir,'/intermediate.scripts/S01A/', thisgwas,'/',  paste0(thismodelID, '.sh')))
  
  # export the command
  writeLines(run.scripts.sh, path_to_run.scripts.sh)
  
  # make the command executable
  system(paste0('chmod +x ', path_to_run.scripts.sh))
  
  return(path_to_run.scripts.sh)
}
