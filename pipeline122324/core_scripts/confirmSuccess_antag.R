# wrapper_main


# Script
# parse args..
libs <- .libPaths()
libs[3] <- "/sc/arion/projects/roussp01a/sanan/Rlibs/231221_R_4.2.0_MultiWAS"
libs[4] <- "/sc/arion/projects/va-biobank/software/Georgios_dev/240702_R_4.2.0_MultiWAS_Antagonist/"
.libPaths(libs)
# Load MultiWAS
library(MultiWAS)
library(antagonist)

option_list = list(
  make_option(c("-r", "--recipe.file"),
              type="character",
              default=NULL,
              help = "Point to recipe file for this project",
              metavar = "character"),
  make_option(c("-g", "--thisgwas"),
              type="character",
              default=NULL,
              help = "Point to recipe file for this project",
              metavar = "character"),
  make_option(c("-m", "--thismodelID"),
              type="character",
              default=NULL,
              help = "Use a number to only run a subset of the signatures, e.g. 2 works well",
              metavar = "character"),
  make_option(c("-d", "--results.dir"),
              type="character",
              default=NULL,
              help = "Use a number to only run a subset of the signatures, e.g. 2 works well",
              metavar = "character"),
  make_option(c("-s", "--results.subdir"),
              type="character",
              default=NULL,
              help = "Use a number to only run a subset of the signatures, e.g. 2 works well",
              metavar = "character"),
  make_option(c("-p", "--parent.subdir"),
              type="character",
              default=NULL,
              help = "Use a number to only run a subset of the signatures, e.g. 2 works well",
              metavar = "character"),
  make_option(c("-c", "--cmapfile"),
              type="character",
              default=NULL,
              help = "Use a number to only run a subset of the signatures, e.g. 2 works well",
              metavar = "character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)



# DEBUG
#opt = list()
#opt$recipe.file = '/sc/arion/projects/va-biobank/PROJECTS/marios_temp/parallel_antagonism/recipes/test_run_recipes/bash_test/S4_test_run_V4.csv'
#opt$thisgwas = 'AD'
#opt$thismodelID = 'Microglia_Genes_.._Microglia_FANS_.._Genes_.._PrediXcan_.._EUR'
#opt$results.dir = 'results/GTP_CDR'
#opt$results.subdir = '/intermediate.files/avgRank'
#opt$parent.subdir = '/intermediate.files/five.rank'
#bash_output = c('results/GTP_CDR/intermediate.files/five.rank/logs/S01A/AD/Microglia_Genes_.._Microglia_FANS_.._Genes_.._PrediXcan_.._EUR/script_2.out', 'results/GTP_CDR/intermediate.files/five.rank/logs/S01A/AD/Microglia_Genes_.._Microglia_FANS_.._Genes_.._PrediXcan_.._EUR/script_1.out')

# more debug
#opt = list()
#opt$recipe.file = '/sc/arion/projects/va-biobank/PROJECTS/marios_temp/parallel_antagonism/recipes/GPNMB/GPNMB.csv'
#opt$thisgwas = 'AD'
#opt$thismodelID = 'ADAM_GPNMB'
#opt$results.dir = 'results/GTP_CDR'
##opt$results.subdir is not needed
#opt$parent.subdir = '/intermediate.files/five.rank/ensureSuccess/V2'
#bash_output = 'results/GTP_CDR/intermediate.files/five.rank/logs/S01A/AD/ADAM_GPNMB/script_190.out'

# DEBUG: bash jobs from postAntag jobs
#jobs_bsub = 'bsub -w "ended(148976380)" -P acc_va-biobank -q express -n 1 -W 11:00 -R span[hosts=1] -R rusage[mem=500] -oo results/GTP_CDR/intermediate.files/five.rank/ensureSuccess/V2/logs/S01A/AD/ADAM_GPNMB.out -eo results/GTP_CDR/intermediate.files/five.rank/ensureSuccess/V2/logs/S01A/AD/ADAM_GPNMB.err -L /bin/bash < results/GTP_CDR/intermediate.files/five.rank/ensureSuccess/V2/scripts/S01A/AD/ADAM_GPNMB.sh'

recipe.file = opt$recipe.file
thisgwas = opt$thisgwas
thismodelID = opt$thismodelID
results.dir = opt$results.dir
results.subdir = opt$results.subdir
parent.subdir = opt$parent.subdir


working.directory = fread(recipe.file)[variable == 'working.directory']$value

setwd(working.directory)

# END OF NON-SPECIFIC SECTION  #
################################

Sys.sleep(5)
if(basename(parent.subdir) %in% c('five.rank', 'VAE')){
  
  message('My parent subdir is ', parent.subdir, ' and because either "five.rank" or "VAE" was found, I will act accordingly.')
  
  tryCatch({
    logs_to_check = ma_paste0(file.path(results.dir, parent.subdir, '/logs/S01A/', thisgwas, thismodelID))
    if(!file.exists(logs_to_check)) stop()
  }, error = function(e){
    logs_to_check = ma_paste0(file.path(results.dir, parent.subdir, '/logs/', thisgwas, thismodelID))
    if(!file.exists(logs_to_check)) stop(paste0('Both of the following directories do NOT exist;\n',
                                                ma_paste0(file.path(results.dir, parent.subdir, '/logs/S01A/', thisgwas, thismodelID)), '\n',
                                                ma_paste0(file.path(results.dir, parent.subdir, '/logs/', thisgwas, thismodelID))))
  })
  
  # you should have -iL and NOT -il
  check.files = paste0('find ', logs_to_check, ' -type f -name "*.out" -exec grep -iL "successfully completed" {} +')
  
  while(TRUE){
        message('I will check if all the files in ', logs_to_check, ' contain the "successfully completed" string in their *.out files in /logs dir.')
        bash_output = system(check.files, intern = TRUE)
        
        if(length(bash_output) == 0){
          message('All the "*.out" files contain the "successfully completed" string. Good job!')
          break
        }else{
          message('Not all files have completed, I will be resubmitting jobs.')
          message('The following files were not completed.')
          for(i in seq_along(bash_output)){
            message('bash_output ', i, ' ', bash_output[i])
          }
          
          # the failed scripts get resubmitted, their job_ids get redeemed.
          job_ids = character()
          message('I queried system ')
          for(i in seq_along(bash_output)){
                
                out_file_path = bash_output[i]
                
                # sh_script is the path to the 'executed cell/script' (not submitted command) that failed (e.g. /path/to/script_2.sh)
                sh_script = gsub('.out', '.sh', out_file_path)
                sh_script = gsub('logs', 'scripts', sh_script)
                
                # this is the number of the script that failed e.g. 2
                int = str_extract(basename(sh_script), '[0-9]+') 
                
                # this is the txt with path to bsub command for the whole array
                jobs_bsub_path = gsub('logs', 'jobs', out_file_path)
                jobs_bsub_path = gsub(paste0('/', basename(jobs_bsub_path)), '.txt', jobs_bsub_path)
                
                # extract the bsub command and modify it so that it submits the script outside of the lsf Array
                jobs_bsub = readLines(jobs_bsub_path)
                jobs_bsub = gsub('intermediate.scripts', 'scripts', jobs_bsub)
                jobs_bsub = gsub('\\[1-[0-9]*\\]', '', jobs_bsub)
                jobs_bsub = gsub('\\.sh', paste0('/', basename(sh_script)), jobs_bsub)
                jobs_bsub = gsub('%I', as.character(int), jobs_bsub)
                
                
                # make script executable (again)
                exec_script = gsub('.*bin/bash|<', '', jobs_bsub)
                exec_script = gsub(' ', '', exec_script)
                message('Now making the script to be executed, executalbe\n')
                system(paste0('chmod +x ', exec_script))
                
                # resubmit job
                temp = system(jobs_bsub, intern = TRUE)
                
                # this is needed so that the lsf has time to process the submissions
                Sys.sleep(3)
                message('\nThe BELOW job has failed. After resubmitting, bash returned: ', temp,
                        '\n', jobs_bsub, '\n')
                # extract job id
                job_ids[i] = str_extract(temp, '[0-9]+')
                
          }
          
          message('Here is the list of the job_ids:')
          for(i in seq_along(job_ids)){
            message(job_ids[i])
          }
          
          # job_status contains the 'status' of the submitted jobs
          job_status = character()
          Sys.sleep(3)
          for(i in seq_along(job_ids)){
                this_job = job_ids[i]
                
                # obtain the job status
                tryCatch({
                  job_status[i] = system(paste0('bjobs -noheader -o stat ', this_job), intern = TRUE)
                }, error = function(e){
                  message('For job id: ', this_job, ' you received the following error message:\n', e, '\n')
                })
                
                # change job priority status as to avoid idle dependencies
                tryCatch({
                  temp = system(paste0('bmod -sp 100 ', this_job), intern = TRUE)
                  message(temp)
                }, error = function(e){
                  message('I caught this error\n', e)
                })
                
                message('Resubmitted Job with ID ', this_job, ' has status ', job_status[i])
          }
          message('Job statuses before entering the while loop are: ', paste(job_status, collapse = ', '))
          message('The while loop condition is ', all(job_status %!in% c('DONE', 'EXIT')))
          
          # the STATUS of the submitted jobs is redeemed every 20 minutes.
          # the loop breaks when their exit code is either 'DONE' or 'EXIT'
          while(all(job_status %!in% c('DONE', 'EXIT'))){
                message('There are ', sum(job_status %!in% c('DONE', 'EXIT')), ' jobs whose status is NOT "DONE" or "EXIT", I will now sleep and check in 10.')
                
                Sys.sleep(60)
                job_status = character()
                for(i in seq_along(job_ids)){
                  this_job = job_ids[i]
                  
                  # collect job status
                  tryCatch({
                    job_status[i] = system(paste0('bjobs -noheader -o stat ', this_job), intern = TRUE)
                  }, error = function(e){
                    message('From nested while loop\n: For job id: ', this_job, ' you received the following error message:\n', e, '\n')
                  })
                  
                  # change job priority status as to avoid idle dependencies
                  tryCatch({
                    temp = system(paste0('bmod -sp 100 ', job_ids[i]), intern = TRUE)
                    message(temp)
                  }, error = function(e){
                    message('I caught this error\n', e)
                  })
                  
                  
                  message('From nested while loop: The job ID ', this_job, ' has status ', job_status[i])
                }
                gc()
          }
        }
      }
} else{
  
  message('My parent subdir is ', parent.subdir, ' and this means these were NOT "five.rank" or "VAE" JOBS.')
          
  tryCatch({
    logs_to_check = ma_paste0(file.path(results.dir, parent.subdir, '/logs/S01A/', thisgwas))
    if(!file.exists(logs_to_check)) stop()
  }, error = function(e){
    logs_to_check = ma_paste0(file.path(results.dir, parent.subdir, '/logs/', thisgwas))
    if(!file.exists(logs_to_check)) stop(paste0('Both of the following directories do NOT exist;\n',
                                                ma_paste0(file.path(results.dir, parent.subdir, '/logs/S01A/', thisgwas)), '\n',
                                                ma_paste0(file.path(results.dir, parent.subdir, '/logs/', thisgwas))))
  })
  
  # you should have -iL and NOT -il
  check.files = paste0('find ', logs_to_check, ' -type f -name "*.out" -exec grep -iL "successfully completed" {} \\;')
  
  while(TRUE){
    
        bash_output = system(check.files, intern = TRUE)
        
        if(length(bash_output) == 0){
              message('All the *.out files contain the "successfully completed" string. Good job!')
              break
        }else{
              message('Not all files have completed, I will be resubmitting jobs.')
              
              job_ids = character()
              for(i in seq_along(bash_output)){
                        
                        out_file_path = bash_output[i]
                        
                        # sh_script is the path to the 'executed cell/script' (not submitted command) that failed (e.g. /path/to/script_2.sh)
                        sh_script = gsub('.out', '.sh', out_file_path)
                        sh_script = gsub('logs', 'scripts', sh_script)
                        
                        # this is the txt with path to bsub command for the whole array
                        jobs_bsub_path = gsub('logs', 'jobs', out_file_path)
                        jobs_bsub_path = gsub('.out', '.sh', jobs_bsub_path)
                        
                        # extract the bsub command and modify it so that it submits the script outside of the lsf Array
                        jobs_bsub = readLines(jobs_bsub_path)
                        jobs_bsub = gsub('-w "done\\([0-9]+\\)" ', '', jobs_bsub)
                        
                        # make script executable (again)
                        exec_script = gsub('.*bin/bash|<', '', jobs_bsub)
                        exec_script = gsub(' ', '', exec_script)
                        message('Now making the script to be executed, executalbe\n')
                        system(paste0('chmod +x ', exec_script))
                        
                        
                        temp = system(jobs_bsub, intern = TRUE)
                        
                        # this is needed so that the lsf has time to process the submissions
                        Sys.sleep(3)
                        
                        message('\nThe BELOW job has failed. After resubmitting, bash returned: ', temp,
                                '\n', jobs_bsub, '\n')
                        job_ids[i] = str_extract(temp, '[0-9]+')
              }
              
              # job_status contains the 'status' of the submitted jobs
              job_status = character()
              
              Sys.sleep(3)
              for(i in seq_along(job_ids)){
                this_job = job_ids[i]
                
                # obtain the job status
                tryCatch({
                  job_status[i] = system(paste0('bjobs -noheader -o stat ', this_job), intern = TRUE)
                }, error = function(e){
                  message('For job id: ', this_job, ' you received the following error message:\n', e, '\n')
                })
                
                # change job priority status as to avoid idle dependencies
                tryCatch({
                  temp = system(paste0('bmod -sp 100 ', this_job), intern = TRUE)
                  message(temp)
                }, error = function(e){
                  message('I caught this error\n', e)
                })
                
                message('Resubmitted Job with ID ', this_job, ' has status ', job_status[i])
              }
              
              
              message('Job statuses before entering the while loop are: ', paste(job_status, collapse = ', '))
              message('The while loop condition is ', all(job_status %!in% c('DONE', 'EXIT')))
              
              # the STATUS of the submitted jobs is redeemed every 10 minutes.
              # the loop breaks when their exit code is either 'DONE' or 'EXIT'
              while(all(job_status %!in% c('DONE', 'EXIT'))){
                    message('There are ', sum(job_status %!in% c('DONE', 'EXIT')), ' jobs whose status is NOT "DONE" or "EXIT", I will now sleep and check in 10.')
                    
                    Sys.sleep(60)
                    job_status = character()
                    for(i in seq_along(job_ids)){
                          this_job = job_ids[i]
                          tryCatch({
                            job_status[i] = system(paste0('bjobs -noheader -o stat ', this_job), intern = TRUE)
                          }, error = function(e){
                            message('From nested while loop\n: For job id: ', this_job, ' you received the following error message:\n', e, '\n')
                          })
                          
                          # change job priority status as to avoid idle dependencies
                          tryCatch({
                            temp = system(paste0('bmod -sp 100 ', job_ids[i]), intern = TRUE)
                            
                          }, error = function(e){
                            message('I caught this error ', e)
                          })
                          
                          message('From nested while loop: The job ID ', this_job, ' has status ', job_status[i])
                    }
                    
                    gc()
                    
              }
        }  
  }
}



#bash_output = c("results/GTP_CDR/intermediate.files/avgRank/logs/S01A/AD/Microglia_Genes_.._Microglia_FANS_.._Genes_.._PrediXcan_.._EUR/script_2.out","results/GTP_CDR/intermediate.files/avgRank/logs/S01A/AD/Microglia_Genes_.._Microglia_FANS_.._Genes_.._PrediXcan_.._EUR/script_1.out")



#


#message('I received this output:\n', bash_output[1], '\n', bash_output[2],
#        ' the class of the output was ', class(bash_output),
#        ' it has length ', length(bash_output))
