# Load MultiWAS

libs <- .libPaths()
libs[3] <- "/sc/arion/projects/roussp01a/sanan/Rlibs/231221_R_4.2.0_MultiWAS"
.libPaths(libs)

# Load MultiWAS
suppressMessages(library(MultiWAS))


#######################
# Parse the recipe file
# DEBUG



#####################################
#####################################
#               STARTUP

path_to_recipes = '/sc/arion/projects/va-biobank/PROJECTS/marios_temp/parallel_antagonism/recipes'

#recipe.file <- ma_paste0(file.path(path_to_recipes, 'bash_test/S4_test_run_V3.csv'))
#recipe.file <- ma_paste0(file.path(path_to_recipes, '/psychAD/V1.csv'))
recipe.file <- ma_paste0(file.path(path_to_recipes, 'GPNMB/GPNMB.csv'))
#recipe.file <- ma_paste0(file.path(path_to_recipes, '/test_synthetic/synthetic_5.csv'))
#recipe.file <- ma_paste0(file.path(path_to_recipes, '/psychAD/V2.csv'))
#recipe.file <- ma_paste0(file.path(path_to_recipes, '/GPNMB/GPNMB_max.csv'))

####################################
# TO-DO: in class Antag make it iterate over the individuals;
# each iteration will run accross all the present steps (up untill ensureSuccess_V3)
# at the end of each iteration you will writeLines() the job_ID of the last ensureSuccess_V3 to an external file
# then when the iteration of individuals end;
# you will submit a 'gatherer lsf array' for each of the job_IDs (or in other words each cell will correspond to an individual)
# you will use the $LSB_JOBINDEX to extract each of the job_IDs from the reference file
# using the job_ID, each job's status will be monitored in a while loop (this could be in bash)
# then the postAntag class will also iterate over participants and 
# their jobs will also be submitted as part of an lsf array that will be waiting cell for cell of the parent array.
# if the parent array is of class antag then it will be the 'gatherer lsf array'.
# Important note: create an ensureSuccess subclass that inherits from postAntag;
# it will have a slot 'local_env' that will state inside which class was the ensureSuccess object created (ie antag or postAntag)
# if it runs inside antag then it will submit independent jobs for each of the lsf-arrays (which in antag correspond to the compound chunks-array for each individual)
# if it run inside postAntag it will callNextMethod() and submitJob() (it will act like a normal postAntag object that expects an array as parentJob)


# unload recipe variables from a recipe-dataframe to a recipe-list.
recipe <- parse_recipe(recipe.file)

# create the working directory (input from recipe)
MultiWAS::gv_dir.create(recipe$working.directory)
setwd(recipe$working.directory)



############################################
############################################
#              CREATE OBJECTS

# Contains disease signatures info
TWASList <- list(path = recipe$df,
                 columnNames =         recipe$columnNames,
                 grep.sig.pattern =    recipe$grep.sig.pattern,
                 model.banlist =       recipe$model.banlist.grep
)

# Contains drug singatures info
SigList <- list(signature.dir = recipe$signature.dir,
                grep.sig.pattern = recipe$grep.sig.pattern,
                prototyping = recipe$prototyping,
                overwrite.intermediate = recipe$overwrite.intermediate
)

# Contains supplementary info (submissions etc). Accompanies most Job-objects (eg. fiverankJob) inside functions.
objRecipe <- new('recipeClass',
                 recipe.file =         recipe.file,
                 account.name =        recipe$account.name,          #'recipe$parameter'
                 queue =               recipe$queue,                 #'recipe$parameter'
                 working.directory =   recipe$working.directory,  #'recipe$parameter'
                 results.dir =         recipe$results.dir,        #'recipe$parameter'
                 dryrun =              recipe$dryrun,             # recipe$dryrun,
                 compoundsList =       list(),                     # handled internally
                 signature.dir =       recipe$signature.dir,
                 path.to.gene.anno.file = recipe$gene.anno.file
)

# create objects environment to create dynamic links between objects
objects.env <- new.env()

# Submits jobs for 5-rank method
fiveRankJob <- new('antag',
                   name = 'fiveRankJob',
                   results.subdir =    '/intermediate.files/five.rank',
                   recipe
)

# Submits jobs for AvgRank calculation (child of fiveRankJob)
avgRankJob <-  new('postAntag',
                   name = 'avgRankJob',
                   results.subdir =    '/intermediate.files/avgRank',
                   recipe
)

# Submits jobs for aggregate_and_prioritize.core.R (child of avgRankJob and VAEJob)
wilcoxRankJob <-  new('postAntag',
                      name = 'wilcoxRankJob',
                      results.subdir =    '/wilcoxRank',
                      recipe
)


robustRankAggJob <- new('collector',
                        name = 'robustRankAggJob',
                        results.subdir = '/robustRankAgg',
                        recipe
                        )

############################################
############################################
#              RUN PIPELINE

# df contains appropriate/filtered TWAS structure
df = prepareTWAS.main(TWASList)

# Filters, Returns and Exports the paths of the .RDS signature files
slot(objRecipe, 'compoundsList') <- ReturnExportSignatures.main(SigList, objRecipe)

gwas.model.combo = unique(df[,c('gwas', 'model_ID')])

injectCollectors.uni(gwas.model.combo)

#
for(i in seq(nrow(gwas.model.combo))){
  
  thisgwas = gwas.model.combo[i,]$gwas
  thismodelID = gwas.model.combo[i,]$model_ID
  
  # exports df.shape. It is used as reference by 'core functions'
  create_df.shape.main(objRecipe, df, thisgwas, thismodelID)
  
  thisgwas = MultiWAS::make_java_safe(thisgwas)
  thismodelID = MultiWAS::make_java_safe(thismodelID)
  
  # passes thisgwas, thismodelID inside classes of type_filter
  injectGwasModelID.uni(thisgwas, thismodelID, type_filter = list('antag', 'postAntag', 'recipeClass'))
  
  
  # fiveRankJob is class antag
  if(fiveRankJob@run.me) submitJob.class(fiveRankJob, objRecipe)
  
  # avgRankJob is class postAntag
  if(avgRankJob@run.me) submitJob.class(avgRankJob, objRecipe)
  
  # wilcoxRankJob is class postAntag
  if(wilcoxRankJob@run.me) submitJob.class(wilcoxRankJob, objRecipe)
  
  # robustRankAggJob is class collector
  if(robustRankAggJob@run.me) collect_SubmitJob.class(robustRankAggJob, objRecipe)
}

