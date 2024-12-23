setClass(
  "recipeClass",
  slots = list(
    
    # The following are defined by the recipe
    account.name =        'character',  #'recipe$parameter'
    queue =               'character',  #'recipe$parameter'
    working.directory =   'character',  #'recipe$parameter'
    results.dir =         'character',  #'recipe$parameter'
    signature.dir =       'character',  #'recipe$parameter'
    recipe.file =         'character',
    recipe =              'list',
    dryrun =              'logical',
    compoundsList =        'ANY',
    thisgwas =            'character',
    thismodelID =         'character',
    path.to.gene.anno.file = 'character'
    
    
  ),
  prototype = list(
    # The following are defined by the recipe
    account.name =        'acc_va-biobank',  #'recipe$parameter'
    queue =               'express',  #'recipe$parameter'
    working.directory =   'recipe$working.directory',  #'recipe$parameter'
    results.dir =         'recipe$results.dir',  #'recipe$parameter'
    signature.dir =       'recipe$signature.dir',  #'recipe$parameter'
    recipe.file =         'opt$recipe',
    recipe =              list(),
    dryrun =              FALSE,
    compoundsList =       list(),
    thisgwas =            'AD',
    thismodelID =         'Microglia',
    path.to.gene.anno.file = NULL
  ),validity = function(object){
    
    errors = character()
    
    if(!file.exists(object@working.directory)) stop('From recipeClass validity(): For some reason the working directory has not been created.')
    if(length(list.files(object@signature.dir)) == 0) stop('From recipeClass validity(): The provided signature directory looks empty')
    if(length(list.files(object@signature.dir)) == 0) stop('From recipeClass validity(): The provided signature directory looks empty')
    
    path.to.gene.anno.file = object@path.to.gene.anno.file
    
   
      if(!file.exists(object@path.to.gene.anno.file)){ 
        message('From recipeClass: The provided gene.anno.file path: ', object@path.to.gene.anno.file, ' does NOT exist as an ABSOLUT path. 
                I will look in the system.file("extdata", package = "MultiWAS")')
        path.to.gene.anno.file = ma_paste0(file.path(system.file('extdata', package = 'MultiWAS'), object@path.to.gene.anno.file))
        if(!file.exists(path.to.gene.anno.file)){
          stop('From recipeClass validity(): You provided the following gene.anno.file variable
                                          in the recipeClass, ', object@path.to.gene.anno.file, ' which does not exist, after looking inside
                                          system.file("extdata", package = "MultiWAS") I still could not find it, please make sure that the
                                          following path ', path.to.gene.anno.file, ' , is either absolut path or relative path to MultiWAS extdata: ')
        } else message('Path ', path.to.gene.anno.file, ' found.')
      }
      
    
    
    
    
  }
)




