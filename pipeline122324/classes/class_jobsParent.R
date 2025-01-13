setClass('jobsParent',
         slots = list(
           # The following are defined when object is created
           name =            'character',
           results.subdir =  'character',  #'/dirname',
           core.function =   'character',  #'/path/to/core/function.R',
           code =            'character',  #'R|python',
           
           # The following are defined by the recipe
           run.me =              'logical',    #'recipe$parameter'
           walltime =            'character',    #'recipe$parameter'
           n.threads =           'character',    #'recipe$parameter'
           mem =                 'character',    #'recipe$parameter'
           
           # The following parameters are handled inside the function
           thisgwas =            'character',  #'handled_internally'
           thismodelID =         'character',  #'handled_internally'
           submission.command =  'character',   #'handled_internally'
           job_id =            'character'
         ),
         prototype = list(
           name =            'object name',
           # The following are defined when object is created
           results.subdir =  '/results.subdir',  #'/dirname',
           core.function =   '/path/to/core/function.R',  #'/path/to/core/function.R',
           code =            'NULL',  #'R|python',
           
           # The following are defined by the recipe
           run.me =              FALSE,
           walltime =            '10:10',    #'recipe$parameter'
           n.threads =           '10',    #'recipe$parameter'
           mem =                 '1500',
           
           # The following parameters are handled inside the function
           thisgwas =                'thisgwas',  #'handled_internally'
           thismodelID =             'thismodelID',  #'handled_internally'
           submission.command =      'bsub ...',   #'handled_internally'
           job_id =                '"done(12345)"'   #'handled_internally'
         ))
