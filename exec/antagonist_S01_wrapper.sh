#!/bin/bash
#BSUB -J Ant_S01_mothership
#BSUB -P acc_va-biobank
#BSUB -q premium
#BSUB -n 20
#BSUB -R span[hosts=1]
#BSUB -R rusage[mem=3000]
#BSUB -W 1440
#BSUB -oo logs/S01mothership.out
#BSUB -eo logs/S01mothership.err
#BSUB -L /bin/bash
# https://medium.com/@wujido20/handling-flags-in-bash-scripts-4b06b4d0ed04
# https://www.baeldung.com/linux/use-command-line-arguments-in-bash-script

# This is just a helper script to load R and then execute Rscript while passing all the arguments

library=/sc/arion/projects/roussp01a/sanan/Rlibs/230919_R_4.2.0_MultiWAS_Antagonist
ml R
Rscript $library/antagonist/exec/antagonist_S01_wrapper.R "$*"

exit 0
