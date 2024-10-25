#!/bin/sh

#$ -N JOBNAME
#$ -o JOBNAME.out 
#$ -e JOBNAME.err
#$ -cwd			    
#$ -l data
#$ -V 			    

# Guix
export PATH="/fast/AG_Landthaler/scripts/GRAND-SLAM_2.0.7b:$PATH"
export PATH="$HOME/.guix-conflicts/bin:$HOME/.guix-profile/bin:$PATH"
export GUIX_PROFILE=$HOME/.guix-profile
source $GUIX_PROFILE/etc/profile

# Main Job
snakemake -s /fast/AG_Landthaler/Pipelines/GRAND-SLAM/v1/Pipeline.py --configfile config.yaml --cluster "qsub {params.qsub} -V -l data -o $PWD/snakemake_output/ -e $PWD/snakemake_output/" --jobs 32 
