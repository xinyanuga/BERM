#PBS -S /bin/bash
#PBS -q batch
#PBS -N testjob
#PBS -l nodes=1:ppn=4
#PBS -l walltime=20:00:00
#PBS -l mem=128gb
#PBS -M jloconn@uga.edu 
#PBS -m abe

## Grab the job id from an environment variable and create a directory for the output
export CUSTOM_JOBID=`echo "$PBS_JOBID" | cut -d"[" -f1`
mkdir ${PBS_O_WORKDIR}/${CUSTOM_JOBID}

cd ${PBS_O_WORKDIR}/${CUSTOM_JOBID}

ml load R/3.5.0-foss-2018a-X11-20180131-GACRC  

R CMD BATCH --slave --quiet --no-save /work/demlab/jloconn/seagrant/analysis/process_all_gce_pixels.r
