#PBS -S /bin/bash
#PBS -q batch
#PBS -N testjob
#PBS -l nodes=1:ppn=1
#PBS -l walltime=480:00:00
#PBS -l mem=64gb
#PBS -M jloconn@uga.edu 
#PBS -m abe

cd $PBS_O_WORKDIR

ml R/3.5.0-foss-2018a-X11-20180131-GACRC  

R -slave -q /work/demlab/jloconn/seagrant/analysis/process_all_gce_pixels.r
