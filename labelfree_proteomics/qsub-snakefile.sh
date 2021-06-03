#PBS -N ms-snakemake
#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=28
#PBS -j oe
#PBS -o snake_$PBS_JOBID.out
#PBS -e snake_$PBS_JOBID.err
#PBS -m abe

echo ------------------------------------------------------
echo -n 'Job is running on node '; cat $PBS_NODEFILE
echo ------------------------------------------------------
echo PBS: qsub is running on $PBS_O_HOST
echo PBS: originating queue is $PBS_O_QUEUE
echo PBS: executing queue is $PBS_QUEUE
echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: execution mode is $PBS_ENVIRONMENT
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: node file is $PBS_NODEFILE
echo PBS: current home directory is $PBS_O_HOME
echo PBS: PATH = $PBS_O_PATH
echo ------------------------------------------------------

#export SINGULARITYENV_OPENMS_DATA_PATH=/usr/local/share/OpenMS

module load python
module load singularity
source activate snakemake

cd $PBS_O_WORKDIR

snakemake --snakefile Snakefile.py --unlock
snakemake --snakefile Snakefile.py --notemp --cores 28 --resources mem_mb=128000 --use-singularity 
