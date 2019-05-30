#!/bin/sh
#ACCOUNT
#PBS -A ONRDC29722307

#QUEUE
##PBS -q background
##PBS -q debug
#PBS -q standard

#WALL TIME
##PBS -l walltime=00:30:00
#PBS -l walltime=60:00:00

#NUMBER OF PROCESSORS
#PBS -l select=10:ncpus=48:mpiprocs=48+1:ncpus=32:mpiprocs=32

#STANDARD ERROR
#PBS -e unstratified.err

#STANDARD OUT
#PBS -o unstratified.out

#JOB NAME
#PBS -N unstratified

#RUN
cd ${WORKDIR}/ekman_bl/rundir_unst/


LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${HOME}/lib/diablo_lib/fftw/2.1.5.3/lib/
export LD_LIBRARY_PATH

mpiexec_mpt -n 480 ./unstratified : -n 32 ./unstratified  >log

