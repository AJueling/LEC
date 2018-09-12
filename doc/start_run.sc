#!/bin/csh -f

# This script starts the run 'run_henk_mixedbc_extravars_viebahn'
# It is a run on 1280 cores. 
# compiler flags are given in: /home/klipdccp/models/pop/code/build/bull.gnu 

# Make sure that:
# 1.The following files are in this directory:
#   - start_run.sc
#   - pop.ll_smt
#   - pop_in
#   - POP_DomainSizeMod.F90 
#   - domain_size.F90 
#   - GNUmakefile
#   - tavg_contents
#   - movie_contents
#   - transport_file_141lines

# For other options:
# 1. change $case below
# 2. change POP_DomainSizeMod.F90 with which you can change the way the domain is divided over the cores
# 3. change @node, @tasks_per_node, @job_name, @wallclock_limit in pop.ll_smt
# 4. change nprocs_clinic and nprocs_tropic in pop_in

# After that this script automatically builds the pop executable and submits the pop.ll_smt job. 
# author: Michael Kliphuis

set case = run_henk_mixedbc_extravars_viebahn

# check if output directories are there. If not make 'm.
set outputdir_base = /projects/0/samoc/pop/tx0.1/output
#if !(-e $outputdir_base/$case) then
#  mkdir -p $outputdir_base/$case
#  mkdir -p $outputdir_base/${case}/movie
#  mkdir -p $outputdir_base/${case}/restart
#  mkdir -p $outputdir_base/${case}/tavg
#endif

# load needed modules
source /hpc/sw/modules/module/init/csh
source ~/.modules.bull

set workdir = /home/klipdccp/models/pop/scripts/samoc/$case

# In case of new build rm object,mod and f90 files:
rm -r $workdir/compile/*

# determine logfile for build
if !($?LID) setenv LID "`date +%y%m%d-%H%M%S`"
set logfile = build.log.pop.$case.$LID

echo ''
echo 'now building pop executable for case '$case
gmake POPEXEDIR=$workdir POPDIR=/home/klipdccp/models/pop/code POPARCH=bull >>& $logfile 

echo 'now submitting job for case '$case

# Submit job to batch queue
#sbatch pop.slurm




