#!/bin/bash
#PBS -l nodes=1:ppn=1,mem=16gb,walltime=23:59:00

source activate scbase
cd $PBS_O_WORKDIR

scbase run_mcmc -m ${ASE_MODEL} ${TGX_MODEL} \
                --hapcode ${MAT_HAPCODE} ${PAT_HAPCODE} \
                -o ${OUTFILE} \
                ${LOOMFILE}

source deactivate
