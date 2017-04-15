#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --nodes=1

#SBATCH --account=soc-kp
#SBATCH --partition=soc-kp

#SBATCH --mail-type=END,FAIL

#SBATCH -o 128_2.out
#SBATCH -e 128_2.err

cd $HOME/FinalProject

srun cat 128x128_2.in | xargs ../factor_seq -t

