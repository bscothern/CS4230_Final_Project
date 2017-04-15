#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --nodes=1

#SBATCH --account=soc-kp
#SBATCH --partition=soc-kp

#SBATCH --mail-type=END,FAIL

#SBATCH -o 96_1.out
#SBATCH -e 96_1.err

cd $HOME/FinalProject

srun cat 96x96_1.in | xargs ../factor_seq -t

