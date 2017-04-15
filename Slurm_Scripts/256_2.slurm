#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --nodes=1

#SBATCH --account=soc-kp
#SBATCH --partition=soc-kp

#SBATCH --mail-type=END,FAIL

#SBATCH -o 256_2.out
#SBATCH -e 256_2.err

cd $HOME/FinalProject

cat 256x256_2.in | xargs ../factor_seq -t

