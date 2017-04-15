#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --nodes=1

#SBATCH --account=soc-kp
#SBATCH --partition=soc-kp

#SBATCH --mail-type=END,FAIL

#SBATCH -o 144_1.out
#SBATCH -e 144_1.err

cd $HOME/FinalProject

cat 144x144_1.in | xargs ../factor_seq -t

