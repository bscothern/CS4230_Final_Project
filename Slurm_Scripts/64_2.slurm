#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --nodes=1

#SBATCH --account=soc-kp
#SBATCH --partition=soc-kp

#SBATCH --mail-type=END,FAIL

#SBATCH -o 64_2.out
#SBATCH -e 64_2.err

cd $HOME/FinalProject/Test_Cases

cat 64x64_2.in | xargs ../factor_seq -t

