#!/bin/bash

#SBATCH --job-name=TestJob
#SBATCH --output=TestJob.out
#SBATCH --time=1-00:10:00
#SBATCH -p compute1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

eval "$(conda shell.bash hook)"

conda activate Space

python Python_test.py
