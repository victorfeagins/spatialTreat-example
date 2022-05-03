#!/bin/bash

#SBATCH --job-name=TestJob
#SBATCH --output=TestJob.out
#SBATCH --time=1-00:10:00
#SBATCH -p compute1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

conda init Space

python Python_test.py
