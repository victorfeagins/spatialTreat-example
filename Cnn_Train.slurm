#!/bin/bash

#SBATCH --job-name=spatialTreat_example_cnn_train
#SBATCH --output=spatialTreat_example_cnn_train.out
#SBATCH --time=1-05:00:00
#SBATCH -p compute1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1


eval "$(conda shell.bash hook)"

conda activate Space

python 02_cnn_train.py



