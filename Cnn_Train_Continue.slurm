#!/bin/bash

#SBATCH --job-name=spatialTreat_example_cnn_train
#SBATCH --output=cnn_train_continue.out
#SBATCH --time=1-05:00:00
#SBATCH -p gpu1v100
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1


eval "$(conda shell.bash hook)"

conda activate Space

python 02_cnn_train_saved_model.py



