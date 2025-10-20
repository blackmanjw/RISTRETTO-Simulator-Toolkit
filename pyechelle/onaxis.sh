#!/bin/sh 
#SBATCH --job-name ristretto 
#SBATCH --error logs/log.e%j 
#SBATCH --output logs/log.o%j  
#SBATCH --gres=gpu:rtx3090:1
#SBATCH --mem-per-gpu=8G
#SBATCH --partition=gpu
#SBATCH --qos=job_gpu
#SBATCH --time 0-00:25:00 

module load GCCcore/12.3.0 CUDA/12.1.1 Python/3.11.3-GCCcore-12.3.0
srun ~/pyechelle/bin/python onaxis.py

