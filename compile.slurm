#!/bin/bash
#SBATCH --job-name=compilation
#SBATCH --cluster=gpu
#SBATCH --partition=a100
#SBATCH --nodes=1      
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --time=00:30:00          
#SBATCH --error=output/job.%J.err
#SBATCH --output=output/job.%J.out
#SBATCH --mail-user=hjl28@pitt.edu
#SBATCH --mail-type=END,FAIL

module purge
export MODULEPATH=/software/rhel9/spack-2026/modules/Core
module load openmpi/4.1.8-cuda
module load hdf5/1.14.6

cd /ihome/eschneider/hjl28/code/cholla
make clean
make -j 4 CUDA_ARCH=sm_80
