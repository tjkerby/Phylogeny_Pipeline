#!/bin/bash

#SBATCH --time=168:00:00   # walltime
#SBATCH --ntasks=32   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=4G   # memory per CPU core
#SBATCH -J "Iqtree" #job name
#SBATCH --mail-user=tjkerby@gmail.com #email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --qos=pws

source ~/.bashrc
conda activate iqtree

iqtree -s /lustre/scratch/usr/tkerby2/pws672/Berlandieri_Phylogeny/6_concatenated_tree/nt/FcC_smatrix_nt.fas -spp /lustre/scratch/usr/tkerby2/pws672/Berlandieri_Phylogeny/6_concatenated_tree/nt/part_def_nt.txt -nt $SLURM_NPROCS -safe -pre conda_test_nt -m TESTMERGEONLY -mset GTR+G
iqtree -s /lustre/scratch/usr/tkerby2/pws672/Berlandieri_Phylogeny/6_concatenated_tree/aa/FcC_smatrix_aa.fas -spp /lustre/scratch/usr/tkerby2/pws672/Berlandieri_Phylogeny/6_concatenated_tree/aa/part_def_aa.txt -nt $SLURM_NPROCS -safe -pre conda_test_aa -m TESTMERGEONLY -mset LG+G
