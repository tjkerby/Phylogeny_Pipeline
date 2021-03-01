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

iqtree -s /lustre/scratch/usr/tkerby2/pws672/Berlandieri_Phylogeny/6_concatenated_tree/nt/FcC_smatrix_nt.fas -spp /lustre/scratch/usr/tkerby2/pws672/Berlandieri_Phylogeny/6_concatenated_tree/nt/part_def_nt.txt -nt $SLURM_NPROCS -safe -pre berlandieri_95_nt -m TESTMERGE -cptime 5000 -bb 1000 -mset GTR+G

