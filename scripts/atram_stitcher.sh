#!/bin/sh

file_path=/lustre/scratch/usr/tkerby2/pws672/Berlandieri_Phylogeny
mkdir -p ${file_path}/atram_scripts/atram_stitcher
for read in ../Reads/last_5/*_1P.fq.gz
do
file=`basename ${read}`
name=`echo ${file} | sed 's/_1P.fq.gz//'`
mkdir -p ${file_path}/atram_stitcher_output/${name}/
cat > ${file_path}/atram_scripts/atram_stitcher/${name}.sh <<EOF
#!/bin/bash
#SBATCH --time=168:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=64G  # memory per CPU core
#SBATCH -J "atram_stitcher"   # job name
#SBATCH --mail-user=tjkerby@gmail.com   # email address
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --qos=pws
module purge
source ~/.bashrc
conda activate aTRAM

${file_path}/aTRAM/atram_stitcher.py \
--assemblies-dir=${file_path}/atram_output/${name} \
--reference-genes=${file_path}/BUSCO/run_ordered_berlandieri_arrow2x_pilon_PGA_assembly_renamed_long_newbusco_arabidopsis/single_copy_busco_sequences/all.faa \
--taxa=${file_path}/BUSCO/run_ordered_berlandieri_arrow2x_pilon_PGA_assembly_renamed_long_newbusco_arabidopsis/single_copy_busco_sequences/taxa.txt \
--output-prefix=${file_path}/atram_stitcher_output/${name}/${name}

EOF

sbatch ${file_path}/atram_scripts/atram_stitcher/${name}.sh

done
