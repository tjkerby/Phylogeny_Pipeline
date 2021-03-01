#!/bin/sh

file_path=/lustre/scratch/usr/tkerby2/pws672/Berlandieri_Phylogeny
mkdir -p ${file_path}/atram_scripts/atram
for read in ../Reads/last_5/*_1P.fq.gz
do
file=`basename ${read}`
name=`echo ${file} | sed 's/_1P.fq.gz//'`
mkdir -p ${file_path}/atram_output/${name}/
cat > ${file_path}/atram_scripts/atram/${name}.sh <<EOF
#!/bin/bash
#SBATCH --time=168:00:00   # walltime
#SBATCH --ntasks=4   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=32G  # memory per CPU core
#SBATCH -J "atram"   # job name
#SBATCH --mail-user=tjkerby@gmail.com   # email address
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --qos=pws
module purge
module load blast/2.7
module load velvet/1.2
module load bwa/0.7.17
module load bowtie2/2.3
module load exonerate/2.2
${file_path}/aTRAM/atram.py \
  --blast-db=${file_path}/atram_db/${name}/${name} \
  --query-split=${file_path}/BUSCO/run_ordered_berlandieri_arrow2x_pilon_PGA_assembly_renamed_long_newbusco_arabidopsis/single_copy_busco_sequences/all.fna \
  --assembler=velvet \
  --output-prefix=${file_path}/atram_output/${name}/${name} \
  --iterations=3 \
  --max-processes=4
EOF

sbatch ${file_path}/atram_scripts/atram/${name}.sh

done
