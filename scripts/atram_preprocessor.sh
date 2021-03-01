#!/bin/sh

file_path=/lustre/scratch/usr/tkerby2/pws672/Berlandieri_Phylogeny
mkdir -p ${file_path}/atram_scripts/preprocessor
for read in ${file_path}/Reads/last_5/*_1P.fq.gz
do
file=`basename ${read}`
name=`echo ${file} | sed 's/_1P.fq.gz//'`
cat > ${file_path}/atram_scripts/preprocessor/${name}.sh <<EOF
#!/bin/bash
#SBATCH --time=72:00:00   # walltime
#SBATCH --ntasks=4   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=16G   # memory per CPU core
#SBATCH -J "atram_preprocessing"   # job name
#SBATCH --mail-user=tjkerby@gmail.com   # email address
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
module purge
module load blast/2.7 
mkdir -p ${file_path}/atram_db/${name}
${file_path}/aTRAM/atram_preprocessor.py \
  --blast-db=${file_path}/atram_db/${name}/${name} \
  --end-1=${file_path}/Reads/last_5/${name}_1P.fq.gz \
  --end-1=${file_path}/Reads/last_5/${name}_2P.fq.gz \
  --max-processes=4 \
  --gzip \
  --fastq
EOF

sbatch ${file_path}/atram_scripts/preprocessor/${name}.sh

done
