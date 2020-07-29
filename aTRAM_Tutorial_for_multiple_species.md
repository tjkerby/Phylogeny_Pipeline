# aTRAM Tutorial for Multiple Species

## Getting Set Up
First thing to do is to get a project directory up and running. Once that is completed the next step is to set up some directories inside of that. You can do so by entering the project directory and running:

```bash
mkdir Jobs atram_scripts Reads BUSCO_Output Reference
```

Also make sure that you have all of the necessary programs installed beforehand. You will need:
* BUSCO
* Blast
* Bowtie2
* BWA
* Exonerate
* Python
* 1 of the 4 (abyss, trinity, velvet, or spades)
  * In this tutorial I use velvet, but any of the four will work

Last but not least, you need to gitclone aTRAM while in the project directory by running:

```bash
git clone https://github.com/juliema/aTRAM.git
```

One other thing to note is that the scripts that I run below where created for the paired end reads and reference genome that I chose and had available. It is possible, and even likely, that if yours are bigger or smaller that you will need more or less time and computational resources. So think of the scripts as a general idea for what you need and adjust as needed. My reference genome was 25 Megabyte .fasta file and the largerst paired-end read was 11 Gigabyte compressed .fastq file.

# Choose a Reference Genome
The reason you are using aTRAM in the first place is because you lack a properly assembled genome but do have low coverage reads. Although we don't have a reference genome there are probably high quality reference genomes already published that are closely related to your taxa of interest that would suffice for our purposes. So in effect, we use a related reference genome as a guide for our paired end reads to map onto. The more closely related our chosen reference genome is to our species of interest the better.

Once you have selected your reference genome you will need to have your reads for the samples of interest that you want to create genome assemblies for. We will be mapping these on to one another using the reference genome as a guide. Once those are obtained move the reads into the Reads directory and the reference into the Reference directory. The next step is to use BUSCO to compare the reads with the reference to find matching genes that can be used to create COGs. 

Below is the script for running BUSCO on a HPC cluster that uses SLURM. Save it in the Jobs directory and then run it for both proteins and genome uisng `sbatch` followed by whatever you decided to name the script. 

```bash
#!/bin/bash

#SBATCH --time=72:00:00   # walltime
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --mem-per-cpu=64G   # memory per CPU core

module purge
module load conda-pws
module load python/3.7
module load busco

assembly=/path_to_project/Reference/reference_genome.fasta
output=BUSCO_protein
species=insert_species_name
m=proteins
l=/path_to_dataset/embryophyta_odb10

run_BUSCO.py -m ${m} -i ${assembly} -o ${output}_${species} -l ${l} -c 1 -f
```

# FIXME

## Does this output to the correct directory?

First wait for the script above to finish. Later on for the aTRAM stitcher part you will need to create a taxa list. This can be created while changing to the single_copy_busco_sequences directory that was created when we ran BUSCO and contains all of the .faa and .fna files. Once that is your working directory run the script below. It essentially takes the names from each of the amino acid sequences (.faa files) by printing all of the .faa files and stripping the .faa at the end and then puts them in a text file named taxa.txt.

```bash
ls *.faa | sed 's/.faa//' > taxa.txt
```

The next step is to concatenate all of the .fna and .faa files.

```bash
cat *.fna > all.fna
cat *.faa > all.faa
```

Next we begin using aTRAM. The first step is the preprocessing step.

# aTRAM Preprocessing

The preprocessing step takes paired-end reads as input and from them creates BLAST database shards as well as a sqlite database for rapid read retrieval. The bigger your reads the more shards you will have. The script below will run atram_preprocessor.py on all of the paired-end reads in the Reads directory.  

```bash
#!/bin/sh

mkdir -p /lustre/scratch/usr/tkerby2/pws672/Oat_Phylogeny/atram_scripts/preprocessor
for read in /path_to_project/Reads/*_1P.fq.gz
do
file=`basename ${read}`
name=`echo ${file} | sed 's/_1P.fq.gz//'`
cat > /path_to_project/atram_scripts/preprocessor/${name}.sh <<EOF
#!/bin/bash
#SBATCH --time=72:00:00   # walltime
#SBATCH --ntasks=4   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=16G   # memory per CPU core
#SBATCH -J "atram_preprocessing"   # job name
#SBATCH --qos=pws

module purge
module load blast/2.7
mkdir -p path_to_project/atram_db/${name}
path_to_project/aTRAM/atram_preprocessor.py \
  --blast-db=/path_to_project/atram_db/${name}/${name} \
  --end-1=/path_to_project/Reads/${name}_1P.fq.gz \
  --end-1=/path_to_project/Reads/${name}_2P.fq.gz \
  --max-processes=4 \
  --gzip \
  --fastq
EOF

sbatch /path_to_project/atram_scripts/preprocessor/${name}.sh

done
```

# aTRAM.py

The next step is the assembler part. It uses the databases built in the preprocessor step and also uses the COGs generated by the BUSCO output to assemble the different loci. Here is the code for that script.

```bash
#!/bin/sh

for read in /path_to_project/Reads/*_1P.fq.gz
do
file=`basename ${read}`
name=`echo ${file} | sed 's/_1P.fq.gz//'`
mkdir -p /path_to_project/atram_output/${name}/
cat > /path_to_project/atram_scripts/atram/${name}.sh <<EOF
#!/bin/bash
#SBATCH --time=168:00:00   # walltime
#SBATCH --ntasks=4   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=32G  # memory per CPU core
#SBATCH -J "atram"   # job name
#SBATCH --qos=pws

module purge
module load blast/2.7
module load velvet/1.2
module load bwa/0.7.17
module load bowtie2/2.3
module load exonerate/2.2
/path_to_project/aTRAM/atram.py \
  --blast-db=/path_to_project/atram_db/${name}/${name} \
  --query-split=/path_to_project/BUSCO_Output/run_BUSCO_genome_oat_avena_atlantica/single_copy_busco_sequences/all.fna \
  --assembler=velvet \
  --output-prefix=/path_to_project/atram_output/${name}/${name} \
  --iterations=3 \
  --max-processes=4
EOF

sbatch /path_to_project/atram_scripts/atram/${name}.sh

done
```

Some things to note from the script is:
* --iterations default is 5. I changed it to 3 since it was taking too long and the decrease in quality was acceptable for my uses.
* --assembler is an argument that lets you choose which program you would like to use for assembling. I chose velvet, but there are three others you can use as well. Just make sure that you have the assembler downloaded and loaded into your environment.
* --query-split is the COGs which is just the concatenated .fna files that were created earlier.

# aTRAM Stitcher

The next step is to do the aTRAM stitcher. Before we run the script we first need to move all of the filtered files from the output since they are redundant and will bring up an error if left unattended. The script below creates a folder and moves all of the filtered ones inside for each species being assembled.

```bash
#!/bin/sh

for directory in /path_to_project/atram_output/*
do
cd ${directory}
mkdir filtered
mv *.filtered* filtered
cd ../
done
```

Now that we have that taken care of we can run the atram_stitcher.py script. It takes the output assemblies generated by atram.py and the reference amino acid targets and stitches them together using an iterative process.

```bash
#!/bin/sh

mkdir -p /path_to_project/atram_scripts/atram_stitcher
for read in path_to_project/Reads/*_1P.fq.gz
do
file=`basename ${read}`
name=`echo ${file} | sed 's/_1P.fq.gz//'`
mkdir -p /path_to_project/atram_stitcher_output/${name}/
cat > /path_to_project/atram_scripts/atram_stitcher/${name}.sh <<EOF
#!/bin/bash
#SBATCH --time=336:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=32G  # memory per CPU core
#SBATCH -J "atram_stitcher"   # job name
#SBATCH --qos=pws

module purge
module load blast/2.7
module load velvet/1.2
module load bwa/0.7.17
module load bowtie2/2.3
module load exonerate/2.2

/path_to_project/aTRAM/atram_stitcher.py \
--assemblies-dir=/path_to_project/atram_output/${name} \
--reference-genes=/path_to_project/BUSCO_Output/run_BUSCO_genome_oat_avena_atlantica/single_copy_busco_sequences/all.faa \
--taxa=/path_to_project/BUSCO_Output/run_BUSCO_genome_oat_avena_atlantica/single_copy_busco_sequences/taxa.txt \
--output-prefix=/path_to_project/atram_stitcher_output/${name}/${name}

EOF

sbatch /path_to_project/atram_scripts/atram_stitcher/${name}.sh

done
```

# You're Finished!!!

Now you should have some genome assemblies for each of the taxa that you sampled. I will be using the output from this to create a phylogeny in the next tutorial.

 
