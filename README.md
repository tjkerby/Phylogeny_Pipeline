# Orphan Crops Lab
I joined the Orphan Crops lab my senior year during my undergrad while at BYU after taking two graduate level Biology/Bioinformatics courses. The lab's website is https://pws.byu.edu/orphaned-crops in case you are interested in learning more about their mission and more recent work.

In the lab I did some comparitive genomic analyses between the genomes of C. Berlandieri and several of its close ancestors. I presented the lab's findings as a speaker at the Plant and Animal Genome Conference in 2020. The presentation slides that I used are in the repository.

![Picture of me at PAG 2020](https://github.com/tjkerby/Phylogeny_Pipeline/blob/master/PAG_2020.jpg)

Apart from the comparitive genomic analyses, I also spent a large amount of time developing a reproducible pipeline for taking raw low coverage reads (pacbio long reads with illumina for polishing) and turning them into a phylogeny. Below is the documentation for how to recreate it. Be warned, using BYU's supercomputer, the process below was able to create phylogenies for 95 different taxa in a little over 2 to 4 weeks of computing time. The amount of computational resources requested is shown in the scripts below.

Here is an image of the final product.

![Astral Amino Acid Phylogeny for 95 relatd taxa to C. Berlandieri](https://github.com/tjkerby/Phylogeny_Pipeline/blob/master/astral_95_aa.tre.jpg)

Although the work I did with the pipeline has not been published, (due mainly to me graduating and there not being someone else to pick up) the project gave me a ton of experience with bash scripting, high performance computing, trouble shooting, and communicating complex results to non-technical audiences.

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

One other thing to note is that the scripts that I run below were created for the paired end reads and reference genome that I chose and had available. It is possible, and even likely, that if yours are bigger or smaller that you will need more or less time and computational resources. So think of the scripts as a general idea for what you need and adjust as needed. My reference genome was 25 Megabyte .fasta file and the largerst paired-end read was 11 Gigabyte compressed .fastq file.

Make sure that the names of the files in the Reads directory are formatted how you want them. The name of the file without the .fq.gz will be the prefix that is used in the remainder of steps and will actas a key for connecting between steps. 

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

First wait for the script above to finish. Later on for the aTRAM stitcher part you will need to create a taxa list. This can be created while changing to the single_copy_busco_sequences directory that was created when we ran BUSCO and contains all of the .faa and .fna files. Once that is your working directory run the script below. It essentially takes the names from each of the amino acid sequences (.faa files) by printing all of the .faa files and stripping the .faa at the end and then puts them in a text file named taxa.txt.

```bash
ls *.faa | sed 's/.faa//' > taxa.txt
```

The next step is to concatenate all of the .fna and .faa files.

```bash
cat *.fna > all.fna
cat *.faa > all.faa
```
Make sure that the headers in both files are not too long. They begin with > and for the files I used have the file path to the reference in there. In one of the runs I did the file path was so long that later portions of the pipeline that use the header in naming files wasn't able to store all of the text and threw an error. If necessary you may need to run a command similar to the one below to shorten it.

```bash
cp all.fna copied_all.fna
sed 's/.fslhome.pjm43.fsl_groups.fslg_lifesciences.pjm43.berlandieri.HiC_assembly.ordered_berlandieri_arrow2x_pilon_PGA_assembly_renamed.fasta/berlandieri.fasta/g' copied_all.fna
mv all.fna old_all.fna
sed -i 's/.fslhome.pjm43.fsl_groups.fslg_lifesciences.pjm43.berlandieri.HiC_assembly.ordered_berlandieri_arrow2x_pilon_PGA_assembly_renamed.fasta/berlandieri.fasta/g' copied_all.fna
mv copied_all.fna all.fna
```

The same steps can be used for the all.faa file as well. The first step is to copy it to make sure that you don't mess it up. Then use the first sed command to get the headers in a format you like. Then rename the original file with the mv command. Once you get what you want then use the 2nd sed command with the -i option to replace the file. Then rename that file  

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

The next step is the assembler part. It uses the databases built in the preprocessor step and also uses the COGs generated by the BUSCO output to assemble the different loci. Here is the code for that script. To this moment I'm still not sure why this is the case but when I run it once some work and some don't. However, when it is run a second time all of them work. I'm not sure if it has something to do with the database not being able to accessed, but it works fine after running it twice.

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

Now that we have that taken care of we can run the atram_stitcher.py script. It takes the output assemblies generated by atram.py and the reference amino acid targets and stitches them together using an iterative process. When I run this with the modules it often has an error right before it finishes because apparently the versions of python and sqlite3 are not compatible. To fix this you can make a conda environment and then use it instead of the module system.

FIXME: Show how to make the conda environment aTRAM and what packages are in it so that there aren't problems when running it.
FIXME: update script with the conda environment.

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

# FIXME
I might consider showing an example of the finished gene file as well as other examples of files at various stages of the aTRAM process.

# Phylogeny

After the above script has finished you should have genes assembled for each of the BUSCO genes. In the phylogeny we will create we will use these BUSCO genes as COGs to build a genetic phylogeny.
 
To format the COGs as well as translate the nucleotide files run the following script:

```bash
#!/bin/sh

module purge
module load emboss
mkdir -p /lustre/scratch/usr/tkerby2/pws672/Berlandieri_Phylogeny/output_summarized/nt
mkdir /lustre/scratch/usr/tkerby2/pws672/Berlandieri_Phylogeny/output_summarized/aa
mkdir /lustre/scratch/usr/tkerby2/pws672/Berlandieri_Phylogeny/output_summarized/temp
for directory in /lustre/scratch/usr/tkerby2/pws672/Berlandieri_Phylogeny/atram_stitcher_output/*/ ;
do
        cd ${directory}
        INDEX=1
        name=`basename "${directory}"`
        for f in *.fasta;
        do
                cog_id=`echo $(ls *.fasta | awk '{print($0)}' | sed 's/^[^\.]*\.//g' | sed 's/_.*//' | tail -n+${INDEX} | head -n1)`
                sequence_id=`echo $(ls *.fasta | awk '{print($0)}' | sed 's/.stitch.*//' | sed 's/.*.fasta_//g' | tail -n+${INDEX} | head -n1)`
                species_name=`echo $(ls *.fasta | awk '{print($0)}' | sed 's/\..*//' | tail -n+${INDEX} | head -n1)`
                nt_sequence=`awk '{if(NR>=2) print$0}' ${f}`
                #length=`awk '{if(NR>=2) print$0}' ${f} | wc -m`
                #header=`echo ">"${cog_id}"|"${species_name}"|"${sequence_id}"|1-"${length}"|.|."`
                header=`echo ">"${cog_id}"|"${species_name}"|"${sequence_id}`
                # Handling if a certain species does not have the gene
                if [ -s ${f} ]
                then
                        # Convert from nucleotides to amino acids
                        `transeq ${f} /lustre/scratch/usr/tkerby2/pws672/Berlandieri_Phylogeny/output_summarized/temp/${f} -clean`
                        aa_sequence=`echo "$(awk '{if(NR>=2) printf $0;}' /lustre/scratch/usr/tkerby2/pws672/Berlandieri_Phylogeny/output_summarized/temp/${f})"`
                        aa_header=`echo ">"${cog_id}"|"${species_name}"|"${sequence_id}`
                        # Append the header and the amino acid sequence for each species in the same cog file
                        echo "${aa_header}" >> /lustre/scratch/usr/tkerby2/pws672/Berlandieri_Phylogeny/output_summarized/aa/${cog_id}.aa.fa
                        echo "${aa_sequence}" >> /lustre/scratch/usr/tkerby2/pws672/Berlandieri_Phylogeny/output_summarized/aa/${cog_id}.aa.fa
                        # Append the header and the nucleotide sequence for each species in the same cog file
                        echo "${header}" >> /lustre/scratch/usr/tkerby2/pws672/Berlandieri_Phylogeny/output_summarized/nt/${cog_id}.nt.fa
                        echo "${nt_sequence}" >> /lustre/scratch/usr/tkerby2/pws672/Berlandieri_Phylogeny/output_summarized/nt/${cog_id}.nt.fa
                else
                        #aa_header=`echo ">"${cog_id}"|"${species_name}"|"${sequence_id}"|0-0|.|."`
                        #aa_header=`echo ">"${cog_id}"|"${species_name}"|"${sequence_id}`
                        #echo "${aa_header}" >> /lustre/scratch/usr/tkerby2/pws672/Oat_Phylogeny/output_summarized/aa/${cog_id}.aa.fa
                        #echo "" >> /lustre/scratch/usr/tkerby2/pws672/Oat_Phylogeny/output_summarized/aa/${cog_id}.aa.fa
                        echo "${f}"
                fi
                # Append the header and the nucleotide sequence for each species in the same cog file
                ((INDEX = INDEX + 1))
        done
cd ../
done
# Remove the temp folder
`rm -r /lustre/scratch/usr/tkerby2/pws672/Berlandieri_Phylogeny/output_summarized/temp`
```

# FIXME
Consider in each for loop submitting a separate job using sbatch to speed up the process.

Now that we have our COG's formed we need to align them so that they can be effectively compared. First run:

```bash
mkdir -p 1_alignment/aa
mkdir 1_alignment/nt 1_alignment/aa_aligned
cp output_summarized/aa/* 1_alignment/aa/
cp output_summarized/nt/* 1_alignment/nt/
rename fa fas 1_alignment/aa/*fa
rename fa fas 1_alignment/nt/*fa
cd 1_alignment/aa
```

mafft.job content:
```bash
#!/bin/bash

#SBATCH --time=01:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=1G   # memory per CPU core
#SBATCH -J "Orthograph" #job name
#SBATCH --mail-user=tjkerby@gmail.com #email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

source ~/.bashrc
conda activate mafft

mafft-linsi $1 > ../aa_aligned/`basename $1 summarized.aa.fas`linsi.aa.fas
```

Then run the following for-loop:
```bash
for i in *fas; do sbatch mafft.job $i; done
```

After waiting for all of the jobs to finish you can then move on and copy the perl script below into the aa_aligned directory.

```perl
#!/usr/bin/perl

# Prints FASTA files with sequences in a single line each

#Change into the alignment directory and run the script. NOTE: The opuput is written on the Screen, so aou have to parse it into a file using >
#$ cd PATH_TO/aa_aligned
#$ perl fastasingleline.pl FASTAFILE.fas > FASTAFILE_noninterleaved.fas
#If you want to run more than one file at once, you have to run it using a loop, for example:
#$ for i in *.linsi.aa.fas; do perl fastasingleline.pl $i > $i.noninterleaved; done
#OUTPUT: #The output files have the additional suffix .fas.noninterleaved (non-interleaved), e.g., EOG50003G.aa.linsi.fas.noninterleaved

use strict;
use warnings;

die "Usage: $0 FASTAFILE\n" unless (scalar @ARGV);

foreach (@ARGV) {
        my $fh = Seqload::Fasta->open($_);
        while (my ($h, $s) = $fh->next_seq()) {
                printf(">%s\n%s\n", $h, $s);
        }
        undef $fh;
}

#
# Module for easy reading of fasta files
#
package Seqload::Fasta;
use strict;
use warnings;
use Carp;
require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw( fasta2csv check_if_fasta );

# Constructor. Returns a sequence database object.
sub open {
  my ($class, $filename) = @_;
  open (my $fh, '<', $filename)
    or confess "Fatal: Could not open $filename\: $!\n";
  my $self = {
    'filename' => $filename,
    'fh'       => $fh
  };
  bless($self, $class);
  return $self;
}

# Returns the next sequence as an array (hdr, seq).
# Useful for looping through a seq database.
sub next_seq {
  my $self = shift;
  my $fh = $self->{'fh'};
        # this is the trick that makes this work
  local $/ = "\n>"; # change the line separator
  return unless defined(my $item = readline($fh));  # read the line(s)
  chomp $item;

  if ($. == 1 and $item !~ /^>/) {  # first line is not a header
    croak "Fatal: " . $self->{'filename'} . " is not a FASTA file: Missing descriptor line\n";
  }

        # remove the '>'
  $item =~ s/^>//;

        # split to a maximum of two items (header, sequence)
  my ($hdr, $seq) = split(/\n/, $item, 2);
        $hdr =~ s/\s+$//;       # remove all trailing whitespace
  $seq =~ s/>//g if defined $seq;
  $seq =~ s/\s+//g if defined $seq; # remove all whitespace, including newlines

  return($hdr, $seq);
}

# Closes the file and undefs the database object.
sub close {
  my $self = shift;
  my $fh = $self->{'fh'};
  my $filename = $self->{'filename'};
  close($fh) or carp("Warning: Could not close $filename\: $!\n");
  undef($self);
}

# Destructor. This is called when you undef() an object
sub DESTROY {
  my $self = shift;
  $self->close;
}
```

You can then run that script by using the following for loop:
```bash
for i in *linsi.aa.fas; do perl fastasingleline.pl $i > $i.noninterleaved; done
```

Then rename them with:
```bash
rename fas.noninterleaved fas *noninterleaved
```

Then make a new folder by running:
```bash
mkdir 2_remove_gaps
```
 
Then copy the alligned amino acid files from step one into the step two folder by running:
```bash
cp 1_alignment/aa_aligned/*fas 2_remove_gaps/
```

Then copy the selectSites.pl script to the 2_remove_gaps directory. (The file is not displayed here since it is so long). Then run the perl script by executing the following for loop in the 2_remove_gaps directory:

```bash
for i in *fas; do perl selectSites.pl -s '1-' -x 1 $i > $i.gapsremoved; done
```

After the for loop finishes then copy the fastasingline.pl file to the directory and run it in it by executing:
```bash
for i in *gapsremoved; do perl fastasingleline.pl $i > $i.noninterleaved; done
```

Then to help keep things tidy delete and replace the extra files by running:
```bash
rm *gapsremoved
rename fas.gapsremoved.noninterleaved fas *noninterleaved
```

Next in the main directory make a new folder by running the following command:
```bash
mkdir 3_pal2nal
cp 2_remove_gaps/*fas 3_pal2nal/
cp 1_alignment/nt/*fas 3_pal2nal/
rename .aa.faslinsi.aa.fas .linsi.aa.fas *.aa.faslinsi.aa.fas
```

Next copy the pal2nal script into the directory as well as the execute.sh file. (The execute.sh file is only necessary because I couldn't get the for loop to work unless it was inside of a bash script). After this is completed run:

```bash
sh execute.sh
```

# FIXME
Consider submitting as a job to slurm

Note that it may take some time to run depending on the number of species and genes you are looking at.

After it is finished look at all of the error files and make sure that they are empty. You can check easily by running:
```bash
cat *error
```

If they are not empty something went wrong earlier on. 

Next create directories to help organize all of the output by running:
```bash
mkdir aa_aligned_non_interleaved error_files nt_aligned_interleaved nt_not_aligned
mv *.error error_files
mv *linsi.aa.fas aa_aligned_non_interleaved
mv *linsi.nt.fas nt_aligned_interleaved
mv *nt.fas nt_not_aligned
```

Copy the script fastasingleline.pl to nt_aligned_interleaved and then change into the directory and then run:
```bash
for i in *linsi.nt.fas; do perl fastasingleline.pl $i > $i.noninterleaved; done
```

Then clean up by running:
```bash
rename fas.noninterleaved fas *noninterleaved
```

Then create a new directory in the main project folder and copy the necessary files into it by running:
```bash
mkdir -p 4_rewrite_headers/aa_headers_renamed
mkdir 4_rewrite_headers/nt_headers_renamed
cp 3_pal2nal/aa_aligned_non_interleaved/*linsi.aa.fas 4_rewrite_headers/aa_headers_renamed/
cp 3_pal2nal/nt_aligned_interleaved/*linsi.nt.fas 4_rewrite_headers/nt_headers_renamed/
```

Copy the following perl script into nt_headers_renamed and aa_headers_renamed directories.
```perl
#!/usr/bin/perl
use strict;
use warnings;
use autodie;

while (my $file = <*.fas>) {                            #definiert "file" als Variable nur in dieser Schleife defines "file" as a variable only in this loop

        open my $ifh, '<', $file;                       #ordnet dann jeweiligen "Namen" zu fuer jedes Gen then assigns each "name" for each gene to
        open my $ofh, '>', "$file.fas.cleaned";

        while (my $line = <$ifh>){
                chomp $line;
                if ($line=~ /^>/){
                        my @fields = split /\|/, $line;
                        if (scalar @fields == 3) {
                                $line = $fields[1];
                        }
                        elsif (scalar @fields == 4) {
                                $line = $fields[2];
                        }
                        $line =~ s/[^A-Za-z0-9_]/_/g;
                        $line = '>' . $line;
                }
                print $ofh $line, "\n";
        }
        close $ifh;
        close $ofh;
}
```

Create a perl-autodie conda environment by running:

```bash
source ~/.bashrc
conda create -n perl-autodie -c bioconda perl-autodie
conda activate perl-autodie
```

Now switch into the aa_headers_renamed directory and run the rewrite_headers.pl script. Then clean it up by using the rename command.
```bash
perl rewrite_headers.pl
rename fas.fas.cleaned fas *cleaned
```

Next make another directory in the main project folder and copy in the necessary files by running:
```bash
mkdir -p 5_fasconcat/nt_concatenated
mkdir 5_fasconcat/aa_concatenated
cp 4_rewrite_headers/nt_headers_renamed/*linsi.nt.fas 5_fasconcat/nt_concatenated/
cp 4_rewrite_headers/aa_headers_renamed/*linsi.aa.fas 5_fasconcat/aa_concatenated/
```

Next copy the FASconCAT_v1.0.pl perl script into each of the subdirectories as well as the extract_part_def_nt.py and extract_part_def_aa.py into their respective directories. Then move into the nt_concatenated directory and run the FASconCAT_v1.0.pl script and then press i enter and then s enter.
```bash
perl FASconCAT_v1.0.pl
```

You will be prompted with a screen that requests the type of analysis that you want to do. Press i to ensure that you print out the information file and then
press s to start the program.
When you are finished, there should be two new files in your directory: the FASTA file containing the concatenated supermatrix, FcC_smatrix.fas , and
a file containing information about what you concatenated, FcC_info.xls . You will be using the information file to make new, special file, called a partition
file. To do this run the script extract_part_def_nt.py on the FcC_info.xls file:

```bash
python extract_part_def_nt.py FcC_info.xls
```

This will create a new file called part_def_nt.txt . Now rename the FcC_smatrix.fas file to FcC_smatrix_nt.fas to ensure you don't get
confused about the file downstream.

```bash
mv FcC_smatrix.fas FcC_smatrix_nt.fas
```

Now, you will be performing the same steps on the amino acids.

```bash
cd ../aa_concatenated
perl FASconCAT_v1.0.pl
python extract_part_def_aa.py FcC_info.xls
mv FcC_smatrix.fas FcC_smatrix_aa.fas
```
