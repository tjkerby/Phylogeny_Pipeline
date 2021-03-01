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
		header=`echo ">"${cog_id}"|"${species_name}"|"${sequence_id}`		
		# Handling if a certain species does not have the gene
		if [ -s ${f} ]
		then 			
			# Convert from nucleotides to amino acids
			`transeq ${f} /lustre/scratch/usr/tkerby2/pws672/Berlandieri_Phylogeny/output_summarized/temp/${f} -clean`
			aa_sequence=`echo "$(awk '{if(NR>=2) printf $0;}' /lustre/scratch/usr/tkerby2/pws672/Berlandieri_Phylogeny/output_summarized/temp/${f})"`
			# Append the header and the amino acid sequence for each species in the same cog file
			echo "${header}" >> /lustre/scratch/usr/tkerby2/pws672/Berlandieri_Phylogeny/output_summarized/aa/${cog_id}.aa.fa
			echo "${aa_sequence}" >> /lustre/scratch/usr/tkerby2/pws672/Berlandieri_Phylogeny/output_summarized/aa/${cog_id}.aa.fa
			# Append the header and the nucleotide sequence for each species in the same cog file
        	        echo "${header}" >> /lustre/scratch/usr/tkerby2/pws672/Berlandieri_Phylogeny/output_summarized/nt/${cog_id}.nt.fa
	                echo "${nt_sequence}" >> /lustre/scratch/usr/tkerby2/pws672/Berlandieri_Phylogeny/output_summarized/nt/${cog_id}.nt.fa
		else
			echo "${f}"
		fi
		((INDEX = INDEX + 1))
        done
cd ../
done
# Remove the temp folder
`rm -r /lustre/scratch/usr/tkerby2/pws672/Berlandieri_Phylogeny/output_summarized/temp`
