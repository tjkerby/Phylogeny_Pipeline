#!/bin/sh

for dir in /lustre/scratch/usr/tkerby2/pws672/Berlandieri_Phylogeny/atram_output/*/ ;
do
	cd ${dir}
	name=`basename "${dir}"`
	mkdir filtered
	mv *.filtered_contigs.fasta filtered/
	cd ..
done
