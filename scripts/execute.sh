#!/bin/sh

for i in *.linsi.aa.fas; do
	a=`basename $i .linsi.aa.fas`
	perl pal2nal.mod.pl $a.linsi.aa.fas $a.nt.fas -output fasta 1>$a.linsi.nt.fas 2>$a.error
done
