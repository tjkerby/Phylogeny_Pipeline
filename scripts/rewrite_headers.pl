#!/usr/bin/perl
use strict;
use warnings;
use autodie;

while (my $file = <*.fas>) {				#definiert "file" als Variable nur in dieser Schleife defines "file" as a variable only in this loop
	
	open my $ifh, '<', $file; 			#ordnet dann jeweiligen "Namen" zu fuer jedes Gen then assigns each "name" for each gene to
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