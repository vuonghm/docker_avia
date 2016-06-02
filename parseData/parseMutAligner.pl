#!/usr/bin/perl
use strict;
my $usage=qq($0 <Mutation Aligner input file> [<output filename>]
	#This script produces a file that can be read by AVIA/ANNOVAR wrapper script
	BRAF:V600->domain->#mutations->#otherRelatedGenes
);

if ($#ARGV==-1){
	print $usage;exit;
}
my $file=$ARGV[0];
my $out="mutAligner.txt";
if ($#ARGV>=1){
	$out=$ARGV[1];
}
#see mut_aligner_snippet.txt in folder for example for this script
open (FH,"<$file") or die "Cannot open $file\n";
open (OUT,">$out") or die "Cannot open $out\n";
open (OUT2,">$out.2") or die "Cannot open $out\n";
print "writing to $out??\n";
while (<FH>){
	next if ($_=~/^\s*$/);#empty line or lines with spaces only
	my ($id,$domain,$numGenes,$posAln,$pVal,$entropyScore,$mutCount,$topGenes,$topCancer)=split("\t",$_);
	chomp $topCancer;
	$topGenes=~s/^\d+\s+//;
	$topGenes=~s/;\s+\d+\s/;/g;
	$topGenes=~s/\s/:/g;
	my @topGenes=split(";",$topGenes);
	$topCancer=~s/\s*\d+\s+//g;
	$topCancer=~s/;/,/g;
	foreach my $topGene (@topGenes){
		$topGene=~s/\s/:/;
		my $related=$topGenes;
		$related=~s/\s+\d+\s//g;
		$related=~s/(;$topGene|$topGene;)//;
		$related=~s/\s/:/g;
		$related=~s/;/,/g;
		print OUT join ("\t",$topGene,$mutCount,$domain,$topCancer,$related)."\n";
		print OUT2 join ("\t",$topGene,"Domain=$domain;CancerTypes=$topCancer;RelatedDomain=$related")."\n";

	}
}
close FH;close OUT;