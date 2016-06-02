#!/usr/bin/perl
use strict;
use Data::Dumper;
my $file=$ARGV[0];
die "Please specify a file\n" if (!$file);
open (FILE,"<$file" ) or die "Cannot open file ($file) for reading\n";
my (%pop,%freq);
while (my $line=<FILE>){
	chomp $line;
	if ($line=~/##INFO=<ID=(\w+),/ ){
		my $stuff=$1;
		next if ($stuff!~/(_|^(CSQ|DP|HWP|HaplotypeScore|)$)/);
		my ($type,$population)=split("_",$stuff);
		if ($population){
			$pop{$population}{$type}="";
		}else{
			$pop{$type}='';
		}
	}elsif ($line!~/#/){
		my ($chr,$start,$id,$ref,$alt,$qual,$filter,$info)=split(/\t/,$line);
		my @infoArr=split(/;/,$info);
		foreach my $allele (@infoArr){
			last if ($allele=~/\|/);
			print "working on $allele\n";
			
		}
		exit;

	}
}