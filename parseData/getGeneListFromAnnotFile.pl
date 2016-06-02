#!/usr/bin/perl
use strict;
my $annotFile=$ARGV[0];
die "You must specify an annot file (FDI or AVIA formatted)\n" if ($#ARGV==-1);
if (! -e $annotFile ){
	die "cannot read your annot file ($annotFile)\n";
}

my $output=($#ARGV>0)?$ARGV[1]:"genelist";
my $col=-1;
# if ($annotFile=~/annovar_wrpr.output/){
	my @headers=split("\t",`head -n1 $annotFile`);chomp $headers[$#headers];
	for (my $i=0;$i<=$#headers;$i++){
		if ($headers[$i]=~/^(Ensembl_){0,1}Gene$/){
			$col=($i+1);#because we are using cut, we have to add one to offset the 0 in perl
			print "Found column #$col for genes\n";
		}
		last if ($col>-1);
	}
# }else{
# 	$col=4;
# }

# print ("grep -e '\\bexonic\\b' -e intronic -e splic -e UTR  $annotFile | cut -f$col |grep -ve NONE |cut -f1 -d '(' | sort -u > $output\n");
my @arr=split("\n",`grep -e '\\bexonic\\b' -e intronic -e splic -e UTR  $annotFile | cut -f$col |grep -ve NONE -ve Gene |cut -f1 -d '(' | sort -u `);
open (OUT,">$output") or die "cannot open $output for writing in $0\n";
my $ensembl=0;
foreach my $gene (@arr){
	chomp $gene;
	if ($gene=~/ENSG\d+/){
		$ensembl++;
	}
	my @list;my $delim;
	if ($gene=~/([;,:])/){
		$delim=$1;
		@list=split("$delim",$gene);
	}else{
		push(@list,$gene);
	}
	foreach my $id (@list){
		next if (!$id);
		print OUT "$id\n";
	}
}
close OUT;
if ($ensembl){
	system ("cp $output $output\_ensembl\n");
}
