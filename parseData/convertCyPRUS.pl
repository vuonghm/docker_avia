#!/usr/bin/perl
use strict;
use Data::Dumper;
die "You must specify an input file\n" if ($#ARGV==-1);
my $input=$ARGV[0];
my %featCount;
open (FILE, "<$input") or die "Could not open file $input\n";
my %uniprotRes;my $input=my $feat='';
while (<FILE>){
	if ($_=~/Input":\s+"(\S+)"/){
		$input=$1;
		$feat='common';
		# $uniprotRes{$input}{$feat}{'gene'}=$$input;
	}elsif($_=~/"(.*?)":\s+\"(.*)\"/){
		my ($key,$value)=($1,$2);
		if (!$value){$value='-';}
		if ($key=~/Feature/){
			# $feat=$2;
			$uniprotRes{$input}{$value}={};
			$featCount{$value}=0 if (!exists $featCount{$value});
			$feat=$value;
		# }elsif($key=~/(Source urls)/ && $feat=~/sequence variant/){
			# $uniprotRes{$input}{$feat}{$key}." "
		}elsif ($key=~/Description/i){
			$uniprotRes{$input}{$feat}{$key}=$value;
			$featCount{$feat}++ if ($value!='-');
			# print "Adding $value to  $feat->$key for $input\n";
		}else{
			$uniprotRes{$input}{$feat}{$key}.="$value," if ($uniprotRes{$input}{$feat}{$key}!~/$value,/);
		}
	}
}
print Dumper (\%featCount);exit;
my $fileinfo= "{\n";
foreach my $key (keys %uniprotRes){
	chop($uniprotRes{$key}{'common'}{'Uniprot accession'});
	$fileinfo.= "\"$key\":{";
	$fileinfo.= "\"accession\":[\"" . $uniprotRes{$key}{'common'}{'Uniprot accession'} . "\"],";
	foreach my $feature (keys(%{$uniprotRes{$key}})){
		next if $feature=~/common/;
		next if ($featCount{$feature}==0);
		$fileinfo.= "\"$feature\":{\"description\":[\"" . 
			$uniprotRes{$key}{$feature}{'Description'} . "\"],\"pubmed\":[\"".
		 	$uniprotRes{$key}{$feature}{'pubmeds'} . "\"]},";

	}
	chop $fileinfo;
	$fileinfo.="},";

}
chop $fileinfo;
$fileinfo.= "}";
print $fileinfo;