#!/usr/bin/perl
use strict;
use Data::Dumper;
umask (000);
die "Usage:  You must supply an input!!\n$0 <ANNOT_FN> [<WORKING DIR>] [1|0|-1]
	where the third argument specifies whether or not to look for header 'Gene'\n" if ($#ARGV==-1);
open (FILE,"<$ARGV[0]") or die "Cannot open $ARGV[0] for reading\n";
my @lines=<FILE>;
my $gene_idx=3;
close FILE;
my $cwd;
my $bin=`dirname $0`;chomp $bin;
my $GNF_base="$bin/GNF";
my %gene_arr;
my $psa=0;
if ($#ARGV==2 && $ARGV[2]==-1){
	#genelist is provided, just read in genelist
	$cwd=$ARGV[1];
	print "$cwd is the cwd\n";
	foreach my $line(@lines){
		$line=~s/[\r\n\s]//g;
		push(@{$gene_arr{$line}},"$line\tN/A");
	}
}elsif ($ARGV<2 || ($#ARGV==2 && $ARGV[2]!=-1)){#die "$#ARGV and $ARGV[2]\n";
	if ($ARGV[0]=~/^\//){
		$cwd=`dirname $ARGV[0]`;chomp $cwd;
	}elsif ($#ARGV>1){
		$cwd=$ARGV[1];
	}else{
		$cwd=`pwd`;chomp $cwd;
	}
	chdir("$cwd");
	#this flag represents that you have to lookup the gene columns!
	my @headers=split("\t",$lines[0]);
	my $feat=-1;my $protpos=-1;
	$gene_idx=-1;
	for (my $i=0;$i<=$#headers;$i++){
		if ($headers[$i]=~/^#{0,1}(Ensembl){0,1}[\s\_]{0,1}Gene$/){
			$gene_idx=$i;
			print "$headers[$i] is found\n";
		}elsif ($headers[$i]=~/(ANNOVAR annot)/i){
			$feat=$i;
		}elsif($headers[$i]=~/(sift|pp2|Provean|ma$|mt$|fathmm|Mutation Assessor|Mutation Taster|Polyphen)/i){
			$psa++;
		}elsif ($headers[$i]=~/ProtPos/){
			$protpos=$i;
		}
	}
	die "Could not find the header 'Gene' in your file $ARGV[0]\n" if ($gene_idx==-1 || $feat==-1);
	$protpos=($feat-1) if ($protpos<0);
	print "Gene index:$gene_idx...feat idx=$feat\n";
	#first aggregate all the mutations per gene
	foreach my $line (@lines){
		my @arr=split("\t",$line);my @trash;
		($arr[$gene_idx])=split('\(',$arr[$gene_idx]);
		chop $arr[$gene_idx] if ($arr[$gene_idx]=~/[;:,]$/);
		my $count=0;
		
		print "$arr[$feat] should be Feat\n";
		next if ($arr[$feat]=~/ncRNA/i);
		if ($arr[$feat]=~/(\bexonic\b)/i){
		 	while ($line =~ /(damaging|disease.causing|High|Med|Deleterious)/gi) { $count++; }
			push(@{$gene_arr{$arr[$gene_idx]}},"$arr[($protpos)]\t$count");
		}elsif ($line=~/(\bexonic\b|intronic|splic|UTR|upstream)/i){
			push(@{$gene_arr{$arr[$gene_idx]}},"$arr[($protpos)]\t$count");
		}
	}
}else{
	die "what?$#ARGV";
}
# print Dumper (\%gene_arr);exit;
my %mapper;
if ($ARGV[2] ne '1' && -e $ARGV[2]){
	print "Reading $ARGV[2]\n";
	my $err=0;
	open (MAPPER,"<$ARGV[2]") or $err=1;
	if (!$err){
		while (<MAPPER>){
			my @arr=split("\t",$_);chomp $arr[$#arr];
			$mapper{$arr[0]}=$arr[1];
			$mapper{"_$arr[1]"}=$arr[0];
		}
	}
}
my $legend;
foreach my $key (sort keys %gene_arr){
	# print "working on $key...?\n";
	my $content= "&nbsp;";
	my $total=0;
	foreach my $foo (@{$gene_arr{$key}}){
		# print $foo."\n";
		my @arr=split("\t",$foo);
		my $addon='';
		$total++;
		next if $arr[1]==0;
		# print 'adding content'. $arr[0]."\n";<STDIN>;
		$content.= "<span title=\"$arr[0]\" $addon>$arr[1]/$psa</span>,";
	}
	chop $content if ($content=~/,$/);
	$gene_arr{$key}=(exists $mapper{$key})?"$mapper{$key} ($key)":$key;
	$gene_arr{$key}.="\t$total\t$content";
}
$gene_arr{'Gene'}="Gene\tNo. of Mutations\tPred. Damaging Counts*";
# print Dumper (\%gene_arr);
my $colcount=0;
foreach my $tissue ("Gland","Brain","Germ","Nerve","Muscle","Cancer","Immune","Digestive","Pineal","Retina","Other"){
	open (FILE2,"<$GNF_base/$tissue.txt") or die "cannot open $GNF_base/$tissue.txt\n";
	open (OUTPUT,">$cwd/GNF_$tissue.txt") or die "Cannot open $cwd/GNF_$tissue.txt";
	while (<FILE2>){
		my @filearr=split("\t",$_);chomp $filearr[$#filearr];
		my $gene=shift(@filearr);
		my $xgene=$gene;my $addon='';
		if (exists $mapper{"_$gene"}){
			$xgene=$mapper{"_$gene"};
			# die "exists in mapper\n";
		}else{
			# print "does not exist in mapper($xgene)\n";<STDIN>;
		}
		if (exists $gene_arr{$xgene}){
			print OUTPUT "$gene_arr{$xgene}\t".join ("\t",@filearr)."\n";
			$gene_arr{$xgene}.="\t".join ("\t",@filearr);
		# }else{
		# 	print "$xgene does not exist in gene_arr!\n";
		}
	}
	# print OUTPUT "Legend: $legend";
	close OUTPUT;
}
open (ALL,">$cwd/GNF_All_Tissues.txt") or die "Cannot open All.txt";
print ALL "$gene_arr{'Gene'}||\n";
my $content='';
foreach my $gene (sort keys %gene_arr){
	next if ($gene=~/Gene/);
	my $addon='';
	print ALL "$gene_arr{$gene}\t";
	# print "($gene_arr{$gene})";<STDIN>;
	if ($gene_arr{$gene}=~/(\&nbsp\;|span\>)$/){
		for (my $i=0;$i<=83;$i++){
			print ALL "NA\t";
		}
	}
	print ALL "\n";
}
close ALL;
