#!/usr/bin/perl
use strict;
use FindBin;
use Getopt::Std;
use Data::Dumper;
use vars qw($opt_o $opt_b $opt_a $opt_f $opt_m $opt_M $opt_h $opt_v);
getopts("a:f:o:b:m:M:v:h");
umask(0000);
my $bin=`dirname $0`;chomp $bin;
my $usage=qq(
	$0
		-a Annotation file from avia looks for exonic mutations only
		-f biodata output file
	[OPTIONAL]
		-o <UNiprot output name>
		-m <mapping filename>
		-M <ANNOVAR output name for aggregation>
		-h <print a header>
	);
die "You must enter a valid file " if (!defined ($opt_f) && !defined $opt_a);
my $in_fn;
if (defined $opt_f){
	$in_fn=$opt_f;
}elsif (defined $opt_a){
	$in_fn=$opt_a;	
}

open (ANNOT,"<$in_fn") or die "Cannot open input file $in_fn $?\n";
$opt_o||= "$in_fn.prot";
my %mapping;
if ($opt_m){
	open (MAP,"<$opt_m") or die "Cannot open mapping $opt_m\n";
	while (<MAP>){
		chomp;
		next if ($_=~/^\s+$/);
		my @arr=split ("\t",$_);
		$mapping{$arr[0]}=$arr[1];
	}
}
$opt_M||=$in_fn.".hg19_ProtPos";
$opt_v||='9606';
my $nsCount=0;my %bl_mapping;
open (VAR,">$opt_M") or die "Cannot open $in_fn.uniprot.map\n";
open (OUTPUT,">$opt_o") or die "Cannot open output file $opt_o\n";
if (defined $opt_a){
	my $count=my $cdsCount=0;my %prothash;
	while (my $line=<ANNOT>){
		my $myline;
		if ($line=~/((chr){0,1}[\dXYMT]{1,}[\t:]\d+[\t:]\d+[:\t][ATCG][\t:][ATCG])/){
			$myline=$1;$myline=~s/^chr//;$myline=~s/\t/:/g;
			# print VAR "$myline\t";
			$myline='' if (!$opt_m);$count++;
		}elsif($line=~/((chr){0,1}[\dXYMT]{1,}[\t:]\d+[\t:]\d+[:\t][ATCG\.\-]*[\t:][\.\-ATCG]*)/){
			$myline=$1;$myline=~s/^chr//;$myline=~s/\t/:/g;$count++;
			# print VAR "$myline\t";
			# print "found $myline...\n";
		}elsif ($line=~/#/){
			print VAR "#ProtPos\n" if ($opt_h);
			next;
		}else{
			# print "Don't know how to parse this input type: in $0:$line\n";next;
		}
		if ($line=~/(stopgain.SNV|stoploss.SNV|synonymous SNV):(\S+)/){
			my @annot=split(",",$2);$cdsCount++;
			
			foreach my $feat (@annot){
				if ($feat=~/(([\w\-\*]+):(NM\_\d+|ENST\d+):.*:p\.([\w\*]+)):?/){
					my ($gene,$nm,$protpos)=($2,$3,$4);
					push (@{$prothash{"$gene:$protpos"}},$line) and $nsCount++ if ($line=~/(nonsynonymous)/);
					push (@{$prothash{"$nm:$protpos"}},$feat);
					print VAR "$nm:$protpos," ;$bl_mapping{$nm}{"$nm:protpos"}=1;
					# print "what to do with $foo??\n";
					# print MAP "$1:$2\t$myline\n" if ($opt_m);
				}else{
					warn "Couldnt parse SNP($feat)...".substr($line,0,100)." in $0\n";$cdsCount--;
				}
			}
			print VAR "\n";
		}elsif($line=~/(frameshift[^:]*):(\S+)/){
			my @annot=split(",",$2);$cdsCount++;
			foreach my $feat (@annot){
				if ($feat=~/([\w\-\*]+):(NM\_\d+|ENST\d+):.*:p\.([^:]+):?/){
					my ($gene,$nm,$protpos)=($1,$2,$3);
					# print "$gene,$nm and $protpos\n";<STDIN>;
					if ($line=~/nonsynonymous/){
						push (@{$prothash{"$gene:$protpos"}},$line) ;
						push (@{$prothash{"$nm:$protpos"}},$feat);
					}
					print VAR "$nm:$protpos,";$bl_mapping{$nm}{"$nm:protpos"}=1;
					# print "\ypushed $myline |$1:$2\n";
					# print MAP "$1:$2\t$myline\n" if ($opt_m);
				}else{
					# print "[WARN]Couldnt parse INDEL($feat)...".substr($line,0,100)." in $0\n";
				}
			}
			print VAR "\n";
		}else{
			print VAR "-\n";
		}
	}
	# print Dumper (\%prothash);die;
	# print keys(\%prothash)."\n";exit;
	my %bl;
	my $avia_lite='/bioinfoA/scripts/biodbnet/easyConverter_avia.php';
	if (-e $avia_lite && keys (%bl_mapping) > 0 ){
		my ($dboutType,$trns);
		my $trns= join (",",keys(%bl_mapping));
		if ($trns=~/ENST/){
			$dboutType='Ensembl Transcript ID';
		}else{
			$dboutType='RefSeq mRNA Accession';
		}
		%bl = map { split(/\n/, $_) } split(/\t/, `php $avia_lite -F $trns -c 0 -d ':' -i '$dboutType' -o 'Ensembl Protein ID' -v $opt_v` );
		$opt_m.=2;
		open(BL,">$opt_o-Ensembl") or die "Cannot open file $opt_M-Ensembl";
		print BL "#$dboutType\tEnsembl Protein ID\n" if ($opt_h);#header
	}
	foreach my $id (keys %prothash){
		my ($gene,$prot_pos)=split(":",$id);
		if ($opt_m=~/\w+/  ){# it is a gene and the opt_m was defined at startup
			if( $gene!~/^(NM_|ENST)/i && (exists $mapping{$gene})){
				print OUTPUT "$mapping{$gene}:$prot_pos\n";
			}elsif ($gene!~/^(NM_|ENST)/){
				print OUTPUT "$gene:$prot_pos\n";
			}
		}else{
			print OUTPUT "$id\n";
		}
		if ($gene=~/^(NM_|ENST)/  && $opt_m=~/2/){# for the extra file...
			if (exists $bl{$gene}){
				my @ids=split(";",$bl{$gene});
				foreach my $idx (@ids){
					next if (!$idx || $idx!~/ENSP/);
					my $foo=$id;
					$foo=~s/$gene/$idx/g;
					print BL "$foo\t$gene->$idx\n";
				}
			}
		}
	}
	if ($nsCount>1){
		print "Total non-synonymous SNPs: <b>$nsCount</b>; Total mutations, including synonymous and indels in protein coding regions: <b>$cdsCount</b>; Total mutations in the input file: <b>$count</b>\n";
	}
	close OUTPUT;close ANNOT;
	if ($cdsCount==0){
		unlink ($opt_o);
	}
	
}else{
	print "converting json from biodata\n";
	while (/JSON:/ .. /\n/){
		print $_;
	}
}
