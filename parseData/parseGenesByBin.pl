#!/usr/local/bin/perl
use strict;
use FindBin;
use Data::Dumper;
use lib "$FindBin::Bin/";
use ABCC_PARSE;
use Getopt::Std;
 umask(0000);
use vars qw($opt_i $opt_g $opt_c $opt_k $opt_x $opt_l $opt_z $opt_I $opt_o $opt_E);
getopts("i:g:c:l:I:o:kxzEa");#opt_z means it is highlighted
die "You must specify an input file $0 -i/I <annovar or inputfile> [-l <genelist fn> -c <cutoff> -k <keyword> -o <basename for output files> -E <use Ensembl Gene>]" if (!defined $opt_i && !defined $opt_I);
my $fn;
if ($opt_i){
	$fn=$opt_i;
}else{
	$fn=$opt_I;
}
if (defined ($opt_o)){
	$opt_o=~s/\.txt//g;
	$opt_o.="." if ($opt_o!~/\.$/);
}else{
	$opt_o="$fn.";
}
my $gene_fn_ct=$opt_g;
my $cutoff=(defined $opt_c)?$opt_c:0;
if ($cutoff!~/^[\<\>\=]/){
	$cutoff="<=$cutoff";
}
my (%genelist);
if (defined $opt_l && -e $opt_l){
	print "Reading genelist($opt_l)\n";
	readGenes($opt_l);
	open (FILTERED,">filtered_genelist_annot.txt") or die "Cannot open filtered_genelist_annot\n";
}
my $keyword=$opt_k;
$gene_fn_ct ||=0;
die "Your file is empty or does not exist\n" if (! -e $fn || -z $fn);
my $addon;
if ($opt_i=~/genes_/ || $opt_I=~/genes_/ 	){
	$addon=" -d ' ' ";
}
die "Your file does not exists\n" if (!-e $fn);
print STDERR "[$0]cut -f1 $addon $fn | sort -u | grep -ve '_' -ve '#'\n";
my @chroms=split("\n",`cut -f1 $addon $fn | sort -u | grep -ve '_' -ve '#'`);
my (%h10k,%h100k,%h1000k);
my $genehead=`head -n1 $fn`;chomp $genehead;
my $genestr;
if ($genehead=~/^(.*\tGene\t)/){
	$genestr=$1;
}elsif ($opt_I){
	$genestr=3;
}
print "$genehead\n" if (tell(FILTERED) != -1);#print header to the new filtered file
my $geneidx=0;
if ($genestr!~/^\d+$/){
	while ($genestr =~ /\t/g) { $geneidx++ }
	$geneidx--;
}elsif ($genestr=~/^\d+$/){	
	$geneidx=$genestr;
}else{
	die "Couldn't find gene idx\n";
}
foreach my $chrom (@chroms){
	chomp $chrom;
	print STDERR "working on CHR($chrom)\n";
	my @arr=split("\n",`grep -P '^$chrom\[\t \s]' $fn`);
	ELEMENT: foreach my $element (@arr){
		chomp $element;#die "workin on $element\n";
		my ($chr,$start,$stop,@rest)=split(/\s/,$element);
		my $id=$rest[$geneidx-4];
		$chr=~s/^chr/hs/;if ($chr!~/^hs/){$chr="hs$chr";}
		next if ($chr!~/hs[XY0-9M]{1,2}/i);
		my $pass=0;
		next if ($element=~/(intergenic)/i);
		if (defined $opt_l){
			my @genes=split(/[;,:]/,$rest[$geneidx-3]);
			if ($#genes>-1){
				foreach my $gene (@genes){
					$gene=uc($gene);
					if (exists $genelist{$gene}){
						$pass++;
						print FILTERED $element."\n";
					}
					last if ($pass);
				}
				next ELEMENT if ($pass==0 && !$opt_z ); #on the gene list and not flagged for highlighting
			}else{
				next ELEMENT;#not in a gene
			}
			# print "pass $rest[$geneidx-3]\n";
		}
		my $annot='';
		if ($id=~/annot=(.*)/){
			$annot=$1;
		}elsif ($opt_c && $id=~/^[\d\.\-]*/ && $opt_x){#opt_x is histograms
			if (eval ($id.$cutoff)){
				$annot='damaging' ;
			}elsif ($id eq '-'){
				if ($rest[0]=~/splic/){
					$annot="damaging";
				}elsif ($rest[0]=~/intronic/){
					$annot="intronic";
				}else{ #syn
					$annot='nondamaging';
				}
			}else{ #nonsyn not above threshold
				$annot='nondamaging';
			}#works
		}elsif ($opt_c && $id=~/^[\d\.\-]*/ && !$opt_x){#not histogram
			if (eval ($id.$cutoff)){
				$annot='passed' ;
			}#works
		}elsif ($opt_k && $id=~/$keyword/){#only count gplexes for example
			$annot='passed';
		}
		for my $binsize (10_000 ,100_000 ,1_000_000 ) {
			my $bin=sprintf("%d",($start/$binsize));
			my $bin2=sprintf("%d",($stop/$binsize));
			my $href;
			if ($binsize==10_000){
				$href=\%h10k;
			}elsif ($binsize==100_000){
				$href=\%h100k;
			}elsif ($binsize==1_000_000){
				$href=\%h1000k;
			}else{
				next;
			}
			$$href{"$chr\t$bin"}{'ct'}++;
			$$href{"$chr\t$bin2"}{'ct'}++ if ($bin!=$bin2);
			if ($annot){
				$$href{"$chr\t$bin"}{$annot}++;
			}
			if ($pass){#if it was in the list
				$$href{"$chr\t$bin"}{'flag'}=1;
			}
		}
	}
}
close FILTERED ;
my (%totals,$color,$trash);
($color,$trash)=color();
for my $binsize (10_000 ,100_000,1_000_000  ) {#10_000 ,,1_000_000
	my $href;my $outname=($binsize/1000);
	if ($gene_fn_ct){
		open (FILE,"<$gene_fn_ct\.$outname\k.txt") or die "[ERROR] LN".__LINE__.":$!-does not exist $gene_fn_ct\.$outname\k.txt\n";
		while (<FILE>){
			my ($bin,$start,$stop,$total)=split(/[\s\t]/,$_);chomp $total;
			$totals{"$bin:$start-$stop"}=$total;
		}
		close FILE;
	}
	
	open (FILE,">$opt_o$outname\k.txt") or die "Cannot open $opt_o$outname\k.txt\n";
	if ($binsize==10_000){
		$href=\%h10k;
	}elsif ($binsize==100_000){
		$href=\%h100k;
	}elsif ($binsize==1_000_000){
		$href=\%h1000k;
	}else{
		next;
	}
	foreach my $key (sort {$a cmp $b}keys %{$href}){
		my ($chr,$bin)=split("\t",$key);
		my $start=($bin*$binsize);
		my $end=(($bin*$binsize)+$binsize-1);
		my $fill_color='';
		if ($$href{$key}{'flag'} ){
			$fill_color= "fill_color=dred";
			my $name;
			if ($opt_i=~/subtractive\.(\S+)\.txt(\.tmp){0,1}/){
				$name=$1;
			}
			if ($name){
				my $type=`grep $name circos.input|cut -f3`;chomp $type;
				if ($type=~/(text|tile|line)/i){
					$fill_color="color=dred";
				}
			}
		}
		if ($gene_fn_ct){
			my $id=(exists $totals{"$chr:$start-$end"})?$totals{"$chr:$start-$end"}:0;
			die "error there are no genes here \n" if (!$id);
			if (scalar (keys %{$href->{$key}})>1){
				foreach my $annot (keys %{$href->{$key}}){
					next if $annot eq 'ct';
					my $value=sprintf("%.03f",($$href{$key}{$annot}/$id));
					print FILE join (" ",$chr,$start,$end,$value,"$fill_color")."\n" if ($value>0);
				}
			}else{
				my $value=sprintf("%.03f",($$href{$key}{'ct'}/$id));
				print FILE join (" ",$chr,$start,$end,$value,$color)."\n" if ($value>0);
			}
		}elsif ($opt_x){#histograms; need multiple values
			if (scalar (keys %{$href->{$key}})>1){
#20120620			print FILE join (" ",$chr,$start,$end,sprintf("%.03f",($$href{$key}{'passed'}/$$href{$key}{'ct'})))."\n";
				$$href{$key}{'damaging'}||=0;
				$$href{$key}{'nondamaging'}||=0;
				$$href{$key}{'intronic'}||=0;
				my $perc_damage=sprintf("%0.10f",($$href{$key}{'damaging'}/$$href{$key}{'ct'}));
				my $perc_nondam=sprintf("%0.10f",($$href{$key}{'nondamaging'}/$$href{$key}{'ct'}));
				my $perc_intronic=sprintf("%0.10f",($$href{$key}{'intronic'}/$$href{$key}{'ct'}));
				print FILE join (" ",$chr,$start,$end,"$perc_damage,$perc_nondam,$perc_intronic",$fill_color)."\n";
			}
		}elsif ($opt_c){
			if ($href->{$key}->{'intronic'}!= $$href{$key}{'ct'} ){
				$$href{$key}{'passed'}||=0;
#delta 20120620 next block
#				$$href{$key}{'damaging'}||=0;
#				$$href{$key}{'nondamaging'}||=0;
#				$$href{$key}{'intronic'}||=0;
#				my $perc_damage=sprintf("%0.10f",($$href{$key}{'damaging'}/$$href{$key}{'ct'}));
#				my $perc_nondam=sprintf("%0.10f",($$href{$key}{'nondamaging'}/$$href{$key}{'ct'}));
#				my $perc_intronic=sprintf("%0.10f",($$href{$key}{'intronic'}/$$href{$key}{'ct'}));
#				print FILE join (" ",$chr,$start,$end,"$perc_damage,$perc_nondam,$perc_intronic")."\n";
				my $perc_passed=sprintf("%0.10f",($$href{$key}{'passed'}/$$href{$key}{'ct'}));
				print FILE join (" ",$chr,$start,$end,$perc_passed,$color)."\n" if ($perc_passed > 0);
			}
		}elsif ($opt_k){#counts
			print FILE join (" ",$chr,$start,$end,$$href{$key}{'ct'},$fill_color)."\n" if ($$href{$key}{'ct'}>0);
		}else{
			my $value=sprintf("%0.10f",($$href{$key}{'ct'}/$binsize));
			print FILE join (" ",$chr,$start,$end,$value,$fill_color)."\n" if ($value>0);
		}
	}
	close FILE;
}


sub readGenes{
	my $rgfn=shift;
	my $annot_fn="/SeqIdx/annovardb/humandb/hg19_refGene.txt" ;
	if ($opt_E){
		$annot_fn="/SeqIdx/annovardb/humandb/hg19_ensGene.txt";
	}
	die "$annot_fn does not exist\n" if (!-e "$annot_fn" && $opt_I);
	open (GENES,"<$rgfn" )  or print "Could not open $rgfn" and return;
	while (<GENES>){
		next if ($_=~/^[\s\r]{0,}$/);#skip blanks
		my @arr=split("\t",$_);chomp $arr[$#arr];
		if ($#arr<1){
			$arr[0]=~s/[\r\s ]//g;
			if (defined $opt_I){
				my @gene_arr=split("\n",`grep -P '\t$arr[0]\t' -i $annot_fn`);
				next if ($#gene_arr==-1);
				my $min= my $max=0;
				foreach my $line (@gene_arr){
					my @features=split("\t",$line);
					if ($min==0 || $min>$features[4]){
						$min=$features[4];
					}
					if ($max==0 || $max<$features[5]){
						$max=$features[5];
					}
				}
				
			}
			$genelist{uc($arr[0])}=1;
		}
	}
	print STDERR Dumper (\%genelist);
	close GENES;
	return;
}