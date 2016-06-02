#!/usr/local/bin/perl
=head
This is the longer version of identify_mark_genes_for_circos.pl
=cut
use strict;
use FindBin;
use lib "$FindBin::Bin/";
use ABCC_PARSE;
use Data::Dumper;
use Getopt::Std;
use vars qw($opt_i $opt_t $opt_z $opt_o $opt_g $opt_b $opt_l $opt_h $opt_e $opt_k $opt_a $opt_x $opt_z $opt_L);#opt_k is to run silently (no STDIN)
getopts("i:o:t:g:b:l:h:xzeakzL");
 umask(0000);
#add this variable for opt_t and every tag "#WHATTODO
my $bin=`dirname $0`;chomp $bin;
my $usage=qq(
	$0 
	[REQUIRED]
		-i <annovar input file>
		-t <s|b|c|a> ... what to do
			(a) - generates all the gene lists below 
			(s) - generate the gene list for all variations with sift scores <=0.05
			(c) - generates the gene list for variations where there are no hits to the first database, or vice versa
				 e.g. Compare CGI and JPT cancer lists
				 use -h to specify the population database (header)
			(b) - generate the gene list for all variations in a nonB motif
			(d) - generates damaging lists for sift and polyphen
	[OPTIONAL]
		-g <gene mappings fn>
			formatted like:
			chr -> start -> stop -> gene name
		-l <gene list>
			gene list, one per line, only these genes will be included in the data
		-b <bin size>
		-h <header name>
		-e < exons only>
		-a <if specified, then all features are included, not just genic>
		-x <histogram flag; check compatibility>
		-o <output filename>
		-z <if specified, then the specified genelist is highlighted instead of filtered>
		-L <for large datasets>
);
my ($color,$nonb);my $ensembl=0;
my $filename=$opt_i;
if (!-e "$filename" || -z "$filename") {
	die "You must supply a variations list($filename does not exist)\n$usage";
}
$opt_b ||= 10_000;
my %genelist;
die "$usage" if (!defined $opt_t);
die "You must specify -p option with -t option to determine the population database you want to use\n$usage\n" if ($opt_t=~/[oc]/i && !defined $opt_h);
#find the corresponding information based on input's header
my @headers=split("\t",`head -n1 $filename`) ;chomp $headers[$#headers];my %header;
print "Found". join ("|",@headers)."\n";
for (my $i=0;$i<=$#headers;$i++){#WHATTODO
	if ($headers[$i]=~/sift/i){
		$header{'s'}=$i;$header{'c'}=$i if ($opt_t=~/c/i && $opt_h=~/sift/i);
	}elsif ($headers[$i]=~/nonb/i){
		$header{'b'}=$i;$header{'c'}=$i if ($opt_t=~/c/i && $opt_h=~/nonb/i);
	}elsif ($headers[$i]=~/^(Ensembl_){0,1}Gene$/i){
		$header{'gene'}=$i;
	}elsif ($opt_t=~/(c|d)/ && defined $opt_h && $headers[$i]=~/$opt_h/i){
		my $tmp=$opt_t;
		if (exists $header{$tmp}){warn "You have multiple databases with '$opt_h' in the header...using last column found\n";}
		$header{$tmp}=$i;
	}elsif ($opt_t=~/o/ && defined $opt_h && $headers[$i]=~/$opt_h/i){
		if (exists $header{'o'}){warn "You have multiple databases with '$opt_h' in the header...using last column found\n";}
		$header{'o'}=$i;
	}elsif ($headers[$i]=~/(Rel pos to feature|Annot Feat)/i){
		$header{'impact'}=$i;
	}elsif ($headers[$i]=~/ANNOVAR annot/i){
		$header{'annot'}=$i;
	}
}
die "There were missing headers in your input file\n" .Dumper (\%header) if ( ($opt_t=~/c/i && !exists $header{'c'}) || (!exists $header{'s'} && $opt_t=~/s/i) || (!exists $header{'b'} && $opt_t=~/b/i) || !exists $header{'gene'} || !exists $header{'annot'} || !exists $header{'impact'});#WHATTODO
my %annot;my $target_population_file;
print Dumper (\%header);
my $parseBinaddon='';
if($opt_x){
	$parseBinaddon="-x ";
}
my $whattodo=$opt_t;
if ($whattodo=~/s/i || $opt_h=~/(avsift)/){
	$opt_h="avsift";
	$parseBinaddon.= "-c '<0.05' ";
}elsif ($opt_h=~/ljb2\d+_pp2/i){
	$parseBinaddon.="-c '<=0.5' ";#benign all others are probably damaging or worse.
}elsif ($whattodo=~/s/i || $opt_h=~/(ljb_(mt|sift))/){
	$opt_h=$1;
	$parseBinaddon.="-c '<0.85' ";
}elsif ($whattodo=~/b/i){
	$opt_h="nonB";
	$parseBinaddon.= "-k 'G_Quadruplex' ";
}elsif (!$opt_h){
	die "LINE ONE HUNDRED\n";
} 
if ($opt_l){
	$parseBinaddon.=" -l $opt_l ";
	$parseBinaddon.=" -z " if ($opt_z);
}
my $outputfilename=(defined $opt_o)?$opt_o:"$whattodo.txt";
$outputfilename.=".txt" if ($outputfilename!~/.txt$/);
next if ($whattodo=~/[ka]/i);#WHATTODO
my $grep_addon="  grep -ve '^-' -ve '#' -ve 'ERR' ";
my $ids;
if ($whattodo=~/c/i){
	$grep_addon= " grep -e '^-'  ";#for subtraction studies you want to get the variations not in the database population
}
$ids=join (" \"\\t\" \$",($header{$whattodo}+1),($header{'annot'}+1),($header{'gene'}+1),($header{'impact'}+1));#add one because unix starts at 1 and perl starts at 0
if ($opt_a && ($whattodo=~/s/ || $opt_h=~/(avsift|ljb2\d{0,1}_mt|ljb2\d{0,1}_pp|sift)/)){##REINSTATE if using $opt_a below ###CHILICOOKOFF
	$grep_addon= "  grep -e splic -e exonic -e intronic -e $opt_h -i| grep -ve ncRNA -i ";
	$ids=join (",",($header{$whattodo}+1),($header{'annot'}+1) ."-".($header{'annot'}+7));
}elsif ($opt_a){
	$ids=join (",",($header{$whattodo}+1),($header{'annot'}+1) ."-".($header{'annot'}+7));
}elsif ($opt_t=~/d/i){
	if ($opt_h=~/s/){
		$ids=substr($opt_o,0,1);
	}else{
		$ids='d';
	}
	
	$ids=($header{$ids}+1);
	system  ("head -n1 $filename| cut -f1-".@headers."> $opt_o.txt;cut -f $ids $filename >1$opt_o.txt;cut -f1-".@headers." $filename >2$opt_o.txt;paste 1$opt_o.txt 2$opt_o.txt | grep -P '^(Probably ){0,1}Damaging' -i | cut -f 2-100>>$opt_o.txt \n");
	 system  ("rm 1$opt_o.txt 2$opt_o.txt\n"); 
	
	# print ("cut -f1-3," . ($header{$ids}+1)." $filename  | grep -P '^(Probably ){0,1}Damaging' -i  > $opt_o.txt");
	exit;
}
# die "Could not parse out data from headers ($ids)\n" if ($ids!~/\s/);chomp $ids;
my $ans;
if ( -e "$outputfilename" && ! -z "$outputfilename" && ! $opt_k){
	until ($ans=~/[ynx]/i){
		print "$outputfilename exists.  Do you wish to re-run and overwrite?? (y/n or 'X' to exit)\n";
		$ans=<STDIN>;
	}
}else{
	$ans="yes";
}
if ($ans=~/y/i){
	print "Generating parsing file for ($whattodo)\nawk -F'\\t' -v OFS=\"\\t\" '{print \$$ids}' $filename |$grep_addon > $outputfilename\n";
	# system ("cut -f $ids $filename |$grep_addon > $outputfilename\n");
	system ("awk -F'\\t' -v OFS=\"\\t\" '{print \$$ids}' $filename |$grep_addon > $outputfilename\n");
}elsif ($ans=~/n/i){#do nothing; reuse
}else{
	exit;
}
if ($opt_a){#places all of the same type into bins (all included) ##CHILICOOKOFF
	my $outputname=`basename $opt_i`;chomp $outputname;
	if ($outputname=~/(.txt)(.annovar_wrpr.output){0,1}/){
		$outputname=~s/(.txt)(.annovar_wrpr.output){0,1}/.$opt_h.txt.tmp/g;
	}elsif($outputname=~/.out$/){
		$outputname=~s/(.out)/.$opt_h.txt.tmp/g;
	}
	if ($opt_i eq "$outputname" || $outputname=~/^\./){
		$outputname=~s/.annovar_wrpr.output/.$opt_h.txt.tmp/g;
		if  ($opt_i eq "$outputname" || $outputname=~/^\./){
			$outputname.="_1";
		}
	}
	_system ("cut -f1-5 $outputfilename >.$outputname.1;cut -f6-12 $outputfilename>.$outputname.2");
	_system ("paste .$outputname.2 .$outputname.1> $outputname");
	_system ("rm $outputfilename .$outputname* -f");
	if ($opt_L){
		_system ("perl $bin/parseBins.pl -i $outputname -H -p chr");
	}else{
		_system ("perl $bin/parseGenesByBin.pl -i $outputname $parseBinaddon -o $outputfilename");
	}
}else{
	open (FILE,"<$outputfilename") or die "Cannot open the generated file $outputfilename\n";
	print "reading out($outputfilename)\n";
	LINE: while (my $line=<FILE>){
		chomp $line;
		next if ($line=~/ERR/);
		my $curr=$line;
		my @arr=split("\t",$line);
		my @gene;chop $arr[2] if ($arr[2]=~/[;:]$/);#do not change $arr[NUM]; relative to grep
		if ($arr[2]=~/^([\d\w\-]*)\(/){#do not change $arr[NUM]; relative to grep
			push (@gene,$1);
		}else{
			@gene=split(/[,:;]/,$arr[2]);
		}
		if ($whattodo=~/[s]/i && !$opt_a){
			my ($pred,$score,$stderr)=split(/[:(]/,$arr[0]);
			if ( $score && $score<=0.05){#do not change $arr[NUM]; relative to grep
				my $lastgene;#do not change $arr[NUM]; relative to grep
				GENE: foreach my $gene (@gene){
					$ensembl=1 if ($gene=~/ENSG\d+/);
					if ($lastgene eq '' || $lastgene ne "$gene"){
						$genelist{$gene}{$whattodo}++;
					}
					$lastgene=$gene;
				}
				print "matches cutoff in $0 $score\n";
			}else{
				print "does not match cutoff($arr[0])\n";
			}
		}elsif ($whattodo=~/[b|c|s]/i){
			#split by type
			##LEAVE $ARR[NUM] alone as now they are relative to what we grep'd
			#do not change $arr[NUM]; relative to grep
			my @id=split(";",$arr[0]);#Mirror_Repeat;chr1_1777264_1777310_MR;Direct_Repeat;chr1_1777264_1777311_DR;Short_Tandem_Repeat;chr1_1777264_1777311_STR;Z_DNA_Motif;chr1_1777264_1777313_ZDNA 
			# print "working on $arr[1]? Should be feat\n";<STDIN>;
			if ($arr[1]=~/(intergenic|downstream)/i){#do not change $arr[NUM]; relative to grep
				next LINE;
			}elsif ($arr[1]=~/(utr|splic|exon|intron|upstream)/i){#do not change $arr[NUM]; relative to grep
			}elsif($line=~/$opt_h/){#header
			}else{
				die "what's this? $arr[1]($line)\n";
			}
			GENE: foreach my $gene (@gene){
					$ensembl=1 if ($gene=~/ENSG\d+/);
					if ($gene=~/\(/ && $arr[1]!~/splic/){
						die "Don't know what to do with this line $curr\n";
					}elsif ($gene=~/^(\w+)\(/ && $arr[1]=~/splic/){
						$gene=$1;
					}
				$genelist{$gene}{$whattodo}++;
				# print "accepting $gene in my $whattodo genelist";
				for (my $k=0;$k<=$#id;$k+=2){
					if (($id[$k]=~/quad/i && $whattodo=~/b/i) || ($whattodo=~/c/i && $id[$k] ne '-')) {print "skipping $id[$k]\n" ;next;}
					if ($arr[1]=~/exonic/ && $arr[1]!~/ncRNA/i){
						#evaluate whether it is syn or nonsyn
						if ($arr[3]=~/p.(.*)?:/ && $1=~/(fs|del|ins|frameshift|X)/i){
							$arr[1]="fs_exonic";#print "working on $1\n";
						}else{
							my $partial_annot=$arr[3];
							my $satisfied=0;
							until ($partial_annot!~/(:p\..*?:)/ || $satisfied){
								my $curr_annot=$1;$partial_annot=$';
								#print "working on $curr_annot\n";
								if ($curr_annot=~/:p\.(\w)\d+(\w):/ && $1 ne $2){
									$satisfied=1;$arr[1]="ns_exonic";
								}elsif ($curr_annot=~/:p\.(\w)\d+(\w):/ && $1 eq $2){
									$arr[1]="syn_exonic";
								}else{
									die "can't parse this($curr_annot)$arr[3]($satisfied)\n";
								}
							}	
						}
					}
					# print "Adding to annothash:$gene,$id[$k],$arr[1]...\n";
					if (exists $annot{$gene}{$id[$k]}){#$annot{UTS2}{Direct_Repeat}=genefeature
						#precendence of importance
						# ns_exon/frameshift >> splice >>  utr5 >> ncrna(NR) >> syn_exon/intron/utr3/upstream >>intergenic
						if ($annot{$gene}{$id[$k]}=~/(.s_exonic|splic)/){#do not remove if it is ns coding#fs or ns exonic
						}elsif ($arr[1]=~/(.s_exonic|splic)/){#do not change $arr[1]; relative to grep
							$annot{$gene}{$id[$k]}=$arr[1];#do not change $arr[1]; relative to grep
						}elsif ($arr[1]=~/utr5/){
							$annot{$gene}{$id[$k]}=$arr[1];#do not change $arr[1]; relative to grep
						}elsif ($annot{$gene}{$id[$k]}=~/(ncRNA)/i){
							$annot{$gene}{$id[$k]}=$1;
						}elsif ($annot{$gene}{$id[$k]}!~/exon/ && $arr[1]=~/exon/){
							$annot{$gene}{$id[$k]}=$arr[1];
						}
					}else{
						$annot{$gene}{$id[$k]}=$arr[1];#do not change $arr[1]; relative to grep
					}
				}
			}
		}elsif ($whattodo=~/[o]/i){
			GENE: foreach my $gene (@gene){					
				$ensembl=1 if ($gene=~/ENSG\d+/);
				$genelist{$gene}{$whattodo}++;
			}
		}
		# print Dumper (\%genelist);die;
	}
}
exit if ($opt_a);
#find the genic coordinates in annovar gene and output for circos
#write to files
my $annot_file="/SeqIdx/circosdb/data/human/genes_labels.human.hg19.txt";
if ($ensembl){
	$annot_file="/SeqIdx/circosdb/data/human/genes_labels_ensembl.human.hg19.txt";
}
die "Cannot locate $annot_file \n" if (! -e $annot_file );
my $namer;
if ($opt_t=~/[sa]/i){#WHATTODO
	open (DAMAGING,">$outputfilename.genes.txt") or die "cannot open the output file for writing\n";
	open (LABELS1,">$outputfilename.genes_label.txt") or die "cannot open labels file for writing\n";
	$namer='avsift';
}
if ($opt_t=~/[ab]/i){
	$namer='nonb';
	open (NONB,">$outputfilename.genes.txt") or die "cannot open nonBgenes for writing\n";
	open (LABELS2,">$outputfilename.genes_label.txt") or die "cannot open labels2 file for writing\n";
}
if ($opt_t=~/(c|d|o)/i){
	$namer=$opt_h;
	open (COMPARE,">$outputfilename.genes.txt") or die "Cannot open $filename.compare.genes for writing\n";
	open (LABELS3,">$outputfilename.genes_label.txt") or die "Can't open $!\n";
}
my (%genecoords,%includelist);
if (defined $opt_g){
	print "Reading $opt_g and opt_z($opt_z):controls highlight or filter gene\n";
	readGenes($opt_g);
}
if (defined $opt_l){#this is a genelist #20120622
	print "Reading $opt_l and opt_z($opt_z):controls highlight or filter gene\n";
	readGenes($opt_l);
}
if ($opt_l || $opt_g){
	print "done reading genes\n" ;
}else{
	readGenes($annot_file);
	print "done reading $annot_file\n";
}
if (keys %{$nonb}==0 || keys %{$color}==0){($color,$nonb)=color();}
my ($type,$fill_color_base);
if ($opt_z ){
	$fill_color_base= "fill_color=dred";
}
if ($namer){
	if ( -e "circos.input"){
		$type=`grep $namer circos.input|cut -f3`;chomp $type;
	}else{
		$type="scatter";
	}
	if ($type=~/(text|tile|line)/i){
		$fill_color_base="color=dred";
	}
}else{
	die "this field ($namer) should be in your circos.input file\n";
}
# print Dumper(\%genelist,\%genecoords);exit;
GENE: foreach my $gene (keys %genelist){
	print "$gene:\n". Dumper $genelist{$gene} if ($gene!~/^[\-\w\.]+$/);
	my ($mychr,$mystart,$myend,$type);
	my $fill_color=$fill_color_base;#lkup once above
	if (exists $genecoords{uc($gene)}){#this  means it was in the genelist ; highlight in red
		($mychr,$mystart,$myend)=split(" ",$genecoords{uc($gene)});
	}elsif ($opt_z){#not in gene list but should be in the circos picture with color neutral
		$fill_color=~s/dred/dgrey/;
		($mychr,$mystart,$myend)=_getMinMax($gene);#($data[2],$data[4],$data[5]);
		next if ($mychr eq '0');
	}elsif (!defined $opt_g && !defined $opt_l){#no gene list;everything in genelist passes
		print __LINE__."($gene) not in genecoords?\n";
	}else{#means not in the gene list 
		print "not in genelist...";next GENE;
	}
	my $regex= ($opt_e)?"([A-Z]*_exonic|utr5|splicing|upstream|intron|ncRNA|utr3)":"(exonic)";
	
	if (exists $genelist{$gene}{'s'} && $genelist{$gene}{'s'}){#WHATTODO
		print DAMAGING join (" ",$mychr,$mystart,$myend,$gene,"score=$genelist{$gene}{'s'}")."\n";
		print LABELS1 join (" ",$mychr,$mystart,$myend,"$gene","$fill_color")."\n";
	}elsif (exists $genelist{$gene}{'b'} && $genelist{$gene}{'b'}){
		print NONB join (" ",$mychr,$mystart,$myend,$gene,"score=$genelist{$gene}{'b'}")."\n";
		foreach my $nonBtype (keys %{$annot{$gene}}){
			print LABELS2 join (" ",$mychr,$mystart,$myend,"$gene","$fill_color")."\n" if ($annot{$gene}{$nonBtype}=~/$regex/i);
		}
	}elsif (exists $genelist{$gene}{'c'} && $genelist{$gene}{'c'}){
		if ($annot{$gene}{'-'}=~/$regex/i){
			my $annot=lc($1);
			print COMPARE join (" ",$mychr,$mystart,$myend,"$gene","annot=$annot")."\n" ;
			if ($annot=~/exonic/){
				print LABELS3 join (" ",$mychr,$mystart,$myend,"$gene",$fill_color)."\n";
			}else{
				print $annot;die;
			}
		}else{
			print "$annot{$gene}{'-'}($gene and $opt_t) doesn't match $regex\n";
		}
	}elsif (exists $genelist{$gene}{'o'} && $genelist{$gene}{'o'}){
		my $annot;
		print COMPARE join (" ",$mychr,$mystart,$myend,"$gene","annot=$annot")."\n" ;
		print LABELS3 join (" ",$mychr,$mystart,$myend,"$gene",$fill_color)."\n";
	}else{
		die "don't know what I'm doing with this\n";
	}
}
#WHATTODO
if ($opt_t=~/[sa]/i){
	close DAMAGING;close LABELS1;
}
if ($opt_t=~/[ba]/i){
	close NONB;close LABELS2;
}
if ($opt_t=~/[c|o]/i){
	close COMPARE;close LABELS3;
}
#WHATTODO
if ($opt_t=~/k/i){
	#separate each nonB subtype into coding, nc, splice
	my $fn="$outputfilename.genes_label.txt";
	die "cannot open $outputfilename.genes_label.txt for reading...\n" if (! -e "$fn");
	print "Reading $outputfilename.genes_label.txt \n";
	foreach my $type ( keys %{$nonb}){
		my @arr=`grep '$type' $fn`;print "working on $type\n";
		if ($#arr>-1){
			open (FILE,">$filename.$type.out" ) or die "cannot open $filename.$type.out\n";#print "opening $filename.$type.out\n";
			open (LABELS,">$filename.$type.labels") or die "Cannot open $filename.$type.labels\n";#print "\topening $filename.$type.labels\n";
			LINE: foreach my $line (@arr){chomp $line;
				my @line=split(/\t/,$line);
				next if ($line[0]!~/hs[\dXYMT]{1,2}$/);
				if ($line=~/name=(\w+).*annot=(.*)/){
					my ($name,$annot)=($1,$2);
#					print "searching for |$name|\t";
#					my $foo=`grep $name -i $opt_g`;chomp $foo;print "($foo)\n";next LINE;
					next LINE if (defined $opt_l && !exists $includelist{uc($name)}  );
					my $rank=(exists $$color{$annot})?$$color{$annot}[0]:'0';
					print FILE join (" ",@line[0..2],"$rank","stroke_color=$$nonb{$type}[0];fill_color=$$nonb{$type}[0];glyph=$$nonb{$type}[1]")."\n";
					print LABELS join (" ",@line[0..2],"$name","color=$$color{annot}[1],id=$annot")."\n";
				}else{
					die "what's this line $line\n";
				}
			}
			close FILE;close LABELS;
		}
	}
}else{#binning
	# print "Not running binning for text \n" and exit (0) if ($type=~/text/i);
	_system ("perl $bin/parseGenesByBin_test.pl -I $outputfilename.genes_label.txt -o $outputfilename\n");
}
exit(0);
sub _godo{

}
sub _system{
	my $cmd=shift;
	eval {
		print "Running $cmd\n";
		system ("$cmd\n");
	};
	if ($@){
		die "Cannot perform $cmd\n";
	}
	print "success running $cmd\n";
	return;
}
sub _getMinMax{
	my $xgenename=shift;
	my $annot_fn="/SeqIdx/annovardb/humandb/hg19_refGene.txt";
	my @gene_arr=split("\n",`grep -P '\t$xgenename\t' -i $annot_fn`);
	return 0 if ($#gene_arr==-1);
	my $min= my $max=0; my $chrom;
	foreach my $line (@gene_arr){
		my @features=split("\t",$line);
		$chrom=$features[2] if (!$chrom);
		if ($min==0 || $min>$features[4]){
			$min=$features[4];
		}
		if ($max==0 || $max<$features[5]){
			$max=$features[5];
		}
	}
	$chrom=~s/chr/hs/;
	return ($chrom,$min,$max);
}
sub readGenes{
	my $rgfn=shift;
	open (GENES,"<$rgfn" )  or print "Could not open $rgfn" and return;
	open (FN,">genelistwithcoords") or die "could not open genelistwithcoords\n";
	
	GENE: while (<GENES>){
		next if ($_=~/^[\s\r]{0,}$/);#skip blanks
		my @arr=split(/\s/,$_);chomp $arr[$#arr];
		if ($#arr<1){
			$arr[0]=~s/[\r\s ]//g;
			$arr[0]=uc($arr[0]);
			my ($c,$m,$ma)=_getMinMax($arr[0]);
			next GENE if ($c eq '0');
			$genecoords{$arr[0]}="$c $m $ma";
			print FN "$genecoords{$arr[0]} $arr[0]\n";
		}else{
			$genecoords{$arr[3]}=join (" ",@arr[0..2]);
		}
	}
	close GENES;close FN;
	# unlink ("genelistwithcoords") if (-z "genelistwithcoords");
	return;
}
