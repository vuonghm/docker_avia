#!/usr/local/perl
use strict;
use Data::Dumper;
use List::MoreUtils 'pairwise';#to add arrays of same size together by element
umask(0000);
open (STDOUT, "|tee summarize.log") or die "Cannot open summarize.log\n";
open (STDERR,">&STDOUT") or die "Cannot redirect STDOUT> STDERR";
my $root_name||=$ARGV[0];
my $renamer=($#ARGV==1)?$ARGV[1]:$root_name;
my $bin=`dirname $0`;chomp $bin;
die "You must specify a root filename (before extensions) \n" if (!$root_name);
opendir ( DIR, './' ) || die "Error in opening dir ./\n";
my @base_arr = qw/0 0 0 0 0 0 0 0 /;
my %genic;
my $header;
#getting the stats from each of the annovar batches and adding them together
while( (my $filename = readdir(DIR))){
	next if ($filename !~/$root_name.*stats$/);
    open ( STATS , "<$filename") or die "Cannot open $filename for reading\n";
    while (my $line=<STATS>){
    	chomp $line;$line=~s/[\n\r]//g;
    	$header=$line and chomp $header and $header=~s/\t/,/g if ($line=~/Gene:/ && !$header);
    	next if ($line=~/Type:/);
    	my @arr=split("\t",$line);
    	@base_arr=pairwise{$a + $b } @base_arr,@arr;
    	if (exists $genic{$arr[0]}){
    		my @arr2= pairwise { $a  + $b } @{$genic{$arr[0]}}, @arr;
    		$genic{$arr[0]}=\@arr2;
	    }else{
	    	$genic{$arr[0]}=\@arr;#print "Adding $arr[0] to hash\n";#print Dumper (\%genic);
	    }
    }
    close STATS;
}
closedir(DIR);
if (!-e "$root_name.STATS"){
	open (STATS," | tee $root_name.STATS") or die "Cannot open $root_name.STATS\n";
	print STATS "$header\n";
	foreach my $gene (sort keys %genic){
		next if $gene=~/^Gene:/;
		shift @{$genic{$gene}};
		print STATS "$gene\t". join ("\t",@{$genic{$gene}})."\n";
	}
	shift @base_arr;print STATS join ("\t","Total:",@base_arr)."\n";
	close STATS;
}
my $continue=1;
if (-e "viz/subtractive.out" ){
	$continue=0;
	exit;
}else{
	#generate subtractive.out
	mkdir ("viz") if (!-e "viz");
	if (! -e "$root_name.annovar_wrpr.output"){
		print "$root_name.annovar_wrpr.output does not exist\n";
	}
	print "Going to make subtractive out\n";
	system ("grep -e 'exonic' -e 'splic' -e 'UTR' -e intron -e '^#' -e 'Variant ID' -i $root_name.annot.txt |grep -ve 'ncRNA_' -i  >> viz/subtractive.out\n");
	print "[INFO] Done\n";
	if (`grep siftv63 searchTheseDBs.txt | wc -l`>0){
		print "Running: perl $bin/generate_genelists_byANNOVAR_fdi.pl -i viz/subtractive.out -t d -h siftv63 -o sift\n";
		system("perl $bin/generate_genelists_byANNOVAR_fdi.pl -i viz/subtractive.out -t d -h siftv63 -o sift\n");
	}else{
		print "sift not found\n";$continue=0;
	}
	my $pp2_name=`grep -P 'ljb\\d{0,2}_pp2(hvar){0,1}' searchTheseDBs.txt `;
	if ( $pp2_name ne ''){
		if ($pp2_name=~/(hg\d{2}.*txt)/){
			$pp2_name=$1;$pp2_name=~s/(hg\d{2}\_|.txt)//g;
		}else{
			if (`which basename`){
				$pp2_name=`basename $pp2_name`;chomp $pp2_name;
			}else{
				exit;
			}
		}
		system("perl $bin/generate_genelists_byANNOVAR_fdi.pl -i viz/subtractive.out -t d -h $pp2_name -o pp2\n");
	}else{
		print "pp2 could not be found\ngrep -P 'ljb\d{0,2}_pp2(hvar){0,1}' searchTheseDBs.txt ";$continue=0;
	}
	if (! -e "DD.txt" && $continue){
		if (-e "$root_name.annot.txt"){
			open (DD,">DD.txt") or die "cannot open DD.txt for writing\n";
			print DD `head -n1 $root_name.annot.txt`;
			system ("grep -i -P 'damaging.*damaging' $root_name.annot.txt >> DD.txt\n");
		}
	}
	# if ( -e "$root_name" && $continue==1){
	# 	system ("head -n 1 viz/subtractive.out |sed 's/#//g' >subtractive.out;grep -ve '^#' viz/subtractive.out >>subtractive.out\n");
		
	# }

	#make graphs for visualization
	if ($pp2_name){
		system ("perl $bin/../R_scripts/make_pie_charts_fdi.pl -f $root_name -p $renamer -k $pp2_name\n");
	}
}

