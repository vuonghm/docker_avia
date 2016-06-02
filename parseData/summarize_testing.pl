#!/usr/local/perl
use strict;
use Data::Dumper;
use List::MoreUtils 'pairwise';#to add arrays of same size together by element
umask(0000);
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
# print Dumper (\%genic);<STDIN>;
open (STATS," | tee $root_name.STATS") or die "Cannot open $root_name.STATS\n";
print STATS "$header\n";
foreach my $gene (sort keys %genic){
	next if $gene=~/^Gene:/;
	shift @{$genic{$gene}};
	print STATS "$gene\t". join ("\t",@{$genic{$gene}})."\n";
}
shift @base_arr;print STATS join ("\t","Total:",@base_arr)."\n";
close STATS;
if (! -e "viz/subtractive.out" ){
	#generate subtractive.out
	mkdir ("viz") if (!-e "viz");
	if (! -e "$root_name.annovar_wrpr.output"){
		die "$root_name.annovar_wrpr.output does not exist\n";
	}
	system ("grep -e 'exonic' -e 'splic' -e 'UTR' -e intron -e '^#' -i $root_name.annovar_wrpr.output |grep -ve 'ncRNA_' -i  >> viz/subtractive.out\n");
}
my $continue=1;
if ( -e "searchTheseDBs.txt" ){
	if ( `grep siftv63 searchTheseDBs.txt | wc -l`>0){
		# print STDERR "About to run perl $bin/generate_genelists_byANNOVAR.pl -i viz/subtractive.out -t d -h siftv63 -o sift\n";<STDIN>;
		system("perl $bin/generate_genelists_byANNOVAR.pl -i viz/subtractive.out -t d -h siftv63 -o sift\n");
	}else{
		print "sift not found\n";$continue=0;
	}
	# print "looking for pp2";<STDIN>;
	my $pp2_name=`grep -P 'ljb2{0,1}_pp2(hvar){0,1}' searchTheseDBs.txt |head -n 1 `;
	print STDERR "Found $pp2_name\n";
	if ( $pp2_name ne ''){
		$pp2_name=`basename $pp2_name`;chomp$pp2_name;$pp2_name=~s/(hg\d{2}_|.txt)//g;
		# print "About to run perl $bin/generate_genelists_byANNOVAR.pl -i viz/subtractive.out -t d -h $pp2_name -o pp2\n";<STDIN>;
		system("perl $bin/generate_genelists_byANNOVAR.pl -i viz/subtractive.out -t d -h $pp2_name -o pp2\n");
	}else{
		print "pp2 could not be found\n";$continue=0;
	}
	# print "Done!\n";
	if (! -e "DD.txt" && $continue){
		my @header2=split("\t",`head -n1 $root_name.annovar_wrpr.output`);
		if ($#header2==-1){die "error at LN".__LINE__."\n";}
		my $start= my $annovar=-1;
		for (my $i=0;$i<=$#header2;$i++){
			next if ($header2[$i]!~/(Chr$|\#*ANNOVAR\sannot)/i);
			if ($header2[$i]=~/Chr/){
				$start=$i;
			}elsif ($header2[$i]=~/ANNOVAR/i){
				$annovar=$i;
			}
			last if ($start>0 && $annovar>0);
		}
		my $err;
		die "Couldn't find your index for DD.txt($start,$annovar)\n" and $err++ if ($start==-1 || $annovar==-1);
		if ($err){
			open (DD,">DD.txt") or die "cannot open DD.txt for writing\n";
			print DD join ("\t",@header2);close DD;
			system ("grep -i -P 'damaging.*damaging' viz/subtractive.out >> DD.txt\n")
		}else{
			$start++;$annovar++;#off by one 0based
			#print Header
			open (DD,">.DD.txt") or die "cannot open DD.txt for writing\n";
			print DD join ("\t",@header2);close DD;
			system ("grep -i -P 'damaging.*damaging' viz/subtractive.out >> .DD.txt\n");#double damaging as predicted by sift and polyphen
			#reorganize
			system ("cut -f$start-".($start+4)." .DD.txt >1; cut -f$annovar-". ($annovar+ 2)." .DD.txt >2;cut -f1-".($annovar-1)." .DD.txt >3;cut -f".($start+6) . "-100 .DD.txt> 4;paste 1 2 3 4 > DD.txt\n");
			system ("rm -f 1 2 3 4 .DD.txt\n");
			system ("rm -f DD.txt\n") and $continue=0 if (`wc -l DD.txt`==1);
		}
	}

	if ( -e "$root_name" && $continue==1){
		system ("head -n 1 viz/subtractive.out |sed 's/#//g' >subtractive.out;grep -ve '^#' viz/subtractive.out >>subtractive.out\n");
		
	}

	system ("perl $bin/../R_scripts/make_pie_charts.pl -f $root_name -p $renamer -k $pp2_name\n");
}
#make graphs for visualization
