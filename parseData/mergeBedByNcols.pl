#!/usr/bin/perl 

use strict;
use Getopt::Std;
use Data::Dumper;

##############################Set up vars and options#################################
umask(0000);
#use vars qw(opt_f opt_s $opt_o $opt_r);
use vars qw($opt_f $opt_s $opt_o $opt_r $opt_t $opt_u $opt_n $opt_c $opt_a $opt_i $opt_h $opt_F $opt_A $opt_k $opt_N);
getopts("f:s:o:n:c:h:i:rtuaFAkN");
my $usage=qq(
	 $0 
	 -f <first file> #this should be the smaller of the two files
	 -s <second file>  ## this should be the main file, all rows are reported from this file!!

	OPTIONAL:
	[-c <column to insert new data> #default 6]  # no zero start!!
	[-n <number of columns to compare with> DEFAULT: 5] 
	[-u <if specified, do not ignore case> DEFAULT: ignore case]
	[-o <outputfile>  #default STDOUT]
	[-h <look for this column header and insert after this column header> ]
		-cannot be empty
		-assumes header is in the first row of the main file
	[-r <record separator> DEFAULT = "\\n";]
	[-i <what to insert into the rows where the key does not exist in the first file>]
		-DEFAULT : '-'
		-should use quotes around this field if specified
	[-a <append to col (-c) instead of inserting new column>]
		-for VCF files where the info column is always populated
	[-F/A <if specified, then outputs in FDI/AVIA format>]
	[-N <force no header>]
 );
my %primary;
die "You must specify the following: \n$usage" if (! (defined $opt_f || defined ($opt_s) ));
$opt_n=(defined $opt_n)?($opt_n-1):4;# offset off by one 0based
$opt_c=(defined $opt_c)?($opt_c-1):5;
$opt_i||='-';
die "$opt_f or $opt_s does not exist\n" if ( ! -e $opt_f || ! -e $opt_s);
my $f_fn=$opt_f;
my $s_fn=$opt_s;
my (%first,%second);#hash to hold differences
if (defined $opt_o){
	if ($opt_A){
		open (OUT,">$opt_o") or die "Cannot open your specified input $opt_o\n";
	}else{
		open (OUT,"|tee $opt_o") or die "Cannot open your specified input $opt_o\n";
	}
}else{
	open (OUT,">&STDERR") or die "what\n";
}
if ($opt_h){
	$opt_c=-1;
	my @header=split("\n",`head -n 1 $opt_s`);
	for (my $i=0;$i=$#header;$i++){
		if ($header[$i] eq "$opt_h"){
			$opt_c=($i+1);
		}
		last if ($opt_c>-1);
	}
	die "Could not find your header $opt_h in row one of $opt_s\n" if ($opt_c==-1);
	print "inserting in $opt_c\n";
}
readFile($f_fn,\%first);
#integrate
open (FILE, "<$opt_s")|| die "Cannot open $opt_s (opt_s) in $0\n";
my ($newdata,$delimiter);
if ($opt_F){
	$delimiter=":";
}else{
	$delimiter="\t";
}
if ($opt_A){
	my $dbname=$opt_o;$dbname=~s/(\.txt$|hg19\_|$opt_s\.)//g;
	if ($dbname!~/funseq/i){
		print OUT "#$dbname\n";
	}elsif (!$opt_N){
		my $head=`head -n1 $opt_f|cut -f5`;
		$head=~s/;/;FunSeq_/g;
		print OUT "#FunSEQ_$head";
	}
}
while (my $line=<FILE>){
	chomp $line;
	$line=uc($line) unless ($opt_u);
	my @arr=split("\t",$line);
	my $idx=join ("\t",@arr[0..$opt_n]);
	if ($opt_F && $idx=~/(^#|start)/i){next;}elsif ($idx=~/(^#|start)/i){
		#do not strip chr off the header
	}else{
		$idx=~s/^chr//i;#strip both sets of chr
	}
	$idx=uc($idx) unless ($opt_u);
	$idx=~s/ //g;
	if ($opt_F){
		$arr[0]=~s/chr//;
	}
	if (exists $first{$idx}){
		$newdata=$first{$idx};
	}else{
		# print "($idx) not found...";<STDIN>;
		if ($line=~/#/ && $opt_A){
			next;
		}elsif ($line=~/#/ && !$opt_a){
			$newdata="#$f_fn;";
		}elsif ($line=~/#/ && $opt_a){
			$newdata="";
		}elsif ($opt_a){
			$newdata=";$opt_i";
		}else{
			$newdata="$opt_i";
		}
	}
	# print "Start\n$line\n";#for debugging
	if ($opt_k){
		$newdata=~s/\t/:/g;
		print OUT join ("\t",@arr[0..$opt_n],$newdata)."\n";
	}elsif ($opt_A){
		print OUT "$newdata\n";
	}elsif ($opt_a){#note; we do not append for the FDI module
		print OUT join ("\t",@arr[0..$opt_n]);
		if ($opt_n+1<$opt_c){
			print OUT "\t". join ("\t",@arr[$opt_n+1..$opt_c-1],"$arr[$opt_c];$newdata",@arr[$opt_c+1..$#arr])."\n";
		}elsif ($line!~/^#/){
			print OUT "\t$arr[$opt_c]$newdata\t". join ("\t",@arr[$opt_c+1..$#arr])."\n";
		}else{
			print OUT "\t$arr[$opt_c]\t". join ("\t",@arr[$opt_c+1..$#arr])."\n";
		}
	}else{
		print OUT join ("$delimiter",@arr[0..$opt_c-1])."\t".join("\t",$newdata,@arr[$opt_c..$#arr])."\n";
	}
	# <STDIN>;#for debugging
}
# if ($type=~/s/i){
# 	foreach my $key (sort keys %second){
# 		if (!exists $first{$key}){
# 			print "$key\n";
# 			#push(@second_only,$key);
# 		}
# 	}
# }
close OUT;
system ("chmod 777 $opt_o \n" ) if (defined $opt_o);
exit;
sub readFile{
	my $fn=shift;
	my $hash_ref=shift;
	open (FILE, "<$fn")|| die "Cannot open $fn at LN".__LINE__." in $0\n";
	my $count=0;
	while (my $line=<FILE>){
		next if ($line eq "");
		chomp $line;
		
		my @arr=split("\t",$line);
		$arr[0]=~s/^chr//i;#strip both sets of chr
		
		my $idx=join ("\t",@arr[0..$opt_n]);
		$idx=uc($idx) unless ($opt_u);
		if (!exists $$hash_ref{$idx}){
			$$hash_ref{$idx}=join("\t",@arr[$opt_n+1..$#arr]);
		}
	}
	close FILE;
	return;
}

