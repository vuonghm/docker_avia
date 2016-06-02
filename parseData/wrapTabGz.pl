#!/usr/bin/perl
use strict;
use Getopt::Std;
use Data::Dumper;
umask(0000);
use vars qw( $opt_i $opt_d $opt_D $opt_v $opt_o $opt_C);
getopts("i:D:d:o:v:C");
my $org=(defined $opt_v)?$opt_v:'hg19';
my $input=$opt_i;
my $db_dir=(defined $opt_D)?$opt_D:'/SeqIdx/tabix/humandb/';
my $dbname=(defined $opt_d)?$opt_d:'cadd';
my $output=(defined $opt_o)?$opt_o:"$opt_i.$org\_$dbname";
my $usage=qq(
	$0 <input file> <database directory>
);
my %dbCols=('cadd'=> {
		'col'=> 6,#1 based position
		'dbformat'=>"$org\_chrMOO_$dbname.txt.gz"
	}
);
if (!$input || !-e $input){
	die "Please specify an input file\n";
}else{
	open (INPUT,"<$input") or die "Cannot open input file\n";
}
if (!$db_dir || !-e $db_dir){
	die "Please specify an database directory \n";
}
my $tabx_cmd="tabix";
my ($ip,$ref,$var,$chr,$addon,$run,$db);
if ($opt_o){
	open (LOG,">$output" ) or die "Cannot open $opt_o\n";
}else{
	open (LOG,">&STDOUT");
}
my $count=0;my $printout='';
while (my $line=<INPUT>){
	$count++;
	chomp $line;next if ($line=~/^#/);
	# print STDERR "working on $line\n";
	$run=1;
	if ($line=~/^(chr)*([XYMT\d]+)\t(\d+)\t(\d+)\t([ATCG\-\.]+)\t([ATCG\.\-]+)/){#bed file with alleles SNPs only
		$chr=$2;
		$ip="$2:$3-$4";$ref=$5;$var=$6;
	}elsif ($line=~/^(chr)*(([XYMT\d]+):\d+-\d+)\s*/){ #UCSC format
		$ip=$2;$chr=$3;
		print "Found $ip and $chr\n";
	}else{
		#skipping
		# die "Skipping for now $line\n";
	}
	#validate inputs
	if ($ip && $ip!~/[\dXYMT]+:\d+-\d+/){
		$run=0;
	}else{
		$db=$dbCols{$dbname}{'dbformat'};
		$db=~s/MOO/$chr/;
		if (!-e "$db_dir/$db"){
			$run=0;
		}
		$addon = "| cut -f $dbCols{$dbname}{'col'} ";
		
	}
	my $res='';
	if ($var=~/^[ATCG\-\.]+/){
		if ($var=~/^[ATCG]$/){
			$addon=" | grep -P '[ATCG]\t$var\t' " . $addon;
			$res=`$tabx_cmd $db_dir/$db $ip $addon`;chomp $res;
		}else{
			$res='-';
		}
		$ip.=":$ref:$var";
	}else{
		my @res=split("\n",`$tabx_cmd $db_dir/$db $ip $addon`);
		chomp $res[$#res];
		$res=join (",",@res);
	}
	# my @res=split("\n",`$tabx_cmd $db $ip $addon`);
	if ($opt_C){
		$ip=~s/-/:/;
		$printout.= "$ip\t$res\n";
	}else{
		$printout.= "$res\n";
	}
	if ($count%10000==0){
		print LOG $printout;
		$printout='';
	}
}
print LOG $printout;
close LOG;