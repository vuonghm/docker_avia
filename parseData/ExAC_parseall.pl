#!/usr/bin/perl -w
# This is a script to read through the ExAC v2 vcf file and parse it for
# Chrom Pos Ref Alt Gene Conseq and calucate/extract AC_Hom AllMAF NFEMAF HighestMAF 
# Written 5 Jan 2015 by Michele Mehaffey

use strict;
use Cwd;
use warnings;
use File::Path;
use File::Basename;
use File::Spec;
use Data::Dumper;
use Getopt::Long qw (:config no_ignore_case no_auto_abbrev);
use Time::localtime;
use List::Util qw[min max];

my ($indir, $DIR,@allfiles, $filename, $outfile);
$outfile = "test.parseall";

# #open dir
# $indir = cwd();




$indir=$ARGV[0];
my $target_file='ExAC.r0.3.sites.vep.vcf';
if ($#ARGV==1){
	$target_file=$ARGV[1];
}
print STDOUT "Current working directory is: $indir\n";
chdir $indir;
open $DIR, $indir;
@allfiles = <*>;
# close $DIR;
# print "looking for $target_file\n";exit;
my ($AFeas,$AFfin,$AFamr,$AFnfe,$AFsas,$AFafr, $AFOth, $max, $i);
my $count = 0;
foreach $filename (@allfiles) {
    if ($filename =~ /^$target_file$/) {
	print "Processing ExAC vcf for all non-zero variants.\n";
	open (OUT, ">>$outfile") || die "Couldn't open output file: $outfile\n";
	open (FIL, "$filename") || die "Couldn't open input file: $filename\n";
	# open (DEL,">dels") || die "Couldn't open output file dels\n";
	print (OUT "#Chr\tStart\tStop\tRef\tAlt" . join (":", "AC_Hom","ExAC_All_MAF","ExAC_NFE_MAF","ExAC_Highest_MAF","ExAC_Ethnicity"). "\n");
	# print "Just printed first line, check it now.\n";
	# sleep (10);
	while (<FIL>) {
	    chomp $_;
	    $count++;
	    if ($_ !~ /^#/) {
		my @fields = split(/\t/, $_);
		my %freq;
		if ($fields[7] !~ /^AC/) {
		    print "$_\n$count\n";
		    exit;
		}
		my $test = scalar split(/[=;]/, $fields[7]);
		if ($test % 2 == 1) {
		    $fields[7] = $fields[7].";dummyfix";
		    %freq = split(/[=;]/, $fields[7]);
		} else {%freq = split(/[=;]/, $fields[7])};
		my @alts = split(/,/, $fields[4]);
		my $altct = scalar @alts;
		my %maxeth;
		my @AC = split(/,/,$freq{AC_Adj});
		my @AC_HOM = split(/,/,$freq{AC_Hom});
		my @AC_EAS = split(/,/,$freq{AC_EAS});
		my @AC_FIN = split(/,/,$freq{AC_FIN});
		my @AC_AMR = split(/,/,$freq{AC_AMR});
		my @AC_NFE = split(/,/,$freq{AC_NFE});
		my @AC_SAS = split(/,/,$freq{AC_SAS});
		my @AC_AFR = split(/,/,$freq{AC_AFR});
		my @AC_OTH = split(/,/,$freq{AC_OTH});

		for ($i=0;$i<$altct;$i++) {
#		    print STDOUT "Original field for AN_EAS is: $freq{AN_EAS},AN_FIN=$freq{AN_FIN},AN_AMR=$freq{AN_AMR},AN_NFE=$freq{AN_NFE},AN_SAS=$freq{AN_SAS},AN_AFR=$freq{AN_AFR}\n";
		    if ($freq{AN_EAS} != 0) {
			$maxeth{AFeas} = $AC_EAS[$i]/$freq{AN_EAS};
		    }else {$maxeth{AFeas} = 0};
		    if ($freq{AN_FIN} != 0) {
			$maxeth{AFfin} = $AC_FIN[$i]/$freq{AN_FIN};
		    }else {$maxeth{AFfin} = 0};
		    if ($freq{AN_AMR} != 0) {
			$maxeth{AFamr} = $AC_AMR[$i]/$freq{AN_AMR};
		    }else {$maxeth{AFamr} = 0};
		    if ($freq{AN_NFE} != 0) {
			$maxeth{AFnfe} = $AC_NFE[$i]/$freq{AN_NFE};
		    }else {$maxeth{AFnfe} = 0};
		    if ($freq{AN_SAS} != 0) {
			$maxeth{AFsas} = $AC_SAS[$i]/$freq{AN_SAS};
		    }else {$maxeth{AFsas} = 0};
		    if ($freq{AN_AFR} != 0) {
			$maxeth{AFafr} = $AC_AFR[$i]/$freq{AN_AFR};
		    }else{$maxeth{AFafr} = 0};
                    if ($freq{AN_OTH} != 0) {
                        $maxeth{AFoth} = $AC_OTH[$i]/$freq{AN_OTH};
                    }else{$maxeth{AFoth} = 0};
		    $max =  max($maxeth{AFeas},$maxeth{AFfin},$maxeth{AFamr},$maxeth{AFnfe},$maxeth{AFsas},$maxeth{AFafr},$maxeth{AFoth});
		    my %rhash = reverse %maxeth;
#		my $eth = grep {$max ~~ $maxeth{$_}}keys %maxeth;
		    my $eth = $rhash{$max};
		    if (($freq{AN_Adj}>0) &&($AC[$i]/$freq{AN_Adj})>0) {
				my $format = sprintf("%d:%.8f:%.8f:%.8f:%s", $AC_HOM[$i],$AC[$i]/$freq{AN_Adj},$maxeth{AFnfe},$max,$eth);
				my ($length1,$length2)=(length($fields[3]),length($alts[$i]));
				if ($length1>1 && $length2>1){##block substitution  ## these can be multinucleotide repeats; we don't care, we will treat them with block substitutions
					if ($length1==$length2){#1       138829  138829  GC      TC due to multi-alleles
						my $k=1;
						until (substr($fields[3],-($length1-$k)) eq substr($alts[$i],-($length1-$k)) || $length1-$k==0){
							$k++;
						}
						if ($length1-$k==0){
							print  "$fields[0]\t$fields[1]\t$fields[1]\t$fields[3]\t$alts[$i]\t$format\n";
							print "??UNCHECKED ERROR".__LINE__."??$fields[0]\t$fields[1]\t".($fields[1]+$length1)."\t$fields[3]\t$alts[$i]\t$format\n";
						}else{
							print  OUT join ("\t", $fields[0],$fields[1]+$k,($fields[1]+$k),substr($fields[3],0,$k),substr($alts[$i],0,$k),$format)."\n";
						}
					}elsif ($length1>$length2){
						if ($fields[3]=~/$alts[$i]/){
							my $idx=index($fields[3],$alts[$i]);
							print OUT join ("\t", $fields[0],$fields[1]+$length2,($fields[1]+$length2),substr($fields[3],$length2,($length1-$length2)),"-",$format)."\n";
						}elsif($alts[$i]=~/$fields[3]/){
							print "??UNCHECKED ERROR".__LINE__."??$length1,$length2)$fields[0]\t$fields[1]\t$fields[1]\t$fields[3]\t$alts[$i]\t$format\n";
						}else{
							#(13,10)1        880453  880453  GTCCTCCTCGCCC   GTCCTCGCCC 
							print OUT "$fields[0]\t".($fields[1])."\t". ($fields[1]+$length1-1)."\t$fields[3]\t$alts[$i]\t$format\n";
						}
					}elsif ($length1<$length2){
						if ($alts[$i]=~/$fields[3]/){
							my $idx=index($alts[$i],$fields[3]);
							print OUT join ("\t", $fields[0],$fields[1]+$length1,($fields[1]+$length1),"-",substr($alts[$i],$length1,($length2-$length1)+1),$format)."\n";
						}else{
							
							print OUT  "$fields[0]\t".($fields[1])."\t". ($fields[1]+$length1-1)."\t$fields[3]\t$alts[$i]\t$format\n";
						}
					}else{
						print "??UNCHECKED ERROR".__LINE__."??$fields[0]\t$fields[1]\t$fields[1]\t$fields[3]\t$alts[$i]\t$format\n";
					}
					
				}elsif (length($fields[3])>1 ){
					print  OUT join ("\t", $fields[0], ($fields[1]+1), ($fields[1]+length($fields[3])-1), substr($fields[3],1,length($fields[3])),'-',$format)."\n"; 
				}elsif ( length($alts[$i])>1){
					print  OUT "$fields[0]\t". ($fields[1]+1). "\t" . ($fields[1]+length($alts[$i])-1). "\t-\t". substr($alts[$i],1,length($alts[$i])). "\t$format\n";
				}else{
					print OUT "$fields[0]\t$fields[1]\t$fields[1]\t$fields[3]\t$alts[$i]\t$format\n";
			    }
			}
		}
	    }
	}
    }
}
close OUT;
close FIL;
