#!/usr/local/bin/perl
=head 
Takes in a genelist and reports values for NetPath
=cut
use strict;
use FindBin;
use Getopt::Std;
use Data::Dumper;
umask(0000);
use vars qw( $opt_g $opt_i $opt_m $opt_o $opt_h $opt_s);
getopts("g:i:m:o:hs");
##BioGrid setup
my $access_key="beb4cd8abbbf56464b21a3baff2d1f6d";#access key for BioGRID
my $url="http://webservice.thebiogrid.org/interactions?searchNames=true&geneList=FOOBAR&includeInteractors=true&includeInteractorInteractions=false&taxId=9606&accesskey=$access_key&includeHeader=true";
##script setup
my $bin=`dirname $0`;chomp $bin;
die "You must specify a gene list!\nUSAGE: $0 -g <genelist> \n\t[-h <if specified, prints header>]\n\t[-o <output filename>]\n\t[-s <if specified, then separate files will be made for each gene>\n\n" if (!defined $opt_g);
my $dir=($opt_g=~/^\//)?`dirname $opt_g`:'./';
chomp $dir;
$opt_o||="$dir/$opt_g.paths";
if ($opt_o!~/^\//){
	$opt_o=$dir.'/'.$opt_o;
}
my $NetPathMapperFile=($opt_m)?$opt_m:"$bin/NetPath/mapping2";
my $NetPathInfoFile=($opt_i)?$opt_i:"$bin/NetPath/NetPath_Gene_regulation_all.txt";
die "One of your reference files does not exist ($NetPathInfoFile or $NetPathMapperFile)\n" if (!-e $NetPathInfoFile || !-e $NetPathMapperFile);
my %pathway;my %mapper;
open (MAPPER,"<$NetPathMapperFile") or die "Cannot open $NetPathMapperFile\n";
while (my $line=<MAPPER>){
	chomp $line;
	my @arr=split("\t",$line);
	$mapper{$arr[0]}=$arr[1];
}
close MAPPER;
open (FILEOUT,">$opt_o") or die "cannot open $opt_o\n";
open (PATHWAY,"<$NetPathInfoFile") or die "Cannot open $NetPathInfoFile\n";
my @report_headers=('pathway_fullname','regulation','expt','pathway_id','pubmed_id','entrez_id');
my @headers=('gene_reg_id','pathway_name','pathway_id','gene_name','entrez_id','regulation','expt','pubmed_id');
while (my $line=<PATHWAY>){
	chomp $line;
	next if ($line=~/Gene regulation/);#header
	my ($gene_reg_id,$pathway_name,$pathway_id,$gene_name,$entrez_id,$regulation,$expt,$pubmed_id)=split("\t",$line);
	my @path=split("\t",$line);
	for (my $i=0;$i<=$#headers;$i++){
		$pathway{$path[3]}{$headers[$i]}=$path[$i];
	}
	$pathway{$path[3]}{'pathway_fullname'}=$mapper{$path[2]};
}
# print Dumper (\%pathway);
close PATHWAY;
open (GENELIST,"<$opt_g") or die "Cannot open $opt_g\n";
print FILEOUT "Source: NetPath Resource for Signal Transduction Pathway and Expression (http://netpath.org); AVIA download date Jan 2014\n". join ("\t","Gene",@report_headers)."\n" if ($opt_h);
my @genelist=<GENELIST>;
close GENELIST;
#get in batches so we don't overwhelm the server
my $batchsize=50;
for (my $i=0;$i<=$#genelist;$i+=$batchsize){
	my $end=(($i+$batchsize)<$#genelist)?($i+$batchsize-1):$#genelist;
	my $geneList=join("|",@genelist[$i..$end]);
	$geneList=~s/\n//g;
	my $getInterxn_cmd ="wget 'http://webservice.thebiogrid.org/interactions?searchNames=true&geneList=$geneList&includeInteractors=true&includeInteractorInteractions=false&taxId=9606&accesskey=$access_key&includeHeader=true'";
	my $outname="BioGrid_All_".sprintf("%03d",($i/$batchsize)).".interxn.txt";
	# system ( "$getInterxn_cmd -O $outname >>BioGrid.log\n");//No longer using this method anymore
	
	for (my $j=$i;$j<=$end;$j++){
		$genelist[$j]=~s/[\s\n\r]//g;
		my $gene=$genelist[$j];
		# system ("grep -e 'Gene' -e '\\b$genelist[$j]\\b' $outname > BioGrid_$genelist[$j].interxn.txt\n") if (!$opt_s);
		if (exists $pathway{$gene}){
			print FILEOUT "$gene\t";
			for (my $i=0;$i<=$#report_headers;$i++){
				print FILEOUT "$pathway{$gene}{$report_headers[$i]}\t";
			}
			print FILEOUT "\n";
		}
	}
	# system ("unlink $outname\n") if (!$opt_s);
	# system ("$getInterxn_cmd -O BioGrid_All.interxn.txt 2>>BioGrid.log \n");
}
close FILEOUT;
# exit;
# foreach my $line (@genelist){
# 	$line=~/^(\S+)/;
# 	my $gene=$1;
# 	if (exists $pathway{$gene}){
# 		print FILEOUT "$gene\t";
# 		for (my $i=0;$i<=$#report_headers;$i++){
# 			print FILEOUT "$pathway{$gene}{$report_headers[$i]}\t";
# 		}
# 		print FILEOUT "\n";
# 	}
	# my $getInterxn_cmd ="wget 'http://webservice.thebiogrid.org/interactions?searchNames=true&geneList=$gene&includeInteractors=true&includeInteractorInteractions=false&taxId=9606&accesskey=$access_key&includeHeader=true'";
	# system ( "echo $getInterxn_cmd -O BioGrid_$gene.interxn.txt". ">>BioGrid.log\n");
	# system ("$getInterxn_cmd -O BioGrid_$gene.interxn.txt 2>>BioGrid.log \n");
# }
