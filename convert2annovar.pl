#!/usr/bin/perl
=head
ABCC_v2.0 2/18/12
	-fixed bug in multiple allele calling in vcf4 files
ABCC v3.0 4/5/12
	-fixed bug in multiple allele calling in vcf4 files using sub smallest_substring
	-previous fix did not take into account the following 
		chrN	pos1	A	ATCG,C
		chrM	pos2	ATC	AGT
		output should be:
		 chrN	pos1 pos1 A C
		 chrN pos1+1 pos1+1 
		 chrM pos2+1 pos2+2 TC GT		#original output was chrM pos2 pos2+2 ATC AGT which is technically incorrect!!! 
ABCC v3.1 3/28/13
	-fixed the multiple sample allele frequency issue  bak saved in convert2annovar_bak.pl
ABCC v3.2 5/16/13 for Mocha
	-added the varscan2 and mutector formats
ABCC v3.3 8/5/13 for Xiqiang Li
	-added Ion Torrent Variant Caller
	-hotspots files
	-targetmutations
ABCC v4.0 includes protein and cdna position conversions!  8/13/13
	-notgenomic flag must be specified
	-input file must be hgvs cdna or hgvs prot
		SAMD11:c.C238T or SAMD11:p.T109T
	-output will have one ore more mutation depending on how it maps
=cut
use warnings FATAL=> 'uninitialized';
use strict;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/";
use ABCC;
use Carp();
$SIG{__DIE__} = \&Carp::confess;
umask(0000);
our $dbtype1;
our $VERSION = 			'$Revision: 496 $';
our $LAST_CHANGED_DATE =	'$LastChangedDate: 2011-11-20 17:27:36 -0800 (Sun, 20 Nov 2011) $';
our $input="$0 ".join (" ",@ARGV);
our ($verbose, $help, $man);
our ($variantfile,%VarAbbr);
our ($outfile, $format, $includeinfo, $snpqual, $snppvalue, $coverage, $maxcoverage, $chr, $chrmt, $altcov, $fraction, $species, $filterword, $confraction, $allallele, $bare, $nofail,
	$withzyg,$af,$logfile,$append,$header,$ensembl,$debug,$allowmismatch);
our ($cmdlineargs)= join (" ",@ARGV);
our %iupac = (R=>'AG', Y=>'CT', S=>'CG', W=>'AT', K=>'GT', M=>'AC', A=>'AA', C=>'CC', G=>'GG', T=>'TT', B=>'CGT', D=>'AGT', H=>'ACT', V=>'ACG', N=>'ACGT', '.'=>'-', '-'=>'-'); ### <<< FOR 5500SOLiD LifeScope ( S=>'GC' is replaced by S=>'CG')
our %iupacrev = reverse %iupac; ### <<< FOR 5500SOLiD LifeScope
my  (%gene_ori,%matches,%cdn2prot,%seqhash);#ABCC
my $sequence_href;
GetOptions('verbose'=>\$verbose, 'help|h'=>\$help, 'man'=>\$man, 'outfile=s'=>\$outfile, 'format=s'=>\$format, 'includeinfo'=>\$includeinfo,
	'snpqual=f'=>\$snpqual, 'snppvalue=f'=>\$snppvalue, 'coverage=i'=>\$coverage, 'maxcoverage=i'=>\$maxcoverage, 'chr=s'=>\$chr, 'chrmt=s'=>\$chrmt, 
	'fraction=f'=>\$fraction, 'altcov=i'=>\$altcov,'af'=>\$af, 'logfile=s'=>\$logfile, 'append'=>\$append, 'header'=>\$header, 
	'species'=>\$species, 'filter=s'=>\$filterword, 'confraction=f'=>\$confraction, 'allallele'=>\$allallele, 'withzyg'=>\$withzyg ,'bare'=>\$bare ,'ensembl'=>\$ensembl, 'debug'=>\$debug, 'allowmismatch=s'=>\$allowmismatch,'nofail'=>\$nofail) or pod2usage ();
my $dbSNPFile;
$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
@ARGV or pod2usage (-verbose=>0, -exitval=>1, -output=>\*STDOUT);
@ARGV == 1 or pod2usage ("Syntax error");

($variantfile) = @ARGV;
$chrmt ||= 'M';
$logfile||='cvrt2anvr.stderr.log';
if (defined $outfile) {
	if ($verbose){
		$verbose="|tee ";
		$verbose.=" -a " if ($append);
	}else{
		$verbose=">";
		$verbose.="> " if ($append);
	}
	if ($append){
		if (-e $outfile){$header=0;}
		open (STDOUT, "$verbose $outfile") or die "Error: cannot append to output file $outfile: $!\n";
		open (STDERR, "$verbose $logfile" ) or die "Cannot open STDERR $logfile\n";
	}else{
		open (STDOUT, "$verbose$outfile") or die "Error: cannot write to output file $outfile: $!\n";
		open (STDERR, "$verbose$logfile" ) or die "Cannot open STDERR $logfile\n";
	}
}
print STDERR "$input\n";
if (not $format ) {
	$format = 'pileup';
	print STDERR "NOTICE: the default --format argument is set as 'pileup'\n";
}
defined $snpqual and $format eq 'pileup' || $format eq 'vcf4' || pod2usage ("Error in argument: the --snpqual is supported only for the 'pileup' or 'vcf4' format");
defined $snppvalue and $format eq 'gff3-solid' || pod2usage ("Error in argument: the --snppvalue is supported only for the 'gff3-solid' format");
if (not defined $snpqual and $format eq 'pileup') {
	$snpqual = 20;
	print STDERR "NOTICE: the default --snpqual argument for pileup format is set as 20\n";
}

if (not defined $snppvalue) {
	$snppvalue = 1;		#by default, do not use any of the P-value cutoff in filtering out GFF3-SOLID files (this is differnt from handling pileup files)
}

if (not defined $coverage) {
	$coverage = 0;
}

if (defined $fraction) {
	$format eq 'pileup' or $format eq 'vcf4' or pod2usage ("Error in argument: the '--fraction' argument is supported for the pileup or vcf4 format only");
	$format eq 'vcf4' and print STDERR "NOTICE: the --fraction argument works ONLY on indels for vcf4 format\n";
	$fraction >= 0 and $fraction <=1 or pod2usage ("Error in argument: the --fraction argument must be between 0 and 1 inclusive");
} else {
	$fraction = 0;
}

if (defined $confraction) {
	$format eq 'vcf4' and print STDERR "NOTICE: the --confraction argument works ONLY on indels for vcf4 format\n";
	$confraction >= 0 and $fraction <=1 or pod2usage ("Error in argument: the --confraction argument must be between 0 and 1 inclusive");
} else {
	$confraction = 0;
}

if (defined $altcov) {
	$format eq 'pileup' or pod2usage ("Error in argument: the '--altcov' argument is supported for the '--format pileup' only");
	$altcov < $coverage or pod2usage ("Error in argument: the --altcov argument must be less than --coverage");
	$altcov > 0 or pod2usage ("Error in argument: the --altcov argument must be a positive integer");
}

if (defined $species) {
	$format eq 'gff3-solid' or pod2usage ("Error in argument: the '--species' argument is only necessary for the '--format gff3-solid'");
}

if ($allallele) {
	$format=~/(vcf4|casava|varscan2|mutector)/ or warn ("Error in argument: the '--allallele' argument is only supported for the '--format vcf4|casava|varscan2|mutector'");
}

if ($format eq 'pileup') {
	convertPileup ($variantfile);
} elsif ($format eq 'cg') {
	convertCG ($variantfile);
} elsif ($format eq 'gff3-solid') {
	convertGFF3SolidSNP ($variantfile);
} elsif ($format eq 'soap') {
	print STDERR "WARNING: the support for '--format soap' is not well developed yet and may contain bugs for indel analysis.\n";
	convertSOAP ($variantfile);
} elsif ($format eq 'maq') {
	print STDERR "WARNING: the support for '--format maq' is not well developed yet and may contain bugs.\n";
	convertMAQSNP ($variantfile);
} elsif ($format eq 'casava') {
	if (not defined $chr) {
		convertCASAVA_noChr($variantfile);
#		pod2usage ("Error in argument: please supply --chr argument for the '--format casava'");
	}else{
		convertCASAVA ($variantfile, $chr);
	}
} elsif ($format eq 'vcf4') {
	convertVCF4 ($variantfile);
} elsif ($format eq 'hgvs') {
	convertFromHGVS($variantfile);
	# convertHGVS_wGpos ($variantfile);
} elsif ($format eq 'modhgvs') {
	convertHGVS_wGpos ($variantfile,1);
} elsif ($format eq 'clcbio') {
	convertCLC ($variantfile);
}elsif ($format eq 'onePos'){
	convertOnePos($variantfile);
}elsif ($format eq 'varscan2'){#ABCC
	convertVS2($variantfile);	
}elsif ($format eq 'mutector'){#ABCC
	convertMutector($variantfile);
}elsif ($format eq 'bambino'){
	convertBambino($variantfile);
}elsif($format eq 'tvc'){#ABCC ion torrent hgvs modified}
	convertTVC($variantfile);
}elsif ($format=~/(bed|anvr)/){
	#we do this because people really don't know what this is!
	convertBED($variantfile);
}elsif ($format=~/notgenomic/){
	convertFromHGVS($variantfile);
}elsif($format=~/dbsnp/){
	$dbSNPFile="/SeqIdx/annovardb/humandb/hg19_avsnp144.txt";
	convertFromDBSNP($variantfile);
} else {
	pod2usage ("[ERR] Error in argument: the --format $format is not currently supported. Please contact ANNOVAR developer for adding the support");
}
print STDERR "Done converting\n";
close STDERR;
#
#sub convertPileup {
#	my ($variantfile) = @_;
#	my ($countline, $countvar, $counthom, $counthet, $countindel, $countsnp, $countti, $counttv) = qw/0 0 0 0 0 0 0 0/;
#	
#	if ($variantfile eq 'stdin') {
#		*VAR = *STDIN;
#	} else {
#		open (VAR, $variantfile) or die "Error: cannot read from variant file $variantfile: $!\n";
#	}
#	print STDERR "NOTICE: Column 6-9 in output are heterozygosity status, SNP quality, total reads, reads with mutation\n";
#
#	while (<VAR>) {
#		s/[\r\n]+$//;
#		$countline++;
#		my $hom = 'hom';
#		my @field = split (/\t/, $_);
#		@field >= 10 or die "Error: invalid record found in pileupfile $variantfile (at least 10 fields expected): <$_>\n";
#		my ($chr, $pos, $wt, $call, @other) = @field;
#		my ($cons_qual, $snp_quality, $readcount, $readallele) = @field[4,5,7,8];
#		$chr =~ s/^chr//;
#		$wt = uc $wt;					#wt may or may not in upper case, it depends on the input FASTA file
#		$call = uc $call;				#usually call is in upper case
#		$readallele = uc $readallele;			#lower case indicate the opposite strand
#		
#		$includeinfo or @other = ();			#unless -includeinfo is set, the other will not be printed
#		
#		$snp_quality >= $snpqual or next;		#quality of the variant call
#		$readcount >= $coverage or next;		#read count of the variant call
#		$maxcoverage and $readcount <= $maxcoverage || next;	#maximum coverage of the variant call
#		
#		if ($wt eq '*') {				#indel
#			#example:
#			#1       970271  *       +C/+C   39      106     44      5       +C      *       1       4       0       0       0
#			#1       1548977 *       */+CCG  29      29      42      3       *       +CCG    2       1       0       0       0
#			#1       1674810 *       */+C    24      24      42      6       *       +C      5       1       0       0       0
#			#1       968466  *       -CT/-CT 53      339     55      5       -CT     *       5       0       0       0       0
#			#1       1093600 *       -GAAA/* 29      29      53      3       -GAAA   *       1       2       0       0       0
#			#1       1110101 *       */-A    41      41      17      6       *       -A      5       1       0       0       0
#			#1       1215395 *       */-TC   26      26      32      4       *       -TC     3       1       0       0       0
#			my @obs = split (/\//, $call);		#example: '+AG/+AG' as homozygotes, '*/-TA' or '*/+T' as heterozygotes
#			@obs == 2 or die "Error: pileup record contains more than two alternative alleles: <$_>\n";
#			my ($end, $ref, $alt);
#			my ($indelreadcount);			#number of reads supporting the indel
#			
#			
#			if ($obs[0] eq $obs[1]) {
#				#something weird may occur in SamTools output: 22      15231121        *       */*     360     0       32      158     *       +GA     156     2       0       0       0
#				$obs[0] eq '*' and next;	
#	
#				#for deletions, SAMTOOLS represent deletion using a location before the first deleted base in the reference sequence coordinate system
#				#for example, a deletion in Samtools output is "1       109266688       *       */-CTT  1429    1429    58      43      *       -CTT    24      19      0       0       0"
#				#the correct ANNOVAR input (for rs35029887) should be "1       109266689       109266691       CTT     -       het     1429"
#				#insertions are fine without change; for example, check rs5745466 in Genome Browser; samtools report "1       76119508        *       +AT/+AT"
#				#for this insertion, ANNOVAR input file (for rs5745466) becomes "1       76119508        76119508        -       AT      hom     1601"
#
#				if ($obs[0] =~ m/^\-/) {
#					$pos++;			#add 1 to position in deletion
#				}
#				
#				$indelreadcount = calculateIndelReadCount ($obs[0], \@field);
#				$indelreadcount/$readcount >= $fraction or next;		#do not meet minimum alternative allele fraction threshold
#				defined $altcov and $indelreadcount >= $altcov || next;
#					
#				($end, $ref, $alt) = recalculateEndRefObs ($pos, $wt, $obs[0]);
#				print STDOUT join ("\t", $chr, $pos, $end, $ref, $alt, $hom, $snp_quality, $readcount, $indelreadcount, @other), "\n";
#				$counthom++;
#			} else {
#				$hom = 'het';
#				if ($obs[0] =~ m/^[\-\+]/) {
#					$obs[0] =~ m/^\-/ and $pos++;
#					($end, $ref, $alt) = recalculateEndRefObs ($pos, $wt, $obs[0]);
#					$indelreadcount = calculateIndelReadCount ($obs[0], \@field);
#					$indelreadcount/$readcount >= $fraction or next;		#do not meet minimum alternative allele fraction threshold
#					defined $altcov and $indelreadcount >= $altcov || next;
#					
#					print STDOUT join ("\t", $chr, $pos, $end, $ref, $alt, $hom, $snp_quality, $readcount, $indelreadcount, @other), "\n";
#					$counthet++;
#				}
#				if ($obs[1] =~ m/^[\-\+]/) {
#					$obs[1] =~ m/^\-/ and $pos++;
#					($end, $ref, $alt) = recalculateEndRefObs ($pos, $wt, $obs[1]);
#					$indelreadcount = calculateIndelReadCount ($obs[1], \@field);
#					$indelreadcount/$readcount >= $fraction or next;		#do not meet minimum alternative allele fraction threshold
#					defined $altcov and $indelreadcount >= $altcov || next;
#					
#					print STDOUT join ("\t", $chr, $pos, $end, $ref, $alt, $hom, $snp_quality, $readcount, $indelreadcount, @other), "\n";
#					$counthet++;
#				}
#			}
#			$countindel++;
#		} else {
#			#1       798494  G       A       36      36      58      3       AAA     bbb
#			#1       798677  T       K       33      33      52      26      ,$.,,G.GG,.,......,..G,,...     b`bbBaJIbFbZWaTNQbb_VZcbbb
#			#1       856182  G       A       42      42      50      5       AAAAA   B\bbb
#			#1       861034  A       M       48      48      49      14      ,$,.,..,cc.c.,c bbBbb`]BFbHbBB
#			#1       864289  T       K       22      22      56      6       .g,,g,  BH^_BB
#			
#			$wt eq $call and next;			#this is not a SNP
#			my $obs = $iupac{$call} or die "Error: invalid best call ($call) in <$_>\n";
#			my @obs = split (//, $obs);
#			@obs == 2 or die "Error: observed IUPAC allele $call should correspond to two nucleotide alleles: <$_>\n";
#			if ($obs[0] ne $obs[1]) {
#				$hom = 'het';
#			}
#				
#			
#			if ($obs[0] eq $wt) {			#obs[0] is guaranteed to be an alternative allele
#				@obs = @obs[1,0];
#			}
#			if ($wt eq 'A' and $obs[0] eq 'G' or $wt eq 'G' and $obs[0] eq 'A' or $wt eq 'C' and $obs[0] eq 'T' or $wt eq 'T' and $obs[0] eq 'C') {
#				unless ($wt ne $obs[0] and $wt ne $obs[1] and $obs[0] ne $obs[1]) {
#					$countti++;
#				}
#				
#			} else {
#				unless ($wt ne $obs[0] and $wt ne $obs[1] and $obs[0] ne $obs[1]) {
#					$counttv++;
#				}
#			}
#			
#			my $mutallelecount;
#			
#			if ($obs[1] eq $wt) {			#het SNP
#				if ($chr eq $chrmt) {
#					$hom = calculateAllelicFraction ($obs[0], $field[8], $readcount);
#				}
#				$mutallelecount = calculateMutAlleleCount ($obs[0], $readallele);
#				$mutallelecount/$readcount >= $fraction or next;		#do not meet minimum alternative allele fraction threshold
#				defined $altcov and $mutallelecount >= $altcov || next;
#				
#				print STDOUT join ("\t", $chr, $pos, $pos, $wt, $obs[0], $hom, $snp_quality, $readcount, $mutallelecount, @other), "\n";
#				$counthet++;
#			} elsif ($obs[1] ne $obs[0]) {		#het SNP but both differ from reference allele
#				if ($chr eq $chrmt) {
#					$hom = calculateAllelicFraction ($obs[1], $field[8], $readcount);
#				}
#				$mutallelecount = calculateMutAlleleCount ($obs[1], $readallele);
#				$mutallelecount/$readcount >= $fraction or next;		#do not meet minimum alternative allele fraction threshold
#				defined $altcov and $mutallelecount >= $altcov || next;
#				
#				print STDOUT join ("\t", $chr, $pos, $pos, $wt, $obs[1], $hom, $snp_quality, $readcount, $mutallelecount, @other), "\n";
#				if ($chr eq $chrmt) {
#					$hom = calculateAllelicFraction ($obs[0], $field[8], $readcount);
#				}
#				$mutallelecount = calculateMutAlleleCount ($obs[0], $readallele);
#				$mutallelecount/$readcount >= $fraction or next;		#do not meet minimum alternative allele fraction threshold
#				defined $altcov and $mutallelecount >= $altcov || next;
#				
#				print STDOUT join ("\t", $chr, $pos, $pos, $wt, $obs[0], $hom, $snp_quality, $readcount, $mutallelecount, @other), "\n";
#				$counthet++;
#				$counthet++;
#			} else {				#homo SNP
#				if ($chr eq $chrmt) {
#					$hom = calculateAllelicFraction ($obs[0], $field[8], $readcount);
#				}
#				$mutallelecount = calculateMutAlleleCount ($obs[0], $readallele);
#				$mutallelecount/$readcount >= $fraction or next;		#do not meet minimum alternative allele fraction threshold
#				defined $altcov and $mutallelecount >= $altcov || next;
#				
#				print STDOUT join ("\t", $chr, $pos, $pos, $wt, $obs[0], $hom, $snp_quality, $readcount, $mutallelecount, @other), "\n";
#				$counthom++;
#			}
#			$countsnp++;
#		}
#		$countvar++;
#	}
#	my $triallelic = $countsnp-$countti-$counttv;
#	print STDERR "NOTICE: Read $countline lines and wrote ${\($counthet+$counthom)} different variants at $countvar genomic positions ($countsnp SNPs and $countindel indels)\n";
#	print STDERR "NOTICE: Among ${\($counthet+$counthom)} different variants at $countvar positions, $counthet are heterozygotes, $counthom are homozygotes\n";
#	print STDERR "NOTICE: Among $countsnp SNPs, $countti are transitions, $counttv are transversions", $triallelic?", $triallelic have more than 2 alleles\n":"\n";
#}
#
sub calculateIndelReadCount {
	my ($obs, $field) = @_;
	#make sure to use upper case in the comparison, for example:
	#chr10   130133  *       */-ca   189     189     59      31      *       -ca     27      4       0       0       0
	if ($obs eq uc $field->[8]) {
		return $field->[10];
	} elsif ($obs eq uc $field->[9]) {
		return $field->[11];
	} else {
		die "Error: invalid record in pileup file (indel counts cannot be inferred): <$obs> vs <@$field>\n";
	}
}

sub calculateMutAlleleCount {
	my ($allele, $string) = @_;	#they should be already both in upper case
	$string =~ s/\^.//g;		#^ is followed by mapping quality
	$string =~ s/\$//g;
	$string =~ s/[+-]1[^\d]//g;	#1 followed by a non-number
	$string =~ s/[+-]2..//g;
	$string =~ s/[+-]3...//g;
	$string =~ s/[+-]4....//g;
	$string =~ s/[+-]5.....//g;
	$string =~ s/[+-]6......//g;
	$string =~ s/[+-]7.......//g;
	$string =~ s/[+-]8........//g;
	$string =~ s/[+-]9.........//g;
	$string =~ s/[+-]10..........//g;
	
	#make sure to use upper case letter
	my @string = split (//, uc $string);
	my $count = 0;
	for my $i (0 .. @string-1) {
		$allele eq $string[$i] and $count++;
	}
	return $count;
}

sub calculateAllelicFraction {
	my ($obs, $readbase, $readcount) = @_;
	my @readbase = split (//, $readbase);
	my $count=0;
	for my $i (0 .. @readbase-1) {
		uc $obs eq uc $readbase[$i] and $count++;
	}
	my $hom = $count/$readcount;
	length ($hom) > 5 and $hom > 0.001 and $hom = sprintf ("%.3f", $hom);
	return $hom;
}

sub recalculateEndRefObs {		#recalculate end position, reference allele and observed allele
	my ($end, $ref, $obs) = @_;
	if ($obs =~ m/^\-(\w+)/) {	#deletion
		$end += (length ($1)-1);
		$ref = $1;
		$obs = '-';
	} elsif ($obs =~ m/^\+(\w+)/) {	#insertion
		$ref = '-';
		$obs = $1;
	} else {
		die "Error: cannot handle $end, $ref, $obs\n";
	}
	return ($end, $ref, $obs);
}

sub convertCG {
	my ($variantfile) = @_;
	
	my ($foundheader, $countline, @field);
	my ($prechr, $prestart, $preend, $prevartype, $preref, $preobs, $prescore, $prexref) = qw/0 0 0 0 0 0 0 0/;
	open (VAR, $variantfile) or die "Error: cannot read from variant file $variantfile: $!\n";
	print STDERR "NOTICE: Converting variants from $variantfile\n";
	while (<VAR>) {
		s/[\r\n]+$//;
		$countline++;
		if (m/^>locus/) {
			$foundheader++;
		}
		if (not $foundheader) {
			$countline > 50 and die "Error: invalid CG-var file format for $variantfile (>locus record is not found within the first 50 lines)\n";
			next;
		}
		my ($locus, $ploidy, $haplo, $chr, $start, $end, $vartype, $ref, $obs, $score, $haplolink, $xref) = split (/\t/, $_);
		$chr =~ s/^chr//;
		$vartype eq 'ins' or $start++;		#CG use zero-start, half-open coordinate. Insertion does not need to be processed (example, "16476   2       2       chr1    751820  751820  ins             T       49              dbsnp:rs59038458")
		$obs eq '' and $obs = '-';
		$ref eq '' and $ref = '-';

		if ($vartype =~ m/^snp|ins|del|delins|sub$/) {		#in new versions of the files, "sub" is used instead of "delins".
			#$chr eq 'M' and next;			#ignore chrM markers as they are not diploid
			if ($chr eq $prechr and $start eq $prestart and $end eq $preend and $obs eq $preobs) {		#homozygous mutation
				print $chr, "\t", $start, "\t", $end, "\t", $ref, "\t", $obs, "\t", $vartype, "\t", ($score+$prescore)/2, "\t", "hom\t", $xref, "\n";
				($prechr, $prestart, $preend, $prevartype, $preref, $preobs, $prescore, $prexref) = qw/0 0 0 0 0 0 0 0/;
			} else {
				if ($prestart and $preend) {
					print $prechr, "\t", $prestart, "\t", $preend, "\t", $preref, "\t", $preobs, "\t", $prevartype, "\t", $prescore, "\thet\t", $prexref, "\n";
				}
				($prechr, $prestart, $preend, $prevartype, $preref, $preobs, $prescore, $prexref) = ($chr, $start, $end, $vartype, $ref, $obs, $score, $xref);
			}
		}
	}
	if ($prestart and $preend) {
		print $prechr, "\t", $prestart, "\t", $preend, "\t", $preref, "\t", $preobs, "\t", $prevartype, "\t", $prescore, "\thet\t", $prexref, "\n";
	}
	print STDERR "NOTICE: Done with $countline lines\n";
}

#sub convertGFF3SolidSNP {
#	my ($variantfile) = @_;
#	my ($countline, $countvar, $countallvar, @other) = (0, 0, 0);
#	my ($unknown_count);		#count of variants with 'unknown' variation type
#	
#	open (VAR, $variantfile) or die "Error: cannot read from variant file $variantfile: $!\n";
#	$_ = <VAR>;
#	s/[\r\n]+$//;
#	m/^##gff-version\s+3/ or die "Error: invalid first line in GFF3 file ('##gff-version 3' expected): <$_>\n";
#	$_ = <VAR>;
#	s/[\r\n]+$//;
#	(m/^##solid-gff-version/ || m/^##source-version/) or print STDERR "WARNING: problematic second line in GFF3-SOLiD file ('##solid-gff-version' or '##source-version' expected): <$_>\n"; ### <<< FOR 5500SOLiD LifeScope
#
#	print STDERR "NOTICE: Column 6-9 in output are heterozygosity status, variant score (P-value), total clipped normal coverage reads, total reads with mutated allele\n";
#	
#	while (<VAR>) {
#		s/[\r\n]+$//;
#		$countline++;
#		m/^##/ and next;		#header of comment lines
#		m/^#/ and next;			#header of results lines
#		
#		my @field = split (/\t/, $_);
#		@field == 9 or die "Error: invalid record found in $variantfile (10 fields expected): <$_>\n";
#		my ($chr, $program, $type, $pos, $end, $score, $attribute) = @field[0,1,2,3,4,5,8];		#score is the P-value for the SNP calls
#		$chr eq 'chr_name' and next;	#header line
#		
#		if ($score ne '.') {
#			$score >=0 and $score <=1 or die "Error: invalid score record found in file (0-1 range expected): <$_>\n";
#			$score <= $snppvalue or next;
#		}
#		
#		if ($species and $species eq 'human') {
#			$chr eq '23' and $chr = 'X';
#			$chr eq '24' and $chr = 'Y';
#			$chr eq '25' and $chr = 'M';
#		}
#
#		$includeinfo and @other = ($attribute);			#unless -includeinfo is set, the other will not be printed
#
#		my ($readcount, $mutallelecount) = ('.', '.');		#total coverage, coverage for mutated alleles
#		
#		if ($type eq 'unknown') {
#			#SOLiD GDD3 may have unknown variation types
#			#chr1    AB_SOLiD Small Indel Tool       unknown 3833062 3833153 1       .       .       ID=5483;len=no_call;allele-call-pos=3833062;allele-call=/CCAC;allele-pos=3833057;alleles=atccatccacccatc/aTCCATCCACCCACCCATC/NO_CALL;allele-counts=REF,2,2;tight_chrom_pos=none;loose_chrom_pos=3833058-3833069;no_nonred_reads=3;coverage_ratio=8.0000;experimental-zygosity=HEMIZYGOUS;experimental-zygosity-score=1.0000;run_names=L1_1_50_10_r,L1_1_50_15_r,L1_1_50_15_r,L1_1_50_12_r;bead_ids=1018_196_970,699_1263_465,220_513_923,2022_1532_1071;overall_qvs=4,6,2,50;no_mismatches=5,4,2,0;read_pos=27,29,31,13;from_end_pos=23,21,19,37;strands=+,+,+,+;tags=R3,F3,F3,F3;indel_sizes=-92,-112,4,4;non_indel_no_mismatches=0,0,8,0;unmatched-lengths=50,50,50,50;ave-unmatched=50.0000;anchor-match-lengths=48,49,49,49;ave-anchor-length=48.7500;read_seqs=G23223321322112233223100132013201320110011001322332,T33223321322112233223100132013201320110013021322332,T33223321322112233223100132013201320110011001322332,T31001320132013201100110013223322113030332233113032;base_qvs=;non_indel_seqs=T21322332211221121322332230321212121223322332233221,G12020202202020012001200213022002130012332310122030,G12020202202020012001000210022012110312331331122030,G22111012101031010100002002321020002202121121313021;non_indel_qvs=
#			$unknown_count++;
#			next;		#do not count this one!
#		}
#		
#		if ($program eq 'SOLiD_diBayes' or $program eq 'AB_SOLiD SNP caller') {		#SNP variants
#			#detailed description can be found at http://solidsoftwaretools.com/gf/download/docmanfileversion/208/866/DiBayes_Documentation_v1.2.pdf
#			#chr1    SOLiD_diBayes   SNP     559817  559817  0.094413        .       .       genotype=Y;reference=T;coverage=9;refAlleleCounts=5;refAlleleStarts=4;refAlleleMeanQV=23;novelAlleleCounts=2;novelAlleleStarts=2;novelAlleleMeanQV=14;diColor1=11;diColor2=33;het=1;flag= 
#			#chr1    SOLiD_diBayes   SNP     714068  714068  0.000000        .       .       genotype=M;reference=C;coverage=13;refAlleleCounts=7;refAlleleStarts=6;refAlleleMeanQV=25;novelAlleleCounts=6;novelAlleleStarts=4;novelAlleleMeanQV=22;diColor1=00;diColor2=11;het=1;flag= 
#			#chr1    SOLiD_diBayes   SNP     714835  714835  0.041579        .       .       genotype=R;reference=A;coverage=5;refAlleleCounts=3;refAlleleStarts=3;refAlleleMeanQV=18;novelAlleleCounts=2;novelAlleleStarts=2;novelAlleleMeanQV=20;diColor1=02;diColor2=20;het=1;flag= 
#
#			$pos == $end or die "Error: invalid record found in GFF3-SOLiD file: start and end discordant: <$_>\n";
#	
#			my ($wt, $call);
#			my ($hit); ### <<< FOR 5500SOLiD LifeScope
#
#			if ($attribute =~ m/ref_base=(\w)/) {
#				$wt = $1;
#			} elsif ($attribute =~ m/reference=(\w)/) {
#				$wt = $1;
#			} else {
#				die "Error: invalid record found in GFF3-SOLiD file (ref_base/reference was not found): <$_>\n";
#			}
#			
#			if ($attribute =~ m/consen_base=(\w)/) {
#				$call = $1;
#			} elsif ($attribute =~ m/genotype=(\w)/) {
#				$call = $1;
#			} elsif ($attribute =~ m/allele-call=([\w\/]+)/) { ### <<< FOR 5500SOLiD LifeScope
#			        $hit = $1;
#			        if ($hit =~ m/\//) {
#			            $call = $iupacrev{join("",sort(split(/\//,$hit)))}; 
#			        } else {
#				    $call = $hit;
#			        }
#			} else {
#				die "Error: invalid record found in GFF3-SOLiD file (consen_base was not found): <$_>\n";
#			}
#						
#			if ($attribute =~ m/coverage=(\d+)/) {
#				$readcount = $1;
#				$readcount >= $coverage or next;		#read count of the variant call
#				$maxcoverage and $readcount <= $maxcoverage || next;
#			}
#			if ($attribute =~ m/novelAlleleCounts=(\d+)/) {
#				$mutallelecount = $1;
#				$mutallelecount/$readcount >= $fraction or next;		#do not meet minimum alternative allele fraction threshold
#				defined $altcov and $mutallelecount >= $altcov || next;
#			}
#			
#			my $obs = $iupac{$call} or die "Error: invalid best call in <$_>\n";
#			my @obs = split (//, $obs);
#			@obs == 2 or die "Error: observed IUPAC allele $call should correspond to two nucleotide alleles: <$_>\n";
#			if ($obs[0] eq $wt and $obs[1] eq $wt) {
#				die "Error: reference alleles are identical to observed alleles: <$_>\n";
#			} elsif ($obs[0] eq $wt) {
#				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[1], "\t", "het\t", "$score\t$readcount\t$mutallelecount\t", join ("\t", @other), "\n";
#			} elsif ($obs[1] eq $wt) {
#				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[0], "\t", "het\t", "$score\t$readcount\t$mutallelecount\t", join ("\t", @other), "\n";
#			} elsif ($obs[1] ne $obs[0]) {
#				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[0], "\t", "het\t", "$score\t$readcount\t$mutallelecount\t", join ("\t", @other), "\n";
#				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[1], "\t", "het\t", "$score\t$readcount\t$mutallelecount\t", join ("\t", @other), "\n";
#				$countallvar++;
#			} else {
#				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[0], "\t", "hom\t", "$score\t$readcount\t$mutallelecount\t", join ("\t", @other), "\n";
#			}
#		} elsif ($program eq 'AB_CNV_PIPELINE') {	#CNV
#			if ($attribute =~ m/copynum=(\d+)/ or $attribute =~ m/copynumber=(\d+)/) {
#				if ($1 < 2) {
#					print $chr, "\t", $pos, "\t", $end, "\t", 0, "\t", '-', "\t", "unk\t", "$score\t.\t.\t", join ("\t", @other), "\n";
#				} elsif ($1 > 2) {
#					print $chr, "\t", $end, "\t", $end, "\t", '-', "\t", 0, "\t", "unk\t", "$score\t.\t.\t", join ("\t", @other), "\n";
#				}
#			} else {
#				print $chr, "\t", $end, "\t", $end, "\t", '-', "\t", 0, "\t", "unk\t", "$score\t.\t.\t", join ("\t", @other), "\n";
#			}
#		} elsif ($program eq 'AB_SOLiD Large Indel Tool') {	#CNV
#			#http://solidsoftwaretools.com/gf/download/docmanfileversion/182/780/Large_Indel_Documentation_v1.0.0.pdf
#			## [FIELDS] (1) chromosome (2) version (3) indel type (4) breakpoint start (5) breakpoint end (6) p-value (7) NA (8) NA (9) attributes
#			#chr10   AB_SOLiD Large Indel Tool       insertion       151910  151969  2.77548e-11     .       .       dev=-71;avgDev=-1.63884;zygosity=HOMOZYGOUS;nRef=0;nDev=14;refDev=0;devDev=-1.60924;refVar=0;devVar=0.0159438;beadIds=1750_720_1641,649_1680_794,503_1756_608,1726_174_1362,1208_1772_353,872_594_1604,1924_990_858,1899_961_1848,901_1226_378,323_1750_1017,1185_180_1237,1519_490_1074,1291_94_324,285_758_922,1135_95_1594,1055_218_1279,
#			#chr10   AB_SOLiD Large Indel Tool       insertion       154109  154729  2.1559e-11      .       .       dev=-66;avgDev=-1.51253;zygosity=HOMOZYGOUS;nRef=0;nDev=15;refDev=0;devDev=-1.02864;refVar=0;devVar=0.133236;beadIds=1728_1671_1739,1250_231_25,811_783_1090,1035_908_491,649_1680_794,503_1756_608,1726_174_1362,1208_1772_353,872_594_1604,1924_990_858,1899_961_1848,901_1226_378,323_1750_1017,1185_180_1237,1519_490_1074,1291_94_324,285_758_922,1135_95_1594,1055_218_1279,
#			my ($call, @call, $zygosity);
#			if ($attribute =~ m#zygosity=HEMIZYGOUS#) {
#				$zygosity = 'het';
#			} elsif ($attribute =~ m#zygosity=HOMOZYGOUS#) {
#				$zygosity = 'hom';
#			} else {
#				$zygosity = 'unk';
#			}
#			if ($type eq 'insertion') {
#				#the true boundary is unknown (start is different from end) so we cannot use "-" to represent reference allele.
#				print $chr, "\t", $pos, "\t", $end, "\t", 0, "\t", 0, "\t", $zygosity, "\t", "$score\t.\t.\t", join ("\t", @other), "\n";
#			} elsif ($type eq 'deletion') {
#				print $chr, "\t", $pos, "\t", $end, "\t", 0, "\t", '-', "\t", $zygosity, "\t", "$score\t.\t.\t", join ("\t", @other), "\n";
#			}
#		} elsif ($program eq 'AB_SOLiD Small Indel Tool') {		#small indels
#			#typical simple insertion and deletions
#			#chr1    AB_SOLiD Small Indel Tool       deletion        1352612 1352614 1       .       .       ID=1290;del_len=3;allele-call-pos=1352612;allele-call=cct/;allele-pos=1352610;alleles=cccctccat/cCCCAT;allele-counts=REF,2;tight_chrom_pos=1352612-1352614;loose_chrom_pos=1352612-1352614;no_nonred_reads=2;coverage_ratio=11.5000;experimental-zygosity=HEMIZYGOUS;experimental-zygosity-score=1.0000;run_names=L1_1_25_3_r,L1_1_25_8_r;bead_ids=1470_2000_506,822_1710_1767;overall_qvs=18,19;no_mismatches=3,3;read_pos=6,13;from_end_pos=19,12;strands=-,+;tags=R3,R3;indel_sizes=-3,-3;non_indel_no_mismatches=1,-1;unmatched-lengths=25,25;ave-unmatched=25.0000;anchor-match-lengths=24,99;ave-anchor-length=61.5000;read_seqs=G0112310001100003120031200,G0300213000011000132110021;base_qvs=;non_indel_seqs=T2120033002022200220000002,;non_indel_qvs=
#			#chr1    AB_SOLiD Small Indel Tool       insertion_site  1311162 1311162 1       .       .       ID=1249;ins_len=1;allele-call-pos=1311162;allele-call=/G;allele-pos=1311161;alleles=gaggggggg/GAGGGGGGGG/NO_CALL;allele-counts=REF,3,1;tight_chrom_pos=none;loose_chrom_pos=1311160-1311169;no_nonred_reads=3;coverage_ratio=4.6667;experimental-zygosity=HEMIZYGOUS;experimental-zygosity-score=1.0000;run_names=L1_1_25_6_r,L1_1_50_10_r,L1_1_25_2_r,L1_1_25_3_r;bead_ids=850_837_429,1160_878_181,404_1050_1881,1084_64_1343;overall_qvs=20,56,25,25;no_mismatches=3,2,2,1;read_pos=11,22,11,11;from_end_pos=14,28,14,14;strands=+,-,-,-;tags=R3,F3,F3,F3;indel_sizes=1,1,1,1;non_indel_no_mismatches=-1,1,0,1;unmatched-lengths=25,50,25,25;ave-unmatched=31.2500;anchor-match-lengths=99,49,24,24;ave-anchor-length=49.0000;read_seqs=G1020001130221020000000020,T03223323210110021000000022122030100020221222222122,T0102210000000221223301000,T0102210000000221220301000;base_qvs=;non_indel_seqs=,G21202030032202013220021321131212021000122300013132,G1331133120001221220120120,G1331133120001221220120220;non_indel_qvs=
#			
#			#sometimes, allele-call is ambiguous that requires a "block substitution" representation (although they were annotated as insertion or deletion by SOLiD, they should be treated as block substitution by ANNOVAR)
#			#sometimes, mutiple allele calls may be present at the same site
#			#chr1    AB_SOLiD Small Indel Tool       deletion        1209357 1209360 1       .       .       ID=1101;del_len=4;allele-call-pos=1209357;allele-call=ggtggg/TT;allele-pos=1209355;alleles=ggggtgggggggtt/gGTTGGGGTT/gGTGTTTTGCCTT/NO_CALL;allele-counts=REF,3,1,1;tight_chrom_pos=none;loose_chrom_pos=1209357-1209363;no_nonred_reads=4;coverage_ratio=3.0000;experimental-zygosity=HEMIZYGOUS;experimental-zygosity-score=0.9888;run_names=L1_1_25_1_r,L1_1_25_2_r,L1_1_25_4_r,L1_1_25_3_r,L1_1_25_7_r;bead_ids=1017_1024_53,1493_1896_615,1794_647_1473,307_1116_687,773_1492_1671;overall_qvs=24,24,28,24,8;no_mismatches=2,3,2,3,2;read_pos=14,9,14,9,15;from_end_pos=11,16,11,16,10;strands=-,+,-,+,+;tags=F3,R3,F3,F3,F3;indel_sizes=-4,-4,-4,-4,3;non_indel_no_mismatches=1,0,0,0,0;unmatched-lengths=25,25,25,25,25;ave-unmatched=25.0000;anchor-match-lengths=24,24,24,24,24;ave-anchor-length=24.0000;read_seqs=T2221100101000101000221100,G0001122000100000101001020,T2221100101000101000221100,T1112200010100010100112000,T1011220000111000130200001;base_qvs=;non_indel_seqs=G0312033221312111022200300,T0111113210210112100001130,G0312133221312111022200300,G0231003132222112000012020,G3121331033101113122312020;non_indel_qvs=
#			#chr1    AB_SOLiD Small Indel Tool       deletion        1209436 1209436 1       .       .       ID=1103;del_len=1;allele-call-pos=1209436;allele-call=ag/A/G;allele-pos=1209434;alleles=tgagggggtt/tGAGGGGTT/tGGGGGGTT;allele-counts=REF,1,1;tight_chrom_pos=none;loose_chrom_pos=1209436-1209441;no_nonred_reads=2;coverage_ratio=5.0000;experimental-zygosity=HEMIZYGOUS;experimental-zygosity-score=1.0000;run_names=L1_1_25_6_r,L1_1_25_2_r;bead_ids=1315_1584_2005,1706_194_437;overall_qvs=28,21;no_mismatches=0,3;read_pos=9,7;from_end_pos=16,18;strands=-,-;tags=R3,R3;indel_sizes=-1,-1;non_indel_no_mismatches=-1,0;unmatched-lengths=25,25;ave-unmatched=25.0000;anchor-match-lengths=99,24;ave-anchor-length=61.5000;read_seqs=G3001010000011001010000001,G3010100022110010111000110;base_qvs=;non_indel_seqs=,T1112003220020013202122300;non_indel_qvs=
#			#chr1    AB_SOLiD Small Indel Tool       insertion_site  1424540 1424540 1       .       .       ID=1376;ins_len=3;allele-call-pos=1424540;allele-call=tt/CCCAC;allele-pos=1424537;alleles=ttttttg/TTTCCCACTG/NO_CALL;allele-counts=REF,1,1;tight_chrom_pos=none;loose_chrom_pos=1424536-1424543;no_nonred_reads=2;coverage_ratio=11.5000;experimental-zygosity=HEMIZYGOUS;experimental-zygosity-score=1.0000;run_names=L1_1_25_7_r,L1_1_50_16_r;bead_ids=703_437_370,1089_878_1744;overall_qvs=1,9;no_mismatches=3,4;read_pos=5,35;from_end_pos=20,15;strands=-,-;tags=R3,F3;indel_sizes=3,3;non_indel_no_mismatches=2,0;unmatched-lengths=25,50;ave-unmatched=37.5000;anchor-match-lengths=24,47;ave-anchor-length=35.5000;read_seqs=G2032002200200000000000020,T30100102220312202103112130230322210121100200002100;base_qvs=;non_indel_seqs=T2121120003012303000000000,G22213300221101011121030022002222300220322213303102;non_indel_qvs=
#			my ($call, @call, $zygosity);
#			my ($refcall, $gapnonred, %temphash); ### <<< FOR 5500SOLiD LifeScope
#			if ($attribute =~ m#experimental-zygosity=HEMIZYGOUS# ||$attribute =~ m#zygosity=HEMIZYGOUS#) { ### <<< FOR 5500SOLiD LifeScope
#				$zygosity = 'het';
#			} elsif ($attribute =~ m#experimental-zygosity=HOMOZYGOUS# || $attribute =~ m#zygosity=HOMOZYGOUS#) { ### <<< FOR 5500SOLiD LifeScope
#				$zygosity = 'hom';
#			} else {
#				$zygosity = 'unk';
#			}
#			$score = '.';			#by default, score=1 in the output
#			
#			#no_nonred_reads: Number of reads with unique start positions (non-redundant reads).
#			#coverage_ratio: Clipped normal coverage/number of non-redundant reads.Clipped coverage is where the parts of the read that are within 5 bp at either end are not counted as a part of coverage.
#			if ($attribute =~ m/no_nonred_reads=(\d+);coverage_ratio=([\d\.]+)/) {
#				$readcount = int ($1*$2);	
#				$readcount >= $coverage or next;		#clipped normal read count of the variant call (this could even be lower than mut allele count)
#				$maxcoverage and $readcount <= $maxcoverage || next;
#			} elsif ($attribute =~ m/gap-nonred-reads=(\d+)/) { ### <<< FOR 5500SOLiD LifeScope
#				$gapnonred = $1;
#				$attribute =~ m/coverage_ratio=(\d+)/;
#				$readcount = int($gapnonred*$1);
#				$readcount >= $coverage or next;
#				$maxcoverage and $readcount <= $maxcoverage || next;
#			} else {
#				$readcount = '.';
#			}
#			if ($attribute =~ m/allele-counts=REF,(\d+)/) {
#				$mutallelecount = $1;
#			} elsif ($attribute =~ m/context-variant-reads=(\d+)/) { ### <<< FOR 5500SOLiD LifeScope
#			    	$mutallelecount = $1;
#			}
#			if ($attribute =~ m#reference=([\w\-]+)#) { ### <<< FOR 5500SOLiD LifeScope (using "reference" tag for an allel-call in the reference) 
#			       	$refcall = $1;
#				$attribute =~ m#;allele\-call=([\w\-\/]+)#;
#				foreach my $item(split(/\//, $1)) { $temphash{$item}++; } # collecting unique alleles
#				@call = keys %temphash;	      
#
#				if ($1 eq '-/-') { # a simple deletion [";allele-call=-/-"]
#					print $chr, "\t", $pos, "\t", $end, "\t", $refcall, "\t", '-', "\t", $zygosity, "\t", "$score\t$readcount\t$mutallelecount\t", join ("\t", @other), "\n";
#				} elsif ($refcall eq '-') { # a simple insertion (single or multiple allele) ["reference=-"]
#					for my $i (0 .. @call-1) {					    
#					    	next if ($refcall eq $call[$i]);
#						print $chr, "\t", $pos, "\t", $pos, "\t", '-', "\t", $call[$i], "\t", $zygosity, "\t", "$score\t$readcount\t$mutallelecount\t", join ("\t", @other), "\n";
#						$i > 0 and $countallvar++;
#					}
#				} else { # an indel that may have several alleles, or may require a block substitution representation (use "context-reference-seq" and "contexit-variant-seq")
#					for my $i (0 .. @call-1) {
#					    	next if ($refcall eq $call[$i]);
#						print $chr, "\t", $pos, "\t", $pos+length($call[0])-1, "\t", $refcall, "\t", $call[$i], "\t", $zygosity, "\t", "$score\t$readcount\t$mutallelecount\t", join ("\t", @other), "\n";
#						$i > 0 and $countallvar++;
#					}
#				} 
#			} elsif ($attribute =~ m#allele\-call=([\w\/]+)#) {
#				@call = split (/\//, $1);
#				
#				if (@call == 1) {		#a simple deletion
#					print $chr, "\t", $pos, "\t", $end, "\t", $call[0], "\t", '-', "\t", $zygosity, "\t", "$score\t$readcount\t$mutallelecount\t", join ("\t", @other), "\n";
#				} elsif ($call[0] eq '') {	#a simple insertion (single or multiple allele)
#					for my $i (1 .. @call-1) {
#						print $chr, "\t", $pos, "\t", $pos, "\t", '-', "\t", $call[$i], "\t", $zygosity, "\t", "$score\t$readcount\t$mutallelecount\t", join ("\t", @other), "\n";
#						$i > 1 and $countallvar++;
#					}
#				} else {			#an indel that may have several alleles, or may require a block substitution representation
#					for my $i (1 .. @call-1) {
#						print $chr, "\t", $pos, "\t", $pos+length($call[0])-1, "\t", $call[0], "\t", $call[$i], "\t", $zygosity, "\t", "$score\t$readcount\t$mutallelecount\t", join ("\t", @other), "\n";
#						$i > 1 and $countallvar++;
#					}
#				}
#			} else {
#				$call = '0';
#				print $chr, "\t", $pos, "\t", $end, "\t", $call, "\t", '-', "\t", $zygosity, "\t", "$score\t$readcount\t$mutallelecount\t", join ("\t", @other), "\n";
#			}
#		} else {
#			die "Error: unrecognizable genotype calling program encountered (valid types are SOLiD_diBayes, AB_CNV_PIPELINE, AB_SOLiD Large Indel Tool, AB_SOLiD Small Indel Tool): <$_>\n";
#		}
#			
#		$countvar++;		#variation positions
#		$countallvar++;		#all variants (maybe several at one variation position)
#	}
#	print STDERR "NOTICE: Finished processing $variantfile with $countline input lines\n";
#	print STDERR "NOTICE: Wrote variants in $countvar variation positions ($countallvar variants considering multi-allelic ones)\n";
#	$unknown_count and print STDERR "WARNING: $unknown_count variants with 'unknown' variation type were skipped\n";
#}
#
#
#
#sub convertSOAP {
#	my ($variantfile) = @_;
#	my ($countline, $countvar, @other);
#	
#	open (VAR, $variantfile) or die "Error: cannot read from variant file $variantfile: $!\n";
#	
#	while (<VAR>) {
#		s/[\r\n]+$//;
#		$countline++;
#		
#		my @field = split (/\t/, $_);
#		if (@field == 18) {		#snp file
#			my ($chr, $pos, $wt, $call, @other) = @field;
#			$chr =~ s/^chr//;
#	
#			$includeinfo or @other = ();			#unless -includeinfo is set, the other will not be printed
#	
#			my $obs = $iupac{$call} or die "Error: invalid best call in <$_>\n";
#			my @obs = split (//, $obs);
#			@obs == 2 or die "Error: observed IUPAC allele $call should correspond to two nucleotide alleles: <$_>\n";
#			if ($obs[0] eq $wt and $obs[1] eq $wt) {
#				die "Error: reference alleles are identical to observed alleles: <$_>\n";
#			} elsif ($obs[0] eq $wt) {
#				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[1], "\t", "het\t", join ("\t", @other), "\n";
#			} elsif ($obs[1] eq $wt) {
#				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[0], "\t", "het\t", join ("\t", @other), "\n";
#			} elsif ($obs[1] ne $obs[0]) {
#				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[0], "\t", "het\t", join ("\t", @other), "\n";
#				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[1], "\t", "het\t", join ("\t", @other), "\n";
#				$countvar++;
#			} else {
#				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[0], "\t", "hom\t", join ("\t", @other), "\n";
#			}
#		} elsif (@field == 6) {		#indel file
#			my ($chr, $pos, $strand, $indellen, $call, $homo) = @field;
#			$homo eq 'Homo' and $homo = 'hom';
#			$homo eq 'Hete' and $homo = 'het';
#			$chr =~ s/^chr//;
#	
#			$includeinfo or @other = ();			#unless -includeinfo is set, the other will not be printed
#	
#			if ($indellen =~ m/^\+(\d+)$/) {		#insertion
#				length ($call) == $1 or die "Error: invalid record found in SOAPindel file: <$_>\n";
#				print join("\t", $chr, $pos, $pos, '-', $call, $homo), "\n";
#			} elsif ($indellen =~ m/^\-(\d+)$/) {		#deletion
#				length ($call) == $1 or die "Error: invalid record found in SOAPindel file: <$_>\n";
#				print join("\t", $chr, $pos, $pos+$1-1, $call, '-', $homo), "\n";
#			} else {
#				die "Error: invalid record found in SOAPindel file: <$_>\n";
#			}
#		} else {
#			die "Error: invalid record found in $variantfile (18 or 6 fields expected, observed ${\(scalar @field)} fields): <$_>\n";
#		}
#		$countvar++;
#	}
#	print STDERR "NOTICE: Read $countline lines and wrote $countvar variants\n";
#}
#
#sub iub{
#	%iub=("K"=>"GT","S"=>"GC","W"=>"AT","M"=>"AC","Y"=>"CT","R"=>"AG" );
#}
#sub convertANNOVAR {
#	my ($variantfile) = @_;
#	my ($countline, $countvar, $countsnp);
#	my ($countti, $counttv) = (0, 0);
#	
#	open (VAR, $variantfile) or die "Error: cannot read from variant file $variantfile: $!\n";
#	open (OUT, ">$outfile") or die "Error could not write to $opt_o:$!\n";
#	open (INVALID , ">$outfile.invalid") or die "Error could not write to $opt_o.invalid:$!\n";
#	while (<VAR>) {
#		$countline++;
#		if ($_=~/^#/){#only want top header
#			print OUT "$_" if ($comment);
#			$comment=0;
#			next;
#		}
#		my @field = split (/\t/, $_);
#		if ($#field<4){
#			#try to split by space 
#			@field =split(/\s/,$_,5);
#			die "Cannot convert $variantfile \n" if ($#field<4);
#		}
#		my ($chr, $start, $end, $ref, $obs, @rest) = @field;
#		if ($chr=~/(chr){0,1}[\dXYMT]{1,2}/ && $start=~/^\d+$/ && $end=~/^\d+$/ && $ref=~/^[ATCG\.\-0]{1,}/i ){
#			$ref=uc($ref);$obs=uc($obs);
#			if ( $obs=~/^[ATCG\.\-0]{1,}$/i){
#				#ok
#			}elsif ($ref=~/^[ATCG]$/ && $obs=~/^[KSWMRY]$/){#snp
#				$obs=~s/$ref//;
#				print INVALID "$_" and next if (!$obs);
#			}elsif ($obs=~/^[KSWMRY]$/){
#				print INVALID "$_" and next;
#			}
#			
#		}
#		print OUT join ("\t",$chr,$start,$end,$ref,$obs,@rest)."\n";
#		if ($ref =~ m/^[ACGT]$/ and $obs =~ m/^[ACGT]$/) {
#			if ($ref eq 'A' and $obs eq 'G' or $ref eq 'G' or $obs eq 'A' or $ref eq 'C' and $obs eq 'T' or $ref eq 'T' and $obs eq 'C') {
#				$countti++;
#			} else {
#				$counttv++;
#			}
#			$countsnp++;
#		}
#		
#		print;
#		$countvar++;
#	}
#	close OUT;close INVALID;close VAR;
#	print STDERR "NOTICE: Read $countline lines and wrote $countvar variants\n";
#	print STDERR "NOTICE: Among $countsnp SNPs, $countti are transitions, $counttv are transversions\n";
#}

#sub convertMAQSNP {
#	my ($variantfile) = @_;
#	my ($countline, $countvar, @other);
#	
#	open (VAR, $variantfile) or die "Error: cannot read from variant file $variantfile: $!\n";
#	
#	while (<VAR>) {
#		s/[\r\n]+$//;
#		$countline++;
#		
#		my @field = split (/\t/, $_);
#		my @other = ();
#		if (@field == 12) {					#SNPs
#			my ($chr, $pos, $wt, $call, @other) = @field;
#			$chr =~ s/^chr//;
#	
#			$includeinfo and @other = @field;			#unless -includeinfo is set, the other will not be printed
#	
#			my $obs = $iupac{$call} or die "Error: invalid best call in <$_>\n";
#			my @obs = split (//, $obs);
#			@obs == 2 or die "Error: observed IUPAC allele $call should correspond to two nucleotide alleles: <$_>\n";
#			if ($obs[0] eq $wt and $obs[1] eq $wt) {
#				die "Error: reference alleles are identical to observed alleles: <$_>\n";
#			} elsif ($obs[0] eq $wt) {
#				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[1], "\t", "het\t", join ("\t", @other), "\n";
#			} elsif ($obs[1] eq $wt) {
#				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[0], "\t", "het\t", join ("\t", @other), "\n";
#			} elsif ($obs[1] ne $obs[0]) {
#				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[0], "\t", "het\t", join ("\t", @other), "\n";
#				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[1], "\t", "het\t", join ("\t", @other), "\n";
#				$countvar++;
#			} else {
#				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[0], "\t", "hom\t", join ("\t", @other), "\n";
#			}
#			$countvar++;
#		} elsif (@field == 13) {				#indels; the deletion start site do not need changes; the duplication start site need additional investigation by ANNOVAR developers
#			my ($chr, $pos, $type, $numread, $call, @other) = @field;
#			$chr =~ s/^chr//;
#	
#			$includeinfo and @other = @field;			#unless -includeinfo is set, the other will not be printed
#	
#			my @obs = split (/:/, $call);
#			@obs == 2 or die "Error: observed IUPAC allele $call should correspond to two nucleotide alleles: <$_>\n";
#			if ($obs[0] =~ m/^\-(\d+)/) {		#deletion
#				my $len = $1;
#				print $chr, "\t", $pos, "\t", $pos+$len-1, "\t", $obs[1], "\t", '-', "\t", "het\t", join ("\t", @other), "\n";
#			} elsif ($obs[0] =~ m/^(\d+)/) {	#insertion
#				my $len = $1;
#				print $chr, "\t", $pos-1, "\t", $pos-1, "\t", '-', "\t", $obs[1], "\t", "het\t", join ("\t", @other), "\n";	#2011jul12: changed pos to pos-1 for insertions
#			}
#			$countvar++;
#		} else {
#			die "Error: invalid record found in $variantfile (12 or 13 fields expected, observed ${\(scalar @field)} fields): <$_>\n";
#		}
#	}
#	print STDERR "NOTICE: Read $countline lines and wrote $countvar variants\n";
#}
#
sub convertCASAVA_noChr {
	my ($variantfile) = @_;
	my ($countline, $countvar, @other);
	
	my ($intype,%bcindex);
	my ($pos_index, $call_index, $reference_index, $type_index, $score_index, $total_index, $used_index,$chr_index);
	my ($ref_indel_index, $quality_index, $maxgtype_index, $maxplygt_index,$bp1_reads_index, $ref_reads_index, $indel_reads_index, $other_reads_index);
	open (VAR, $variantfile) or die "Error: cannot read from variant file $variantfile: $!\n";
	if ($af){
		if (defined $outfile){
			if ($append){
				open (AF,">>$outfile.hg19_af") or die "Cannot open $outfile.hg19_af\n";
			}else{
				open (AF,">$outfile.hg19_af") or die "Cannot open $outfile.hg19_af\n";
			}
		}else{
			if ($append){
				open (AF,">>ANNOVAR.input.hg19_af") or die "Cannot open AF\n";
			}else{
				open (AF,">ANNOVAR.input.hg19_af") or die "Cannot open AF\n";
			}
		}
		print AF "#AlleleFreq\t#Depth\t#Variant Type\n" if ($includeinfo || $header);
	}
	while (<VAR>) {
		s/[\r\n]+$//;
		$countline++;
		my @field;

		if ($countline==1) {
			if (s/^\$\sCOLUMNS\s//) {
				@field = split (/\s+/, $_);
			} else {
				@field = split (/\t/, $_);
			}
			if (m/\bposition\b/ or m/\bpos\b/) {
				for my $i (0 .. @field-1) {
					$field[$i]=~s/\s+$//g;
					if ($field[$i] eq 'position' or $field[$i] eq 'pos') {
						$pos_index = $i;
					}elsif ( $field[$i] eq 'seq_name' || $field[$i]=~/^chr/i){
						$chr_index=$i;
					} elsif ($field[$i] eq 'modified_call') {#not used
						$intype = 'snp';
						print STDERR "NOTICE: Automatically detected input type as $intype\n";
						$call_index = $i;
					} elsif ($field[$i] eq 'reference' || $field[$i] eq 'ref') {#snps only
						$reference_index = $i;
					} elsif ($field[$i] eq 'type') {#indels only
						$type_index = $i;
					} elsif ($field[$i] eq 'score') {#not used
						$score_index = $i;
					} elsif ($field[$i] eq 'depth') {#indels only hv
						$total_index = $i;
					} elsif ($field[$i] eq 'bcalls_used') {#snps only #deltahv
						$used_index = $i;
					} elsif ($field[$i] eq 'ref/indel') {#indels only
						$intype = 'indel';
						print STDERR "NOTICE: Automatically detected input type as $intype\n";
						$ref_indel_index = $i;
					} elsif ($field[$i] eq 'Q(indel)') {#indels only
						$quality_index = $i;
					} elsif ($field[$i] eq 'max_gtype') {#snps only
						$maxgtype_index = $i;
					}elsif ($field[$i] eq 'max_gt|poly_site'){#snps only hv
						$maxplygt_index=$i;						
					} elsif ($field[$i] eq 'bp1_reads') {#not used
						$bp1_reads_index = $i;
					} elsif ($field[$i] eq 'ref_reads') {#not used
						$ref_reads_index = $i;
					} elsif ($field[$i] eq 'indel_reads') {#indels only
						$indel_reads_index = $i;
					} elsif ($field[$i] eq 'other_reads') {#indels only
						$other_reads_index = $i;
					} elsif ($field[$i]=~/([ATCG])\_used/){#snps only hv
						$bcindex{$1}=$i;
						$intype='snp';
					}
				}
				print "#zygosity" if ($includeinfo || $header);
				print "\t#$_" if ($includeinfo);
				print "\n" if ($includeinfo || $header);
			}
			next;
		}
		
		$intype or die "Error: unable to recognize the correct type of the input file (make sure that header line is present in $variantfile)\n";
		@field = split (/\t/, $_);
		
		if ($intype eq 'snp') {					#SNPs
			defined $chr_index and defined $pos_index and defined $reference_index and defined $used_index and defined $maxplygt_index and keys(%bcindex)==4 or die "Error: unable to find the position, reference and modified_call column header in $variantfile\n";
			my ($chr,$pos, $wt,$obs,$total) = @field[$chr_index,$pos_index, $reference_index, $maxplygt_index,$used_index];
			my (@other);
			defined $pos and defined $wt and defined $obs or die "Could not find one or mmore of your columns correctly at file line $countline\n";
			$includeinfo and @other = @field;
			length ($obs) == 1 and $obs .= $obs;#what does this do?hv
			my @obs = split (//, uc($obs));
			if (!$allallele){
				my $myobs;
				@obs == 2 or die "Error: observed allele $obs should correspond to two nucleotide alleles: <$_>\n";
				if ($obs[0] eq $wt and $obs[1] eq $wt) {#do not print out
					next;# "Error: reference alleles are identical to observed alleles: <$_>\n";
				} elsif ($obs[0] eq $wt) {
					$myobs=$obs[1];
					print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[1], "\t", "het\t", join ("\t", @other), "\n";
					
				} elsif ($obs[1] eq $wt) {
					$myobs=$obs[0];
					print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[0], "\t", "het\t", join ("\t", @other), "\n";
					
				} elsif ($obs[1] ne $obs[0]) {
					print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[0], "\t", "het\t", join ("\t", @other), "\n";
					print AF sprintf("%0.03f",($field[$bcindex{$obs[0]}]/$total))."\t$total\n" if ($af);
					$myobs=$obs[1];
					print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[1], "\t", "het\t", join ("\t", @other), "\n";
					$countvar++;
					
				} else {
					$myobs=$obs[0];
					print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[0], "\t", "hom\t", join ("\t", @other), "\n";
					
				}
				print AF sprintf("%0.03f",($field[$bcindex{$myobs}]/$total))."\t$total\t$intype\n" if ($af);
				$countvar++;
			}elsif ($allallele){
				foreach my $allele (keys %bcindex){
					next if $allele eq $wt;
					next if ($field[$bcindex{$allele}]==0);
					if ($obs!~/$allele/){#not a major allele
						print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $allele, "\t", "minor_allele\t", join ("\t", @other), "\n";
					}elsif ($obs[0] ne $obs[1]){
						print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $allele, "\t", "het\t", join ("\t", @other), "\n";
					}else{
						print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $allele, "\t", "hom\t", join ("\t", @other), "\n";
					}
					print AF sprintf("%0.03f",($field[$bcindex{$allele}]/$total))."\t$total\t$intype\n" if ($af);
				}
			}
		} elsif ($intype eq 'indel') {		#indels
			defined $total_index and defined $indel_reads_index and defined $chr_index and defined $pos_index and defined $ref_indel_index and defined $maxgtype_index or die "Error: unable to find the $chr_index,$pos_index, $ref_indel_index, $maxgtype_index(chr,pos, ref_indel and max_gtype) column header in $variantfile\n";
			my ($chr ,$pos, $call, $total, $indel_count, $hom,@other) = @field[$chr_index,$pos_index, $ref_indel_index, $total_index,$indel_reads_index, $maxgtype_index];
			$includeinfo and @other = @field;
			#hg19 coordinate below; insertion needs position adjustment!!! deletion is fine
			#948847  1I      CCTCAGGCTT      -/A     ATAATAGGGC      969     hom     47      het     22      0       16      6       A       1       2
			#978604  2D      CACTGAGCCC      CT/--   GTGTCCTTCC      251     hom     20      het     8       0       4       4       CT      1       0
			#1276974 4I      CCTCATGCAG      ----/ACAC       ACACATGCAC      838     hom     39      het     18      0       14      4       AC      2       4
			#1289368 2D      AGCCCGGGAC      TG/--   GGAGCCGCGC      1376    hom     83      het     33      0       25      9       TG      1       0
			#185137455     11I10M2I        TATGTGTCCT      -----------TTTTTTATTT--/AAATGATAGACTTTTTTTTTTAA ATTTCAGAAA      1126    het     988     hom    45       20      24      7       N/A     0       0
			#1276931 2D41M4I CACACACATG      CACACACACGCACACACGTGCAATGTGAAAACACCTCATGCAG----/--CACACACGCACACACGTGCAATGTGAAAACACCTCATGCAGACAC ACACATGCAC      548     hom     16      het     8       0       11      11      N/A     0       0
			
			my @obs = split (/\//, $call);
			@obs == 2 or die "Error: observed indel allele $call should correspond to two alleles: <$_>\n";
			if ($obs[0] =~ m/^\-+$/) {		#insertion
				my $len = length ($obs[0]);
				print $chr, "\t", $pos-1, "\t", $pos-1, "\t", '-', "\t", $obs[1], "\t", $hom, "\t", join ("\t", @other), "\n";
			} elsif ($obs[1] =~ m/^\-+$/) {		#deletion
				my $len = length ($obs[0]);
				print $chr, "\t", $pos, "\t", $pos+$len-1, "\t", $obs[0], "\t", '-', "\t", $hom, "\t", join ("\t", @other), "\n";
			} elsif (length ($obs[0]) eq length ($obs[1])) {	#block substitution
				$obs[0] =~ s/\-//g;
				$obs[1] =~ s/\-//g;
				print STDERR $chr, "\t", $pos, "\t", $pos+length($obs[0])-1, "\t", $obs[0], "\t", $obs[1], "\t", $hom, "\t", join ("\t", @other), "\n";
			} else {
				die "Error: invalid record found in indel line: $obs[0],$obs[1]\n<$_>\n";
			}
			print AF sprintf("%0.03f",($indel_count/$total))."\t$total\t$intype\n" if ($af);
			$countvar++;
		} else {
			die "Error: invalid record found in $variantfile (11 or 15 fields expected, observed ${\(scalar @field)} fields): <$_>\n";
		}
	}
	print STDERR "NOTICE: Read $countline lines and wrote $countvar variants\n";
}
sub convertCASAVA {
	my ($variantfile, $chr) = @_;
	my ($countline, $countvar, @other);
	
	my ($intype);
	my ($pos_index, $call_index, $reference_index, $type_index, $score_index, $total_index, $used_index);
	my ($ref_indel_index, $quality_index, $maxgtype_index, $bp1_reads_index, $ref_reads_index, $indel_reads_index, $other_reads_index);
	open (VAR, $variantfile) or die "Error: cannot read from variant file $variantfile: $!\n";
	
	while (<VAR>) {
		s/[\r\n]+$//;
		$countline++;
		my @field;

		if (m/^#/) {
			s/^#//;
			if (s/^\$\sCOLUMNS\s//) {
				@field = split (/\s+/, $_);
			} else {
				@field = split (/\t/, $_);
			}
			if (m/\bposition\b/ or m/\bpos\b/) {
				for my $i (0 .. @field-1) {
					if ($field[$i] eq 'position' or $field[$i] eq 'pos') {
						$pos_index = $i;
					} elsif ($field[$i] eq 'modified_call') {
						$intype = 'snp';
						print STDERR "NOTICE: Automatically detected input type as $intype\n";
						$call_index = $i;
					} elsif ($field[$i] eq 'reference') {
						$reference_index = $i;
					} elsif ($field[$i] eq 'type') {
						$type_index = $i;
					} elsif ($field[$i] eq 'score') {
						$score_index = $i;
					} elsif ($field[$i] eq 'total') {
						$total_index = $i;
					} elsif ($field[$i] eq 'used') {
						$used_index = $i;
					} elsif ($field[$i] eq 'ref/indel') {
						$intype = 'indel';
						print STDERR "NOTICE: Automatically detected input type as $intype\n";
						$ref_indel_index = $i;
					} elsif ($field[$i] eq 'Q(indel)') {
						$quality_index = $i;
					} elsif ($field[$i] eq 'max_gtype') {
						$maxgtype_index = $i;
					} elsif ($field[$i] eq 'bp1_reads') {
						$bp1_reads_index = $i;
					} elsif ($field[$i] eq 'ref_reads') {
						$ref_reads_index = $i;
					} elsif ($field[$i] eq 'indel_reads') {
						$indel_reads_index = $i;
					} elsif ($field[$i] eq 'other_reads') {
						$other_reads_index = $i;
					}
				}
			}
			next;
		}
		
		$intype or die "Error: unable to recognize the correct type of the input file (make sure that header line is present in $variantfile)\n";
		@field = split (/\t/, $_);
		
		if ($intype eq 'snp') {					#SNPs
			defined $pos_index and defined $reference_index and defined $call_index or die "Error: unalbe to find the position, reference and modified_call column header in $variantfile\n";
			my ($pos, $wt, $obs) = @field[$pos_index, $reference_index, $call_index];
			my (@other);
			defined $pos and defined $wt and defined $obs or die;
			$includeinfo and @other = @field;
			
			length ($obs) == 1 and $obs .= $obs;
			my @obs = split (//, $obs);
			@obs == 2 or die "Error: observed allele $obs should correspond to two nucleotide alleles: <$_>\n";
			if ($obs[0] eq $wt and $obs[1] eq $wt) {
				die "Error: reference alleles are identical to observed alleles: <$_>\n";
			} elsif ($obs[0] eq $wt) {
				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[1], "\t", "het\t", join ("\t", @other), "\n";
			} elsif ($obs[1] eq $wt) {
				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[0], "\t", "het\t", join ("\t", @other), "\n";
			} elsif ($obs[1] ne $obs[0]) {
				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[0], "\t", "het\t", join ("\t", @other), "\n";
				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[1], "\t", "het\t", join ("\t", @other), "\n";
				$countvar++;
			} else {
				print $chr, "\t", $pos, "\t", $pos, "\t", $wt, "\t", $obs[0], "\t", "hom\t", join ("\t", @other), "\n";
			}
			$countvar++;
		} elsif ($intype eq 'indel') {				#indels
			defined $pos_index and defined $ref_indel_index and defined $maxgtype_index or die "Error: unable to find the pos, ref_indel and max_gtype column header in $variantfile\n";
			my ($pos, $call, $hom, @other) = @field[$pos_index, $ref_indel_index, $maxgtype_index];
			$includeinfo and @other = @field;

			#hg19 coordinate below; insertion needs position adjustment!!! deletion is fine
			#948847  1I      CCTCAGGCTT      -/A     ATAATAGGGC      969     hom     47      het     22      0       16      6       A       1       2
			#978604  2D      CACTGAGCCC      CT/--   GTGTCCTTCC      251     hom     20      het     8       0       4       4       CT      1       0
			#1276974 4I      CCTCATGCAG      ----/ACAC       ACACATGCAC      838     hom     39      het     18      0       14      4       AC      2       4
			#1289368 2D      AGCCCGGGAC      TG/--   GGAGCCGCGC      1376    hom     83      het     33      0       25      9       TG      1       0
			#185137455     11I10M2I        TATGTGTCCT      -----------TTTTTTATTT--/AAATGATAGACTTTTTTTTTTAA ATTTCAGAAA      1126    het     988     hom    45       20      24      7       N/A     0       0
			#1276931 2D41M4I CACACACATG      CACACACACGCACACACGTGCAATGTGAAAACACCTCATGCAG----/--CACACACGCACACACGTGCAATGTGAAAACACCTCATGCAGACAC ACACATGCAC      548     hom     16      het     8       0       11      11      N/A     0       0
			
			my @obs = split (/\//, $call);
			@obs == 2 or die "Error: observed indel allele $call should correspond to two alleles: <$_>\n";
			if ($obs[0] =~ m/^\-+$/) {		#insertion
				my $len = length ($obs[0]);
				print $chr, "\t", $pos-1, "\t", $pos-1, "\t", '-', "\t", $obs[1], "\t", $hom, "\t", join ("\t", @other), "\n";
			} elsif ($obs[1] =~ m/^\-+$/) {		#deletion
				my $len = length ($obs[0]);
				print $chr, "\t", $pos, "\t", $pos+$len-1, "\t", $obs[0], "\t", '-', "\t", $hom, "\t", join ("\t", @other), "\n";
			} elsif (length ($obs[0]) eq length ($obs[1])) {	#block substitution
				$obs[0] =~ s/\-//g;
				$obs[1] =~ s/\-//g;
				print $chr, "\t", $pos, "\t", $pos+length($obs[0])-1, "\t", $obs[0], "\t", $obs[1], "\t", $hom, "\t", join ("\t", @other), "\n";
			} else {
				die "Error: invalid record found in indel line: <$_>\n";
			}
			$countvar++;
		} else {
			die "Error: invalid record found in $variantfile (11 or 15 fields expected, observed ${\(scalar @field)} fields): <$_>\n";
		}
	}
	print STDERR "NOTICE: Read $countline lines and wrote $countvar variants\n";
}
#
sub convertCLC{#added loop for jason lih
	my ($variantfile) = @_;
	my ($countline, $countvar, $counthom, $counthet, $countunknown, $countindel, $countsnp, $countti, $counttv) = qw/0 0 0 0 0 0 0 0 0/;#ANNOVAR default
	my %stats;#hv added 2/10/12
	#CLC Bio requires the following header
	#not counting allele freq or count freq
	#Mapping	Reference Position	Consensus Position	Variation Type	Length	Reference	Variants	Allele Variations	Frequencies	Counts	Coverage	Variant #1	Frequency of #1	Count of #1	Variant #2	Frequency of #2	Count of #2
	open (VAR, $variantfile) or die "Error: cannot read from variant file $variantfile: $!\n";
	VAR: while (my $line=<VAR>) {
		$countline++;
		$line=~s/[\n\r]//g;
		next if ($line =~/^[\n\r\s]*$/);
			
		my @field=split(/\t/,$line);chomp $field[$#field];
		@field >= 14 or die "Error: invalid record found in CLCBio file (at least 14 tab-delimited fields expected): <$_> at line $countline\n";
		my $otherinfo="";
		if ($includeinfo){
			$otherinfo=join ("\t",@field);
		}
		my $end;
		my ($chr,$start,$ref,$mut)=@field[0,1,5,7];$ref=uc($ref);$mut=uc($mut);
		if ($countline==1){
			die "Expecting the following Fields: Mapping,Reference Position,Reference,Allele Variations.  Found:$chr,$start,$ref,$mut\n" if (join (",",$chr,$start,$ref,$mut)!~/Mapping,Reference Position,Reference,Allele Variations/i);
			print "#" if ($field[0]!~/#/);
			print join("\t","ABCC_RefId",@field)."\n";
		
		}else{
			#this is your variant lines
			$chr=~s/\smapping//g;
			$end=$start;
			my $note='';
			my @mutations=split("/",$mut);
			foreach my $var (@mutations){$note='';
				next if ($var eq "$ref" || $var=~/N/);
				my ($myref,$myvar,$mystart,$end)=($ref,$var,$start,$end);#we do this for complex mutations
				if ($field[3]=~/^(complex\s){0,1}snp/i){
				}elsif ($ref eq '-' || $ref=~/^\-/){#This is an insertion
					#start==end and the insertion ALLELE,$note
					$note=",ins$var";	
				}elsif ($var eq '-' || $var=~/^\-/){
					$end+=length($ref)-1;
					$note=",del$ref";
				}elsif (length($ref) == length($var)){	
					if (my $head=smallest_substring($ref,$var)){
						chomp $head;
						if ($head eq $var) {#deletion
								$mystart=$mystart+length($head);
								$myref=substr ($ref, length($head),length($myref)-length($head));
								$end=$start+length($myref)-1;
								$myvar="-";
								$note=",del$head*";
						}elsif ($head eq $myref){#insertion
								$myref='-';
								$myvar=~s/$head//;
								$note=",ins$head*";
						}else{#snp or multiallelic snp
								$myref=~s/$head//;
								$myvar=~s/$head//;
								$mystart+=length($head);
								$end=$mystart+length($myref)-1;
								if (length($myref)==1){
									$note=",snp*" ;
								}else{
									$note=",mnp*";
								}
						}	
					}else{
						$note=",mnp*" if (length($myref)>1);
						$end+=length($myref)-1;
		#				print "($head)".join (",",@field)."\n";<STDIN>;
					}			
				}else{
					die "[ERROR]($myvar and $myref)\n".join ("..",@field)."\n";
				}
				print join ("\t",$chr,$mystart,$end,$myref,$myvar,"ABCC_id=$countline$note",$otherinfo)."\n";
			}
		}
	}
	close VAR;
	print STDERR "NOTICE: Completed at ". `date`;
}
sub convertHGVS_wGpos {#added loop for jason lih
	my ($variantfile,$hgvs_flag) = @_;
	my ($countline, $countvar, $counthom, $counthet, $countunknown, $countindel, $countsnp, $countti, $counttv) = qw/0 0 0 0 0 0 0 0 0/;#ANNOVAR default
	my %stats;#hv added 2/10/12
	my ($source_program, $gtpos);		#the program that generated the VCF4 file; the GT position within FORMAT record
	open (VAR, $variantfile) or die "Error: cannot read from variant file $variantfile: $!\n";
	#this requires knowledge about the strand of the gene.
	open (STRAND,"</SeqIdx/annovardb/humandb/hg19_refGene2Strand.txt") or die 'Your file for genic strands does not exist!';
	while (<STRAND>){
		my ($gene,$ori)=split("\t",$_);chomp $ori;
		$gene_ori{uc($gene)}=$ori;
	}
	close STRAND;
	VAR: while (my $line=<VAR>) {
		$countline++;
		$line=~s/[\n\r]//g;
		my @field=split(/\t/,$line);chomp $field[$#field];
		@field >= 4 or die "Error: invalid record found in HGVS file (at least 5 tab-delimited fields expected): <$_>\n";
		my $otherinfo;
		if ($includeinfo){
			$otherinfo = join("\t",@field);	#this is the complete line (when -includeinfo is set, the entire line will be included in output file)
		}
		my ($chr, $start, $stop,$gene,$hgvs,$hgvs2) = @field;
		$gene=uc($gene);
		$chr='chr'.$chr if ($chr!~/^chr/);
		if ($countline==1 ){#first line must have header info
			die "Your inputs are not correctly formatted by Chr\tStart\tStop\tGene\tHGVS Mutation..<$_>\n" if ($line!~/Chr\tStart\t(Stop|End)\tGene/i);
			print "#" if ($field[0]!~/#/);
			print join("\t",@field)."\n";
		}else{
			my ($strand,$err_flag);
			if (exists $gene_ori{$gene}){ 
				$strand=$gene_ori{$gene};
			}else{
				$strand = '+';
			}
			my ($ref,$mut);
			my $nts="[ATCG]+";
			if (!$hgvs_flag){#default nothing
				if ($hgvs=~/($nts)>($nts)/){
					$ref=$1;$mut=$2;
					if ( ($stop-$start) != (length($ref)-1)){
						$stop=$start+length($ref)-1;
					}
				}elsif ($hgvs=~/del($nts)/ && $hgvs!~/ins/){
					$ref=$1,$mut='-';
				}elsif ($hgvs=~/ins($nts)/ && $hgvs!~/del/){
					$ref='-';$mut=$1;
					if ($strand eq '-'){
						$start=$stop;
					}else{
						$stop=$start;
					}
				}elsif ($hgvs=~/intronic/i){
					$ref=0;$mut=0;
				}else{
					$ref=0;$mut=0;
				}
			}else{
				$ref=$hgvs;
				$mut=$hgvs2;
			}
			if ($strand eq '-'){
				$ref=revcompl($ref);
				$mut=revcompl($mut);
			}
			print join ("\t",$chr,$start,$stop,$ref,$mut);
			if ($includeinfo){
				print "\t$otherinfo";
			}
			print "\n";
						
		}
		
		$countvar++;
	}
	my $triallelic = $countsnp-$countti-$counttv;
	print STDERR "NOTICE: Completed at ". `date`;
}
sub getMemory{
	my $size=`ps -p $$ -o size | grep -ve SIZE`;chomp $size;
	print STDERR "MEMORY USAGE:$size\n";
}
sub readAnnot{
	my $file=shift;
	open (STRAND,"<$file") or die 'Your file for genic strands does not exist!';
	while (<STRAND>){
		my (@info)=split("\t",$_);chomp $info[$#info];
		next if ($info[1]=~/^NR/ && !$ensembl);
		if ($info[2]=~/^(chr){0,1}[\dXYMT]{1,2}$/ ){
			# $matches{$info[1]}++ if (exists $gene_ori{uc($info[12])}{$info[1]});
			$gene_ori{uc($info[12])}{$info[1]}.=join("\t",@info[1..3],@info[6..7],@info[9..10],0,0)."\n";
			$gene_ori{$info[1]}=$gene_ori{uc($info[12])}{$info[1]};
			# print Dumper (\%gene_ori) and die if ($info[12] eq 'SETD2');
		}
	}
	close STRAND;
	my %alias;
	my $xreffile=$file;$xreffile=~s/\.txt/Xref.txt/;
	if (-e $xreffile){
		open (XREF,"<$xreffile") or die "Cannot open Xref synonym file\n";
		while (<XREF>){
			my ($alias,$gene)=split("\t",$_);chomp $gene;
			$alias{$alias}=$gene;
		}
		close XREF;
	}

	return (\%alias);
}
sub convertFromHGVS {#added loop for jack
	my %found;#keeps track of which lines have not been found
	my ($variantfile) = @_;
	my ($countline, $countvar, $counthom, $counthet, $countunknown, $countindel, $countsnp, $countti, $counttv) = qw/0 0 0 0 0 0 0 0 0/;#ANNOVAR default
	open (VAR, $variantfile) or die "Error: cannot read from variant file $variantfile: $!\n";
	open (FAIL,">notconverted.txt") or die "Error: cannot write to $variantfile.notcoverted\n";
	#this requires knowledge about the strand of the gene.
	my $annot;my $dbfa;
	if (!$ensembl){
		$annot='/SeqIdx/annovardb/humandb/hg19_refGene.txt';
		$dbfa='refGene';

	}else{
		$annot='/SeqIdx/annovardb/humandb/hg19_ensGene.txt';
		$dbfa="ensGene";
	}
	my $script_dir=`dirname $0`;chomp $script_dir;
	my $prot2trx_exe="php $script_dir/../easyConverter_avia.php ";
	readSeqFromFASTADB($dbfa);
	my $alias_ref=readAnnot($annot);#by default
	if ($annot=~/hg19_(\w+).txt/){$dbtype1=$1;}
	# my $num=keys (%matches);print $num;die;
	my $failedresponse='';# keeps track of the 3 genomic positions if the aa1 => aa2 conversion fails
	my %protConverter;
	print join ("\t","#Chr","Start","End","RefAllele","MutAllele","Comment")."\n";
	VAR: while (my $line=<VAR>) {
		$countline++;
		# print "working on $countline and $line" if ($verbose);
		$line=~s/[\n\r]//g;
		next if ($line=~/^\s*$/);
		next if ($line=~/(^#|exon)/);
		my @field=split(/\t/,$line);chomp $field[$#field];
		@field>0 or die "Error: invalid record found in HGVS file (at least 5 tab-delimited fields expected): <$_>\n";
		my ($stop,$start);
		my $otherinfo;
		if ($includeinfo){
			$otherinfo = join("\t",@field);	#this is the complete line (when -includeinfo is set, the entire line will be included in output file)
		}
		my ($input,@other) = @field;#ignore the rest
		my ($genetrx,$hgvs)=split(":",$input);#genetrx can either be a gene name or transcript
		if ($genetrx=~/^(NM_\d+)[\_\.]\d+/){
			$genetrx=$1;
		}else{
			if ($genetrx=~/^nm(\d+)/i){
				$genetrx='NM_'.$1;
			}elsif ($genetrx=~/^([A-Z0-9\-]*)\_HUMAN/i){
				$genetrx=$1;
				$genetrx=~s/\-//;
			}elsif ($genetrx=~/(ENSG\d+)/){
				$genetrx=$1;
			}elsif($genetrx=~/(ENSP\d+|NP\_*\d+)/i){
				#convert
				my $prot=uc($1);my ($ip,$op);
				if ($prot=~/NP\_*(\d+)/){
					$prot="NP_$1";$ip="RefSeq Protein Accession";$op="RefSeq mRNA Accession";
				}else{
					$ip="Ensembl Protein ID";$op="Ensembl Transcript ID";
				}
				if (keys %protConverter==0){#not loaded into memory
					print STDERR "[INFO] not loaded into memory yet...does it exist already??";
					if ($prot && $ip && $op && !-e "easyConverter.lkup"){##file not previously run, creating
						print STDERR "[INFO] This doesn't exist yet!... Running..\n";
						system ("grep -ve '#' $variantfile|cut -f 1 -d :|cut -f1 -d '.' |sort -u > .$variantfile");
						system ("$prot2trx_exe -f .$variantfile -i '$ip' -o '$op' > easyConverter.lkup\n" );
					}
					print STDERR "[INFO] OK, now did it actually run or exist??\n";
					if(-e "easyConverter.lkup" && !-z "easyConverter.lkup"){#biodbnet run and was successful, read contents into memory
						print STDERR "Yay! it exists...Reading into memory\n";
						open (CONVERTER,"<easyConverter.lkup") or die "Cannot open easyConverter.lkup file\n";
						while (<CONVERTER>){	
							my ($p1,$t1)=split("\t",$_);chomp $t1;
							next if ($t1 eq "-");
							$protConverter{$p1}=$t1;
						}
					}else{
						print STDERR "[ERR] Tried run bioDBnet and it can't connect...Exiting";
						print STDERR "[ERR] $prot could not be converted because bioDBnet took too long to respond...\n";exit;
						next;
					}
				}
				if (exists $protConverter{$prot}){
					$genetrx=$protConverter{$prot};
				}else{
					print STDERR "[ERR] $prot could not be converted to a transcript identifier using bioDBnet\n";next;
				}
			}elsif ($genetrx=~/^(\w+):/i){
				$genetrx=$1;
			}else{
			}
			$genetrx=uc($genetrx);
		}
		# print "$genetrx is it\n";
		if ($hgvs){
			$hgvs=~s/\s//g;
		}else{
			print STDERR "hgvs was empty and $line\n";
		}

		print STDERR "Skipping $line because complex mutation\n" and next if ($line=~/(delins|del.*ins|ins.*del|insdel)/);#skipping indels

		print STDERR "Could not deal with $hgvs\n" and next if (!$hgvs or $hgvs=~/\?/ or $input=~/(fs)/);
		my ($type,$pos,$ref,$var,$cvrtmeth,$indeltype);
		$indeltype='';
		$hgvs=~s/\*/X/g;
		if ($hgvs=~/^([cp])\.([A-Z\*]+)([\d\_]+)([A-Z\*]+)$/i){
			($type,$ref,$pos,$var)=($1,$2,$3,$4);$cvrtmeth=1;
		}elsif ($hgvs=~/^([cp])\.([\d\_\-\+]+\d)([\-A-Z])(>|&gt;)([A-Z\*\-]+)$/i){
			($type,$pos,$ref,$var)=($1,$2,$3,$5);$cvrtmeth=2;
			#deletion notation in an intronic splice site
			#c.285-1->CCCC
		}elsif ($hgvs=~/^([cp])\.([\d\_\-\+]+)([A-Z\*]+)$/i){
			# p.311_312del
			#c.120_121insAA
			($type,$pos,$ref,$var)=($1,$2,'-',$3);$cvrtmeth=3;
			if ($var eq 'del'){
				#can calculate?
				$var='0';
				$ref='-';$cvrtmeth.=".1";
			}elsif($var=~/(ins|del)(\w+)$/){
				$indeltype=$1;$var=$2;
				if ($indeltype=~/del/){
					$ref=$var;
					$var='-';
				}
				$cvrtmeth.=".2";
			}else{
				warn "don't know what this is($hgvs)\n";
				next;
			}
		# }elsif ($hgvs=~/^(p)\.([A-Z]*)(\d+)>([A-Z\*]+)/i){#protein insertion?
		# 	($type,$ref,$pos,$var)=($1,$2,$3,$4);
		# }elsif ($hgvs=~/^([cp])\.([\d\_\-\+]+)((ins|del)[A-Z\*]+)/i){
		# 	($type,$pos,$ref,$var)=($1,$2,'.',$3);
		# 	my $indeltype=$4;
		# 	$var=~s/$indeltype/0/g;
		# 	print join ("\t",$type,$pos,$ref,$var)."\n";
		}elsif ($hgvs=~/^([cp])\.{1,}([\d\+\-]+)([ATCG\*]{1,})(>|&gt;)([ATCG\*]{1,})/){
			($type,$pos,$ref,$var)=($1,$2,$3,$5);$cvrtmeth=4;
		}elsif ($hgvs=~/^([A-Z\*]+)([\d\_]+)([A-Z\*]+)$/i){#no type specified
			($ref,$pos,$var)=($1,$2,$3);$cvrtmeth=1;
			$type='p';$hgvs='p.'.$hgvs;
			# print  "now working on $hgvs instead of $input\n";
		}else{
			next if ($hgvs=~/(\_|intronic|\d\w{0,1}>\w)/i);
			warn "what's this ($hgvs)$line\n";
			next;
		}
		# print "$genetrx and $ref,$pos,$var and $type\n" if ($verbose);
		##finished checking input protein or cdna variant format...	

####################################################################
# Check the gff database to see that the TRX or gene exists in %gene_ori
###################################################################
		my ($err_flag);
		if (!exists $gene_ori{$genetrx}){
			if (!exists $$alias_ref{$genetrx}){
				print STDERR "$genetrx does not exist in ANNOVAR database as $genetrx; try a synonym\n" and next ;
			}else{
				print STDERR "Found an alias for $genetrx ($$alias_ref{$genetrx})\n";
				print STDERR "Not really\n" and next if (!exists $gene_ori{$$alias_ref{$genetrx}});
				$genetrx=$$alias_ref{$genetrx};

			}
		# }else{
		# 	print STDERR Dumper (${gene_ori{$genetrx}})." exists!\n";
		}
####################################################################
# Split into each transcript for processing
###################################################################
		my @singlets;my $origpos=$pos;
		if ($genetrx=~/^(NM\_|ENST)/){
			# print __LINE__.":$genetrx\n" and <STDIN> if ($verbose);
			@singlets=split ("\n",$gene_ori{$genetrx});
		}else{
			print __LINE__.":$genetrx\n" and <STDIN> if ($debug);
			# print __LINE__.":$gene_ori{$genetrx}\n" and <STDIN> if ($debug);
			foreach my $key (keys %{$gene_ori{$genetrx}}){
				# print "adding $key for $genetrx...\n" and <STDIN>;#;
				my @ids=split("\n",$gene_ori{$genetrx}{$key});
				foreach my $id(@ids){
					push (@singlets,$id);
				}
			}
		}
####################################################################
# Declare variables needed for the next section
###################################################################		
		$debug=0;#global variable being set if you wish to debug section
		my $split_cdn=-1;#keeps track if triplet spans 2 exons!
		my @whichEx;#keeps track of whether a span is across exons
		my $nts="[ATCG]+";
		my $prots="[ABCDEFGHIKLMNPWRSTVWXYZ\*]";
####################################################################
# For each transcript for gene protein variant input, otherwise this only runs /once
# Transcript or protein is more specific and therefore better, but we allow for gene as well
#################################################################
		my %found_gpos;
		SINGLET: for (my $singlet_idx=0;$singlet_idx<=$#singlets;$singlet_idx++){
			$pos=$origpos;@whichEx=();#reset variables $origpos is the aapos from user input
			my ($trx,$chr,$strand,$cds_leftmost,$cds_rightmost,$ex_lefts,$ex_rights,$cds_offset,$fiveprimerutr)=split("\t",$singlets[$singlet_idx]);## this is extracted from the ANNOVAR genes file
			print __LINE__. "TESTING ($singlet_idx) $singlets[$singlet_idx]\n" and <STDIN> if ($debug);
			my $addon='';
			if ($pos=~/(\d+)([\-\+]\d+)/){#indels
				$pos=$1;$addon=$2;
			}
			my @lefts_ex=split(",",$ex_lefts);
			my @right_ex=split(",",$ex_rights);
			print __LINE__.":($singlet_idx]working on $genetrx,$trx,$pos,$ref,$var($cds_offset)\n\t$singlets[$singlet_idx]" if ($debug);
			#find cds offset
			my $first_cds_ex=-1;my $tally=0;my $last_cds_ex=-1;
			
			if (!$cds_offset){##this is for the first time we encounter the transcript and the cds offset is not set
				$cds_offset='';#first clear so we can append to it
				if ($strand eq '-'){
					$first_cds_ex=$#right_ex;
	# print "start with $first_cds_ex\n";
					print "$right_ex[$first_cds_ex]<=$cds_rightmost && $lefts_ex[$first_cds_ex]>=$cds_rightmost\n$singlets[$singlet_idx]\n" if $debug;
					until (($lefts_ex[$first_cds_ex]<=$cds_rightmost && $right_ex[$first_cds_ex]>=$cds_rightmost)|| $first_cds_ex==0){
						$fiveprimerutr+=abs($lefts_ex[$first_cds_ex]-$right_ex[$first_cds_ex]);
						$first_cds_ex--;
	# print "Decrement $first_cds_ex\n";<STDIN>;
					}
					# die "here looking for the correct fiveprimerutr for KRAS in viz5213b4dd37988-dev\n";
					$fiveprimerutr+=abs($right_ex[$first_cds_ex]-$cds_rightmost) ;#if ($first_cds_ex==$#right_ex);#do not re-add the one you already have in the until loop
					# print "Adding  $fiveprimerutr\n" and <STDIN> if ($debug);
					$last_cds_ex=0;
					until ($lefts_ex[$last_cds_ex]<=$cds_leftmost && $right_ex[$last_cds_ex]>=$cds_leftmost){
						$last_cds_ex++;
					}
					print "got first coding exon from the right:$first_cds_ex(out of $#right_ex)\n" if $debug;
					die "You cannot exceed this $singlets[$singlet_idx]\n" if ($first_cds_ex<0);
					$tally=0;my $last=$cds_rightmost;
###############################################################
# Now calculate the CDS end for each of the exons (-ve strand)
###############################################################		
					for (my $k=$#lefts_ex;$k>=0;$k--){
						#find all cds offsets
						if ($first_cds_ex<$k || $cds_leftmost >$right_ex[$k]){
							$cds_offset="0,$cds_offset";
						}elsif ($k==$#lefts_ex ||  $cds_leftmost <=$right_ex[$k]){
							print "($k)$last-$lefts_ex[$k]" if $debug;
							my $end=($lefts_ex[$k]<=$cds_leftmost)?$cds_leftmost:$lefts_ex[$k];
							# print "\t$last-$end+$tally\n" if $debug;#calculating cds offset at each exon
							$cds_offset=($last-$end+$tally).",$cds_offset";
							$tally=($last-$end+$tally);
							$last=$right_ex[$k-1];
						}else{
							warn "[COULD NOT PARSE] ($k) $cds_leftmost >$right_ex[$k]\n";
						}
					}
					print __LINE__.":OFFSETs: $cds_offset\n" if $debug;
					$cds_offset=~s/,$//;
				}else{
					#find the cds positions
					$first_cds_ex=0;
					print __LINE__.":($first_cds_ex)$lefts_ex[$first_cds_ex]<$cds_leftmost && $right_ex[$first_cds_ex]>$cds_leftmost\n$singlets[$singlet_idx]" and <STDIN> if $debug;
					until (($lefts_ex[$first_cds_ex]<=$cds_leftmost && $right_ex[$first_cds_ex]>=$cds_leftmost)|| $first_cds_ex>=$#lefts_ex){
						$fiveprimerutr+=abs($lefts_ex[$first_cds_ex]-$right_ex[$first_cds_ex]);
						print "\nAdding 5' $fiveprimerutr $lefts_ex[$first_cds_ex]-$right_ex[$first_cds_ex])\n"  if ($debug);
						print "\t($first_cds_ex)\n" if ($debug);
						$first_cds_ex++;
						print "\tAfter?($first_cds_ex)\n" if ($debug);
					}
					print "$lefts_ex[$first_cds_ex]-$cds_leftmost\n" if ($debug);
					$fiveprimerutr+=abs($lefts_ex[$first_cds_ex]-$cds_leftmost) ;#do not re-add the one you already have in the until loop
					print __LINE__."$first_cds_ex.Final length:5' UTR $fiveprimerutr\n"  if ($debug);
					print __LINE__."got $first_cds_ex\n" and <STDIN> if $debug;
					die "You cannot exceed this $singlets[$singlet_idx]\n" if ($first_cds_ex>$#lefts_ex);
					my $last=$cds_leftmost;$tally=0;
	###############################################################
	# Now calculate the CDS end for each of the exons (+ve strand)
	###############################################################				
					for (my $k=0;$k<=$#lefts_ex;$k++){
						#find all cds offsets
						if ($k<$first_cds_ex ){
							print __LINE__.":($k)$k<$first_cds_ex\n" if $debug;
							$cds_offset.="0,";
						# }elsif ($k==$first_cds_ex){
						# 	print __LINE__.":($k)$k==$first_cds_ex\n" if $debug;
						# 	my $end=($right_ex[$k]-$cds_leftmost+1);
						# 	$cds_offset.="$end,";
						# 	$tally+=$end;
						# 	print __LINE__.":$end($tally)\n" if $debug;
						}elsif ($k==$#lefts_ex || $right_ex[$k]<=$cds_rightmost){
							print __LINE__.":($k)$right_ex[$k]-$last\n" if $debug;
							my $end=($right_ex[$k]>$cds_rightmost)?$cds_rightmost:$right_ex[$k];
							$cds_offset.=($end-$last+$tally).",";
							$tally=($end-$last+$tally);
							$last=$lefts_ex[$k+1];
						}elsif($right_ex[$k]>=$cds_rightmost && $lefts_ex[$k]<=$cds_rightmost){
							print __LINE__.":$right_ex[$k]>=$cds_rightmost && $lefts_ex[$k]<=$cds_rightmost\n" if $debug;
							my $end=($right_ex[$k]>$cds_rightmost)?$cds_rightmost:$right_ex[$k];
							$cds_offset.=($end-$last+$tally).",";
							$tally=($end-$last+$tally);
						}elsif ($lefts_ex[$k]>$cds_rightmost ){
							print __LINE__.":$lefts_ex[$k]>$cds_rightmost\n" if ($debug);
							$cds_offset.="0,";
						}else{
							warn "=================[COULD NOT PARSE] ($k) $right_ex[$k]>$cds_rightmost\n";
						}
					}
				}
			}
			if ($singlets[$singlet_idx]=~/\t0\t0$/){
				$singlets[$singlet_idx]=~s/\t0\t0$//;
				$singlets[$singlet_idx].="\t$cds_offset\t$fiveprimerutr";
				print "New $singlets[$singlet_idx]\n"  if ($debug);
				if ($genetrx=~/^(NM_|ENST)/){
					$gene_ori{$genetrx}=join ("\n",@singlets);
				}else{
					$gene_ori{$genetrx}{$trx}=$singlets[$singlet_idx];
				}
				# print "Just added new singlets:\n". Dumper ($gene_ori{$genetrx}) and <STDIN> if ($debug);
			}else{
				print "Calculated or found $singlets[$singlet_idx]\n" and <STDIN> if ($debug);
			}
			my @offsets=split(",",$cds_offset);
			print "(CVRTMETH$cvrtmeth):$pos,$ref,$var($hgvs) $line\n"  if ($debug);
# print STDERR "Skipping $line(TRX:$trx,TALLY:$tally,POS:$pos) because of tally\n" and next if($tally>0 && $pos>$tally);
			my $prot2cds='';
			if ($type eq 'p'){
				$prot2cds=$pos;
				$pos=($pos*3);my $l=($pos%3);
				if ($l==0){#last in triplet
					$pos=($pos-2)."-".($pos);
				}elsif ($l==1){#middle
					$pos=($pos-1)."-".($pos+1);
				}elsif ($l==2){#first
					$pos="$pos-".($pos+2);
				}else{
					die "what? this should be a triplet!\n";
				}
				_define()  if (keys %matches == 0);
			}
			my $gpos;
			my ($p,$p2)=split("-",$pos);
			$p2||=$p;
			# $pos=$p;
			print __LINE__.":now working on". join (",",$genetrx,$pos,$var,$ref,$hgvs,$prot2cds,"STRAND:$strand")."\n\t($cds_offset)\n"  and <STDIN> if ($debug);
			if ($strand eq '-'){
				if ($first_cds_ex<0){#we do this because it is not always set due to the storage of data
					$first_cds_ex=$#lefts_ex;#print "\n\tstarts with $first_cds_ex\n";
					until ($offsets[$first_cds_ex]!=0){
						$first_cds_ex--;
					}
					$last_cds_ex=0;
					until ($lefts_ex[$last_cds_ex]<=$cds_leftmost && $right_ex[$last_cds_ex]>=$cds_leftmost){
						$last_cds_ex++;
					}
				}else{
					print __LINE__.":FIRST AND LAST CDS EX: $first_cds_ex and $last_cds_ex\n" if ($debug);
				}
				$gpos=0;my $padding=0;
				foreach my $cur_pos ($p, $p2){
					next if (!$cur_pos);
					my $found=0;
					print __LINE__."\t\t\t\tworking on $cur_pos: FOUND $found($first_cds_ex)\n" if ($debug);
					if ($#lefts_ex>0 && $cur_pos>=$offsets[$first_cds_ex]){
						print __LINE__.":$cur_pos>=$offsets[$first_cds_ex]??\n" if ($debug);
						for (my $i=0;$i<=$#lefts_ex-1;$i++){
							$padding=($i==$first_cds_ex)?abs($right_ex[$i]-$cds_rightmost):0;
							if ($cur_pos==$offsets[$first_cds_ex]){
								
								$gpos=$lefts_ex[$first_cds_ex];
								$found++;print __LINE__.":Looking for gpos($gpos)\n" if $debug;
							}elsif ($i+1<=$#offsets && $offsets[$i+1]==0 && $offsets[$i]!=0){
								if ($i==$first_cds_ex){
									$gpos=$cds_rightmost-$cur_pos+2;
								}else{
									$gpos=$lefts_ex[$i]+$padding + ($cur_pos-$offsets[$i-1]);
								}
								print __LINE__.":Looking for gpos\n" if $debug;
								print __LINE__."\tDefault $i th from the left $lefts_ex[$i]+($cur_pos-$offsets[$i-1])\n" if ($debug);
								$found++;
								
							}elsif ($cur_pos<=$offsets[$i] && $cur_pos>$offsets[$i+1]){
								print __LINE__."\t$i th from the left is the exon $right_ex[$i]-(abs($cur_pos-$offsets[$i+1])+1)\n" if ($debug);
								$gpos=$right_ex[$i]+1-(abs($cur_pos-$offsets[$i+1])-1);
								print __LINE__.":(Success!)Found a position for $cur_pos:($gpos)!\n" if $debug;
								$found++;
							}elsif ($cur_pos==$offsets[$i]){
								print __LINE__.":Looking for gpos\n" if $debug;
								print "($i)right on $cur_pos and $offsets[$i+1].....$last_cds_ex\n" if ($debug);
								$gpos=($i+1==$last_cds_ex)?($cds_leftmost+1):$lefts_ex[$i+1]+1;
								$found++;
							}else{
								print __LINE__.":Looking for gpos(none) $cur_pos is not between $offsets[$i] and $offsets[$i+1]\n" if $debug;
							}
							# print  "Adding $i to array($found)\n" if ($found);
							push(@whichEx,"$i") and last if ($found);
						}
					}elsif($offsets[$first_cds_ex]>=$cur_pos){
						print __LINE__.":Looking for gpos\n" if $debug;
						$gpos=$cds_rightmost-$cur_pos+2;
						push(@whichEx,"$first_cds_ex") ;
						print "\tFound first coding exon\n" and <STDIN> if ($debug);
					}else{
						$padding=($cds_leftmost-$lefts_ex[0]);
						$gpos=$cds_rightmost-$cur_pos-1;
						print "\tOnly one exon?$lefts_ex[0]+$padding + $cur_pos;\n" and <STDIN> if ($debug);
						push(@whichEx,0);
					}
				}
				die "ERR" if ($#whichEx>1);
				print STDERR __LINE__."skipping trx $trx b/c it is too short for c.$p/$p2 ($line)|TRX". length($seqhash{$trx})."{". join ("\n\t==>|",@whichEx). "}|\n" and $failedresponse.="\t($trx)\tCDS Offsets ($cds_offset) is too short for calculated mrna position $p-$p2\n" and  next if ($p2 && $#whichEx<1);
#########################################################################
# This lets me know which triplet belongs to which exon for genomic position
# Calculation (-ve strand)
##########################################################################

				print "$pos is what comes out of that last block....\n" if $debug;
				if ($p2 && $whichEx[0]!=$whichEx[1]){
					# print STDERR "Skipping mapping to multiple exons for now  $whichEx[0]!=$whichEx[1])$line\n";next SINGLET;
					print "waiting $p,$p2 and $offsets[$whichEx[0]]\n" and <STDIN> if ($debug);
					if (($p)==$offsets[$whichEx[0]]){
						$split_cdn=1;
						print "The first/second codon is the split: $p==$offsets[$whichEx[0]]\n" and <STDIN> if ($debug);
					}else{#die "Testing split codons at LN".__LINE__;
						print "The second/third codon is the split\n" and <STDIN> if ($debug);
						$split_cdn=2;
					}
					print "keep track of where the split between two exons occur: $split_cdn\n" and <STDIN> if $debug;
				}
				print "SplitCodon?? $split_cdn\n" if ($debug);
				if ($addon){
					print "starting with $gpos\t" and <STDIN> if ($debug);
					$gpos-=$addon;die 'havent coded';
					print "end with $gpos\n" and <STDIN> if ($debug);
				}
				print "The genomic position is ($gpos) [#trx: $#singlets and my current:$singlet_idx] $line\n"  and <STDIN> if ($debug);#and <STDIN> if ($gpos!=$field[2]);#test for later
			}else{#positive strand
				if ($first_cds_ex<0){#we do this because it is not always set due to the storage of data
					$first_cds_ex=0;
					# print "offsets: $first_cds_ex\n " and <STDIN> if ($verbose);
					my $exit_flag=0;
					until ($exit_flag || $offsets[$first_cds_ex]!=0 ){
						$first_cds_ex++;
						# print "\t$first_cds_ex...$offsets[$first_cds_ex]:". $#offsets. "\n" if ($verbose);
						if ($first_cds_ex>=$#offsets){
							$exit_flag=1;
						}
					}
					next if ($exit_flag);
				}
				print __LINE__.":$first_cds_ex,$p,$p2\n" and <STDIN> if ($debug);
				my $padding=0;
#################################################################
# This gets the exon for which the protein position is between
# Makes sure that the SNP is in the same exon so when we calculate 
# t1,t2,t3, we know that it is a matter of adding to the beginning 
# triplet position (+ve)
################################################################
				$gpos=0;
				foreach my $cur_pos ($p, $p2){
					next if (!$p2);
					
					print __LINE__. "\nworking on one of two $cur_pos($p,$p2)\n"  and <STDIN> if ($debug);#.join ("\n\t",@offsets);
					my $found=0;
					if ($#lefts_ex>0 && $cur_pos>$offsets[0]){
						for (my $i=1;$i<=$#lefts_ex;$i++){
							$padding=($i==$first_cds_ex)?abs($lefts_ex[$i]-$cds_leftmost):0;
							if ($cur_pos>=$offsets[$i-1] && $cur_pos<=$offsets[$i]){
								print "\nTRUE: $cur_pos>=$offsets[$i-1] && $cur_pos<=$offsets[$i]\n" if $debug;
								 if (!$gpos){
								 	$gpos=$lefts_ex[$i]+$padding + abs($cur_pos-$offsets[$i-1]) ;
								 	print "$i is the exon $lefts_ex[$i]+$padding + abs($cur_pos-$offsets[$i-1])\n" and <STDIN> if ($debug);
							 	}else{
							 		print __LINE__. ":already have a gpos $gpos!\n" if $debug;
							 	}
							 	$found++;
							}elsif ($i+1<=$#offsets && $offsets[$i+1]==0 && $offsets[$i]!=0){
								print "\nTRUE: $i+1<=$#offsets && $offsets[$i+1]==0 && $offsets[$i]!=0\n" if $debug;
								if (!$gpos){
									print "Default $i $lefts_ex[$i]+($cur_pos-$offsets[$i-1])\n" and <STDIN> if ($debug);
									$gpos=$lefts_ex[$i]+$padding + abs($cur_pos-$offsets[$i-1]) ;
								}else{
									print __LINE__. ":already found gpos! $gpos\n" if $debug;
								}
								$found++;
							}
							push(@whichEx,"$i") and last if ($found);
						}
					}else{# ($cur_pos<=$offsets[0]){
						$padding=($cds_leftmost-$lefts_ex[0]);
						$gpos=$lefts_ex[0]+$padding + $cur_pos if (!$gpos);
						push(@whichEx,0);
						$found++;#fatcat
						print "GPOS: $lefts_ex[0]+$padding + $cur_pos-1 = $gpos\n" if $debug;
						print "CDS smaller than first exon or only one exon!$lefts_ex[0]+$padding + $cur_pos;\n" and <STDIN> if ($debug);
					# }else{
					# 	$padding=($cds_leftmost-$lefts_ex[0]);
					# 	$gpos=$lefts_ex[0]+$padding + $cur_pos;
					# 	push(@whichEx,0);
					# 	print "Only one exon?$lefts_ex[0]+$padding + $cur_pos;\n" and <STDIN> if ($debug);
					}
				}
				print Dumper (\@whichEx) and <STDIN> if $debug;
				print STDERR __LINE__."skipping trx $trx b/c it is too short for p.$p/$p2 ($line)===>" . join ("\n\t==>",@whichEx). "\n" and next SINGLET if ($p2 && $#whichEx<1);
#########################################################################
# This lets me know which triplet belongs to which exon for genomic position
# Calculation (+ve strand)
##########################################################################
				if ($p2 && $whichEx[0]!=$whichEx[1]){
					# print STDERR "Skipping mapping to multiple exons for now  $whichEx[0]!=$whichEx[1]) ($line)\n";next SINGLET;
					#find the codon position that is the odd man out
					if ($p==$offsets[$whichEx[0]]){
						$split_cdn=0;
						print "$p==$offsets[$whichEx[0]]\n" and <STDIN> if ($debug);
					}elsif ($p+1==$offsets[$whichEx[0]]){
						$split_cdn=1;
						print "$p+1==$offsets[$whichEx[0]]\n" and <STDIN> if ($debug);
					}else{die "Testing split codons at LN".__LINE__;
						$split_cdn=2;
					}
					print "keep track of where the split between two exons occur: $split_cdn\n" and <STDIN> if $debug;
				# }elsif (!$p2){
				# 	$gpos.="\t$gpos";
				}else{
				}
				print "SplitCodon?? $split_cdn\n" if ($debug);

				#now if this is protein coordinates; need to figure out alleles
				
				if ($addon){
					print "starting with $gpos\t" and <STDIN> if ($debug);
					$gpos+=$addon;
					print "end with $gpos\n" and <STDIN> if ($debug);
				}
				print __LINE__."The genomic position of the first triplet is:($gpos) [#trx: " . ($#singlets+1)." and my current:$singlet_idx] $line\n"  and <STDIN> if ($debug);
			}#end pos strand
			my @group=();#keeps track of the different variants possible
			if ($type eq 'p'){
				print "Testing...Distance from exon start to 5'UTR:$fiveprimerutr\n" and <STDIN> if ($debug);
				
				my ($refp,$varp)=(uc($ref),uc($var));
				if (length($refp)>2){
					$refp=substr($refp,0,1).lc(substr($refp,1,length($refp)-1));
					$varp=substr($varp,0,1).lc(substr($varp,1,length($varp)-1));
					if (exists $VarAbbr{$refp}){
						$refp=$VarAbbr{$refp};
					}
					if (exists $VarAbbr{$varp}){
						$varp=$VarAbbr{$varp};
					}
				}
				print STDERR "$trx does not exist in seqhash\n" and next SINGLET if (!exists $seqhash{$trx});
				if ( ($p+$fiveprimerutr+2) > length($seqhash{$trx})){ #sometimes you have multiple transcripts for a gene and one of the transcripts is shorter than the protein position being interrogated
					
					if ($ref ne 'X' && $ref ne '*'){
						print "$p+$fiveprimerutr+3 is greater than the length of the $trx...hoping this is multi-transcipt\n" and <STDIN> if $debug;
						next SINGLET ;
					}else{
						print "$p+$fiveprimerutr+3 is greater than the length of the $trx...but this is a stoploss!\n" and <STDIN> if $debug;
					}
				}
				my $ref_triplet=substr($seqhash{$trx},$p+$fiveprimerutr-1,3);
				print "RefTriplet:$ref_triplet... got from $p+$fiveprimerutr-1\n" if $debug;
				# print $refp;exit;
				# print "Reftriplet? $ref_triplet\n" .Dumper $matches{$refp}{orig}."\n" and <STDIN> and next SINGLET 
				next if (!$matches{$refp}{'orig'});
				if ($matches{$refp}{'orig'}!~/$ref_triplet/i){
					$failedresponse.="\tRefTriplet($ref_triplet) in mrna ($trx) did not match any user's Input reference protein codon ($matches{$refp}{'orig'})\n";
					# warn "$matches{$refp}{'orig'} does not match the mrna matched sequence triplet $ref_triplet\n";next;
					next;
				}
				#can A turn into B?
				my $regex;my $cdnnum;
################################################################
# Here we try all permutations by triplet of the possible aa changes
#################################################################
				for (my $i=0;$i<=2;$i++){
					my $t;
					if ($i==0){
						$t='\\w'.substr($ref_triplet,1,2);
					}elsif ($i==1){
						$t=substr($ref_triplet,0,1)."\\w".substr($ref_triplet,2,1);
					}else{
						$t=substr($ref_triplet,0,2)."\\w";
					}
					my @matches;#holds all of the aa1=>aa2 conversions that are feasible e.g AAA cannot get to CCC in s single nt change!!!
					if (!exists $matches{$refp}{$varp}){
						$varp='orig';
					}
					print __LINE__.":ORIGINAL P? $gpos\n" and <STDIN> if $debug;
					if (exists $matches{$refp}{$varp}){
						print __LINE__.":working on $t $matches{$refp}{$varp}\n" if ($debug);
						push (@matches,$1) while($matches{$refp}{$varp}=~/($t)/g);
						if ($#matches>-1){
							print __LINE__.":MATCHED CODONS\n\t". join ("\n\t",@matches)."!!!\n" if $debug;
							foreach my $match(@matches){
								print __LINE__. ":current codon:$match ($i)\n" if $debug;
								my $tmpPos=$gpos;
								my $newvar=substr($match,$i,1);
								my $cur_ref=substr($ref_triplet,$i,1);
								# next if ($newvar eq $cur_ref);#do not print synonymous G-G 
								my $gen;
								if ($strand eq '-'){
									print __LINE__.":working on current $gpos want to use?$lefts_ex[$whichEx[0]] vs $lefts_ex[$whichEx[1]]\n" and <STDIN> if $debug;
									$newvar=revcompl($newvar);
									$cur_ref=revcompl($cur_ref);
									print "____________EVAL:$split_cdn,$i??\n" and <STDIN> if ($debug);#2,0 tested correct
									if ($split_cdn>-1 && $i<$split_cdn){
										if ($i==0 ){
											print __LINE__.":prev $gpos\t"  if ($debug);
											$tmpPos=($lefts_ex[$whichEx[0]]+1);#checked
											# $gpos=$lefts_ex[$whichEx[1]] ;
											print "new $gpos\t" and <STDIN> if ($debug);
										}elsif ($i==1 ){
											print __LINE__.":prev $gpos\t"  if ($debug);##works
											$tmpPos=$lefts_ex[$whichEx[0]] +1;
											print "new $gpos\t" and <STDIN> if ($debug);
										}elsif ($i==2 ){
											print __LINE__.":prev $gpos\t"  if ($debug);
											$tmpPos=$lefts_ex[$whichEx[1]] +2;
											print "new $gpos\t" and <STDIN> if ($debug);
										# }else{
										# 	$gpos=0;
										}else{
											$$tmpPos='';
										}
										# print "changed to $gpos....$lefts_ex[$whichEx[0]]+1,CDNID:$i>SPLIT:$split_cdn\n" and <STDIN> if ($debug);
									}else{
										$tmpPos=($tmpPos+(1-$i));
									}
									# $gen=join ("\t",$chr,($gpos+(1-$i)),($gpos+(1-$i)),$cur_ref,$newvar);
									$gen=join ("\t",$chr,$tmpPos,$tmpPos,$cur_ref,$newvar);
								}else{#+ve
									print "____________\n" and <STDIN> if ($debug);
									if ($split_cdn>-1 && $i>$split_cdn){
										print "prev $gpos\n" and <STDIN> if ($debug);
										if ($i==1 ){
											$tmpPos=($lefts_ex[$whichEx[1]] + 1+$split_cdn);#left coordinate needs to add one
											# $gpos=$lefts_ex[$whichEx[1]] ;
										}elsif ($i==2 && $split_cdn){#checked
											$tmpPos=$lefts_ex[$whichEx[1]] +1;
										}elsif ($i==2 && !$split_cdn){
											$tmpPos=$lefts_ex[$whichEx[1]] +2;
										# }else{
										# 	$gpos=0;
										}
										print "changed to $tmpPos....$lefts_ex[$whichEx[1]]+1,CDNID:$i>SPLIT:$split_cdn\n" and <STDIN> if ($debug);
										$gen=join ("\t",$chr,$tmpPos,$tmpPos,$cur_ref,"$newvar");
									}else{
										print __LINE__. "new pos?($i)$chr,($tmpPos+$i),($tmpPos+$i),$cur_ref,$newvar\n" if ($debug);
									 	$gen=join ("\t",$chr,($tmpPos+$i),($tmpPos+$i),$cur_ref,"$newvar");
									 }
								}
								if (!exists $found{$countline} || (exists $found{$countline} && $found{$countline}!~/$gen/)){
									if ($includeinfo){
										# print "\tLN:$countline\t$otherinfo";
										$gen.= "\t$otherinfo($ref_triplet:$strand)";
									}
									if (!exists $found_gpos{$gen}){
										print $gen;
										print "\n";
										<STDIN> if $debug;
									}
									$found_gpos{$gen}++;
									$found{$countline}.="$gen\n";
								}
							}
						}else{
							$failedresponse.= "\t$trx ## (Gene Orientation:$strand)could not match the protein change using SNP, may be MNP" if ($includeinfo);
							$failedresponse.= "\n";
						}
					}else{
						print __LINE__.":working on $t $refp->$varp\n" if ($debug);
						if (exists $matches{$refp}{'orig'}){
							print __LINE__.":working on $t $matches{$refp}{'orig'}\n" if ($debug);
						}
						print "This does not exist in the aa1=>aa2 hash conversion\n" and <STDIN> if ($debug);
					}
				}
			}else{
				my $gpos2;
				if ($strand eq '-'){
					$ref=revcompl($ref);
					$var=revcompl($var);
				}
				if($indeltype=~/del/){
					my $v=length($ref);
					if ($strand eq '+'){
						$gpos2=($gpos+$v-1);
					}else{
						$gpos2=$gpos;
						$gpos-=$v;
					}
				}else{
					$gpos2=$gpos;
				}
				if ($ref eq '0'|| $var eq '0'){
					# print FAIL "[FAILED]$line--------->".join ("\t","$chr",$gpos,$gpos2,$ref,$var). "\n";
				}else{
					print join ("\t","$chr",$gpos,$gpos2,$ref,$var);
					$found{$countline}++;
					if ($includeinfo){
						print "\tLN#$countline:$otherinfo";
					}
					print "\n";<STDIN> if $debug;
				}
			}	
		}#ends foreach singlet			
	# 	$countvar++;

		print FAIL "[FAILED] $line\tLN:$countline\n" and $failedresponse='' if (!exists $found{$countline});
		print "\t\t_______________NEXT VAR______________________\n" if $debug;
	}#ends foreach var
	close VAR;close FAIL;
	# my $triallelic = $countsnp-$countti-$counttv;
	# print STDERR "NOTICE: Completed at ". `date`;
}
sub readSeqFromFASTADB {
	my ($dbtype1) = @_;
	
	#the four statements below should be condensed in the future (they are identical)
	my $seqdbfile = File::Spec->catfile('/SeqIdx/annovardb/humandb/hg19_'."$dbtype1" . "Mrna.fa");
	open (SEQ,"<$seqdbfile") or die "Cannot find $seqdbfile\n";
	my ($seqid, $curseq) = ('', '');
	while (<SEQ>) {
		if (m/^>(\S+)/ ) {
			$seqid = $1;
			$curseq = '';
			$seqhash{$seqid}='';
		# }elsif ($seqid && m/^>/){
		# 	close SEQ;
		# 	return (\%seqhash);
		} else {
			s/[\r\n]+$//;
			$seqhash{$seqid} .= uc $_;					#only use upper case characters
		}

	}
	# print STDERR "Reading all of FASTA\t:";
	# getMemory() ;
	# die "Could not find $refseq in annotation file\n";
	close (SEQ);
}
sub convertOnePos{#added for JasonLih 2012/12/11
	my ($variantfile) = @_;
	my ($countline, $countvar, $counthom, $counthet, $countunknown, $countindel, $countsnp, $countti, $counttv) = qw/0 0 0 0 0 0 0 0 0/;#ANNOVAR default
	my %stats;#hv added 2/10/12
	my ($source_program, $gtpos);		#the program that generated the VCF4 file; the GT position within FORMAT record
	open (VAR, $variantfile) or die "Error: cannot read from variant file $variantfile: $!\n";
	VAR: while (<VAR>) {
		my $bpcount=0;
		$countline++;
		m/^#/ and next;		#skip comment lines
		s/[\r\n]+$//;		#delete trailing new lines
	
		my @field=split(/\t/,$_);
		my ($chr, $start, $ref_allele, $mut_allele, $quality_score, $ID, @sample) = @field;
		if ($start!~/^\d+/ and $chr=~/Chrom{0,1}$/i){
			print "#$_\n";$af=0;
			next;
		}
		my $otherinfo = join("\t",$quality_score,$ID,@sample);	#this is the complete line (when -includeinfo is set, the entire line will be included in output file)
		my ($end);
		my (@mut_alleles, $zygosity);
		
		#sometimes the alleles are not in the same case
		#chr1    1869771 1869774 actc    aCTctc          43.5    13      INDEL;DP=13;AF1=0.5;CI95=0.5,0.5;DP4=0,4,4,0;MQ=37;PV4=0.029,0.45,1,0.46
		$ref_allele = uc $ref_allele;
		$mut_allele = uc $mut_allele;
		
		
		if ($mut_allele eq '.') {			#no variant call was made at this position
			$stats{'unknown'}++;
			next VAR;
		}elsif ($mut_allele=~/^\d+$/){#deletion; added for ANNOVAR 1000G
			$bpcount=$mut_allele-1;
			$mut_allele='-';
		}elsif ($mut_allele=~/^0([ATCG]*)/i){#insertion; added for ANNOVAR 1000G
			$mut_allele=$1;$ref_allele='-';
		}elsif ($mut_allele=~/^\d+/){
			next VAR;
		}
		#hv changed the next block
		if ($mut_allele =~ m/,/ ) {
			if ($allallele){
				@mut_alleles=split(",",$mut_allele);
			}else{
				$mut_allele=~/^([ACGNT]*),/;
				push (@mut_alleles,$1);
			}
		}else{
			push(@mut_alleles,$mut_allele);
		}
		my $ct=($#mut_alleles+1);
		$stats{"$ct-allelic"}++;	
		#ends hv block
		
		foreach my $mymut (@mut_alleles){
			my ($mystart,$end,$myref);#define these inside loop because it may change per mutation
			
			
			($mystart,$myref)=($start,$ref_allele);$end=$mystart+$bpcount;
			my ($ref_len,$mut_len)=(length($myref),length($mymut));
			if (my $head=smallest_substring($myref,$mymut)){
				chomp $head;
				if ($head eq $mymut) {#deletion
						$mystart=$mystart+length($head);
						$myref=substr ($myref, length($head),length($myref)-length($head));
						$end=$mystart+length($myref)-1;
						$mymut="-";
				}elsif ($head eq $myref){#insertion
						$myref='-';
						$mymut=~s/$head//;
				}else{#snp or multiallelic snp
						$myref=~s/$head//;
						$mymut=~s/$head//;
						$mystart+=length($head);
						$end=$mystart+length($myref)-1;
						
				}
			}
			if( ($ref_len==$mut_len && $myref ne '-')) {  	### output snv
				##chr1    1887091 13018;13019     CG      C,CA    .    #abcc_hv
				if (defined $format) {
					my @format = split (/:/, $format);
					undef $gtpos;
					for my $i (0 .. @format-1) {
						if ($format[$i] eq 'GT') {
							$gtpos = $i;
							last;
						}
					}
				} else {
					$zygosity = 'unknown_gt';
					$countunknown++;
				}
	
				if ($myref=~/A$/ and $mymut=~/G$/ or $myref=~/G$/ and $mymut=~/A$/ or $myref=~/C$/ and $mymut=~/T$/ or $myref=~/T$/ and $mymut=~/C$/) {
					$countti++;
					$stats{'snp_ct_ti'}++;
				} else {
					$counttv++;
					$stats{'snp_ct_tv'}++;
				}
				if ($includeinfo) {
					print $chr, "\t", $mystart, "\t", $end, "\t", $myref, "\t", $mymut, "\t", $otherinfo, "\n";
				} elsif ($bare){
					print $chr, "\t", $mystart, "\t", $end, "\t", $myref, "\t", $mymut,  "\n";
				}else{
					print $chr, "\t", $mystart, "\t", $end, "\t", $myref, "\t", $mymut,   "\t", $quality_score, "\n";
				}
				#print $chr, "\t", $start, "\t", $start, "\t", $ref_allele, "\t", $mut_allele, "\t$zygosity",  "\t", $quality_score, (defined $unfiltered_read_depth)? "\t$unfiltered_read_depth" : '', (defined $MappingQuality) ? "\t$MappingQuality" : '', (defined $QualityByDepth) ? "\t$QualityByDepth" : '', $includeinfo ? "\t$otherinfo" : '', "\n";	#commented Sep 2011
				$countsnp++;
			} else{#}if ($myref eq '-' || $mymut eq '-') {  ### output indel
				$stats{'indel'}++;$stats{"indel_$ct"}++;
				print $chr, "\t", $mystart, "\t", $end, "\t", $myref, "\t", $mymut;
	
				if ($includeinfo) {
					 print $withzyg?"\t$zygosity":"", "\t", $otherinfo;
				}elsif ($bare){
					
				} else {
					print "\t", $zygosity;
					print "\t", $quality_score;
					
					#$includeinfo and print "\t", $otherinfo;	#commented Sep 2011
				}
				print "\n";
				$countindel++;
			}
		}
		$countvar++;
	}
	my $triallelic = $countsnp-$countti-$counttv;
	print STDERR "NOTICE: Read $countline lines and wrote ${\($counthet+$counthom)} different variants at $countvar genomic positions ($countsnp SNPs and $countindel indels)\n";
	print STDERR "NOTICE: Among ${\($counthet+$counthom+$countunknown)} different variants at $countvar positions, $counthet are heterozygotes, $counthom are homozygotes\n";
	print STDERR "NOTICE: Among $countsnp SNPs, $countti are transitions, $counttv are transversions", $triallelic?", $triallelic have more than 2 alleles\n":"\n";
	print STDERR Dumper (\%stats);
	
}



sub convertBambino{#added for JasonLih 2012/12/11
	my ($variantfile) = @_;
	my ($countline, $countvar, $counthom, $counthet, $countunknown, $countindel, $countsnp, $countti, $counttv) = qw/0 0 0 0 0 0 0 0 0/;#ANNOVAR default
	my %stats;#hv added 2/10/12
	my ($source_program, $gtpos);		#the program that generated the VCF4 file; the GT position within FORMAT record
	open (VAR, $variantfile) or die "Error: cannot read from variant file $variantfile: $!\n";
#NormalSample    TumorSample     Name    Chr     Pos     Type    Size    Coverage        Percent_alternative_allele      Chr_Allele      Alternative_Allele      Score   Text    unique_alternative_ids  reference_normal_count  reference_tumor_count   alternative_normal_count        alternative_tumor_count count_ref_normal_fwd    count_ref_normal_rev    count_ref_tumor_fwd     count_ref_tumor_rev     count_var_normal_fwd    count_var_normal_rev    count_var_tumor_fwd     count_var_tumor_rev     alternative_fwd_count   alternative_rev_count   alternative_bidirectional_confirmation  broad_coverage  broad_reference_normal_count    broad_reference_tumor_count     broad_alternative_normal_count  broad_alternative_tumor_count   dbSNP   unique_alt_read_starts  unique_alt_read_starts_fwd      unique_alt_read_starts_rev      avg_mapq_alternative    somatic_or_germline     loh_flag        alt_to_ref_ratio_normal alt_to_ref_ratio_tumor  alt_to_ref_ratio_diff_normal_tumor      strand_skew     annovar_region  annovar_region_gene     annovar_exonic_function annovar_exonic_function_gene    annovar_dbsnp_id        annovar_avsift_score    alt_n_freq      alt_t_freq      tSpvalue        call
#OV_0001_BL_b.trim.bwa.sort.rmdup.bam    OV_0001_FP.trim.bwa.sort.rmdup.bam      chr1.11769373   chr1    11769373        SNP     1       240     0.167   G       A       1.000   CCAGCTCCCTGATGAAGAAG[G/A]CAGAGCTCTCCGAAGCCCAG   31      115     85      0       40      48      67      40      45      0       0       22      18      22      18      1       261     130     88      0       43              34      19      18      60      S               0.000   0.328   0.328   0.591   exonic  C1orf187        nonsynonymous SNV       C1orf187:NM_198545:exon3:c.G493A:p.A165T,               0.46    0       0.32824427480916        3.620556e-15    Somatic
#OV_0001_BL_b.trim.bwa.sort.rmdup.bam    OV_0001_FP.trim.bwa.sort.rmdup.bam      chr1.44432453   chr1    44432453        SNP     1       371     0.197   T       C       1.000   CAAGGCCTCCTGTGGCTTCT[T/C]TGTGAGTCCCATGCTGAACC   59      162     136     0       73      60      102     52      84      0       0       30      43      30      43      1       405     184     139     0       82              47      25      34      59      S               0.000   0.371   0.371   0.529   exonic;splicing IPO13;IPO13     nonsynonymous SNV       IPO13:NM_014652:exon17:c.T2522C:p.F841S,                0       0       0.3710407239819 6.287957e-26    Somatic
#OV_0001_BL_b.trim.bwa.sort.rmdup.bam    OV_0001_FP.trim.bwa.sort.rmdup.bam      chr1.75078365   chr1    75078365        SNP     1       795     0.192   C       G       1.000   AAAGTAGCCTCGTTTGCCTC[C/G]AAGCCTGGAACCTTTCCGAT   112     411     231     0       153     207     204     107     124     0       0       78      75      78      75      1       860     450     249     0       161             72      49      52      59      S               0.000   0.393   0.393   0.517   exonic  C1orf173        nonsynonymous SNV       C1orf173:NM_001002912:exon9:c.G1129C:p.G377R,           0       0       0.392682926829268       1.870061e-61    Somatic

##Columns of interest are 3-5,9-10 Chr	Pos     Type	Chr_Allele      Alternative_Allele
	VAR: while (<VAR>) {
		my $bpcount=0;
		$countline++;
		m/^#/ and next;		#skip comment lines
		s/[\r\n]+$//;		#delete trailing new lines
		my @field=split(/\t/,$_);
		my ($norm,$tumor,$name,$chr, $start,$type,$size,$cvg,$perc_alt, $ref_allele, $mut_allele, @sample) = @field;
		if ($start!~/^\d+/ and $chr=~/Chr$/i){
			print "#$_\n";$af=0;#header
			next;
		}
		my $otherinfo = join("\t",@field);	#this is the complete line (when -includeinfo is set, the entire line will be included in output file)
		my ($end);
		my (@mut_alleles, $zygosity);
		
		#sometimes the alleles are not in the same case
		#chr1    1869771 1869774 actc    aCTctc          43.5    13      INDEL;DP=13;AF1=0.5;CI95=0.5,0.5;DP4=0,4,4,0;MQ=37;PV4=0.029,0.45,1,0.46
		$ref_allele = uc $ref_allele;
		$mut_allele = uc $mut_allele;
		
		
		if ($mut_allele eq '.' ) {			#no variant call was made at this position
			$stats{'unknown'}++;
			next VAR;
		}elsif ($type=~/del/){#deletion
			$bpcount=length($ref_allele)-1;
			$mut_allele='-';
		}elsif ($type=~/ins/){#insertion
			$ref_allele='-';
		}
		
		push(@mut_alleles,$mut_allele);
		my $ct=($#mut_alleles+1);
		$stats{"$ct-allelic"}++;	
		#ends hv block
		
		foreach my $mymut (@mut_alleles){
			my ($mystart,$end,$myref);#define these inside loop because it may change per mutation
			
			
			($mystart,$myref)=($start,$ref_allele);$end=$mystart+$bpcount;
			my ($ref_len,$mut_len)=(length($myref),length($mymut));
			if (my $head=smallest_substring($myref,$mymut)){
				chomp $head;
				if ($head eq $mymut) {#deletion
						$mystart=$mystart+length($head);
						$myref=substr ($myref, length($head),length($myref)-length($head));
						$end=$mystart+length($myref)-1;
						$mymut="-";
				}elsif ($head eq $myref){#insertion
						$myref='-';
						$mymut=~s/$head//;
				}else{#snp or multiallelic snp
						$myref=~s/$head//;
						$mymut=~s/$head//;
						$mystart+=length($head);
						$end=$mystart+length($myref)-1;
						
				}
			}
			if( ($ref_len==$mut_len && $myref ne '-')) {  	### output snv
				##chr1    1887091 13018;13019     CG      C,CA    .    #abcc_hv
				if (defined $format) {
					my @format = split (/:/, $format);
					undef $gtpos;
					for my $i (0 .. @format-1) {
						if ($format[$i] eq 'GT') {
							$gtpos = $i;
							last;
						}
					}
				} else {
					$zygosity = 'unknown_gt';
					$countunknown++;
				}
	
				if ($myref=~/A$/ and $mymut=~/G$/ or $myref=~/G$/ and $mymut=~/A$/ or $myref=~/C$/ and $mymut=~/T$/ or $myref=~/T$/ and $mymut=~/C$/) {
					$countti++;
					$stats{'snp_ct_ti'}++;
				} else {
					$counttv++;
					$stats{'snp_ct_tv'}++;
				}
				if ($includeinfo) {
					print $chr, "\t", $mystart, "\t", $end, "\t", $myref, "\t", $mymut, "\t", $otherinfo, "\n";
				} else{
					print $chr, "\t", $mystart, "\t", $end, "\t", $myref, "\t", $mymut,  "\n";
				}
				#print $chr, "\t", $start, "\t", $start, "\t", $ref_allele, "\t", $mut_allele, "\t$zygosity",  "\t", $quality_score, (defined $unfiltered_read_depth)? "\t$unfiltered_read_depth" : '', (defined $MappingQuality) ? "\t$MappingQuality" : '', (defined $QualityByDepth) ? "\t$QualityByDepth" : '', $includeinfo ? "\t$otherinfo" : '', "\n";	#commented Sep 2011
				$countsnp++;
			} else{#}if ($myref eq '-' || $mymut eq '-') {  ### output indel
				$stats{'indel'}++;$stats{"indel_$ct"}++;
				print $chr, "\t", $mystart, "\t", $end, "\t", $myref, "\t", $mymut;
	
				if ($includeinfo) {
					 print $withzyg?"\t$zygosity":"", "\t", $otherinfo;
				}
				print "\n";
				$countindel++;
			}
		}
		$countvar++;
	}
#	my $triallelic;# = $countsnp-$countti-$counttv;
#	print STDERR "NOTICE: Read $countline lines and wrote ${\($counthet+$counthom)} different variants at $countvar genomic positions ($countsnp SNPs and $countindel indels)\n";
#	print STDERR "NOTICE: Among ${\($counthet+$counthom+$countunknown)} different variants at $countvar positions, $counthet are heterozygotes, $counthom are homozygotes\n";
#	print STDERR "NOTICE: Among $countsnp SNPs, $countti are transitions, $counttv are transversions", $triallelic?", $triallelic have more than 2 alleles\n":"\n";
#	print STDERR Dumper (\%stats);
	
}#end convertBambino

sub convertVS2{#Added for Michelle M.
	my ($variantfile) = @_;
	
	my ($countline, $countvar, $counthom, $counthet, $countunknown, $countindel, $countsnp, $countti, $counttv, $countaf) = qw/0 0 0 0 0 0 0 0 0 0/;#ANNOVAR default
	my %stats;#hv added 2/10/12
	my ($source_program, $gtpos);		#the program that generated the VCF4 file; the GT position within FORMAT record
	open (VAR, $variantfile) or die "Error: cannot read from variant file $variantfile: $!\n";
	if ($af){
		print STDERR "Notice:  --af is not supported for this format --$format\n";
	}
	my $id;
	VAR: while (<VAR>) {
		$countline++;
		m/^##/ and next;		#skip comment lines
		s/[\r\n]+$//;		#delete trailing new lines
		
		my @field=split(/\t/,$_);
		@field >=5 or die "Error: invalid record found in VarScan2 file (at least 5 tab-delimited fields expected): <$_>\n";
		my ($chr, $start, $ref_allele, $mut_allele, @info) = @field;
		if ($start!~/^\d+/ and $chr=~/Chrom{0,1}$/i){
			if ($field[2]=~/ID/i){
				($chr,$start,$id,$ref_allele,$mut_allele,@info)=@field;
			}
			if ($chr=~/#/){print "$_\n";}else{
			print "#$_\n";}$af=0;
			next;
		}
		if ($id){
			($chr,$start,$id,$ref_allele,$mut_allele,@info)=@field;
		}
		if ($includeinfo){
			
		}
		my $af_idx;
		my ($end);
		my (@mut_alleles, $zygosity);
		my $otherinfo = join("\t",@field);	#this is the complete line (when -includeinfo is set, the entire line will be included in output file)
		#sometimes the alleles are not in the same case
		#chr1    1869771 1869774 actc    aCTctc          43.5    13      INDEL;DP=13;AF1=0.5;CI95=0.5,0.5;DP4=0,4,4,0;MQ=37;PV4=0.029,0.45,1,0.46
		$ref_allele = uc $ref_allele;
		$mut_allele = uc $mut_allele;
		
		if ($mut_allele eq '.') {			#no variant call was made at this position
			$stats{'unknown'}++;
			next VAR;
		}
		if ($mut_allele =~ m/,/ ) {
			if ($allallele){
				@mut_alleles=split(",",$mut_allele);
			}else{
				$mut_allele=~/^([ACGNT]*),/;
				push (@mut_alleles,$1);
			}
		}else{
			push(@mut_alleles,$mut_allele);
		}
		my $ct=($#mut_alleles+1);
		$stats{"$ct-allelic"}++;	
		#ends hv block
		
		foreach my $mymut (@mut_alleles){
			my ($mystart,$end,$myref,$mytype);#define these inside loop because it may change per mutation
			($mystart,$myref)=($start,$ref_allele);$end=$mystart;
			if ($mymut eq '-'){#indel ok
				$mytype='del';
			}elsif ($mymut =~/^\+([ATCG]+)/){#insertion explicit
				$mymut=$1;
				$mytype='ins';
			}elsif ($mymut=~/^\-([ATCG]+)/){#deletion explicit
				$mymut=$1;
				$mytype='del';
			}elsif (my $head=smallest_substring($myref,$mymut)){
				chomp $head;
				if ($head eq $mymut) {#deletion
						$mystart=$mystart+length($head);
						$myref=substr ($myref, length($head),length($myref)-length($head));
						$end=$mystart+length($myref)-1;
						$mymut="-";
						$mytype='del';
				}elsif ($head eq $myref){#insertion
						$myref='-';
						$mymut=~s/$head//;
						$mytype='ins';
				}else{#snp or multiallelic snp
						$myref=~s/$head//;
						$mymut=~s/$head//;
						$mystart+=length($head);
						$end=$mystart+length($myref)-1;
						$mytype='snp';
						
				}
				# print "FOUND:($head)".join (",",@field)."\n";<STDIN>;
			}else{
				$mytype='snp';
			}
			my ($ref_len,$mut_len)=(length($myref),length($mymut));
			# if ($mytype eq 'del') {#deletion
			# 		$mystart=$mystart+1;
			# 		$myref=$mymut;
			# 		$end=$mystart+length($myref)-1;
			# 		$mymut="-";
			# }elsif ($mytype eq 'ins'){#insertion
			# 		$mystart++;
			# 		$end++;
			# 		$myref='-';
			# }else{#snp 
					
			# }
			if( ($ref_len==$mut_len && $myref ne '-')) {  	### output snv
				##chr1    1887091 13018;13019     CG      C,CA    .    #abcc_hv
	
				if ($myref=~/A$/ and $mymut=~/G$/ or $myref=~/G$/ and $mymut=~/A$/ or $myref=~/C$/ and $mymut=~/T$/ or $myref=~/T$/ and $mymut=~/C$/) {
					$countti++;
					$stats{'snp_ct_ti'}++;
				} else {
					$counttv++;
					$stats{'snp_ct_tv'}++;
				}
				if ($includeinfo) {
					print $chr, "\t", $mystart, "\t", $end, "\t", $myref, "\t", $mymut,  "\t", $otherinfo, "\n";
				} else{
					print $chr, "\t", $mystart, "\t", $end, "\t", $myref, "\t", $mymut,  "\n";
				}
				#print $chr, "\t", $start, "\t", $start, "\t", $ref_allele, "\t", $mut_allele, "\t$zygosity",  "\t", $quality_score, (defined $unfiltered_read_depth)? "\t$unfiltered_read_depth" : '', (defined $MappingQuality) ? "\t$MappingQuality" : '', (defined $QualityByDepth) ? "\t$QualityByDepth" : '', $includeinfo ? "\t$otherinfo" : '', "\n";	#commented Sep 2011
				$countsnp++;
			} else{
				if ($includeinfo) {
					print $chr, "\t", $mystart, "\t", $end, "\t", $myref, "\t", $mymut,  "\t", $otherinfo, "\n";
				} else{
					print $chr, "\t", $mystart, "\t", $end, "\t", $myref, "\t", $mymut,  "\n";
				}
			}
		}
	
		$countvar++;
	}
}


sub _define{
	%cdn2prot=(
		  "AAA","K","AAC","N","AAG","K","AAT","N",
          "ACA","T","ACC","T","ACG","T","ACT","T",
          "AGA","R","AGC","S","AGG","R","AGT","S",
          "ATA","I","ATC","I","ATG","M","ATT","I",
          "CAA","Q","CAC","H","CAG","Q","CAT","H",
          "CCA","P","CCC","P","CCG","P","CCT","P",
          "CGA","R","CGC","R","CGG","R","CGT","R",
          "CTA","L","CTC","L","CTG","L","CTT","L",
          "GAA","E","GAC","D","GAG","E","GAT","D",
          "GCA","A","GCC","A","GCG","A","GCT","A",
          "GGA","G","GGC","G","GGG","G","GGT","G",
          "GTA","V","GTC","V","GTG","V","GTT","V",
          "TAA","X","TAC","Y","TAG","X","TAT","Y",
          "TCA","S","TCC","S","TCG","S","TCT","S",
          "TGA","X","TGC","C","TGG","W","TGT","C",
          "TTA","L","TTC","F","TTG","L","TTT","F",
          );
	%VarAbbr=(
		"Ala"=>'A',	"Cys"=>'C',	"Asp"=>'D',
		"Glu"=>"E",	"Phe"=>'F',	'Gly'=>'G',
		'His'=>'H',	'Ile'=>'I',	'Lys'=>'K',
		'Leu'=>'L',	'Met'=>'M',	'Asn'=>'N',
		'Pro'=>'P',	'Gln'=>'Q',	'Arg'=>'R',
		'Ser'=>"S",	'Thr'=>"T",	"Val"=>'V',
		'Trp'=>'W',	'Tyr'=>'Y','Stop'=>'X',
		'*'=>'X'
		);
	foreach my $triplet (keys %cdn2prot){#SNPs
	  my @triple=split("",$triplet);
	  for (my $i=0;$i<=2;$i++){
	     my $newcdn;
	     foreach my $letter ('A', 'T', 'C', 'G'){
	          next if ($letter eq "$triple[$i]");
	          if ($i==0){
	               $newcdn=$letter.substr($triplet,1,2);
	          }elsif ($i==1){
	               $newcdn=$triple[0].$letter.$triple[2];
	          }elsif ($i==2){
	               $newcdn=$triple[0].$triple[1].$letter;
	          }
	          $matches{$cdn2prot{$triplet}}{'orig'}='' if (!exists $matches{$cdn2prot{$triplet}}{'orig'});
	          $matches{$cdn2prot{$triplet}}{'orig'}.="$triplet," if ($matches{$cdn2prot{$triplet}}{'orig'}!~/$triplet,/);
	          $matches{$cdn2prot{$triplet}}{$cdn2prot{$newcdn}}.="$newcdn," ;
	     }
	  }
	}
	my $helperseq='ATGCATG';#000,011,110,111,101 
	#block substitution
	# print Dumper (\%matches);exit;
}






sub convertVCF4 {
	my ($variantfile) = @_;
	$debug=0;
	my ($countline, $countvar, $counthom, $counthet, $countunknown, $countindel, $countsnp, $countti, $counttv, $countaf) = qw/0 0 0 0 0 0 0 0 0 0/;#ANNOVAR default
	my %stats;#hv added 2/10/12
	my ($source_program, $gtpos);		#the program that generated the VCF4 file; the GT position within FORMAT record
	open (VAR, $variantfile) or die "Error: cannot read from variant file $variantfile: $!\n";
	if ($af){
		if (defined $outfile){
			if ($append){
				open (AF,">>$outfile.hg19_af") or die "Cannot open $outfile.hg19_af\n";
			}else{
				open (AF,">$outfile.hg19_af") or die "Cannot open $outfile.hg19_af\n";
			}
		}else{
			if ($append){
				open (AF,">>ANNOVAR.input.hg19_af") or die "Cannot open AF\n";
			}else{
				open (AF,">ANNOVAR.input.hg19_af") or die "Cannot open AF\n";
			}
			$outfile="ANNOVAR.input";
		}
	}
	if ($withzyg){
		if (defined $outfile){
			my $zygfile=$variantfile;
			# For AVIA
			$zygfile=~s/(.txt|.vcf)//;
			if ($append){
				open (ZYG,">>$outfile.hg19_zygosity") or die "Cannot open $outfile.hg19_zygosity\n";
			}else{
				open (ZYG,">$outfile.hg19_zygosity") or die "Cannot open $outfile.hg19_zygosity\n";
			}
		}else{
			if ($append){
				open (ZYG,">>ANNOVAR.input.hg19_zygosity") or die "Cannot open ZYG\n";
			}else{
				open (ZYG,">ANNOVAR.input.hg19_zygosity") or die "Cannot open ZYG\n";
			}
			$outfile="ANNOVAR.input";
		}
	}
	VAR: while (<VAR>) {
		$countline++;
		print "-----------------working on variant on line $countline\n" if ($debug);
		if (m/^##fileformat=VCFv(\d+\.)/) {
			$1<4 and print STDERR "ERROR: Input file is not in VCF version 4 format but is $_" and exit;
		}
		if (m/##source=VarScan2/){
			$cmdlineargs=~s/vcf4/varscan2/;
			exec ("$0 $cmdlineargs\n");
			exit;
		}
		if (m/^##UnifiedGenotyper/) {
			$source_program = 'gatksnp';
			print STDERR "NOTICE: Detected that the VCF4 file is generated by GATK UnifiedGenotyper\n";
			$includeinfo or print STDERR "NOTICE: column 6-10 represent heterozygosity status, quality score, read depth, RMS mapping quality, quality by depth\n";
			$fraction and print STDERR "WARNING: the --fraction argument will be ignored for GATK SNP calls!!!\n";
			$confraction and print STDERR "WARNING: the --confraction argument will be ignored for GATK SNP calls!!!\n";
		}
		if (m/^##IndelGenotyper/) {
			$source_program = 'gatkindel';
			print STDERR "NOTICE: Detected that the VCF4 file is generated by GATK IndelGenotyper\n";
			$includeinfo or print STDERR "NOTICE: column 6-10 represent heterozygosity status, quality score, read depth, read count supporting indel call, RMS mapping quality\n";
		}
		
		if (not m/^#/ and not $source_program) {	#finished reading header line but did not detect the source program
			$includeinfo or print STDERR "NOTICE: for SNPs, column 6 and beyond MAY BE heterozygosity status, quality score, read depth, RMS mapping quality, quality by depth, if these information can be recognized automatically\n";
			$includeinfo or print STDERR "NOTICE: for indels, column 6 and beyond MAY BE heterozygosity status, quality score, read depth, read count supporting indel call, RMS mapping quality, if these information can be recognized automatically\n";
			$source_program = 'unknown';
		}
		if (m/^#CHROM/ and ($includeinfo || $header)){#for ABCC pipeline only
			my $ln=$_;$ln=~s/[\n\r]//g;
			print "$ln\n";
			if ($af){
				print AF join ("\t","#AlleleFreq","#Depth")."\n";#header for secondary file;
			}
			if ($withzyg){
				print ZYG "#zygosity\n";
			}
		}
		m/^#/ and next;		#skip comment lines
		s/[\r\n]+$//;		#delete trailing new lines
		next if $_=~/^\s*$/;
		#format description: http://www.1000genomes.org/wiki/Analysis/vcf4.0
		#standard VCF4 should have 8 columns, but some software may produce more columns (for example, for genotype calls). The first 8 columns should follow the specification
		
		#example of VCF4 generated by GATK SNP caller
		#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE
		#1       55      .       T       G       34.82   .       DP=2;Dels=0.00;HRun=0;HaplotypeScore=0.00;MQ=14.16;MQ0=0;QD=17.41;SB=-10.00     GT:DP:GL:GQ     0/1:1:-6.66,-0.30,-0.00:1.76
		#1       2646    .       G       A       40.91   .       DP=4;Dels=0.00;HRun=0;HaplotypeScore=0.00;MQ=7.50;MQ0=3;QD=10.23;SB=-10.00      GT:DP:GL:GQ     0/1:1:-7.27,-0.30,-0.00:1.76
		
		#example of VCF4 generated by GATK indel caller
		#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE
		#1       2525324 .       G       GC      .       PASS    AC=5,5;DP=12;MM=4.8,3.7142856;MQ=29.0,42.285713;NQSBQ=33.0,46.463768;NQSMM=0.24,0.20289855;SC=0,5,1,6  GT       0/1
		#1       3553372 .       GC      G       .       PASS    AC=6,6;DP=6;MM=0.8333333,0.0;MQ=60.0,0.0;NQSBQ=63.533333,0.0;NQSMM=0.0,0.0;SC=0,6,0,0   GT      1/0
		#1       6093011 .       CG      C       .       PASS    AC=31,31;DP=32;MM=0.7096774,2.0;MQ=59.64516,60.0;NQSBQ=64.192184,39.666668;NQSMM=0.0,0.11111111;SC=23,8,0,1     GT      1/0
		
		#example of VCF4 generated by 1000G
		#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
		#1       533     .       G       C       .       PASS    AA=.;AC=6;AN=120;DP=423
		#1       41342   .       T       A       .       PASS    AA=.;AC=29;AN=120;DP=188
		#1       41791   .       G       A       .       PASS    AA=.;AC=5;AN=120;DP=192
		#1       44449   .       T       C       .       PASS    AA=C;AC=2;AN=120;DP=166
		#1       44539   rs2462492       C       T       .       PASS    AA=T;AC=2;AN=120;DP=131    
		
		#example of VCF4 generated by 1000G
		#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
		#1       1000153 .       TCACA   T       100     PASS    AF=0.115095;HP=1;NF=16;NR=13;NS=52;CA=0;DP=615
		#1       1000906 .       CA      C       48      PASS    AF=0.0772696;HP=1;NF=2;NR=9;NS=51;CA=0;DP=281
		#1       1000950 rs60561655;-/G  CG      C       100     PASS    AF=0.447771;HP=5;DB;NF=10;NR=20;NS=50;CA=M;DP=291
		#1       1010786 rs36095298;-/G,mills,venter     A       AG      100     PASS    AF=0.774334;HP=1;DB;NF=21;NR=27;NS=51;CA=0;DP=306
		#1       1026158 .       T       TGGGGG  100     PASS    AF=0.115637;HP=1;NF=5;NR=2;NS=52;CA=0;DP=591
                
                #example of VCF4 generated by SamTools mpileup (Note that GT was not the first field in the FORMAT string)
                ##fileformat=VCFv4.0
		##samtoolsVersion=0.1.16 (r963:234)
		##fileformat=VCFv4.0
		#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  1247MFL0003.NOVO.srt.bam
		#chr1    14574   .       A       G       3.54    .       "DP=3;AF1=0.4999;CI95=0.5,0.5;DP4=1,0,2,0;MQ=21;PV4=1,1,1,1"    PL:GT:GQ        "31,0,34:0/1:32"
		#chr1    14930   .       A       G       37      .       "DP=19;AF1=0.5;CI95=0.5,0.5;DP4=7,5,5,1;MQ=25;PV4=0.6,6.3e-05,1,0.23"   PL:GT:GQ        "67,0,103:0/1:70"
		#chr1    16495   .       G       C       28      .       "DP=4;AF1=0.5;CI95=0.5,0.5;DP4=0,0,4,0;MQ=32"   PL:GT:GQ        "70,12,70:0/1:58"
		#chr1    59040   .       T       C       4.77    .       "DP=4;AF1=0.4999;CI95=0.5,0.5;DP4=0,2,2,0;MQ=22;PV4=0.33,0.21,1,1"      PL:GT:GQ        "33,0,39:0/1:35"
		#chr1    69270   .       A       G       46      .       "DP=20;AF1=0.5;CI95=0.5,0.5;DP4=2,0,18,0;MQ=24;PV4=1,1,1,0.28"  PL:GT:GQ        "94,18,100:0/1:78"
		#chr1    69511   .       A       G       24      .       "DP=5;AF1=0.5;CI95=0.5,0.5;DP4=1,0,2,1;MQ=25;PV4=1,0.46,1,0.039"        PL:GT:GQ        "54,0,57:0/1:55"
                
		#reserved VCF4 sub-fields in the INFO field
		#    * AA ancestral allele
		#    * AC allele count in genotypes, for each ALT allele, in the same order as listed
		#    * AF allele frequency for each ALT allele in the same order as listed: use this when estimated from primary data, not called genotypes
		#    * AN total number of alleles in called genotypes
		#    * BQ RMS base quality at this position
		#    * CIGAR cigar string describing how to align an alternate allele to the reference allele
		#    * DB dbSNP membership
		#    * DP combined depth across samples, e.g. DP=154
		#    * END end position of the variant described in this record (esp. for CNVs)
		#    * H2 membership in hapmap2
		#    * MQ RMS mapping quality, e.g. MQ=52
		#    * MQ0 Number of MAPQ == 0 reads covering this record
		#    * NS Number of samples with data
		#    * SB strand bias at this position
		#    * SOMATIC indicates that the record is a somatic mutation, for cancer genomics
		#    * VALIDATED validated by follow-up experiment


		#SAMtools/BCFtools specific information
		#SAMtools/BCFtools may write the following tags in the INFO field in VCF/BCF.
		#Tag	Description
		#I16	16 integers:
		#1	#reference Q13 bases on the forward strand 	2	#reference Q13 bases on the reverse strand
		#3	#non-ref Q13 bases on the forward strand 	4	#non-ref Q13 bases on the reverse strand
		#5	sum of reference base qualities 	6	sum of squares of reference base qualities
		#7	sum of non-ref base qualities 	8	sum of squares of non-ref base qualities
		#9	sum of ref mapping qualities 	10	sum of squares of ref mapping qualities
		#11	sum of non-ref mapping qualities 	12	sum of squares of non-ref mapping qualities
		#13	sum of tail distance for ref bases 	14	sum of squares of tail distance for ref bases
		#15	sum of tail distance for non-ref bases 	16	sum of squares of tail distance for non-ref
		#INDEL	Indicating the variant is an INDEL.
		#DP	The number of reads covering or bridging POS.
		#DP4	Number of 1) forward ref alleles; 2) reverse ref; 3) forward non-ref; 4) reverse non-ref alleles, used in variant calling. Sum can be smaller than DP because low-quality bases are not counted.
		#PV4	P-values for 1) strand bias (exact test); 2) baseQ bias (t-test); 3) mapQ bias (t); 4) tail distance bias (t)
		#FQ	Consensus quality. If positive, FQ equals the phred-scaled probability of there being two or more different alleles. If negative, FQ equals the minus phred-scaled probability of all chromosomes being identical. Notably, given one sample, FQ is positive at hets and negative at homs.
		#AF1	EM estimate of the site allele frequency of the strongest non-reference allele.
		#CI95	Equal-tail (Bayesian) credible interval of the site allele frequency at the 95% level.
		#PC2	Phred-scaled probability of the alternate allele frequency of group1 samples being larger (,smaller) than of group2 samples.
		#PCHI2	Posterior weighted chi^2 P-value between group1 and group2 samples. This P-value is conservative.
		#QCHI2	Phred-scaled PCHI2
		#RP	Number of permutations yeilding a smaller PCHI2

		#example of triallelic variants generated by mpileup/bcftools
		#1       156706559       .       A       C,G     114     .	DP=20;AF1=1;CI95=1,1;DP4=0,0,1,19;MQ=60;FQ=-63  GT:PL:GQ	1/2:237,126,90,162,0,138:99
		#6       31129642        .       A       G,C     76      .	DP=31;AF1=1;CI95=1,1;DP4=0,0,28,3;MQ=60;FQ=-75  GT:PL:GQ	1/2:255,194,146,164,0,119:99
		#1       11297762        .       T       C,A     98      .	DP=19;AF1=1;CI95=1,1;DP4=0,0,17,1;MQ=60;FQ=-78  GT:PL:GQ	1/1:131,51,0,120,28,117:99
		
		my @field=split(/\t/,$_);
		@field >=5 or warn "Warn: invalid record found in VCF4 file (at least 8 tab-delimited fields expected): <$_>\n";
		my ($chr, $start, $ID, $ref_allele, $mut_allele, $quality_score, $filter, $info, $format, @samples) = @field;
		$quality_score||=0;
		$filter||=0;
		$info||='';
		$format||='';
		my $sample=$samples[0];
		if ($start!~/^\d+/ and $chr=~/Chrom{0,1}$/i){
			print "#$_\n";$af=0;
			close AF and system ("rm $outfile.hg_af\n") if ($af);
			next;
		}
		$chr=~s/\s//g;
		$start=~s/\s//g;
		$ref_allele=~s/\s//g;
		$mut_allele=~s/\s//g;
		if ($includeinfo){
			chop $field[7] if ($field[7]=~/[;,:]$/);
			$field[2]=~s/;/,/g;
			 $field[7].="ABCC_TAG=LINE".$countline.";";## need this for converting back to VCF for AVIA v2.0 --changed 20160104 from : $field[7].="($field[2] eq '.')?";"ABCC_TAG=LINE".$countline.";":";ABCC_TAG=".$field[2].";"; 
			$field[8].=($af && $#samples>0 && $format=~/AD:/)?':SMAF':'';#for Michelle add AF to each of the sample fields instead of a different column		
		}
		my $af_idx;
		if ($af){
			if ($format=~/^(.*:AD:)/ ){#jason lih
				my $fmt=$1;
				$af_idx=$fmt=~tr/://;$af_idx--;#offset by 1 b/c starts at 0
				for (my $i=0;$i<=$#samples;$i++){
					my @sampleinfo=split(":",$samples[$i]);
					if ($#samples==0){
						if ($af_idx<0){
							print AF "-\t-\n";$stats{no_af}++;
						}else{
							my ($ref_ct,$mut_ct)=split(",",$sampleinfo[$af_idx]);
							my $total=($ref_ct+$mut_ct);
							if ($total){
								my $afreq=sprintf("%.03f",$mut_ct/($total));
								print AF "$afreq\t$total\n";$stats{afct}++;
							}else{
								print AF "-\t-\n";$stats{no_af}++;
							}
						}
					}else{
						close AF and system ("rm $outfile.hg19_af\n") if (-e "$outfile.hg19_af");
						if ($af_idx<0 || $samples[$i]=~/\.\/\./){
							$samples[$i].=":0.00" if ($samples[$i]!~/\.\/\./);$stats{no_af}++;
						}else{
							my ($ref_ct,@mut_cts)=split(",",$sampleinfo[$af_idx]);
							my $mut_ct=my $mut_idx=-1;#the counts of the current maf and the index of the maf
							for (my $k=0;$k<=$#mut_cts;$k++){
								if ( $mut_cts[$k]>$mut_ct){
									$mut_ct=$mut_cts[$k];
									$mut_idx=$k;
								}
							}
							my $total=($mut_ct < 0 || $mut_idx<0)?-1:($ref_ct+$mut_cts[$mut_idx]);
							if ($total>0){
								my $afreq=sprintf("%.03f",$mut_cts[$mut_idx]/($total));
								$samples[$i].=":$afreq";$stats{afct}++;
							}else{
								$samples[$i].=":".sprintf("%.03f",0) if ($samples[$i]!~/\.\/\./);$stats{no_af}++;
							}
						}
						$field[9+$i]=$samples[$i];
					}
				}
			}elsif ($info=~/AF=([\d\.\,]+);/){
				my $tmp_af=$1;
				my @mut_allele=split(",",$mut_allele);
				my @afs_allele=split(",",$tmp_af);
				my $depth;
				if ($info=~/NS=(\d+)/){
					$depth=$1;
				}
				for (my $i=0;$i<=$#mut_allele;$i++){
					print AF "$afs_allele[$i]\t$depth\n";$stats{afct}++;
					last if (!$allallele);
				}			
			}else{
				my $mut_count=($mut_allele=~tr/,//);
				for (my $i=0;$i<=$mut_count;$i++){
					print AF "-\t-\n";$stats{no_af}++;
					last if (!$allallele);
				}
			}
		}
		my ($end);
		my (@mut_alleles, $zygosity);
		my $otherinfo = join("\t",@field);	#this is the complete line (when -includeinfo is set, the entire line will be included in output file)
		if ($filterword) {		#ensure that the filter field contains the filterword
			$filter =~ m/\b$filterword\b/i or next;
		}
		$info =~ s/^"//; $info =~ s/"$//;
		#sometimes the alleles are not in the same case
		#chr1    1869771 1869774 actc    aCTctc          43.5    13      INDEL;DP=13;AF1=0.5;CI95=0.5,0.5;DP4=0,4,4,0;MQ=37;PV4=0.029,0.45,1,0.46
		$ref_allele = uc $ref_allele;
		$mut_allele = uc $mut_allele;
		
		#if ($ID eq '.' || $ID =~ /^rs/) {		#per MISHIMA, Hiroyuki suggestion (vcf4's third column (ID column) are not always ".")
		#	$end = $start;				#this block is commented out on 2011feb19
		#}
		
		if ($mut_allele eq '.') {			#no variant call was made at this position
			$stats{'unknown'}++;
			next VAR;
		}
		#hv changed the next block
		if ($mut_allele =~ m/,/ ) {
			if ($allallele){
				@mut_alleles=split(",",$mut_allele);
			}else{
				$mut_allele=~/^([ACGNT]*),/;
				push (@mut_alleles,$1);
			}
#			$mut_allele2 = $2;
		}else{
			push(@mut_alleles,$mut_allele);
		}
		my $ct=($#mut_alleles+1);
		$stats{"$ct-allelic"}++;	
		#ends hv block
		
		foreach my $mymut (@mut_alleles){
			my ($mystart,$end,$myref);#define these inside loop because it may change per mutation
			next if ($mymut eq '*');
			
			($mystart,$myref)=($start,$ref_allele);$end=$mystart;
			# print "working on $mymut\n";
			if (my $head=smallest_substring($myref,$mymut)){
				chomp $head;
				# if ($debug){print "smallest substring...working on $head||$myref,$mymut";<STDIN>;	}
				if ($head eq $mymut) {#deletion
						$mystart=$mystart+length($head);
						$myref=substr ($myref, length($head),length($myref)-length($head));
						$end=$mystart+length($myref)-1;
						$mymut="-";
				}elsif ($head eq $myref){#insertion
						$myref='-';
						$mymut=~s/$head//;
				}else{#snp or multiallelic snp
						$myref=~s/^$head//;
						$mymut=~s/^$head//;
						$mystart+=length($head);
						$end=$mystart+length($myref)-1;
				}
				# print "FOUND:($head)".join (",",@field)."\n";<STDIN>;
			}else{
				$end=$mystart+length($ref_allele)-1;
				# print __LINE__. ":($head)".join (",",@field)."\n";<STDIN>;
			}
			# print "$myref,$mymut...\n";
			my ($ref_len,$mut_len)=(length($myref),length($mymut));
			if( ($ref_len==$mut_len && $ref_len==1 && $myref ne '-')) {  	### output snv
				##chr1    1887091 13018;13019     CG      C,CA    .    #abcc_hv
				my ($unfiltered_read_depth) = $info =~ /DP=(\d+)/;
				my ($MappingQuality) = $info =~ /MQ=([^;]+)/; 
				my ($QualityByDepth) = $info =~ /QD=([^;]+)/;		
				if ($coverage) {
					$stats{'failed_rd_depth'}++;
					# defined $unfiltered_read_depth and $unfiltered_read_depth >= $coverage || next;
					$stats{'failed_rd_depth'}--;
					if ($maxcoverage) {
						$stats{'failed_maxcvg'}++;
						# defined $unfiltered_read_depth and $unfiltered_read_depth <= $maxcoverage || next;
						$stats{'failed_maxcvg'}--;
					}
				}
				if ($snpqual) {
					$stats{'failed_qual'}++;
					# defined $QualityByDepth and $QualityByDepth >= $snpqual || next;		#the QD was used here as quality score
					$stats{'failed_qual'}--;
				}		
				$zygosity='-';	
				if (defined $format) {
					my @format = split (/:/, $format);
					undef $gtpos;
					for my $i (0 .. @format-1) {
						if ($format[$i] eq 'GT') {
							$gtpos = $i;
							last;
						}
					}
					if (defined $sample and defined $gtpos) {
						my @sample = split (/:/, $sample);
						if ($sample[$gtpos] =~ /^(\d+)[\/\|](\d+)$/ &&  $1!=$2) {
							$zygosity = 'het';
							$counthet++;
						} elsif ($sample[$gtpos] =~ /^(\d+)[\/\|](\d+)$/) {
							$zygosity = 'hom';
							$counthom++;
						} else {
							$zygosity = '-';
							$countunknown++;
						}
						print ZYG join (":",$chr, $mystart,$end, $myref,  $mymut). "\t" .$zygosity ."\n" if ($withzyg);
					}
				} else {
					$zygosity = '-';
					$countunknown++;
					print ZYG join (":",$chr, $mystart,$end, $myref,  $mymut). "\t" .$zygosity ."\n" if ($withzyg);
				}
				$stats{$zygosity}++;
	
	#			#the subject is called as homozygous for the first alternative allele (genotype 1/1. i.e. C/C), but since there was one read containing A, samtools still keep both alleles in the VCF file (but gives a very low probabilities for it).
	#			#1       11297762        .       T       C,A     98      . DP=19;AF1=1;CI95=1,1;DP4=0,0,17,1;MQ=60;FQ=-78  GT:PL:GQ 1/1:131,51,0,120,28,117:99			
	#			if ($mut_allele2 and $zygosity eq 'hom') {
	#				$mut_allele2 = '';
	#			}
				if ($myref=~/A$/ and $mymut=~/G$/ or $myref=~/G$/ and $mymut=~/A$/ or $myref=~/C$/ and $mymut=~/T$/ or $myref=~/T$/ and $mymut=~/C$/) {
					$countti++;
					$stats{'snp_ct_ti'}++;
				} else {
					$counttv++;
					$stats{'snp_ct_tv'}++;
				}
				if ($includeinfo) {
					print $chr, "\t", $mystart, "\t", $end, "\t", $myref, "\t", $mymut, $withzyg?"\t$zygosity":"", "\t", $otherinfo, "\n";
				} elsif ($bare){
					print $chr, "\t", $mystart, "\t", $end, "\t", $myref, "\t", $mymut,  "\n";
				}else{
					print $chr, "\t", $mystart, "\t", $end, "\t", $myref, "\t", $mymut, $withzyg?"\t$zygosity":"",  "\t", $quality_score, (defined $unfiltered_read_depth)? "\t$unfiltered_read_depth" : '', (defined $MappingQuality) ? "\t$MappingQuality" : '', (defined $QualityByDepth) ? "\t$QualityByDepth" : '', "\n";
				}
				#print $chr, "\t", $start, "\t", $start, "\t", $ref_allele, "\t", $mut_allele, "\t$zygosity",  "\t", $quality_score, (defined $unfiltered_read_depth)? "\t$unfiltered_read_depth" : '', (defined $MappingQuality) ? "\t$MappingQuality" : '', (defined $QualityByDepth) ? "\t$QualityByDepth" : '', $includeinfo ? "\t$otherinfo" : '', "\n";	#commented Sep 2011
				$countsnp++;
			} else{#}if ($myref eq '-' || $mymut eq '-') {  ### output indel
				if ($debug){print __LINE__.":working on $mystart-$end,$myref,$mymut";<STDIN>;	}
# print "TESTING========================================\n$chr,$mystart,$end,$ref_len,$mut_len,$myref,$mymut\t|$otherinfo\n========================================\n";<STDIN>;
				my ($indel_read_depth1, $indel_read_depth2) = $info =~ /AC=([^,;]+),([^,;]+)/;		#number of reads supporting consensus indel, any indel
				my ($unfiltered_read_depth) = $info =~ /DP=(\d+)/;
				$stats{'indel'}++;$stats{"indel_$ct"}++;
									# if ($coverage) {die "coverage";
									# 	$stats{'failed_rd_depth'}++;
									# 	defined $unfiltered_read_depth and $unfiltered_read_depth >= $coverage || next;
									# 	$stats{'failed_rd_depth'}--;
									# 	if ($maxcoverage) {
									# 		$stats{'failed_max_cvg'}++;
									# 		defined $unfiltered_read_depth and $unfiltered_read_depth <= $maxcoverage || next;
									# 		$stats{'failed_max_cvg'}--;
									# 	}
									# }
									
									# if (defined $indel_read_depth1 and defined $unfiltered_read_depth) {die "read depth\n"
									# 	$stats{'failed_frac_rddepth'}++;
									# 	$indel_read_depth1/$unfiltered_read_depth >= $fraction or next;		#do not meet minimum alternative allele fraction threshold
									# 	$indel_read_depth2 and $indel_read_depth1/$indel_read_depth2 >= $confraction or next;
									# 	$stats{'failed_frac_rddepth'}--;
									# }
									
				my ($MappingQuality) = $info =~ /MQ=([^;]+),/;
			
				#example VCF4 records below:
				#20      2       .       TCG     T       .       PASS    DP=100
				#Chr1    5473    .       AT      ATT     23.5    .       INDEL;DP=16;AF1=0.5;CI95=0.5,0.5;DP4=4,2,3,1;MQ=42;PV4=1,0.41,0.042,0.24
				#Chr1    6498    .       ATTTT   ATTTTT  53.5    .       INDEL;DP=9;AF1=1;CI95=1,1;DP4=0,0,5,3;MQ=28
				#chr1   889158  .       GA      CC   
				my $mycoords;
				# if ($ref_len>=$mut_len && $myref ne '-' && $mymut ne '-'){
				# 	my $diff;
				# 	if ($ref_len==$mut_len && $ref_len>1){
				# 		$diff=$ref_len-1;
				# 	}elsif ($ref_len==$mut_len){
				# 		$diff=0;
				# 	}else{
				# 		$diff=$ref_len-$mut_len+1;
				# 	}
				# 	print $chr, "\t", $mystart, "\t", ($end+$diff), "\t", $myref, "\t", $mymut ; 
				# 	$mycoords=join(":",$chr, $mystart, ($end+$diff), $myref, $mymut );
				# }else{
					print $chr, "\t", $mystart, "\t", $end, "\t", $myref, "\t", $mymut;
					$mycoords=join(":",$chr, $mystart, $end, $myref, $mymut );
				# }
				$zygosity='-';
				if (defined $format) {
					my @format = split (/:/, $format);
					undef $gtpos;
					for my $i (0 .. @format-1) {
						if ($format[$i] eq 'GT') {
							$gtpos = $i;
							last;
						}
					}
					if (defined $sample and defined $gtpos) {
						my @sample = split (/:/, $sample);
						if ($sample[$gtpos] =~ /^(\d+)[\/\|](\d+)$/ &&  $1!=$2) {
							$zygosity = 'het';
							$counthet++;
						} elsif ($sample[$gtpos] =~/^(\d+)[\/\|](\d+)$/) {
							$zygosity = 'hom';
							$counthom++;
						} else {
							$zygosity = '-';
							$countunknown++;
						}
						print ZYG "$mycoords\t$zygosity\n" if ($withzyg);
					}
				} else {
					$zygosity = '-';
					$countunknown++;
					print ZYG "$mycoords\t$zygosity\n" if ($withzyg);
				}
				$stats{$zygosity}++;
				#commented out on 2011May21
				#if (defined $sample) {
				#	if ($sample =~ m#^0/1# or $sample =~ m#^1/0#) {
				#		print "\thet";
				#		$counthet++;
				#	} elsif ($sample =~ m#^1/1#) {
				#		print "\thom";
				#		$counthom++;
				#	} # BEGIN ARQ
				#	elsif ($sample =~ m#^./.#) {
				#		print "\tunknown";
				#		$countunknown++;
				#	} # END ARQ
				#}
				
				if ($includeinfo) {
					 print $withzyg?"\t$zygosity":"", "\t", $otherinfo;
				}elsif ($bare){
					
				} else {
					print "\t", $zygosity;
					print "\t", $quality_score;
					defined $unfiltered_read_depth and print "\t", $unfiltered_read_depth;
					
					defined $indel_read_depth1 and print "\t", $indel_read_depth1;
					defined $MappingQuality and print "\t", $MappingQuality;
					#$includeinfo and print "\t", $otherinfo;	#commented Sep 2011
				}
				print "\n";
				$countindel++;
	
#			}else{
#				print "========================================\n$chr,$mystart,$end,$ref_len,$mut_len,$myref,$mymut\t|$otherinfo\n========================================\n";<STDIN>;
			}
		}
		$countvar++;
	}
	my $triallelic = $countsnp-$countti-$counttv;
	print STDERR "NOTICE: Read $countline lines and wrote ${\($counthet+$counthom)} different variants at $countvar genomic positions ($countsnp SNPs and $countindel indels)\n";
	print STDERR "NOTICE: Among ${\($counthet+$counthom+$countunknown)} different variants at $countvar positions, $counthet are heterozygotes, $counthom are homozygotes\n";
	print STDERR "NOTICE: Among $countsnp SNPs, $countti are transitions, $counttv are transversions", $triallelic?", $triallelic have more than 2 alleles\n":"\n";
	print STDERR Dumper (\%stats);
	if ($stats{no_af} && $stats{afct} && $stats{no_af}>$stats{afct} && -e "$outfile.hg19_af"){
		system ("rm $outfile.hg19_af\n");
	}
}
sub revcompl{
	my $string=shift;
	return '' if (!$string);
	my $rev_string;
	for (my $i=length($string);$i>=0;$i--){
	        $rev_string.=substr($string,$i,1);
	}
	$rev_string=~tr/ACGTRYKMSWacgtrykmsw/TGCAYRMKSWtgcayrmksw/;
	return $rev_string;
}
sub convertTVC{
	my ($variantfile)=@_;
	my $revcmpl_exe='/users/abcc/vuonghm/scripts.dir/util/revCmpl.pl';
	open (VAR,"<$variantfile") or die "Cannot open variant file $variantfile\n";
	my $header_found=0;
	my %Expected_Headers=('Chrom'=>-1,
		'Position'=>-1,
		'Ref'=>-1,
		'Variant'=>-1,
		);
	#Get strand information
	# open (STRAND,"</SeqIdx/annovardb/humandb/hg19_refGene2Strand.txt") or die 'Your file for genic strands does not exist!';
	# my %gene_ori;
	# while (<STRAND>){
	# 	my ($gene,$ori)=split("\t",$_);chomp $ori;
	# 	$gene_ori{uc($gene)}=$ori;
	# }
	# close STRAND;
	my %header_idx;my $count=0;
	while (my $line=<VAR>){
		chomp $line;
		next if ($line=~/^\s{0,}$/);
		my @array=split("\t",$line);
		if ($header_found){
			my ($ori,$end);
			$end=$array[$Expected_Headers{'Position'}];
			if ($array[$Expected_Headers{'Chrom'}]=~/((chr){0,1}[\dXYMT]{1,})$/){
				$chr=$1;
			}
			# $gene_ori{$array[$Expected_Headers{'Gene Sym'}]}=uc($gene_ori{$array[$Expected_Headers{'Gene Sym'}]});
			# if (!exists $gene_ori{$array[$Expected_Headers{'Gene Sym'}]}){
			# 	print STDERR "Skipping $array[$Expected_Headers{'Gene Sym'}]...does not exist\n" and next;
			# }else{
			# 	$ori=$gene_ori{$array[$Expected_Headers{'Gene Sym'}]};
			# }
			if ($array[$Expected_Headers{'Position'}]!~/^\d+$/){
				print STDERR "Skipping $array[$Expected_Headers{'Position'}]...is not an integer\n" and next;
			}
			
			if ($array[$Expected_Headers{'Ref'}]!~/^[ATCG]*$/i){
				print STDERR "Skipping $array[$Expected_Headers{'Ref'}]...is not ATCG\n" and next;
			}
			if ($array[$Expected_Headers{'Variant'}]!~/^[ATCG]*$/i){
				print STDERR "Skipping $array[$Expected_Headers{'Variant'}]...is not ATCG\n" and next;
			}
			# if ($array[$Expected_Headers{'Type'}]=~/DEL/i){
			# 	$end+=length($array[$Expected_Headers{'Ref'}])-1;
			# # }elsif($array[$Expected_Headers{'Type'}]=~/INS/i){
			# }els
			my $ref=$array[$Expected_Headers{'Ref'}];my $myref=$ref;
			my $mymut=$array[$Expected_Headers{'Variant'}];
			if(length($ref)>1 || length($mymut)>1){
				if (my $head=smallest_substring($myref,$mymut)){
					chomp $head;
					if ($head eq $mymut) {#deletion
						$array[$Expected_Headers{'Position'}]=$array[$Expected_Headers{'Position'}]+length($head);
						$myref=substr ($myref, length($head),length($myref)-length($head));
						$end=$array[$Expected_Headers{'Position'}]+length($myref)-1;
						$mymut="-";
					}elsif ($head eq $myref){#insertion
						$myref='-';
						$mymut=~s/$head//;
					}else{#snp or multiallelic snp
						$myref=~s/$head//;
						$mymut=~s/$head//;
						$array[$Expected_Headers{'Position'}]+=length($head);
						$end=$array[$Expected_Headers{'Position'}]+length($myref)-1;
					}
	#				print "FOU
				}else{
					die "Don't know what to do with this $line\n";
				}
			}
			print join ("\t",@array[$Expected_Headers{'Chrom'},$Expected_Headers{'Position'}],$end,$myref,$mymut);
			if ($includeinfo){
				print "\t". join ("\t",$line);
			}
			print "\n";#<STDIN> if ($flag);

		}elsif ($line=~/(\bChrom\b)/i){
			
			for (my $i=0;$i<=$#array;$i++){
				if (exists $Expected_Headers{$array[$i]}){
					$Expected_Headers{$array[$i]}=$i;
					$count++;
				}
			}
			die "Could not find one or more of your headers". Dumper (\%Expected_Headers) if ($count!=keys(%Expected_Headers));
			print "#$line" ;
			print "\n";
			$header_found=1;
		}else{
			die "Could not find header!";
		}
	}
}
sub convertMutector{#Added for Michele 2013-05
	my ($variantfile) = @_;
	
	my ($countline, $countvar, $counthom, $counthet, $countunknown, $countindel, $countsnp, $countti, $counttv, $countaf) = qw/0 0 0 0 0 0 0 0 0 0/;#ANNOVAR default
	my %stats;#hv added 2/10/12
	my ($source_program, $gtpos);		#the program that generated the VCF4 file; the GT position within FORMAT record
	open (VAR, $variantfile) or die "Error: cannot read from variant file $variantfile: $!\n";
	if ($af){
		print STDERR "Notice:  --af is not supported for this format --$format\n";
	}
	VAR: while (<VAR>) {
		$countline++;
		
		m/^#/ and next;		#skip comment lines
		print "#".$_  and next if ($_=~/contig\tpos/ && ($header || $includeinfo));#print header
		s/[\r\n]+$//;		#delete trailing new lines
		
		my @field=split(/\t/,$_);
		@field >29 or die "Error: invalid record found in Mutector file (at least 5 tab-delimited fields expected): <$_>\n";
		my ($chr, $start, $context,$ref_allele, $mut_allele, @info) = @field;
		if ($start!~/^\d+/ and $chr!~/chr[\dXYM]{1,2}$/i){
			print "#$_\n";$af=0;
			next;
		}
		my ($end);
		my (@mut_alleles, $zygosity);
		my $otherinfo = join("\t",@field);	#this is the complete line (when -includeinfo is set, the entire line will be included in output file)
		#sometimes the alleles are not in the same case
		## muTector v1.0.47986
#contig  position        context ref_allele      alt_allele      tumor_name      normal_name     score   dbsnp_site      covered power   tumor_power     normal_power    total_pairs     improper_pairs  map_Q0_reads    t_lod_fstar     tumor_f contaminant_fraction    contaminant_lod t_ref_count     t_alt_count     t_ref_sum       t_alt_sum       t_ref_max_mapq  t_alt_max_mapq  t_ins_count     t_del_count     normal_best_gt  init_n_lod      n_ref_count     n_alt_count     n_ref_sum       n_alt_sum       judgement
#chr1    14976   TGCxCCT G       A       0063_FP_SSv4AllExome    0063_BL_SSv4AllExome    0       DBSNP   UNCOVERED       0       0.043787        0       19      0       0       8.568941        0.6     0.02    -0.061305       2       3       61      86      56      56      0       0       AG      -5.448858       3       2       93      60      REJECT
#chr1    16963   GTCxTTG G       A       0063_FP_SSv4AllExome    0063_BL_SSv4AllExome    0       NOVEL   UNCOVERED       0       0.024074        0       14      0       1       5.092379        0.4     0.02    -0.043852       3       2       91      56      3       1       0       0       GG      1.805146        6       0       180     0       REJECT
		$ref_allele = uc $ref_allele;
		$mut_allele = uc $mut_allele;
		
		if ($mut_allele eq '.') {			#no variant call was made at this position
			$stats{'unknown'}++;
			next VAR;
		}
		if ($mut_allele =~ m/,/ ) {
			if ($allallele){
				@mut_alleles=split(",",$mut_allele);
			}else{
				$mut_allele=~/^([ACGNT]*),/;
				push (@mut_alleles,$1);
			}
		}else{
			push(@mut_alleles,$mut_allele);
		}
		my $ct=($#mut_alleles+1);
		$stats{"$ct-allelic"}++;	
		#ends hv block
		
		foreach my $mymut (@mut_alleles){
			my ($mystart,$end,$myref,$mytype);#define these inside loop because it may change per mutation
			($mystart,$myref)=($start,$ref_allele);$end=$mystart;
			if ($mymut eq '-'){#indel ok
				$mytype='del';
			}elsif ($mymut =~/^\+([ATCG]+)/){#insertion explicit
				$mymut=$1;
				$mytype='ins';
			}elsif ($mymut=~/^\-([ATCG]+)/){#deletion explicit
				$mymut=$1;
				$mytype='del';
			}else{
				$mytype='snp';
			}
			my ($ref_len,$mut_len)=(length($myref),length($mymut));
			if ($mytype eq 'del') {#deletion
					$mystart=$mystart+1;
					$myref=$mymut;
					$end=$mystart+length($myref)-1;
					$mymut="-";
			}elsif ($mytype eq 'ins'){#insertion
					$mystart++;
					$end++;
					$myref='-';
			}else{#snp 
					
			}
			if( ($ref_len==$mut_len && $myref ne '-')) {  	### output snv
				##chr1    1887091 13018;13019     CG      C,CA    .    #abcc_hv
	
				if ($myref=~/A$/ and $mymut=~/G$/ or $myref=~/G$/ and $mymut=~/A$/ or $myref=~/C$/ and $mymut=~/T$/ or $myref=~/T$/ and $mymut=~/C$/) {
					$countti++;
					$stats{'snp_ct_ti'}++;
				} else {
					$counttv++;
					$stats{'snp_ct_tv'}++;
				}
				if ($includeinfo) {
					print $chr, "\t", $mystart, "\t", $end, "\t", $myref, "\t", $mymut,  "\t", $otherinfo, "\n";
				} else{
					print $chr, "\t", $mystart, "\t", $end, "\t", $myref, "\t", $mymut,  "\n";
				}
				#print $chr, "\t", $start, "\t", $start, "\t", $ref_allele, "\t", $mut_allele, "\t$zygosity",  "\t", $quality_score, (defined $unfiltered_read_depth)? "\t$unfiltered_read_depth" : '', (defined $MappingQuality) ? "\t$MappingQuality" : '', (defined $QualityByDepth) ? "\t$QualityByDepth" : '', $includeinfo ? "\t$otherinfo" : '', "\n";	#commented Sep 2011
				$countsnp++;
			} else{
				if ($includeinfo) {
					print $chr, "\t", $mystart, "\t", $end, "\t", $myref, "\t", $mymut,  "\t", $otherinfo, "\n";
				} else{
					print $chr, "\t", $mystart, "\t", $end, "\t", $myref, "\t", $mymut,  "\n";
				}
			}
		}
	
		$countvar++;
	}
}
sub convertFromDBSNP{
	my ($variantfile) = @_;
	open (VAR ,$variantfile) or die 'Cannot open ' . $variantfile. "\n";
	if ($header){
		print join ("\t","#Chr","Start","End","RefAllele","MutAllele","dbSNP id")."\n";
	}
	while (<VAR>){

	}
}
sub convertBED {#hv added
	my ($variantfile) = @_;
	$debug=0;
	my ($countline, $countvar, $counthom, $counthet, $countunknown, $countindel, $countsnp, $countti, $counttv, $countaf) = qw/0 0 0 0 0 0 0 0 0 0/;#ANNOVAR default
	my %stats;#hv added 2/10/12
	my ($source_program, $gtpos);		#the program that generated the VCF4 file; the GT position within FORMAT record
	open (VAR, $variantfile) or die "Error: cannot read from variant file $variantfile: $!\n";
	my $header=`head $variantfile -n1 |grep '^#'`;chomp $header;
	if ($header){
		print "$header\n";
	}elsif($nofail){
		print join ("\t","#Chr","Start","End")."\n";
	}else{
		print join ("\t","#Chr","Start","End","RefAllele","MutAllele",$header)."\n";
	}
	VAR: while (<VAR>) {
		m/^#/ and $countline++ and next;		#skip comment lines
		s/[\r\n]+$//;		#delete trailing new lines
		my $line=$_;
		$line=~s/\s{2,}/ /g;
		next if ($_=~/^$/);
		$countline++;
		my @field=split(/\t/,$line);
		my $fail;
		if (@field <=4){
			#try to split by spaces
			next if ($line=~/^\s+$/);
			@field=split(/\s/,$line,6);
			@field>=4 or $fail=1;
			if ($nofail){
				if ($line=~/^\s{0,}$/){next;}
				@field=split(/[\s\t]/,$line);
				print join ("\t",@field[0..2])."\n";next;
			}elsif ($fail){
				die "Error: invalid record found in BED file (at least 4 tab-delimited fields expected): <$_>\n";
			}
		}
		
		my ($chr, $start, $stop, $ref_allele, $mut_allele, @samples) = @field;
		my $sample=$samples[0];
		next if ($start!~/^\d+/ and $stop!~/^\d+/ and $chr!~/(chr){0,1}[\dXYMT]*/i);#ignore anything that doesn't resemble a bed file
		
		my ($end);
		my (@mut_alleles, $zygosity,$filter);
		my $otherinfo = join("\t",@field);	#this is the complete line (when -includeinfo is set, the entire line will be included in output file)
		if ($filterword) {		#ensure that the filter field contains the filterword
			$filter =~ m/\b$filterword\b/i or next;
		}
		#sometimes the alleles are not in the same case
		#chr1    1869771 1869774 actc    aCTctc          43.5    13      INDEL;DP=13;AF1=0.5;CI95=0.5,0.5;DP4=0,4,4,0;MQ=37;PV4=0.029,0.45,1,0.46
		$ref_allele = uc $ref_allele;
		$mut_allele = uc $mut_allele;
# print STDERR __LINE__.":working on $chr, $start, $stop, $ref_allele, $mut_allele\n";<STDIN>;
		#if ($ID eq '.' || $ID =~ /^rs/) {		#per MISHIMA, Hiroyuki suggestion (vcf4's third column (ID column) are not always ".")
		#	$end = $start;				#this block is commented out on 2011feb19
		#}
		
		if ($mut_allele eq '.') {			#no variant call was made at this position
			$stats{'unknown'}++;
			next VAR;
		}elsif ($mut_allele!~/^[\w\-\.]+$/){
			next;
		}
		#hv changed the next block
		if ($mut_allele =~ m/,/ ) {
			if ($allallele){
				@mut_alleles=split(",",$mut_allele);
			}else{
				$mut_allele=~/^([ACGNT]*),/;
				push (@mut_alleles,$1);
			}
#			$mut_allele2 = $2;
		}else{
			push(@mut_alleles,$mut_allele);
		}
		my $ct=($#mut_alleles+1);
		$stats{"$ct-allelic"}++;	
		#ends hv block

		if ($mut_allele ne "$ref_allele"){
			foreach my $mymut (@mut_alleles){
				my ($mystart,$end,$myref);#define these inside loop because it may change per mutation
if ($debug){print STDERR __LINE__.":working on $mymut($ref_allele)\n";<STDIN>;}
				
				($mystart,$myref)=($start,$ref_allele);$end=$mystart;

				if (my $head=smallest_substring($myref,$mymut)){
					chomp $head;
					if ($debug){print __LINE__.":$head\n";}
					if ($head eq $mymut) {#deletion
							$mystart=$mystart+length($head);
							$myref=substr ($myref, length($head),length($myref)-length($head));
							$end=$mystart+length($myref)-1;
							$mymut="-";
							# print "going to use $mystart,$end\n";<STDIN>;
					}elsif ($head eq $myref){#insertion
							$myref='-';
							$mymut=~s/$head//;
					}elsif ($start==$stop && length($ref_allele)>1){
						$end=$mystart+length($ref_allele)+1;
					}else{#snp or multiallelic snp
							$myref=~s/$head//;
							$mymut=~s/$head//;
							$mystart+=length($head);
							$end=$mystart+length($myref)-1;
							
					}
				}elsif ($mymut eq '-' && $myref=~/[ATCG]*/){#hv added 20140604
					# deletion
					$end=$mystart+length($myref)-1;
					# print "$mystart:$end and old:$stop\n";
				}elsif ($start!=$stop){
					if ($debug){print __LINE__. smallest_substring($myref,$mymut). "?\n";<STDIN>;}
					#check that there is a deletion and that it is correct
					# if (($stop-$start+1)==length($ref_allele)){
					# 	$end=$stop;
					# 	$mymut='-';
					# 	$myref=$ref_allele;
					# }else{
						#have not tested
						if ($debug){print __LINE__.":". length($ref_allele) . " and length of mut:" .length($mymut)."\n";<STDIN>;}
						$end=$mystart+length($ref_allele)-1;
					# }
					# print "SMLL".smallest_substring($myref,$mymut);<STDIN>;
				}else{
					$end=$mystart;
				}


					
				my ($ref_len,$mut_len)=(length($myref),length($mymut));
				
				# if ($mymut eq '-') {#deletion
				# 		$end=$mystart+length($myref)-1;
				# 		$mymut="-";
				# }elsif ($myref eq '-'){#insertion
				# 		$myref='-';
				# }elsif ($ref_len>1){
				# 	if ($ref_len==$mut_len){
				# 		$end=($mystart+$ref_len-1);#AA TT
				# 	}elsif ($ref_len>$mut_len){#AA G
				# 		$end=($mystart+$ref_len-1);
				# 	}else{#AG TGA
				# 		$end=($mystart+$ref_len-1);
				# 	}
				# }
		
				if ($includeinfo) {

					print $chr, "\t", $mystart, "\t", $end, "\t", $myref, "\t", $mymut, $withzyg?"\t$zygosity":"", "\t", $otherinfo, "\n";
				} else{
					print $chr, "\t", $mystart, "\t", $end, "\t", $myref, "\t", $mymut,  "\n";
				}
				#print $chr, "\t", $start, "\t", $start, "\t", $ref_allele, "\t", $mut_allele, "\t$zygosity",  "\t", $quality_score, (defined $unfiltered_read_depth)? "\t$unfiltered_read_depth" : '', (defined $MappingQuality) ? "\t$MappingQuality" : '', (defined $QualityByDepth) ? "\t$QualityByDepth" : '', $includeinfo ? "\t$otherinfo" : '', "\n";	#commented Sep 2011
				$countsnp++;
			}
			$countvar++;
		}else{
			if ($includeinfo) {
				print $chr, "\t", $start, "\t", $stop, "\t", $ref_allele, "\t", $mut_allele, $withzyg?"\t$zygosity":"", "\t", $otherinfo, "\n";
			} else{
				print $chr, "\t", $start, "\t", $stop, "\t", $ref_allele, "\t", $mut_allele,  "\n";
			}
		}
	}
	my $triallelic = $countsnp-$countti-$counttv;
	print STDERR "NOTICE: Read $countline lines and wrote ${\($counthet+$counthom)} different variants at $countvar genomic positions ($countsnp SNPs and $countindel indels)\n";
	print STDERR "NOTICE: Among ${\($counthet+$counthom+$countunknown)} different variants at $countvar positions, $counthet are heterozygotes, $counthom are homozygotes\n";
	print STDERR "NOTICE: Among $countsnp SNPs, $countti are transitions, $counttv are transversions", $triallelic?", $triallelic have more than 2 alleles\n":"\n";
	print STDERR Dumper (\%stats);
	if ($stats{no_af} && $stats{afct} && $stats{no_af}>$stats{afct} && -e "$outfile.hg19_af"){
		system ("rm $outfile.hg19_af\n");
	}
}
=head1 SYNOPSIS

 convert2annovar.pl [arguments] <variantfile>

 Optional arguments:
        -h, --help                      print help message
        -m, --man                       print complete documentation
        -v, --verbose                   use verbose output
            --format <string>		input format (default: pileup)
            --outfile <file>		output file name (default: STDOUT)
            --snpqual <float>		quality score threshold in pileup file (default: 20)
            --snppvalue <float>		SNP P-value threshold in GFF3-SOLiD file (default: 1)
            --coverage <int>		read coverage threshold in pileup file (default: 0)
            --maxcoverage <int>		maximum coverage threshold (default: none)
            --includeinfo		include supporting information in output
            --chr <string>		specify the chromosome (for CASAVA format)
            --chrmt <string>		chr identifier for mitochondria (default: M)
            --altcov <int>		alternative allele coverage threshold (for pileup format)
            --fraction <float>		minimum allelic fraction to claim a mutation (for pileup/vcf4_indel format)
            --species <string>		if human, convert chr23/24/25 to X/Y/M (for gff3-solid format)
            --filter <string>		output variants with this filter (case insensitive, for vcf4 format)
            --confraction <float>	minimum consensus indel / all indel fraction (for vcf4 format)
            --allallele			print all alleles when multiple calls are present (for vcf4 format)
            --withzyg			print zygosity when -includeinfo is used (for vcf4 format)

 Function: convert variant call file generated from various software programs 
 into ANNOVAR input format
 
 Example: convert2annovar.pl -format pileup -outfile variant.query variant.pileup
          convert2annovar.pl -format cg -outfile variant.query variant.cg
          convert2annovar.pl -format gff3-solid -outfile variant.query variant.snp.gff
          convert2annovar.pl -format soap variant.snp > variant.avinput
          convert2annovar.pl -format maq variant.snp > variant.avinput
          convert2annovar.pl -format casava -chr 1 variant.snp > variant.avinput
          convert2annovar.pl -format vcf4 variantfile > variant.avinput
          convert2annovar.pl -format vcf4 -filter pass variantfile > variant.avinput

 Version: $LastChangedDate: 2011-11-20 17:27:36 -0800 (Sun, 20 Nov 2011) $

=head1 OPTIONS

=over 8

=item B<--help>

print a brief usage message and detailed explanation of options.

=item B<--man>

print the complete manual of the program.

=item B<--verbose>

use verbose output.

=item B<--format>

the format of the input files.

=item B<--outfile>

specify the output file name. By default, output is written to STDOUT.

=item B<--snpqual>

quality score threshold in the pileup file, such that variant calls with lower 
quality scores will not be printed out in the output file. When VCF4 file is 
used, this argument works on the Quality-by-Depth measure, rather than the raw 
quality measure.

=item B<--coverage>

read coverage threshold in the pileup file, such that variants calls generated 
with lower coverage will not be printed in the output file.

=item B<--includeinfo>

specify that the output should contain additional information in the input line. 
By default, only the chr, start, end, reference allele, observed allele and 
homozygosity status are included in output files.

=item B<--chr>

specify the chromosome for CASAVA format

=item B<--chrmt>

specify the name of mitochondria chromosome (default is MT)

=item B<--altcov>

the minimum coverage of the alternative (mutated) allele to be printed out in 
output

=item B<--fraction>

specify the minimum fraction of alternative allele, to print out the mutation. 
For example, a site has 10 reads, 3 supports alternative allele. A -fraction of 
0.4 will not allow the mutation to be printed out.

=item B<--species>

specify the species from which the sequencing data is obtained. For the GFF3-
SOLiD format, when species is human, the chromosome 23, 24 and 25 will be 
converted to X, Y and M, respectively.

=item B<--filter>

for VCF4 file, only print out variant calls with this filter annotated. For 
example, if using GATK VariantFiltration walker, you will see PASS, 
GATKStandard, HARD_TO_VALIDATE, etc in the filter field. Using 'pass' as a 
filter is recommended in this case.

=item B<--confraction>

consesus indel fraction, calculated as reads supporting consensus indels divided 
by reads supporting any indels

=item B<--allallele>

print all alleles for mutations at a locus, rather than the first allele, if the 
input VCF4 file contains multiple alternative alleles for a mutation. By 
default, this option is off. When it is on, two lines will be printed out in the 
output, and both will have the same quality scores as VCF4 does not provide 
separate quality scores for individual alleles.

=back

=head1 DESCRIPTION

This program is used to convert variant call file generated from various 
software programs into ANNOVAR input format. Currently, the program can handle 
Samtools genotype-calling pileup format, Solid GFF format, Complete Genomics 
variant format, SOAP format. These formats are described below.

=over 8

=item * B<pileup format>

The pileup format can be produced by the Samtools genotyping calling subroutine. 
Note that the phrase pileup format can be used in several instances, and here I 
am only referring to the pileup files that contains the actual genotype calls. 

Using SamTools, given an alignment file in BAM format, a pileup file with 
genotype calls can be produced by the command below:

	samtools pileup -vcf ref.fa aln.bam> raw.pileup
	samtools.pl varFilter raw.pileup > final.pileup

ANNOVAR will automatically filter the pileup file so that only SNPs reaching a 
quality threshold are printed out (default is 20, use --snpqual argument to 
change this). Most likely, users may want to also apply a coverage threshold, 
such that SNPs calls from only a few reads are not considered. This can be 
achieved using the -coverage argument (default value is 0).

An example of pileup files for SNPs is shown below:

	chr1 556674 G G 54 0 60 16 a,.....,...,.... (B%A+%7B;0;%=B<:
	chr1 556675 C C 55 0 60 16 ,,..A..,...,.... CB%%5%,A/+,%....
	chr1 556676 C C 59 0 60 16 g,.....,...,.... .B%%.%.?.=/%...1
	chr1 556677 G G 75 0 60 16 ,$,.....,...,.... .B%%9%5A6?)%;?:<
	chr1 556678 G K 60 60 60 24 ,$.....,...,....^~t^~t^~t^~t^~t^~t^~t^~t^~t B%%B%<A;AA%??<=??;BA%B89
	chr1 556679 C C 61 0 60 23 .....a...a....,,,,,,,,, %%1%&?*:2%*&)(89/1A@B@@
	chr1 556680 G K 88 93 60 23 ..A..,..A,....ttttttttt %%)%7B:B0%55:7=>>A@B?B;
	chr1 556681 C C 102 0 60 25 .$....,...,....,,,,,,,,,^~,^~. %%3%.B*4.%.34.6./B=?@@>5.
	chr1 556682 A A 70 0 60 24 ...C,...,....,,,,,,,,,,. %:%(B:A4%7A?;A><<999=<<
	chr1 556683 G G 99 0 60 24 ....,...,....,,,,,,,,,,. %A%3B@%?%C?AB@BB/./-1A7?

The columns are chromosome, 1-based coordinate, reference base, consensus base, 
consensus quality, SNP quality, maximum mapping quality of the reads covering 
the sites, the number of reads covering the site, read bases and base qualities.

An example of pileup files for indels is shown below:

	seq2  156 *  +AG/+AG  71  252  99  11  +AG  *  3  8  0

ANNOVAR automatically recognizes both SNPs and indels in pileup file, and process them correctly.

=item * B<GFF3-SOLiD format>

The SOLiD provides a GFF3-compatible format for SNPs, indels and structural 
variants. A typical example file is given below:

	##gff-version 3
	##solid-gff-version 0.3
	##source-version 2
	##type DNA
	##date 2009-03-13
	##time 0:0:0
	##feature-ontology http://song.cvs.sourceforge.net/*checkout*/song/ontology/sofa.obo?revision=1.141
	##reference-file 
	##input-files Yoruban_snp_10x.txt
	##run-path 
	chr_name        AB_SOLiD SNP caller     SNP     coord   coord   1       .       .       coverage=# cov;ref_base=ref;ref_score=score;ref_confi=confi;ref_single=Single;ref_paired=Paired;consen_base=consen;consen_score=score;consen_confi=conf;consen_single=Single;consen_paired=Paired;rs_id=rs_id,dbSNP129
	1       AB_SOLiD SNP caller     SNP     997     997     1       .       .       coverage=3;ref_base=A;ref_score=0.3284;ref_confi=0.9142;ref_single=0/0;ref_paired=1/1;consen_base=G;consen_score=0.6716;consen_confi=0.9349;consen_single=0/0;consen_paired=2/2
	1       AB_SOLiD SNP caller     SNP     2061    2061    1       .       .       coverage=2;ref_base=G;ref_score=0.0000;ref_confi=0.0000;ref_single=0/0;ref_paired=0/0;consen_base=C;consen_score=1.0000;consen_confi=0.8985;consen_single=0/0;consen_paired=2/2
	1       AB_SOLiD SNP caller     SNP     4770    4770    1       .       .       coverage=2;ref_base=A;ref_score=0.0000;ref_confi=0.0000;ref_single=0/0;ref_paired=0/0;consen_base=G;consen_score=1.0000;consen_confi=0.8854;consen_single=0/0;consen_paired=2/2
	1       AB_SOLiD SNP caller     SNP     4793    4793    1       .       .       coverage=14;ref_base=A;ref_score=0.0723;ref_confi=0.8746;ref_single=0/0;ref_paired=1/1;consen_base=G;consen_score=0.6549;consen_confi=0.8798;consen_single=0/0;consen_paired=9/9
	1       AB_SOLiD SNP caller     SNP     6241    6241    1       .       .       coverage=2;ref_base=T;ref_score=0.0000;ref_confi=0.0000;ref_single=0/0;ref_paired=0/0;consen_base=C;consen_score=1.0000;consen_confi=0.7839;consen_single=0/0;consen_paired=2/2
	
Newer version of ABI BioScope now use diBayes caller, and the output file is given below:

	##gff-version 3
	##feature-ontology http://song.cvs.sourceforge.net/*checkout*/song/ontology/sofa.obo?revision=1.141
	##List of SNPs. Date Sat Dec 18 10:30:45 2010    Stringency: medium Mate Pair: 1 Read Length: 50 Polymorphism Rate: 0.003000 Bayes Coverage: 60 Bayes_Single_SNP: 1 Filter_Single_SNP: 1 Quick_P_Threshold: 0.997000 Bayes_P_Threshold: 0.040000 Minimum_Allele_Ratio: 0.150000 Minimum_Allele_Ratio_Multiple_of_Dicolor_Error: 100
	##1     chr1
	##2     chr2
	##3     chr3
	##4     chr4
	##5     chr5
	##6     chr6
	##7     chr7
	##8     chr8
	##9     chr9
	##10    chr10
	##11    chr11
	##12    chr12
	##13    chr13
	##14    chr14
	##15    chr15
	##16    chr16
	##17    chr17
	##18    chr18
	##19    chr19
	##20    chr20
	##21    chr21
	##22    chr22
	##23    chrX
	##24    chrY
	##25    chrM
	# source-version SOLiD BioScope diBayes(SNP caller)
	#Chr    Source  Type    Pos_Start       Pos_End Score   Strand  Phase   Attributes
	chr1    SOLiD_diBayes   SNP     221367  221367  0.091151        .       .       genotype=R;reference=G;coverage=3;refAlleleCounts=1;refAlleleStarts=1;refAlleleMeanQV=29;novelAlleleCounts=2;novelAlleleStarts=2;novelAlleleMeanQV=27;diColor1=11;diColor2=33;het=1;flag= 
	chr1    SOLiD_diBayes   SNP     555317  555317  0.095188        .       .       genotype=Y;reference=T;coverage=13;refAlleleCounts=11;refAlleleStarts=10;refAlleleMeanQV=23;novelAlleleCounts=2;novelAlleleStarts=2;novelAlleleMeanQV=29;diColor1=00;diColor2=22;het=1;flag= 
	chr1    SOLiD_diBayes   SNP     555327  555327  0.037582        .       .       genotype=Y;reference=T;coverage=12;refAlleleCounts=6;refAlleleStarts=6;refAlleleMeanQV=19;novelAlleleCounts=2;novelAlleleStarts=2;novelAlleleMeanQV=29;diColor1=12;diColor2=30;het=1;flag= 
	chr1    SOLiD_diBayes   SNP     559817  559817  0.094413        .       .       genotype=Y;reference=T;coverage=9;refAlleleCounts=5;refAlleleStarts=4;refAlleleMeanQV=23;novelAlleleCounts=2;novelAlleleStarts=2;novelAlleleMeanQV=14;diColor1=11;diColor2=33;het=1;flag= 
	chr1    SOLiD_diBayes   SNP     714068  714068  0.000000        .       .       genotype=M;reference=C;coverage=13;refAlleleCounts=7;refAlleleStarts=6;refAlleleMeanQV=25;novelAlleleCounts=6;novelAlleleStarts=4;novelAlleleMeanQV=22;diColor1=00;diColor2=11;het=1;flag= 
	The file conforms to standard GFF3 specifications, but the last column is solid-
	specific and it gives certain parameters for the SNP calls.

An example of the short indel format by GFF3-SOLiD is given below:

	##gff-version 3
	##solid-gff-version 0.3
	##source-version SOLiD Corona Lite v.4.0r2.0, find-small-indels.pl v 1.0.1, process-small-indels v 0.2.2, 2009-01-12 12:28:49
	##type DNA
	##date 2009-01-26
	##time 18:33:20
	##feature-ontology http://song.cvs.sourceforge.net/*checkout*/song/ontology/sofa.obo?revision=1.141
	##reference-file 
	##input-files ../../mp-results/JOAN_20080104_1.pas,../../mp-results/BARB_20071114_1.pas,../../mp-results/BARB_20080227_2.pas
	##run-path /data/results2/Yoruban-frag-indel/try.01.06/mp-w2x25-2x-4x-8x-10x/2x
	##Filter-settings: max-ave-read-pos=none,min-ave-from-end-pos=9.1,max-nonreds-4filt=2,min-insertion-size=none,min-deletion-size=none,max-insertion-size=none,max-deletion-size=none,require-called-indel-size?=T
	chr1    AB_SOLiD Small Indel Tool       deletion        824501  824501  1       .       .       del_len=1;tight_chrom_pos=824501-824502;loose_chrom_pos=824501-824502;no_nonred_reads=2;no_mismatches=1,0;read_pos=4,6;from_end_pos=21,19;strands=+,-;tags=R3,F3;indel_sizes=-1,-1;read_seqs=G3021212231123203300032223,T3321132212120222323222101;dbSNP=rs34941678,chr1:824502-824502(-),EXACT,1,/GG
	chr1    AB_SOLiD Small Indel Tool       insertion_site  1118641 1118641 1       .       .       ins_len=3;tight_chrom_pos=1118641-1118642;loose_chrom_pos=1118641-1118642;no_nonred_reads=2;no_mismatches=0,1;read_pos=17,6;from_end_pos=8,19;strands=+,+;tags=F3,R3;indel_sizes=3,3;read_seqs=T0033001100022331122033112,G3233112203311220000001002

The keyword deletion or insertion_site is used in the fourth column to indicate 
that file format.

An example of the medium CNV format by GFF3-SOLiD is given below:

	##gff-version 3
	##solid-gff-version 0.3
	##source-version SOLiD Corona Lite v.4.0r2.0, find-small-indels.pl v 1.0.1, process-small-indels v 0.2.2, 2009-01-12 12:28:49
	##type DNA
	##date 2009-01-27
	##time 15:54:36
	##feature-ontology http://song.cvs.sourceforge.net/*checkout*/song/ontology/sofa.obo?revision=1.141
	##reference-file 
	##input-files big_d20e5-del12n_up-ConsGrp-2nonred.pas.sum
	##run-path /data/results2/Yoruban-frag-indel/try.01.06/mp-results-lmp-e5/big_d20e5-indel_950_2050
	chr1    AB_SOLiD Small Indel Tool       deletion        3087770 3087831 1       .       .       del_len=62;tight_chrom_pos=none;loose_chrom_pos=3087768-3087773;no_nonred_reads=2;no_mismatches=2,2;read_pos=27,24;from_end_pos=23,26;strands=-,+;tags=F3,F3;indel_sizes=-62,-62;read_seqs=T11113022103331111130221213201111302212132011113022,T02203111102312122031111023121220311111333012203111
	chr1    AB_SOLiD Small Indel Tool       deletion        4104535 4104584 1       .       .       del_len=50;tight_chrom_pos=4104534-4104537;loose_chrom_pos=4104528-4104545;no_nonred_reads=3;no_mismatches=0,4,4;read_pos=19,19,27;from_end_pos=31,31,23;strands=+,+,-;tags=F3,R3,R3;indel_sizes=-50,-50,-50;read_seqs=T31011011013211110130332130332132110110132020312332,G21031011013211112130332130332132110132132020312332,G20321302023001101123123303103303101113231011011011
	chr1    AB_SOLiD Small Indel Tool       insertion_site  2044888 2044888 1       .       .       ins_len=18;tight_chrom_pos=2044887-2044888;loose_chrom_pos=2044887-2044889;no_nonred_reads=2;bead_ids=1217_1811_209,1316_908_1346;no_mismatches=0,2;read_pos=13,15;from_end_pos=37,35;strands=-,-;tags=F3,F3;indel_sizes=18,18;read_seqs=T31002301231011013121000101233323031121002301231011,T11121002301231011013121000101233323031121000101231;non_indel_no_mismatches=3,1;non_indel_seqs=NIL,NIL
	chr1    AB_SOLiD Small Indel Tool       insertion_site  74832565        74832565        1       .       .       ins_len=16;tight_chrom_pos=74832545-74832565;loose_chrom_pos=74832545-74832565;no_nonred_reads=2;bead_ids=1795_181_514,1651_740_519;no_mismatches=0,2;read_pos=13,13;from_end_pos=37,37;strands=-,-;tags=F3,R3;indel_sizes=16,16;read_seqs=T33311111111111111111111111111111111111111111111111,G23311111111111111111111111111111111111111311011111;non_indel_no_mismatches=1,0;non_indel_seqs=NIL,NIL

An example of the large indel format by GFF3-SOLiD is given below:

	##gff-version 3
	##solid-gff-version 0.3
	##source-version ???
	##type DNA
	##date 2009-03-13
	##time 0:0:0
	##feature-ontology http://song.cvs.sourceforge.net/*checkout*/song/ontology/sofa.obo?revision=1.141
	##reference-file 
	##input-files /data/results5/yoruban_strikes_back_large_indels/LMP/five_mm_unique_hits_no_rescue/5_point_6x_del_lib_1/results/NA18507_inter_read_indels_5_point_6x.dat
	##run-path 
	chr1    AB_SOLiD Large Indel Tool       insertion_site  1307279 1307791 1       .       .       deviation=-742;stddev=7.18;ref_clones=-;dev_clones=4
	chr1    AB_SOLiD Large Indel Tool       insertion_site  2042742 2042861 1       .       .       deviation=-933;stddev=8.14;ref_clones=-;dev_clones=3
	chr1    AB_SOLiD Large Indel Tool       insertion_site  2443482 2444342 1       .       .       deviation=-547;stddev=11.36;ref_clones=-;dev_clones=17
	chr1    AB_SOLiD Large Indel Tool       insertion_site  2932046 2932984 1       .       .       deviation=-329;stddev=6.07;ref_clones=-;dev_clones=14
	chr1    AB_SOLiD Large Indel Tool       insertion_site  3166925 3167584 1       .       .       deviation=-752;stddev=13.81;ref_clones=-;dev_clones=14

An example of the CNV format by GFF3-SOLiD if given below:

	##gff-version 3
	##solid-gff-version 0.3
	##source-version ???
	##type DNA
	##date 2009-03-13
	##time 0:0:0
	##feature-ontology http://song.cvs.sourceforge.net/*checkout*/song/ontology/sofa.obo?revision=1.141
	##reference-file 
	##input-files Yoruban_cnv.coords
	##run-path 
	chr1    AB_CNV_PIPELINE repeat_region   1062939 1066829 .       .       .       fraction_mappable=51.400002;logratio=-1.039300;copynum=1;numwindows=1
	chr1    AB_CNV_PIPELINE repeat_region   1073630 1078667 .       .       .       fraction_mappable=81.000000;logratio=-1.409500;copynum=1;numwindows=2
	chr1    AB_CNV_PIPELINE repeat_region   2148325 2150352 .       .       .       fraction_mappable=98.699997;logratio=-1.055000;copynum=1;numwindows=1
	chr1    AB_CNV_PIPELINE repeat_region   2245558 2248109 .       .       .       fraction_mappable=78.400002;logratio=-1.042900;copynum=1;numwindows=1
	chr1    AB_CNV_PIPELINE repeat_region   3489252 3492632 .       .       .       fraction_mappable=59.200001;logratio=-1.119900;copynum=1;numwindows=1
	chr1    AB_CNV_PIPELINE repeat_region   5654415 5657276 .       .       .       fraction_mappable=69.900002;logratio=1.114500;copynum=4;numwindows=1
	chr1    AB_CNV_PIPELINE repeat_region   9516165 9522726 .       .       .       fraction_mappable=65.850006;logratio=-1.316700;numwindows=2
	chr1    AB_CNV_PIPELINE repeat_region   16795117        16841025        .       .       .       fraction_mappable=44.600002;logratio=1.880778;copynum=7;numwindows=9

The keyword repeat_region is used here, although it actually refers to CNVs.

An example of the inversion format by GFF3-SOLiD is given below:

	##gff-version 3
	##solid-gff-version 0.2
	##generated by SOLiD inversion tool
	chr10   AB_SOLiD        inversion       46443107        46479585        268.9   .       .       left=chr10:46443107-46443146;right=chr10:46479583-46479585;leftscore=295.0;rightscore=247.0;count_AAA_further_left=117;count_AAA_left=3;count_AAA_right=3;count_AAA_further_right=97;left_min_count_AAA=chr10:46443107-46443112;count_AAA_min_left=0;count_AAA_max_left=3;right_min_count_AAA=chr10:46479585-46479585;count_AAA_min_right=1;count_AAA_max_right=3;homozygous=UNKNOWN
	chr4    AB_SOLiD        inversion       190822813       190850112       214.7   .       .       left=chr4:190822813-190822922;right=chr4:190850110-190850112;leftscore=140.0;rightscore=460.0;count_AAA_further_left=110;count_AAA_left=78;count_AAA_right=74;count_AAA_further_right=77;left_min_count_AAA=chr4:190822813-190822814;count_AAA_min_left=69;count_AAA_max_left=77;right_min_count_AAA=chr4:190850110-190850112;count_AAA_min_right=74;count_AAA_max_right=74;homozygous=NO
	chr6    AB_SOLiD        inversion       168834969       168837154       175.3   .       .       left=chr6:168834969-168835496;right=chr6:168836643-168837154;leftscore=185.4;rightscore=166.2;count_AAA_further_left=67;count_AAA_left=43;count_AAA_right=40;count_AAA_further_right=59;left_min_count_AAA=chr6:168835058-168835124,chr6:168835143-168835161,chr6:168835176-168835181,chr6:168835231-168835262;count_AAA_min_left=23;count_AAA_max_left=29;right_min_count_AAA=chr6:168836643-168836652;count_AAA_min_right=23;count_AAA_max_right=31;homozygous=NO

The program should be able to recognize all the above GFF3-SOLiD format 
automatically, and handle them accordingly.

=item * B<Complete Genomics format>

This format is provided by the Complete Genomics company to their customers. The 
file var-[ASM-ID].tsv.bz2 includes a description of all loci where the assembled 
genome differs from the reference genome.

An example of the Complete Genomics format is shown below:

	#BUILD  1.5.0.5
	#GENERATED_AT   2009-Nov-03 19:52:21.722927
	#GENERATED_BY   dbsnptool
	#TYPE   VAR-ANNOTATION
	#VAR_ANN_SET    /Proj/Pipeline/Production_Data/REF/HUMAN-F_06-REF/dbSNP.csv
	#VAR_ANN_TYPE   dbSNP
	#VERSION        0.3
	
	>locus  ploidy  haplotype       chromosome      begin   end     varType reference       alleleSeq       totalScore      hapLink xRef
	1       2       all     chr1    0       959     no-call =       ?                       
	2       2       all     chr1    959     972     =       =       =                       
	3       2       all     chr1    972     1001    no-call =       ?                       
	4       2       all     chr1    1001    1008    =       =       =                       
	5       2       all     chr1    1008    1114    no-call =       ?                       
	6       2       all     chr1    1114    1125    =       =       =                       
	7       2       all     chr1    1125    1191    no-call =       ?                       
	8       2       all     chr1    1191    1225    =       =       =                       
	9       2       all     chr1    1225    1258    no-call =       ?                       
	10      2       all     chr1    1258    1267    =       =       =                       
	12      2       all     chr1    1267    1275    no-call =       ?                       
	13      2       all     chr1    1275    1316    =       =       =                       
	14      2       all     chr1    1316    1346    no-call =       ?                       
	15      2       all     chr1    1346    1367    =       =       =                       
	16      2       all     chr1    1367    1374    no-call =       ?                       
	17      2       all     chr1    1374    1388    =       =       =                       
	18      2       all     chr1    1388    1431    no-call =       ?                       
	19      2       all     chr1    1431    1447    =       =       =                       
	20      2       all     chr1    1447    1454    no-call =       ?                       

The following information is provided in documentation from Complete Genomics, that describes the var-ASM format.

	1. locus. Identifier of a particular genomic locus
	2. ploidy. The ploidy of the reference genome at the locus (= 2 for autosomes, 2 for pseudoautosomal regions on the sex chromosomes, 1 for males on the non-pseudoautosomal parts of the sex chromosomes, 1 for mitochondrion, '?' if varType is 'no-ref' or 'PAR-called-in-X'). The reported ploidy is fully determined by gender, chromosome and location, and is not inferred from the sequence data.
	3. haplotype. Identifier for each haplotype at the variation locus. For diploid genomes, 1 or 2. Shorthand of 'all' is allowed where the varType field is one of 'ref', 'no-call', 'no-ref', or 'PAR-called-in-X'. Haplotype numbering does not imply phasing; haplotype 1 in locus 1 is not necessarily in phase with haplotype 1 in locus 2. See hapLink, below, for phasing information.
	4. chromosome. Chromosome name in text: 'chr1','chr2', ... ,'chr22','chrX','chrY'. The mitochondrion is represented as 'chrM'. The pseudoautosomal regions within the sex chromosomes X and Y are reported at their coordinates on chromosome X.
	5. begin. Reference coordinate specifying the start of the variation (not the locus) using the half-open zero-based coordinate system. See section 'Sequence Coordinate System' for more information.
	6. end. Reference coordinate specifying the end of the variation (not the locus) using the half-open zero-based coordinate system. See section 'Sequence Coordinate System' for more information.
	7. varType. Type of variation, currently one of:
		snp: single-nucleotide polymorphism
		ins: insertion
		del: deletion
		sub: Substitution of one or more reference bases with the bases in the allele column
		'ref' : no variation; the sequence is identical to the reference sequence on the indicated haplotype
		no-call-rc: 'no-call reference consistent 'one or more bases are ambiguous, but the allele is potentially consistent with the reference
		no-call-ri: 'no-call reference inconsistent' one or more bases are ambiguous, but the allele is definitely inconsistent with the reference
		no-call: an allele is completely indeterminate in length and composition, i.e. alleleSeq = '?'
		no-ref: the reference sequence is unspecified at this locus.
		PAR-called-in-X: this locus overlaps one of the pseudoautosomal regions on the sex chromosomes. The called sequence is reported as diploid sequence on Chromosome X; on chromosome Y the sequence is reported as varType = 'PAR-called-in-X'.
	8. reference. The reference sequence for the locus of variation. Empty when varType is ins. A value of '=' indicates that the user must consult the reference for the sequence; this shorthand is only used in regions where no haplotype deviates from the reference sequence.
	9. alleleSeq. The observed sequence at the locus of variation. Empty when varType is del. '?' isused to indicate 0 or more unknown bases within the sequence; 'N' is used to indicate exactly one unknown base within the sequence.'=' is used as shorthand to indicate identity to the reference sequence for non-variant sequence, i.e. when varType is 'ref'.
	10. totalScore. A score corresponding to a single variation and haplotype, representing the confidence in the call.
	11. hapLink. Identifier that links a haplotype at one locus to haplotypes at other loci. Currently only populated for very proximate variations that were assembled together. Two calls that share a hapLink identifier are expected to be on the same haplotype,
	12. xRef. Field containing external variation identifiers, currently only populated for variations corroborated directly by dbSNP. Format: dbsnp:[rsID], with multiple entries separated by the semicolon (;).

In older versions of the format specification, the sub keyword used to be insdel 
keyword. ANNOVAR takes care of this.

=item * B<SOAPsnp format>

An example of the SOAP SNP caller format is shown below:

	chr8  35782  A  R  1  A  27  1  2  G  26  1  2  5   0.500000  2.00000  1  5   
	chr8  35787  G  R  0  G  25  4  6  A  17  2  4  10  0.266667  1.60000  0  5   

The following information is provided in documentation from BGI who developed 
SOAP suite. It differs slightly from the description at the SOAPsnp website, and 
presumably the website is outdated.

	Format description:(left to right)
	1. Chromosome name
	2. Position of locus
	3. Nucleotide at corresponding locus of reference sequence
	4. Genotype of sequencing sample
	5. Quality value
	6. nucleotide with the highest probability(first nucleotide)
	7. Quality value of the nucleotide with the highest probability
	8. Number of supported reads that can only be aligned to this locus 
	9. Number of all supported reads that can be aligned to this locus
	10. Nucleotide with higher probability 
	11. Quality value of nucleotide with higher probability 
	12. Number of supported reads that can only be aligned to this locus 
	13. Number of all supported reads that can be aligned to this locus 
	14. Total number of reads that can be aligned to this locus 
	15. Order and quality value
	16. Estimated copy number for this locus 
	17. Presence of this locus in the dbSNP database. 1 refers to presence and 0 refers to inexistence
	18. The distance between this locus and another closest SNP

=item * B<SOAPindel format>

The current version of ANNOVAR handles SoapSNP and SoapIndel automatically via a 
single argument '--format soap'. An example of SOAP indel caller format is shown 
below:

	chr11   44061282        -       +2      CT      Hete
	chr11   45901572        +       +1      C       Hete
	chr11   48242562        *       -3      TTC     Homo
	chr11   57228723        *       +4      CTTT    Homo
	chr11   57228734        *       +4      CTTT    Homo
	chr11   57555685        *       -1      C       Hete
	chr11   61482191        -       +3      TCC     Hete
	chr11   64608031        *       -1      T       Homo
	chr11   64654936        *       +1      C       Homo
	chr11   71188303        +       -1      T       Hete
	chr11   75741034        +       +1      T       Hete
	chr11   76632438        *       +1      A       Hete
	chr11   89578266        *       -2      AG      Homo
	chr11   104383261       *       +1      T       Hete
	chr11   124125940       +       +4      CCCC    Hete
	chr12   7760052 *       +1      T       Homo
	chr12   8266049 *       +3      ACG     Homo

I do not see a documentation describing this format yet as of September 2010.

=item B<--SOAPsv format>

An example is given below:

	Chr2 Deletion 42894 43832 43167 43555 388 0-0-0 FR 41

An explanation of the structural variation format is given below:

	Format description (from left to right)
	1. Chromosome name
	2. Type of structure variation
	3. Minimal value of start position in cluster
	4. Maximal value of end position in cluster
	5. Estimated start position of this structure variation
	6. Estimated end position of this structure variation
	7. Length of SV
	8. Breakpoint of SV (only for insertion)
	9. Unusual matching mode (F refers to align with forward sequence, R refers
	to align with reverse
	sequence)
	10. number of paired-end read which support this structure variation

=item * B<MAQ format>

MAQ can perform alignment and generate genotype calls, including SNP calls and 
indel calls. The format is described below:

For indel header: The output is TAB delimited with each line consisting of chromosome, start 
position, type of the indel, number of reads across the indel, size of the indel 
and inserted/deleted nucleotides (separated by colon), number of indels on the 
reverse strand, number of indels on the forward strand, 5' sequence ahead of the 
indel, 3' sequence following the indel, number of reads aligned without indels 
and three additional columns for filters.

An example is below:

	chr10   110583  -       2       -2:AG   0       1       GCGAGACTCAGTATCAAAAAAAAAAAAAAAAA        AGAAAGAAAGAAAAAGAAAAAAATAGAAAGAA        1       @2,     @72,   @0,
	chr10   120134  -       8       -2:CA   0       1       CTCTTGCCCGCTCACACATGTACACACACGCG        CACACACACACACACACATCAGCTACCTACCT        7       @65,62,61,61,45,22,7,   @9,12,13,13,29,52,67,   @0,0,0,0,0,0,0,
	chr10   129630  -       1       -1:T    1       0       ATGTTGTGACTCTTAATGGATAAGTTCAGTCA        TTTTTTTTTAGCTTTTAACCGGACAAAAAAAG        0       @       @      @
	chr10   150209  -       1       4:TTCC  1       0       GCATATAGGGATGGGCACTTTACCTTTCTTTT        TTCCTTCCTTCCTTCCTTCCCTTTCCTTTCCT        0       @       @      @
	chr10   150244  -       2       -4:TTCT 0       1       CTTCCTTCCTTCCTTCCCTTTCCTTTCCTTTC        TTCTTTCTTTCTTTCTTTCTTTTTTTTTTTTT        0       @       @      @
	chr10   159622  -       1       3:AGG   0       1       GAAGGAGGAAGGACGGAAGGAGGAAGGAAGGA        AGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGA        0       @       @      @
	chr10   206372  -       2       2:GT    1       0       ATAATAGTAACTGTGTATTTGATTATGTGTGC        GTGTGTGTGTGTGTGTGTGTGTGTGCGTGCTT        1       @37,    @37,   @8,
	chr10   245751  -       11      -1:C    0       1       CTCATAAATACAAGTCATAATGAAAGAAATTA        CCACCATTTTCTTATTTTCATTCATTTTTAGT        10      @69,64,53,41,30,25,22,14,5,4,   @5,10,21,33,44,49,52,60,69,70,  @0,0,0,0,0,0,0,0,0,0,
	chr10   253066  -       1       2:TT    0       1       TATTGATGAGGGTGGATTATACTTTAGAACAC        TATTCAAACAGTTCTTCCACATATCTCCCTTT        0       @       @      @
	chr10   253455  -       2       -3:AAA  1       0       GTTGCACTCCAGCCTGGCGAGATTCTGTCTCC        AAAAAAAAAAAAAAAAATTGTTGTGAAATACA        1       @55,    @19,   @4,

For snp output file: Each line consists of chromosome, position, reference base, 
consensus base, Phred-like consensus quality, read depth, the average number of 
hits of reads covering this position, the highest mapping quality of the reads 
covering the position, the minimum consensus quality in the 3bp flanking regions 
at each side of the site (6bp in total), the second best call, log likelihood 
ratio of the second best and the third best call, and the third best call.

An example is below:

	chr10   83603   C       T       28      12      2.81    63      34      Y       26      C
	chr10   83945   G       R       59      61      4.75    63      62      A       47      G
	chr10   83978   G       R       47      40      3.31    63      62      A       21      G
	chr10   84026   G       R       89      22      2.44    63      62      G       49      A
	chr10   84545   C       T       54      9       1.69    63      30      N       135     N
	chr10   85074   G       A       42      5       1.19    63      38      N       108     N
	chr10   85226   A       T       42      5       1.00    63      42      N       107     N
	chr10   85229   C       T       42      5       1.00    63      42      N       112     N
	chr10   87518   A       G       39      4       3.25    63      38      N       9       N
	chr10   116402  T       C       39      4       1.00    63      38      N       76      N


=item * B<CASAVA format>

An example of Illumina CASAVA format is given below:

	#position       A       C       G       T       modified_call   total   used    score           reference       type
	14930   3       0       8       0       GA      11      11      29.10:11.10             A       SNP_het2
	14933   4       0       7       0       GA      11      11      23.50:13.70             G       SNP_het1
	14976   3       0       8       0       GA      11      11      24.09:9.10              G       SNP_het1
	15118   2       1       4       0       GA      8       7       10.84:6.30              A       SNP_het2

An example of the indels is given below:

	# ** CASAVA depth-filtered indel calls **
	#$ CMDLINE /illumina/pipeline/install/CASAVA_v1.7.0/libexec/CASAVA-1.7.0/filterIndelCalls.pl--meanReadDepth=2.60395068970547 --indelsCovCutoff=-1 --chrom=chr1.fa /data/Basecalls/100806_HARMONIAPILOT-H16_0338_A2065HABXX/Data/Intensities/BaseCalls/CASAVA_PE_L2/Parsed_14-08-10/chr1.fa/Indel/varling_indel_calls_0000.txt /data/Basecalls/100806_HARMONIAPILOT-H16_0338_A2065HABXX/Data/Intensities/BaseCalls/CASAVA_PE_L2/Parsed_14-08-10/chr1.fa/Indel/varling_indel_calls_0001.txt /data/Basecalls/100806_HARMONIAPILOT-H16_0338_A2065HABXX/Data/Intensities/BaseCalls/CASAVA_PE_L2/Parsed_14-08-10/chr1.fa/Indel/varling_indel_calls_0002.txt /data/Basecalls/100806_HARMONIAPILOT-H16_0338_A2065HABXX/Data/Intensities/BaseCalls/CASAVA_PE_L2/Parsed_14-08-10/chr1.fa/Indel/varling_indel_calls_0003.txt /data/Basecalls/100806_HARMONIAPILOT-H16_0338_A2065HABXX/Data/Intensities/BaseCalls/CASAVA_PE_L2/Parsed_14-08-10/chr1.fa/Indel/varling_indel_calls_0004.txt
	#$ CHROMOSOME chr1.fa
	#$ MAX_DEPTH undefined
	#
	#$ COLUMNS pos CIGAR ref_upstream ref/indel ref_downstream Q(indel) max_gtype Q(max_gtype) max2_gtype bp1_reads ref_reads indel_reads other_reads repeat_unit ref_repeat_count indel_repeat_count
	948847  1I      CCTCAGGCTT      -/A     ATAATAGGGC      969     hom     47      het     22      0       16      6       A       1       2
	978604  2D      CACTGAGCCC      CT/--   GTGTCCTTCC      251     hom     20      het     8       0       4       4       CT      1       0
	1276974 4I      CCTCATGCAG      ----/ACAC       ACACATGCAC      838     hom     39      het     18      0       14      4       AC      2       4
	1289368 2D      AGCCCGGGAC      TG/--   GGAGCCGCGC      1376    hom     83      het     33      0       25      9       TG      1       0

=item * B<VCF4 format>

VCF4 can be used to describe both population-level variation information, or for 
reads derived from a single individual.

One example of the indel format for one individual is given below:

	##fileformat=VCFv4.0
	##IGv2_bam_file_used=MIAPACA2.alnReAln.bam
	##INFO=<ID=AC,Number=2,Type=Integer,Description="# of reads supporting consensus indel/any indel at the site">
	##INFO=<ID=DP,Number=1,Type=Integer,Description="total coverage at the site">
	##INFO=<ID=MM,Number=2,Type=Float,Description="average # of mismatches per consensus indel-supporting read/per reference-supporting read">
	##INFO=<ID=MQ,Number=2,Type=Float,Description="average mapping quality of consensus indel-supporting reads/reference-supporting reads">
	##INFO=<ID=NQSBQ,Number=2,Type=Float,Description="Within NQS window: average quality of bases from consensus indel-supporting reads/from reference-supporting reads">
	##INFO=<ID=NQSMM,Number=2,Type=Float,Description="Within NQS window: fraction of mismatching bases in consensus indel-supporting reads/in reference-supporting reads">
	##INFO=<ID=SC,Number=4,Type=Integer,Description="strandness: counts of forward-/reverse-aligned indel-supporting reads / forward-/reverse-aligned reference supporting reads">
	##IndelGenotyperV2=""
	##reference=hg18.fa
	##source=IndelGenotyperV2
	#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Miapaca_trimmed_sorted.bam      
	chr1    439     .       AC      A       .       PASS    AC=5,5;DP=7;MM=7.0,3.0;MQ=23.4,1.0;NQSBQ=23.98,25.5;NQSMM=0.04,0.0;SC=2,3,0,2   GT      1/0
	chr1    714048  .       T       TCAAC   .       PASS    AC=3,3;DP=9;MM=3.0,7.1666665;MQ=1.0,10.833333;NQSBQ=23.266666,21.932203;NQSMM=0.0,0.15254237;SC=3,0,3,3 GT      0/1
	chr1    714049  .       G       GC      .       PASS    AC=3,3;DP=9;MM=3.0,7.1666665;MQ=1.0,10.833333;NQSBQ=23.233334,21.83051;NQSMM=0.0,0.15254237;SC=3,0,3,3  GT      0/1
	chr1    813675  .       A       AATAG   .       PASS    AC=5,5;DP=8;MM=0.4,1.0;MQ=5.0,67.0;NQSBQ=25.74,25.166666;NQSMM=0.0,0.033333335;SC=4,1,1,2       GT      0/1
	chr1    813687  .       AGAGAGAGAGAAG   A       .       PASS    AC=5,5;DP=8;MM=0.4,1.0;MQ=5.0,67.0;NQSBQ=24.54,25.2;NQSMM=0.02,0.06666667;SC=4,1,1,2    GT      1/0


=back

The code was written by Dr. Kai Wang and modified by Dr. Germ�n Gast�n Leparc. 
Various users have provided sample input files for many SNP callin software, for 
the development of conversion subroutines. We thank these users for their 
continued support to improve the functionality of the script.

For questions or comments, please contact kai@openbioinformatics.org.

=cut
