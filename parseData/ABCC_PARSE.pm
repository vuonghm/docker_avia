#!/usr/local/bin/perl
 umask(0000);
sub color{
	#establishing the color pallate for circos
	my %color=(
		'ns_exonic'=> ['0.1','red','circle'],
		'syn_exonic'=>['0.1','vlred','circle'],
		'exonic' => ['0.1','red','circle'],
		'splicing' => ['0.2','blue','circle'],
		'ncrna'=>['0.3','green','circle'],
		'upstream'=>['0.4','orange','triangle'],
		'intron'=>['0.6','grey','rectangle'],
		'utr'=>['0.5','yellow','circle'],
		'intergenic'=>['0.7','white','rectangle'],
		'utr5'=>['0.5','black','circle'],
		'utr3'=>['0.5','yellow','circle'],
	);
	#add new nonB subtypes here
	my %nonb=(
		'A_Phased_Repeat'=> ['green','triangle'],
		'Direct_Repeat'=>['blue','circle'],
		'G_Quadruplex_Motif'=>['red','rectangle'],
		'Inverted_Repeat'=>['orange','triangle'],
		'Mirror_Repeat'=>['vvlblue','rectangle'],
		'Short_Tandem_Repeat'=>['lgrey','circle'],
		'Z_DNA_Motif'=>['black','triangle']
	);
	return (\%color,\%nonb);
}
sub readCGI{
	
}
sub reverse_colorcode{
	my $id=shift;#value, color or shape
	my $type=shift;#nonb or color
	if ($id=~/[\d\.]+/){
		$myidtype=0;
	}elsif ($id=~/(triangle|circle|rectangle)/i){
		$myidtype=2;
	}else{
		$myidtype=1;
	}
	my %href;
	if ($type=~/nonb/){
		$myidtype--;
		die "ERR, You cannot use this datatype ($id) for nonB\n";
		return %href;#return empty
	}
	my ($colorhref,$nonbhref)=color();
	my $tmphref=$colorhref;
	if ($type=~/nonb/i){
		$tmphref=$nonbhref;
	}
	foreach my $key (keys %${$tmphref}){
		my $revkey=$$tmphref{$key}[$myidtype];
		$href{$myidtype}=$key;
	}
	return \%href;	
}

sub _IUB{
	%iub=(
		'A'=>'A','AA'=>'A',
		'C'=>'C','CC'=>'C',
		'T'=>'T','TT'=>'T',
		'G'=>'G','GG'=>'G',
		'AG'=>'R', 'GA'=>'R',
		'CT'=>'Y', 'TC'=>'Y',
		'GT'=>'K','TG'=>'K',
		'AC'=>'M','CA'=>'M',
		'AT'=>'W','TA'=>'W',
		'GC'=>'S','CG'=>'S'
	);
	%reviub=(
		'A'=>'A',
		'C'=>'C',
		'T'=>'T',
		'G'=>'G',
		'R'=>'AG',
		'Y'=>'CT',
		'K'=>'GT',
		'M'=>'CA',
		'W'=>'TA',
		'S'=>'CG'
	);
	return;
}
1;
