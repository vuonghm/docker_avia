#!/usr/bin/perl
# use strict;
use FindBin;
use Getopt::Std;
use Data::Dumper;
use FindBin;
use Tie::File;
use lib "$FindBin::Bin/../perl_modules";#config
use lib "$FindBin::Bin/";
use Config::IniFiles;
use lib "/bioinfoC/hue/annovar/current/Parallel-ForkManager-0.7.9/lib";
use Parallel::ForkManager;
use lib "/bioinfoC/AVA/prod_scripts/util/perl_modules";#parse and sendmail
use Mail::Sendmail;
use MIME::Lite;
use parse;
use vars qw($opt_f $opt_c $opt_d $opt_s $opt_F $opt_N);
 getopts("f:o:c:d:F:saN");
 umask(0000);
 #for testing sendMsg


my $bin=`dirname $0`;chomp $bin;
my $preconfigured_path="/SeqIdx/circosdb/conf/predefined/";#preconfigured conf files (if doesn't exist)
my $preconfigured_data="/SeqIdx/circosdb/data/human/"; #preconfigured data paths (if doesn't exist)
my $preconfigured_plots="/SeqIdx/circosdb/plots/human/"; #preconfigured png files 
my $usage=qq{
	$0 
		[REQUIRED]
			-f/F <ANNOVAR filename or UNannotated ANNOVAR input file>
			-d <workdir>  #everything will be relative to this workdir
		[OPTIONAL]
			-c <config file from web>
				DEFAULT: config.ini
			-s <if specified, separate each into it's own circos png file>
			-N <if specified, then notifies the user that the run is completed
};	
my $workdir=$opt_d;
die "$usage\n" if ( !defined $opt_d || (!defined $opt_f && !defined $opt_F ) );
$opt_c||="config.ini";
if ($opt_c!~/^\//){
	$opt_c="$workdir/$opt_c";
}

if (defined $opt_f && $opt_f!~/^\//){
	$opt_f=$workdir . "/$opt_f";
}elsif (defined $opt_F&& $opt_F!~/^\//){
	$opt_F=$workdir."/$opt_F";
}
my %cfg;

if (-e $opt_c){
	if ($opt_c=~/abcc/){#web config
		writeConf();
	}else{
		tie %cfg , 'Config::IniFiles', (-file =>"$opt_c");
	}
}
my $pngid="0";
mkdir ("$workdir/viz",0777) if (! -e "$workdir/viz");
chdir ("$workdir/viz");
my (%dbs,%sub);
open (LOG,"| tee runviz.log") or die "cannot open runviz.log\n";
open (STDERR , '>&LOG') or die "cannot redirect stderr to log\n";
print LOG "[INFO] Started at ". `date` ;
print STDERR "[ERR] Outputing stderr to logfile\n";
open (CIRCOS,"|tee circos.input") or die "Cannot open circos.input for writing\n";
my $addon;
my $usegenelist=`grep abcc_genelists $opt_c | cut -f2 -d '='`;chomp $usegenelist;
my $runData=0;#this is the flag if the user specified only preconfigured viz; runs user data anyway
if (defined $opt_f){
	my %dbmap;
	print LOG "working on $opt_f\n";
	print LOG "[INFO] Read in $opt_c\n";
	print LOG "[INFO] In $workdir/viz\n";
	if ($usegenelist || (exists $cfg{'UserInfo'}{'genelists'} && $cfg{'UserInfo'}{'genelists'}) || (exists $cfg{'UserInfo'}{'filter_file'} && $cfg{'UserInfo'}{'filter_file'})){
		print LOG "[INFO] ($usegenelist)Using $cfg{'UserInfo'}{'genelists'},$cfg{'UserInfo'}{'filter_file'} and $cfg{'UserInfo'}{'abcc_genelists'}\n";
		system ("rm genelist\n") if (-e "genelist");
		if (exists $cfg{'UserInfo'}{'genelists'} &&  -e $cfg{'UserInfo'}{'genelists'}){
			system ("cat $cfg{'UserInfo'}{'genelists'} >> .genelist\n");
		}
		if ($usegenelist && -e "/SeqIdx/circosdb/genelists/$usegenelist.txt"){
			system ("cat /SeqIdx/circosdb/genelists/$usegenelist.txt >> .genelist\n");
		}
		if (exists $cfg{'UserInfo'}{'filter_file'} && -e $cfg{'UserInfo'}{'filter_file'}){
			system ("cat $cfg{'UserInfo'}{'filter_file'} >>.genelist\n");
		}
		system ("sort -u .genelist> genelist;rm .genelist\n");
		$usegenelist=' -l genelist';
		system ("ln -s genelist viz_genelist\n");
	}else{
		$usegenelist='';
	}
	if (-e "$workdir/searchTheseDBs.txt" && !-z "$workdir/searchTheseDBs.txt"){
		print LOG "Reading searchTheseDBs.txt\n";
		open (FILE,"<$workdir/searchTheseDBs.txt") or die "cannot open or read $workdir/searchTheseDBs.txt\n";
		while (<FILE>){
			my $db_dir=`dirname $_`;chomp $db_dir;
			my $name=`basename $_`;chomp $name;
			$name=~s/^(hg1\d|mm\d+)\_//g;
			$name=~s/.txt//g;
			print "working on $name...$db_dir\n";
			if ($name=~/hom_only/){##AVIA enforced;cgi-pop in AVIA
				#subtractive
				$sub{$name}=-1;
			}elsif (exists $cfg{'UserInfo'}{'groupfile'}){
				$sub{$name}=$cfg{'UserInfo'}{'groupfile'};
				die "haven't coded yet. This is using the subtractive method and the user's input file\n";
			}elsif (exists $cfg{'UserInfo'}{"circos_type_$name"} && $cfg{'UserInfo'}{"circos_type_$name"}=~/on/){
				$dbs{$name}=-1;
			}elsif ($db_dir =~/$cfg{UserInfo}{label}/ || $db_dir=~/mnt\/fr-s-abcc-avia\d/){
				#Note to self:these are user custom annotations and  if circos was selected are in the config.ini file as:
				#circos_type_userdefined_db1=on
				#userdefined_annotdb1=20131115145213.userdb1_tassdb.txt

				## in searchTheseDbs they will be renamed as "hg19_userdb1_tassdb.txt"
				#$name= userdb1_tassdb ##in this example above

				##find out if the user requested circos
				my $found=`grep $name $opt_c`;chomp $found;
				if ($found=~/userdefined_annotdb(\d+)/ && exists $cfg{'UserInfo'}{"circos_type_userdefined_db$1"}){
					$dbs{$name}=-100;
					$dbmap{$name}=$1;
				}
			}
		}
	}
	die "Your file does not exists $opt_f\n" if (!-e "$opt_f");
	system ("rm forthelegend.txt\n") if ( -e "forthelegend.txt");
	
	my @header=split("\t",`head -n1 $opt_f`);chomp $header[$#header];
	if (join("",@header)!~/(siftv|ljb2\d)/ && -e "../subtractive_headers.out"){
		@header=split("\t",`head -n1 ../subtractive_headers.out`);chomp $header[$#header];
	}
	##get all information about subtractive analysis
	my $finder="^\w{0,}\t";
	for (my $i=0;$i<=$#header;$i++){
		if ($header[$i]=~/#{0,}(.*)/){
			my $name=$1;
			if (exists $sub{$name}){
				$sub{$name}=$i;
				$finder.="\-\\t";
			}elsif (exists $dbs{$name}){
				$dbs{$name}=$i;
				$finder.="\\S+\\t";
			}else{
				$finder.="\\S+\\t";
			}
		}#ignore else
	}
	if (!-e "../subtractive_headers.out"){
		_system ("head -n1 $opt_f> ../subtractive_headers.out");
	}else{
		print "../subtractive_headers.out exists already!";
	}
	if (keys %sub>0){
		_system ("cat ../subtractive_headers.out> subtractive.out");
		_system ("echo ''>>subtractive.out");
		_system ("grep -P '$finder' $opt_f | grep -e 'exonic' -e 'splic' -e 'UTR' -e intron  -i |grep -ve 'ncRNA_' -i  >> subtractive.out");#all the variants without any mutations in the normal population
	}else{
		_system ("cat ../subtractive_headers.out > subtractive.out");
		_system ("echo ''>>subtractive.out");
		_system ("grep -e 'exonic' -e 'splic' -e 'UTR' -e intron -e  -i $opt_f |grep -ve 'ncRNA_' -i   >> subtractive.out");
	}
	my $exe="perl $bin/generate_genelists_byANNOVAR_fdi.pl -k -i subtractive.out ";#base, other elements added
	if (-e "genelist" && !-z "genelist"){
		$exe.=" -l genelist ";
		if (exists $cfg{'UserInfo'}{'genelist_type'} && $cfg{'UserInfo'}{'genelist_type'}=~/highlight/i){
			$exe.= "-z ";
		}
	}
	
	my ($success,$type);
	my $resolution=$cfg{'UserInfo'}{'resolution'};
	$resolution ||= '10k';
	my $true_pop= $cfg{'UserInfo'}{'cgi_pop'};
	
	foreach my $key (keys (%{$cfg{'UserInfo'}})){
		if ($cfg{'UserInfo'}{$key}=~/(on|1)/ && $key=~/(preconfigure_db_\S+)/  ){
			$dbs{$1}="preconfigure";
		}
	}
	my $MAX_PROCESSES=1;
#OLIVE	my $pm = new Parallel::ForkManager($MAX_PROCESSES); 
#OLIVE	print LOG "About to fork processes in $0($MAX_PROCESSES)\n";
	DB: foreach my $db (keys %dbs){
		print "found $db...$dbs{$db}\n";
		next if ($dbs{$db}<0);#this should means that it does not exist in the header file of the annovar file; annotation failed or could not be aggregated
		my ($xfile,$precfg);
#OLIVE		$pm->start and next;
		$addon="";
		if ($db=~/preconfigure_db_(\S+)/){
			$db=$1;$precfg=1;
			$success=1;#indicates that the executable for non-config files does not need to be run
			$type='preconfigure_scatter';#gene model tiles; all others scatter
			if ($db=~/genes/ ){ #it is the gene model  or it doesn't have a population but requests that population be used
				$type='preconfigure_tile';
				$xfile="$preconfigured_data/genes.human.hg19.txt";
			}elsif ($db=~/(.*Rpt)$/){
				$type='preconfigure_line';
				$xfile="$preconfigured_data/$1.human.10k.hg19.txt";
			}elsif ($cfg{'UserInfo'}{"key_$db"}=~/pop/ && $true_pop ){#requests population and has a population
				$xfile="$preconfigured_data/CGI_Population/$resolution/hg19_$true_pop.$db.$resolution.txt";
			}elsif ($cfg{'UserInfo'}{"key_$db"}=~/all/ || !$true_pop){
				$xfile="$preconfigured_data/CGI_Population/$resolution/hg19_CGI-genes.$db.$resolution.txt";
			}elsif (!exists $cfg{'UserInfo'}{"key_$db"} || !$cfg{'UserInfo'}{"key_$db"}){#if it does not exist or if it is empty
				$success="0";
			}
			$success=0 if (! -e $xfile || !$xfile);#if the file doesn't exist; do not include in the report
		}elsif (exists $cfg{'UserInfo'}{"viz_type_$db"}  && $cfg{'UserInfo'}{"viz_type_$db"}){#if exists and is not null
			$type=$cfg{'UserInfo'}{"viz_type_$db"};
			$xfile="$workdir/viz/COUNT$db.$resolution.txt";
			if ($cfg{'UserInfo'}{'genelist_count'}=~/multi/i){#this is not the one gene counted version; this is the one variant counts version
				$addon.= " -a " ;
			}
			if ($type=~/hist/i){
					$addon.=" -x ";
			}elsif ($type=~/text/i){
				$xfile="$workdir/viz/COUNT$db.txt.genes_label.txt";
			}
		}elsif (exists $cfg{'UserInfo'}{"viz_type_userdefined_db$dbmap{$db}"}){
			$type=$cfg{'UserInfo'}{"viz_type_userdefined_db$dbmap{$db}"};
			$xfile="$workdir/viz/COUNT$db.$resolution.txt";
			if ($cfg{'UserInfo'}{'genelist_count'}=~/multi/i){#this is not the one gene counted version; this is the one variant counts version
				$addon.= " -a " ;
			}
			if ($type=~/hist/i){
					$addon.=" -x ";
			}elsif ($type=~/text/i){
				$xfile="$workdir/viz/COUNT$db.txt.genes_label.txt";
			}

		}else{
			$type="text";#defaults to text for subtractive analysis
			$xfile="$workdir/viz/COUNT$db.txt.genes_label.txt";
		}
		if ($precfg){
		}elsif ($db=~/((s)ift|non(B))/i){
			my $id=lc($2).lc($3);$runData++;
			print "[INFO] $exe  -t $id $addon -o COUNT$db \n";
			$success=_system ("$exe  -t $id $addon -o COUNT$db \n",0);
		}else{$runData++;
			print "[INFO] $exe  -t o -h $db -o COUNT$db $addon \n";
			$success=_system ("$exe  -t o -h $db -o COUNT$db $addon \n",0);
		}
		
		if ($opt_s){
			$pngid++;
		}
		if ($success){
			 $pngid=$db if ($opt_s);
			 $pngid="pre".$db if ($precfg);
			 print CIRCOS "USER_$pngid\t$xfile\t$type\n";
			 if ($cfg{'UserInfo'}{'genelist_mt'}=~/yes/){#this is if it had a gene list
			 	open (LIST,">$workdir/viz/$db.list\n");
			 	print LIST "#generated from AVIA_VIZ pipeline ".`date`;
			 	close LIST;
			 	_system ("cut -f4 $xfile >>$workdir/viz/$db.list\n");
			 }
		}else{
			print LOG "[ERROR] $exe  -t c -h $db -o COUNT$db $addon\n\tUSER_$pngid\t$xfile\t$type failed\n";
		}
#OLIVE		$pm->finish;#end child process
		
	}
#OLIVE	$pm->wait_all_children;
}
# die "For testing at ln".__LINE__."\n";
if (defined $opt_F || (defined $opt_f && (!$runData || $cfg{'UserInfo'}{'modulesid'}==0)) ){#if this is a viz only or subtractive with no databases or only preconfigured databases
	my $resolution=$cfg{'UserInfo'}{'resolution'};
	$resolution||="1000k";
	if (-e "genelist" && !-z "genelist"){
		$addon=" -l genelist ";
		if (exists $cfg{'UserInfo'}{'genelist_type'} && $cfg{'UserInfo'}{'genelist_type'}=~/highlight/i){
			$addon.= "-z ";
		}
	}
	if ($opt_F){
		_system("head -n1 $opt_F > subtractive.out");
		open (FDI_OUTS,"<$opt_F") or die "Cannot open this file $opt_F\n";# ("grep -e 'exonic' -e 'splic' -e 'UTR' -e intron -i $opt_F |grep -ve 'ncRNA_' -i  >> subtractive.out");
		open (NEW,">subtractive.out" ) or die "Cannot open this file subtractive.out\n";
		while (my $line=<FDI_OUTS>){
			chomp $line;
			my @arr=split("\t",$line);
			if ($line=~/Variant\sID\tANNOVAR\sannot/){
				print NEW "#Chr\tQuery Start\tQuery End\t#" . join ("\t#",@arr[3..$#arr],@arr[0..3])."\n";
			}elsif($line=~/(\bexonic\b|splic|UTR|intron)/){
				my @cord=split(":",$arr[0]);
				print NEW join("\t",@cord[0..2])."\t".join("\t",@arr[3..$#arr],@arr[0..3])."\n";
			}
		}
		close NEW; close FDI_OUTS;
	}
	_system ("perl $bin/parseGenesByBin_test.pl -I subtractive.out $addon");
	if (`grep 'includeannot' $opt_c | wc -l` >0){
		print CIRCOS "USER_pregenes\t/SeqIdx/circosdb/data/human/genes.human.hg19.txt\ttile\n" ;
		print CIRCOS "USER_presnp\t/SeqIdx/circosdb/data/human/dataByBin/$bin/hg19.snp135.$resolution.txt\tscatter\n";
	}
	print CIRCOS "USER_data\t$workdir/viz/subtractive.out.$resolution.txt\tscatter\n";#user data binned
}
close CIRCOS;
my $proj_name=`basename $workdir`;chomp $proj_name;
_system ("zip $proj_name.zip *list *png\n");
if ($opt_s){
	$addon=" -I ";
}else{
	$addon='';
}
_system ("perl $bin/../circos_utils/write_circos_conffiles.pl -f circos.input -N $addon");#die 'for testing at '.__LINE__."\n";
_system ("chmod +x $workdir/viz/circos.bat");
_system ("$workdir/viz/circos.bat");
_system ("touch Done.circos\n");
_system ("ls -al >>manifest.txt\n");
if ($opt_N && $cfg{'UserInfo'}{'email'}){
	my $dev='';
	if ($cfg{'UserInfo'}{'file'}=~/avia0/ || $cfg{'UserInfo'}{'label'}=~/\-dev/){
		$dev='dev';
	}
	if (`ls $workdir/viz/*png | wc -l` > 0){
		# if (!$dev){
			#sendMsg("Thank you for using the AVIA software" , "<h2>Your analysis <font color=\"red\">$cfg{'UserInfo'}{'label'}</font> is now complete. </h2> You can directly link to your page by clicking below or by cutting the link below and pasting into any web browser: <br />".
			#"<a href=\"http://avia$dev.abcc.ncifcrf.gov/apps/site/results/?id=$cfg{'UserInfo'}{'label'}\">http://avia$dev.abcc.ncifcrf.gov/apps/site/results/?id=$cfg{'UserInfo'}{'label'}</a> <br /><br />You can also retrieve other submissions by using our data retrieval page at : <br />".
			#"<a href=\"http://avia$dev.abcc.ncifcrf.gov/apps/site/retrieve_a_request\">http://avia$dev.abcc.ncifcrf.gov/apps/site/retrieve_a_request</a> to retrieve your results by providing your id above.  Your results will be stored for 1 week from the date of submission.","$cfg{'UserInfo'}{'email'}");
		# }
#		sendMsg("Thank you for using the AVIA software" ,"Your analysis $cfg{'UserInfo'}{'label'} is now complete.  Please visit your retrieval page at : ".
#			"http://avia$dev.abcc.ncifcrf.gov/apps/site/retrieve_a_request to see your results and input your id and analysis type.  Your results will be stored for 1 week from the date of submission","$cfg{'UserInfo'}{'email'}");
	}else{
		sendMsg("ERROR: $cfg{'Target_Info'}{'label'}", "Please check the directory :$cfg{'Target_Info'}{'label'}.  running circos failed<br />","vuonghm\@mail.nih.gov");
	}
}elsif (!$cfg{'UserInfo'}{'email'}){
	print STDERR "$opt_N:$cfg{'UserInfo'}{'email'} does not exist\n".Dumper (\%cfg);	
}

print LOG "[INFO] Ended on ". `date`;
close STDERR;
close LOG;
#_system ("cp $proj_name.zip /user
sub _system {
	my $cmd=shift;
	my $die=shift;
	print LOG "Running $cmd\n";
	eval{
		if ($cmd!~/\>/){
			system ("$cmd>>runviz.log\n");
		}else{
			system ("$cmd\n");
		}
	};
	if ($?){
		if ($die){
			die "[ERROR] Cannot run command $cmd\n";
		}else{
			warn "[WARN] Could not run $cmd\n";
			return 0;
		}
	}else{
		print LOG "[SUCCESS] Running $cmd\n";
	}
	return 1;
}

# sub sendMsg{ #old email using sendmail
# 	my $subject=shift;
# 	my $msg=shift;
# 	my $email_addy=shift;
# 	$msg.="\n\nPlease do not reply to this message. It goes to an unattended inbox.  Please use our webform at http://avia.abcc.ncifcrf.gov/apps/site/submit_a_question\n";

# 	my  %mail = ( To  =>   "$email_addy",
# 			BCC => "vuonghm\@mail.nih.gov",
#             From    => 'NCI-FrederickAVIA@mail.nih.gov',
#             Subject=> "$subject",
#             Message => "<a href=\"http://avia.abcc.ncifcrf.gov target=\"_blank\"><img src=\"http://avia.abcc.ncifcrf.gov/modules/site/library/images/avia_logo_01.png\" width=\"360\" height=\"82\" alt=\"AVIA logo\" border=\"0\" style=\"margin-bottom:50px;\"></a>\n
# $msg\n"
#            );

#   sendmail(%mail) or die $Mail::Sendmail::error;
  
#  }


sub sendMsg{#new with attachment
	my $subject=shift;
	my $text=shift;
	my $to=shift;
  	my $path_to_image = qq(/bioinfoC/hue/scripts.dir/parseData/avia_logo_01.png);
	my $message = MIME::Lite->new(
    From    => "NCI-FrederickAVIA\@mail.nih.gov",
    To      => $to,
    Bcc      => 'vuonghm@mail.nih.gov',
    Subject => $subject,
    Type    => 'multipart/related',
);
	$msg.="<h6>Please do not reply to this message. It goes to an unattended inbox.  To contact us, please use our webform at <a href=\"http://avia.abcc.ncifcrf.gov/apps/site/submit_a_question\">http://avia.abcc.ncifcrf.gov/apps/site/submit_a_question</a><br />";

# Now, we have to attach the message in HTML. First the HTML
my $html_message = qq(
	<body>
     
     <img src="cid:avia_logo.png">
     <br />
     $text
     $msg
     </body>
);

# Now define the attachment
$message->attach (
    Type => 'text/html',
    Data => $html_message,
);

# Let's not forget to attach the image too!
$message->attach (
    Type => 'image/png',
    Id   => 'avia_logo.png',
    Path => "$path_to_image",
);

$message->send
    or die qq(Message wasn't sent: $!\n);
}

sub writeConf{
	open (FILE,"<$opt_c") or die "Cannot open file";
	while (my $param=<FILE>){
		chomp $param;
		my ($key,$value)=split("=",$param);
		if ($key=~/\.(\S+)/){
			$key=$1;
		}
		$cfg{'UserInfo'}{$key}=$value;
	}
	# print Dumper (\%cfg);
}

