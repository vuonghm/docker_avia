#!/usr/bin/perl
use strict;
# This script will try to find all of the dependencies of each script and set the env var so it's ok to run 
my $env_file="set-vars.sh";
open (ENVLOG,">log.env") or die "Cannot open Log for checkEnv.pl script\n";
#standard ones that are run
my @exe_arr=(
	"cp",
	"mv",
	"perl",
	"basename",
	"dirname",
	"qsub",
	"qstat",
	"grep",
	"head",
	"cat",
	"date",
	"wc",
	"pwd",
	"mysql",
	"ls",
	"echo",
	"ps",
	"file",
	"find",
	"diff",
	"cut",
	"uname",
	);
if (!`which locate`){
	print ENVLOG "ERR: Cannot run this script\n";
	exit;
}
if (!`which which`){
	print ENVLOG "ERR: Cannot run this script\n";
	exit
}
my $env='';
foreach my $cmd (@exe_arr){
	my $path=`which $cmd`;chomp $path;
	print ENVLOG "[Info] Checking PATH $cmd=$path\n";
	if (!$path){
		my $newpath=`locate $cmd 2>/dev/null| /bin/grep $cmd\$ |head -n1`;chomp $newpath;
		print ENVLOG "[Info] Looking for new path?(locate $cmd 2>/dev/null| /bin/grep $cmd\$ |head -n1)?\n";
		if ($newpath){
			chomp $newpath;
			my $currpath='';
			if ($newpath=~/(.*)\/$cmd/){
				print ENVLOG "\t[Info] Adding $1\n";
				$env.="$1:" if ($env!~/$1/);
			}
		}else{
			print ENVLOG "[ERR] $cmd will not run on this node";
		}
	}else{
		# print "Found $path\n";
	}
}
if ($env){
	open (ENV,">$env_file") or die "Cannot open $env_file\n";
	print ENVLOG "Writing env file ($env_file)\n";
	print ENV "export PATH=$env\${PATH}\n";
	print ENV "exec \"\$\@\"\n";
	close ENV;
	chmod oct("0755"),$env_file;
}
close ENVLOG;
