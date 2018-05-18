#!/usr/bin/perl
use strict;
use warnings;
# script to make a new subdirectory directory (and new directory, if necessary) to save a set of camera images.
use Getopt::Long;
use DateTime;

my $datasetdir="IG_acq_lib";
my $subdir="test1";
my $imdatadir="/home/ccs/mvIMPACT_acquire/dat";
my $timestampsubdir;
my $stacked;
my $mosaic;
my $usage=join("\n",
	       "usage::save_camerashot.perl --datasetdir=<datasetdir>",
	       "                            --subdir=<subdir>",
	       "                            --subdir=<subdir>",
	       "                            --timestamp (appends subdir name with timestamp)",
	       "",
	       "Will create the directory <datasetdir>/<subdir> and save 4 most recent",
	       "images into that directory. If --timestamp is specified, the subdir name",
	       "will be appended with a UT timestamp (YYMMDDhhmmss).",
	       "Will terminate if the destination directory exists.");

GetOptions("datasetdir=s" => \$datasetdir,
	   "subdir=s"     => \$subdir,
	   "stacked"      => \$stacked,
	   "mosaic"       => \$mosaic,
	   "timestamp"    => \$timestampsubdir) 
    or die (join("\n","Error in command line arguments.",$usage,"\n"));

my $ts;
if (defined($timestampsubdir)) {
    printf "here we go\n";
    my $dt=DateTime->now;
    $ts=sprintf("%02d%02d%02d%02d%02d%02d",
		$dt->year-2000,  $dt->month,  $dt->day,
		$dt->hour,       $dt->minute, $dt->second);
    printf "timestamp=$ts\n";
    $subdir=$subdir."_".$ts;
}

if (! -d $datasetdir) {
    printf STDERR "target dataset directory $datasetdir doesn't exist. creating..\n";
    mkdir $datasetdir;
} else {
    printf STDERR "target dataset directory $datasetdir exists. using it.\n";
}

$subdir=$datasetdir."/".$subdir;
if (! -d $subdir) {
    printf STDERR "target subdirectory $subdir doesn't exist. creating..\n\n";
    mkdir $subdir ;
} else {
    printf STDERR "target subdirectory $subdir exists. will *NOT* use it because this will potentially overwrite previously acquired files.\n\n";
    exit(1);
}

# find the most recent raw images from each camera to save them in the $subdir

printf STDERR "\ncopying recent files into $subdir.\n\n";

my $files;
my $recent;

foreach my $cix (0..3) {
    $files=`ls -t ${imdatadir}/BF${cix}_r*fits`;
    $recent=(split("\n",$files))[0];
    printf "most recent raw $cix file is: %s\n",$recent;
    `cp $recent $subdir`;
    
    if (defined($stacked)) {
	$files=`ls -t ${imdatadir}/BF${cix}_st*fits`;
	$recent=(split("\n",$files))[0];
	printf "most recent stacked $cix file is: %s\n",$recent;
	`cp $recent $subdir`;
    }
}

if (defined($mosaic)) {
    $files=`ls -t ${imdatadir}/BF_mo*fits`;
    $recent=(split("\n",$files))[0];
    printf "most recent mosaic file is: %s\n",$recent;
    `cp $recent $subdir`;
}

# and gzip them, saves lots of disk space.

printf STDERR "compressing the copied images..\n";

`gzip $subdir/*`;

printf STDERR "finished!\n";

# done
