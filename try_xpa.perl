#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use DateTime;
use POSIX ":sys_wait_h";
use Term::ReadKey;

# script to use XPA to repeatedly a series of fits files
my $orig_dataset="MTR2-install_b13_180516/setup_180516223343";
my $match;

my $usage=join("\n",
	       "usage::try_xpa.perl --orig_dataset=<orig_dataset>",
	       "                    --match=<match>",
	       "",
	       "Will create a ds9 window in which the orig_dataset and currently incoming datasets",
	       "will be overlaid for alignment purposes.");

GetOptions("orig_dataset=s" => \$orig_dataset,
	   "match=s"        => \$match) 
    or die (join("\n","Error in command line arguments.",$usage,"\n"));

die ($usage) if (!defined($match));

my $bkg=`ls -t ${orig_dataset}/*${match}*`;
my @bkg=split("\n",$bkg);
printf "bkg = $bkg[0]\n";
my $fg=`ls -t dat/*${match}*fits`;
my @fg=split("\n",$fg);

printf "fg = $fg[0]\n";
printf "bkg = $bkg[0]\n";

my $cmd;
my $pid=fork();

if ($pid) {
    # this is the parent
    # use xpa to repeatedly send images to the child
    sleep 3;
#    open(F,"|-","bash") || die;
#    select F; $|=1;
    # setup initial RGB display (already done by child)
    my @flist=@fg;
    while (1) {
	my @cmd=();
	# commands to fill green & blue channels, leaving red alone.
	my ($green,$blue)=@flist[0,1];

	push(@cmd,sprintf("xpaset ds9 rgb green < /dev/null"));
	push(@cmd,sprintf("xpaset ds9 fits < %s",$green));
	push(@cmd,sprintf("xpaset ds9 rgb blue < /dev/null"));
	push(@cmd,sprintf("xpaset ds9 fits < %s",$blue));

	$cmd=join("\n",@cmd);
#	printf F "%s\n",$cmd;
	`$cmd`;
	printf STDERR "assigning green=$green and blue=$blue (press any key to pause)\n";

	my $last_frame=$flist[0];
	do {
	    $fg=`ls -t dat/*${match}*fits`;
	    @flist=split("\n",$fg);
	    process_key();
	    select(undef,undef,undef,0.1) if ($flist[0] eq $last_frame);
	    exit if (waitpid(-1,WNOHANG)>0);
	} while ($flist[0] eq $last_frame);
    }
} else {
    # this is the child
    $cmd=sprintf("ds9 -rgb -red %s -green %s -blue %s",$bkg[0],$fg[0],$fg[0]);
    `ds9 -rgb -red $bkg[0] -green $fg[0] -blue $fg[0]`;
}

sub process_key {
    my $key;
    ReadMode 4; # listen for keyboard events
    if (defined($key=ReadKey(-1))) {
	# somebody punched a key
	printf STDERR "got key $key\nPress any key to continue..\n";
	do {
	    select(undef,undef,undef,0.1);
	} while (!defined($key=ReadKey(-1)));
	printf STDERR "got key $key so continuing..\n";
    }
    ReadMode 0; # revert to normal keyboard input mode.
}
