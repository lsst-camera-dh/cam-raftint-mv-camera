#!/usr/bin/perl
use warnings;
use strict;
# script to parse list:
# <camera ID> <reg file> <rotation>
# and output on the command line the appropriate 
# arguments needed by QuadCameraScope, i.e.
# --roi=<camera_ID>,<xmin>:<xmax>,<ymin>:<ymax>,<theta>
my @ret;

while (<>) {
    chomp;
    my ($cix,$filename,$rot)=split(' ',$_);
    push(@ret,extract_roi_lims($cix,$filename,$rot));
}

printf "%s\n",join(' ',@ret);

sub extract_roi_lims {
    my ($cix,$fn,$rot)=@_;
    my $ret="";
    open(F,"<",$fn) || die;
    while (<F>) {
	next if !/^box/;
	$_ =~ tr/,()/\ \ \ /;
	my @list=split(' ');
	$ret .= sprintf("--rot=%s,%s:%s,%s:%s,%s",$cix,
			$list[1]-0.5*$list[3],
			$list[1]+0.5*$list[3],
			$list[2]-0.5*$list[4],
			$list[2]+0.5*$list[4],
			$rot);
    }
    close(F);
    $ret;
}
