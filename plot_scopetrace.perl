#!/usr/bin/perl
use strict;
use warnings;
use POSIX;
use File::Basename;
# script to locate the most recent scopetrace file
# and repeatedly filter it, plotting out trace by trace.
# use camera index, roi index & rotation angle (columns 2,3,4)

while (1) {
    plot_one();
    select(undef,undef,undef,0.25);
}

sub plot_one {
    my $st=[split("\n",`ls -t dat/scopetrace*`)];
    my $file=$st->[0];
    open(F,$file) || die;
    my $dat={};

    while (<F>) {
	chomp;
	my ($ts,$cix,$rix,$theta,$p01,$p05,$p09)=split(' ');
	if (!defined($dat->{$cix,$rix,$theta})) {
	    $dat->{$cix,$rix,$theta}={};
	    $dat->{$cix,$rix,$theta}->{"timestamp"}=[];
	    $dat->{$cix,$rix,$theta}->{"lvl_0.1"}=[];
	    $dat->{$cix,$rix,$theta}->{"lvl_0.5"}=[];
	    $dat->{$cix,$rix,$theta}->{"lvl_0.9"}=[];
	} 
	my $subset=$dat->{$cix,$rix,$theta};
	push(@{$subset->{"timestamp"}},$ts);
	push(@{$subset->{"lvl_0.1"}},$p01);
	push(@{$subset->{"lvl_0.5"}},$p05);
	push(@{$subset->{"lvl_0.9"}},$p09);
    }

    my $by_nom_theta={};
    my $kys=[sort {$a cmp $b} keys %{$dat}];
    foreach my $kix (0..$#{$kys}) {
	my ($cix,$rix,$theta)=split($;,$kys->[$kix]);
	$theta+=360 while ($theta<0); 
	my $nom_theta=(90*floor($theta/90+0.5) % 360);
	$by_nom_theta->{$nom_theta}=[] if (!defined($by_nom_theta->{$nom_theta}));
	push(@{$by_nom_theta->{$nom_theta}},$kys->[$kix]);
    }

    my $mmppx={("0"  =>+0.025,
		"90" =>+0.025,
		"180"=>-0.025,
		"270"=>-0.025)};

    my $outfile=basename($file);
    my $outfile_ps=$outfile;
    $outfile_ps =~ s/qdp/ps/;
    open(F,">",$outfile) || die;
    my @deltas=();
    foreach my $nom_theta ( 90,270,0,180 ) {
	foreach my $array_ix (0..$#{$by_nom_theta->{$nom_theta}}) {
	    my $array=$dat->{$by_nom_theta->{$nom_theta}->[$array_ix]};
	    my $last;
	    my $every=1;
	    my $st=0;
	    my $sdelta=0;
	    my $sn=0;
	    my $startpos=0;
	    for (my $i=0;$i<=$#{$array->{"timestamp"}};$i++) {
		my ($t,$pos)=($array->{"timestamp"}->[$i],$array->{"lvl_0.5"}->[$i]);
		if (($i!=0) && ($#{$array->{"timestamp"}}-$i)%$every == 0) {
		    if ($sn==0) {
			$st=$t;
			$sn=1;
			$sdelta=$pos;
		    }
		    $startpos=$sdelta/$sn if ($startpos == 0);
		    my $delta=($sdelta/$sn-$startpos)*$mmppx->{$nom_theta};
		    printf F ("%f %f\n",$st/$sn,$delta);
		    push(@deltas,$delta);
		    $st=$sdelta=$sn=0;
		    my $stale=$#{$array->{"timestamp"}}-$i;

		    $every = $stale/10;
		    $every=10 if ($every>10);
		    $every=1 if ($every<1);

		} else {
		    $st     += $t;
		    $sdelta += $pos;
		    $sn     += 1;
		}
	    }
	    printf F "no no\n";
	}
	printf F "no no\n";
    }

# printf "keys:\n%s\n",join("\n",keys %{$a});

    my $cmdlist=[];
    if (0) {
	push(@{$cmdlist},"/null");
	push(@{$cmdlist},"skip on");
	
	push(@{$cmdlist},"col off 2");
	push(@{$cmdlist},"col off 4");
	push(@{$cmdlist},"col off 6");
	push(@{$cmdlist},"col off 8");

	push(@{$cmdlist},"plot vert");
	
	push(@{$cmdlist},"win 1\nyplot 2");
	push(@{$cmdlist},"win 3\nyplot 4");
	push(@{$cmdlist},"win 5\nyplot 6");
	push(@{$cmdlist},"win 7\nyplot 8");

	push(@{$cmdlist},"col 1 on 1 2");
	push(@{$cmdlist},"lst 2 on 2");
	push(@{$cmdlist},"col 2 on 3 4");
	push(@{$cmdlist},"lst 2 on 4");
	push(@{$cmdlist},"col 3 on 5 6");
	push(@{$cmdlist},"lst 2 on 6");
	push(@{$cmdlist},"col 4 on 7 8");
	push(@{$cmdlist},"lst 2 on 8");

#	push(@{$cmdlist},"win all");
	push(@{$cmdlist},"lab x elapsed time [msec]");
	push(@{$cmdlist},"lab 10 vpos 0.15 0.85 jus lef col 1 \"\\gdX\\dCCS\\u (BF2 & BF3)\"");
	push(@{$cmdlist},"lab 11 vpos 0.15 0.65 jus lef col 2 \"\\gdX\\dCCS\\u (BF0 & BF1)\"");
	push(@{$cmdlist},"lab 12 vpos 0.15 0.45 jus lef col 3 \"\\gdY\\dCCS\\u (BF1 & BF3)\"");
	push(@{$cmdlist},"lab 13 vpos 0.15 0.25 jus lef col 4 \"\\gdY\\dCCS\\u (BF0 & BF2)\"");
	push(@{$cmdlist},
	     "lab 20 vpos 0.05 0.5 rot 90 \"deviation from starting position [mm]\"");
	push(@{$cmdlist},"font roman");
	push(@{$cmdlist},"cpd /xw");
	push(@{$cmdlist},"plot");
	push(@{$cmdlist},"q");
    } else {
	my @sorted_deltas=sort {$a<=>$b} @deltas;
	my ($llim,$ulim)=@sorted_deltas[0,$#sorted_deltas];
	$llim=($llim>-0.025)?-0.025:$llim*1.1;
	$ulim=($ulim<0.025)?0.025:$ulim*1.1;
	push(@{$cmdlist},"/null");
	push(@{$cmdlist},"skip on");

	push(@{$cmdlist},"col 1 on 1 2");
	push(@{$cmdlist},"col 2 on 3 4");
	push(@{$cmdlist},"col 3 on 5 6");
	push(@{$cmdlist},"col 4 on 7 8");

	push(@{$cmdlist},"lst 1 on 1 3");
	push(@{$cmdlist},"lst 1 on 5 7");

	push(@{$cmdlist},"lst 2 on 2 4");
	push(@{$cmdlist},"lst 2 on 6 8");
	
	push(@{$cmdlist},"font roman");
	push(@{$cmdlist},"r y $llim $ulim");
	push(@{$cmdlist},"lab 10 vpos 0.15 0.85 jus lef col 1 \"\\gdX\\dCCS\\u (BF2 & BF3)\"");
	push(@{$cmdlist},"lab 11 vpos 0.15 0.80 jus lef col 2 \"\\gdX\\dCCS\\u (BF0 & BF1)\"");
	push(@{$cmdlist},"lab 12 vpos 0.15 0.75 jus lef col 3 \"\\gdY\\dCCS\\u (BF1 & BF3)\"");
	push(@{$cmdlist},"lab 13 vpos 0.15 0.70 jus lef col 4 \"\\gdY\\dCCS\\u (BF0 & BF2)\"");
	push(@{$cmdlist},"lab x elapsed time [msec]");
	push(@{$cmdlist},"lab y deviation from starting position [mm]");
	push(@{$cmdlist},"plot over");
	
	push(@{$cmdlist},"cpd /xw");
	push(@{$cmdlist},"plot");
	push(@{$cmdlist},"cpd /cps");
	push(@{$cmdlist},"plot");
	push(@{$cmdlist},"cpd");
	push(@{$cmdlist},"\$mv pgplot.ps $outfile_ps");
	push(@{$cmdlist},"q");
    }

    open(F,"|qdp $outfile") || die;
    foreach my $cmd (@{$cmdlist}) {
	printf F "%s\n",$cmd;
    }
    close(F);
}

