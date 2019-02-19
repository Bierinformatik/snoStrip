#!/usr/bin/env perl

use strict;

&usage() if scalar @ARGV != 2;
my $in = $ARGV[0];
my $out = $ARGV[1];


open(IN, "<".$in) or die "$!"."$in";
my @lines = <IN>;
close(IN);

my $row=0;
my $front=0;
my @alnstr;
my @desc;
my $length=0;

chomp(@lines);
open(OUT, ">".$out);
foreach(@lines){
    if ($_ =~ /^CLUSTAL.*$/){
	#print OUT $_."\n\n\n";
    }
    elsif ($_ =~ /^(\S+)\s+(\S+)$/){
	if ($front==0){
	    $front = length($_) - length($2);
	}
	if ($#alnstr >= $row){
	    $alnstr[$row].=$2;
	}
	else {
	    push(@alnstr, $2);
	    push(@desc, $1);
	}
	$row++;
    }
    else{
	$row = 0;
    }
}

print OUT "# STOCKHOLM 1.0\n\n\n";
for(my $i=0; $i<=$#alnstr; $i++){
    my $tmpalnstr = $alnstr[$i];
    if ($length == 0){
	$length = length($alnstr[$i]);
    }
    $desc[$i] =~ tr/ /_/;
    $alnstr[$i] =~ tr/-/./;
    print OUT $desc[$i];
    print OUT " " x ($front-(length $desc[$i]));
    print OUT $alnstr[$i]."\n";
}
my $tmpstr="#=GC SS_cons";
print OUT $tmpstr;
print OUT " " x ($front-(length $tmpstr));
print OUT "." x $length;
print OUT "\n//\n\n\n";
close(OUT);


sub usage
{

    print STDERR "USAGE:\n
         \t ./stockholm.pl inputfile outpufile
         \n";
    exit(0);

}
