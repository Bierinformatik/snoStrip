#!/usr/bin/env perl

use strict;
use warnings;

my ($filename,$fileoutname) = ($ARGV[0],$ARGV[1]);
open(my $fh, '<:encoding(UTF-8)', $filename)
or die "Could not open file '$filename' $!";

my %overlap_region;
my @original_regions;
my @chains;
my @chain;

my @line;
my ($chr,$start,$end,$strand);

my @blast_line;
my @ov;

#Iterate through blast file and make a hash with one value (array of
#arrays) for each chromosome/strand-pair
#Array of arrays are a list of regions with 'start end index'
while (my $row = <$fh>) {
    unless($row eq "\n"){ #Ignore empty lines
	chomp $row;
	push(@original_regions,$row);

	@line=split(/\s+/,$row);
	($chr,$start,$end)=($line[1],$line[8],$line[9]);
    
	if($start<=$end){
	    $strand="+";
	}else{
	    $start=$line[9];
	    $end=$line[8];
	    $strand="-";
	}
	
	#Last element is an index to identify the original line later
	@blast_line=($start,$end,scalar(@original_regions)-1);
	
	unless($overlap_region{$chr.$strand}){
	    @ov=();
	    push (@{$overlap_region{$chr.$strand}},[@blast_line]);
	    
	}else{
	    push (@{$overlap_region{$chr.$strand}},[@blast_line]);
	}
    }
}

#Find chains of overlapping regions
for my $keys (keys %overlap_region){
    my @rec=@{$overlap_region{$keys}};
    my ($ov_start,$ov_end); 
    
    if(scalar(@rec)>1){
	@rec = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @rec;
	$ov_start=@{$rec[0]}[0];
	$ov_end=@{$rec[0]}[1];
	push(@chain,@{$rec[0]}[2]);

	for (my $i=1; $i <scalar(@rec); $i++){
	    my @t=@{$rec[$i]};
	    
	    if($t[0]<=$ov_end){
		if($t[1]>$ov_end){
		    $ov_end=$t[1];
		}
		push(@chain,$t[2]);
	    }else{
		push(@chains,[@chain]);
		@chain=();
		push(@chain,$t[2]);
		$ov_start=$t[0];
		$ov_end=$t[1];
	    }
	}
	push(@chains,[@chain]);
	@chain=();	
    }else{
	push(@chain,@{$rec[0]}[2]);
	push(@chains,[@chain]);
	@chain=();
    }
}

#Find best scoring line per chain and write them to $fileoutname
open ( OUT, ">>$fileoutname" ) or die $!;
for (my $i=0; $i <scalar(@chains); $i++){
    my @t=@{$chains[$i]};

    my $maxScore=0;
    my $maxID=-1;
    
    for (my $j=0; $j <scalar(@t); $j++){
	
	@line=split(/\s+/,$original_regions[$t[$j]]);
	if($line[11]>$maxScore){
	    $maxScore=$line[11];
	    $maxID=$t[$j];
	}
    }
    print OUT $original_regions[$maxID],"\n";
}
close(OUT);
    

##########################################################################
sub overlap
#Checks if two intervals are overlapping, and if yes, if they second
#interval can extend the first one
# 0 - no overlap
# !0 - array with (possibly) extended positions
    
##########################################################################
{
    my ($s1,$e1,$s2,$e2)=@_;
    my ($sR,$eR);
    my @ret;
    
    if(($s2>=$s1 && $s2<=$e1) || ($e2>=$s1 && $e2<=$e1)){
	
	if($s2<$s1){
	    $sR=$s2;
	}else{
	    $sR=$s1;
	}

	if($e2>$e1){
	    $eR=$e2;
	}else{
	    $eR=$e2;
	}
	@ret=($sR,$eR);
	return @ret;
	
    }else{
    	return 0;
    }
}
    
