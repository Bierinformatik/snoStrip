#!/usr/bin/env perl

use Cwd 'abs_path';


BEGIN{
    my $path = abs_path($0);
    $path =~ s/scripts.*/packages\//;
    push @INC, $path;
}

use strict;
use Getopt::Long;
use CONFIG;
use find_alignment_position;

my ( $opt_h, $opt_f, $opt_o, $opt_k ) = ( "", "", "", "" );
my ( $family, $targets, $prime_targets, $organism );
my ( @file );


my $usage = << "END";

  usage:     perl $0 -f|file STRING -o|output STRING
    
  purpose:   This script should remap the alignment positions
             that can be found in the given modifications  
             file to the alignments located in the alignmnet
             path found in the CONFIG.pm .
             A modification file should look like that:
             CD_22    18S-200,25S-2673    18S-298

  options:
             -h|help      Prints this help message.
                          [OPTIONAL]

             -f|file      Input modification file.
                          [REQUIRED]

             -o|output    Output modification file.
                          [REQUIRED]

             -k|kingdom   Target RNAs whose kingdom shall be mapped.
                          Either 'fun' or 'deu'.
                          [REQUIRED]

END


GetOptions(
    "h|help"      =>  \$opt_h,
    "f|file=s"    =>  \$opt_f,
    "o|output=s"  =>  \$opt_o,
    "k|kingdom=s" =>  \$opt_k,
);

if ( $opt_h ){
    print STDERR $usage, "\n";
    exit(0);  
} 

if ( !$opt_f && !$opt_h && !$opt_o && !$opt_k ){
    print STDERR $usage, "\nERROR: seems like you forgot to specifiy a mandatory argument.\n";
    exit(0);
}
if ( !$opt_f && !$opt_o && !$opt_k ){
    print STDERR $usage, "\n";
    exit(0);
}

##################################################################################################################
## START


## set the correct reference organism
if( $opt_k eq "fun" ){ $organism = "S_cerevisiae"; }
elsif( $opt_k eq "pro" ){$organism = "D_melanogaster"; }
elsif( $opt_k eq "deu" ){$organism = "H_sapiens"; }
else{
    print STDERR $usage, "\nERROR: The kingdom wasn't correctly specified!\n";
   exit(0); 
}

## open the modification file
open FILE, $opt_f or die $!;
@file = <FILE>;
close FILE;


open OUT, ">" . $opt_o  or die $!;
## go through prediction file and replace the 
## old alignment position with the new one
#while ( <FILE> ){ 
foreach ( @file ){

    chomp $_;

    ( $family, $prime_targets, $targets ) = split( /\s+/, $_ );

    ## map targets if this family contains targetsites at their 3' end
    if ( $targets ne "." ){
	$targets = &map_targets( $targets );
    }

    ## map targets it this family contains prime targetsites
    if ( $prime_targets ne "." ){
	$prime_targets = &map_targets( $prime_targets );
    }

    print OUT $family, "\t", $prime_targets, "\t", $targets, "\n";
    

}

close OUT;



##################################################################################################################
sub map_targets
## this method split the list of target sites and maps the yeast positions
## to the alignment
##################################################################################################################
{

    my @targets = split ( /,/, shift) ;
    my ( $position, $target_RNA, $return_list );

    foreach my $target ( @targets ){
	
	( $target_RNA, $position ) = split( /-/, $target );
	
	$return_list .= &map_position( $organism . "_" . $target_RNA, $target_RNA, $position, $organism, $opt_k ) . ",";

    }

    chop $return_list;

    return $return_list;

}
