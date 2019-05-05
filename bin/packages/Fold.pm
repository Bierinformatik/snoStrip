package Fold;

use strict;
use CONFIG;
use vars qw(@ISA @EXPORT);

require Exporter;

@ISA = qw( Exporter );
@EXPORT = qw( getFoldingCD getFoldingHACA getIndependentFoldingHACA get_subopt_folding);

####################################################################################################
## VARIABLES
####################################################################################################

my $RNASUBOPT = $CONFIG::RNASUBOPT;
my $RNAEVAL = $CONFIG::RNAEVAL;


####################################################################################################
sub getFoldingCD
#tries to identify the correct folding of a given sequence with respect to its snoRNA-type
####################################################################################################
{


    #declare some variables and paramaters
    my ($seq, $boxA, $boxB,$tmpPath) = @_;
    
    my ($folding, $mfe, $constraint, $rnafold, $struc, $energy);
    my (@tmp, @line);
    my $up_paired_region = $boxA - 1;
    my $down_paired_region = (length $seq) - $boxB + 1;
    my $last_paired_nt;
    my $p=150;
    
    ## at the beginning one has to create a constraint which forces RNAfold to fold the sequence into such a structure
    ## example ((((((...............................................))))))
    ## parameters for RNAfold -C are: '<','>' for bases which have to pair, and 'x' for bases which do not pair at all
    $constraint = "<" x $up_paired_region . "x" x ($boxB - $boxA) . ">" x $down_paired_region;

    #make a first folding
    #$rnafold = qx(echo -e "$seq\n$constraint" | $RNASUBOPT -C -e 3 -s);
    #@tmp = split(/\n/, $rnafold);

    @tmp=get_subopt_folding($seq,$constraint,$p,$tmpPath);

    
    for(my $i = 1; $i < scalar @tmp; $i++){
	($struc, $mfe) = split(/\s+/, $tmp[$i]);

	if($struc =~ /^\.{0,$up_paired_region}\(+\.{0,1}\(*\.+\)*\.{0,1}\)+\.{0,$down_paired_region}$/){

	    ## analyze how many unpaired nucleotides lie between the last paired one and the C-box
	    $last_paired_nt =  &get_last_paired_nt_position( $struc, $up_paired_region );

	    next if $up_paired_region - $last_paired_nt > 5;

	    ## remove all foldings with less than 3 basepairs
	    next if ( $struc =~ tr/\(/\(/ ) <= 2;

	    $folding = $struc;
	    $energy = $mfe;
	    last;
	}
    }

    if(!$folding){
	($folding, $energy) = split(/\s+/, $tmp[1]);
	$energy = "atypical structure";
    }

    return ($folding, $energy);
    
}



####################################################################################################
sub getIndependentFoldingHACA
#tries to calculate the correct folding for a given sequence of HACA-snoRNAs
####################################################################################################
{
    ##
    ## VARIABLES
    ##
    my ( $seq, $cons, $start_box1, $start_box2 , $tmpPath) = @_;
    my ( $hairpin1, $hairpin2 );
    my ( $folding_hairpin1, $mfe_hairpin1 );
    my ( $folding_hairpin2, $mfe_hairpin2 );
    my ( @RNAfold_hairpin1, @RNAfold_hairpin2 );
    my ( $folding, $mfe );
    my $p=100;

    ##
    ## split sequence into two distinct regions (according to its box motifs)
    ##
    $hairpin1 = substr( $seq, 0, $start_box1 -1 );
    $hairpin2 = substr( $seq, $start_box1 + 5, $start_box2 - $start_box1 - 6 );

    ##
    ## use RNAsubopt to produce all suboptimal foldings of each hairpin within a given energy-range
    ##
    #@RNAfold_hairpin1 = split( /\n/, `echo -e '$hairpin1' | $RNASUBOPT -e 2 -s` );
    #@RNAfold_hairpin2 = split( /\n/, `echo -e '$hairpin2' | $RNASUBOPT -e 2 -s` );
    @RNAfold_hairpin1=get_subopt_folding($seq,0,$p,$tmpPath);
    @RNAfold_hairpin2=get_subopt_folding($seq,0,$p,$tmpPath);
    shift( @RNAfold_hairpin1 );
    shift( @RNAfold_hairpin2 );

    ##
    ## search for hairpins in the RNAsubopt output
    ##
    foreach ( @RNAfold_hairpin1 ){
	if( $_ =~ /^[\.\(]+\(\.{3,}\)[\.\)]+\s+/ ){
	    ( $folding_hairpin1, $mfe_hairpin1 ) = split( /\s+/, $_ );
	    last;
	}
    }
    foreach ( @RNAfold_hairpin2 ){
	if( $_ =~ /^[\.\(]+\(\.{3,}\)[\.\)]+\s+/ ){
	    ( $folding_hairpin2, $mfe_hairpin2 ) = split( /\s+/, $_ );
	    last;
	}
    }

    ##
    ## combine both individual structures
    ##
    if( $folding_hairpin1 && $folding_hairpin2 ){
	$folding = $folding_hairpin1 . "......" . $folding_hairpin2 . "......";
	$mfe = $mfe_hairpin1 . " " . $mfe_hairpin2;
    }


    return ( $folding, $mfe, $mfe_hairpin1, $mfe_hairpin2, length($hairpin1), length($hairpin2) );

}



####################################################################################################
sub getFoldingHACA
#tries to calculate the correct folding for a given sequence of HACA-snoRNAs
####################################################################################################
{
 
    ##
    ## VARIABLES
    ##
    my ($seq, $cons,$tmpPath) = @_;
    my ($fold, $mfe, $struc, $tmp, $folding, $energy);
    my @tmp;
    my ($open, $close, $constmp, $startHairpinB, $hairpins);

    #added by JAN for changed RNAsubopt method (-p)
    my $p=50;

    @tmp=get_subopt_folding($seq,$cons,$p,$tmpPath);

#    print "search for folding...\n";

    for(my $i = 1; $i < scalar @tmp; $i++){
	($struc, $mfe) = split(/\s/, $tmp[$i]);

	#to test if a folding consists of two hairpins
	#($tmp = $struc) =~ s/\.*[^\)]+\.+[^\(]+\.+[^\)]+\.+[^\(]+\.*//;

	($tmp = $struc) =~ s/\.*\([\.\(]+\)[\.\)]+\([\.\(]+\)[\.\)]+//;
	if($tmp eq ""){
	
	    #print $tmp[$i],"\n";
    
	    #test if both hairpins are independent of each other
	    #remove first hairpin and count openning and closing parentheses in second hairpin
	    #if both amounts are equal -> independent
	    ($tmp = $struc) =~ s/^\.*\([\.\(]+\)[\.\)]+//;
	    #print "hairpin 2: $tmp\n";
	    $startHairpinB = (length $cons) - (length $tmp);
	    $open = ($tmp =~ s/\(/\(/g);
	    $close = ($tmp =~ s/\)/\)/g);
	    


	    #maybe one should test if the first hairpin is longer than it should be (does not be longer than till first box)
	    #therefore remove last hairpin and check if folding is longer than constraint-line which is cut off before boxA 
	    ($tmp = $struc) =~ s/\.*\([\.\(]+\)[\.\)]+$//;
	    ($constmp = $cons) =~ s/x.*//;
	    #print "hairpin 1: $tmp\n\n";
	    #print "open: $open,\tclose: $close\n";

	    #print "\nLaengen:\n hairpin1: ", length $tmp, "\nmax hairpin1: ", length $constmp, "\nconstraint: ",length $cons,"\nstartHP2: ", $startHairpinB, "\n\n";

	    if($open == $close && length $tmp <= length $constmp && (length $constmp) + 6 <= $startHairpinB ){ 
		#print "\n\nfolding gefunden!\n";
		#print $tmp , "\n", $constmp ,"\n";
		$folding = $struc;
		$energy = $mfe;
		last; 
	    }
	}
    }
    
    #no characteristical folding was found
    if(!$folding){
	#print "FOLDING: no suitable folding was found!\n\n";
	($folding, $energy) = split(/\s/, $tmp[1]);
	$energy = "atypical structure";
	
	#find out whats wrong with this folding
	#how many hairpins?
	$hairpins = ($folding =~ s/(\([\.\(]+\)[\.\)]+\([\.\(]+\)[\.\)]+)/$1/g);
	
	#first hairpin ends before boxA
	($tmp = $folding) =~ s/\.*\([\.\(]+\)[\.\)]+$//;
	($constmp = $cons) =~ s/x.*//;
	
	#second hairpin starts before boxA
	$startHairpinB = (length $cons) - (length $tmp);


	#print "hairpins: $hairpins\nconstmp: ", length $constmp, "\ntmp: ", length $tmp, "\nhairpinB:  $startHairpinB\n\n";

	if($hairpins > 2){$energy .= " ($hairpins hairpins)";}
	elsif(length $constmp < length $tmp){$energy .= " (hairpin 1 too long)";}
	elsif($startHairpinB <( length $constmp) + 6){$energy .= " (hairpin 2 too long)";}
	else{$energy .= " (strange)";}
    }

    return ($folding, $energy);

}



####################################################################################################
sub get_last_paired_nt_position
## retrieves the position of the last paired nucleotide upstream of box C
####################################################################################################
{

    my $folding = shift;
    my $box_start = shift;
    my $pos = 0;

    my @folding = split( //, $folding );

    for ( my $i = 0; $i < $box_start ; $i++){
	$pos = $i if $folding[$i] =~ /\(/;
    }

    return $pos;

}

####################################################################################################
sub get_subopt_folding
## retrieves the position of the last paired nucleotide upstream of box C
####################################################################################################
{
    my($seq,$cons, $p,$tmpPath) = @_;

    my ($fold, $mfe, $struc);
    my @tmp;
    my ($open, $close, $constmp, $startHairpinB, $hairpins);

    #added by JAN for changed RNAsubopt method (-p)
    my $rnaeval;
    my @tmp2;
    my $subopt_file="$tmpPath/rnasubopt_out-tmp-".int(rand(1000));
    my @sorted;

    #$fold = `echo -e '$seq\n$cons' | $RNASUBOPT -C -e 1 -s`;   #fold with constraint #OLD METHOD from Sebastiaon
    if($cons){
	$fold = `echo -e '$seq\n$cons' | $RNASUBOPT -C -p $p`;   #fold with constraint   #NEW METHOD by JAN
    }else{
	$fold = `echo -e '$seq' | $RNASUBOPT -p $p`;   #fold without constraint   #NEW METHOD by JAN
    }

    @tmp = split(/\n/, $fold);

    open( FOLD_OUT, '>', $subopt_file) or print $!;
    for(my $i = 1; $i < scalar @tmp; $i++){
	print FOLD_OUT "$tmp[0]\n$tmp[$i]\n";
    }
    close FOLD_OUT;
    
    $rnaeval = `cat  '$subopt_file' | $RNAEVAL`;   #fold with constraint
    unlink $subopt_file;
    
    @tmp = split(/\n/, $rnaeval);
    my $line;
    for(my $i = 1; $i < scalar @tmp; $i+=2){
	$tmp[$i]=~/^([.()]*)\s+\(\s*([-]?\d+[.]\d\d)\)$/;
	($struc, $mfe) = ($1,$2);
	push(@tmp2, "$struc $mfe\n");
    }

    @sorted = sort { (split(/\s+/, $a))[1] <=> (split(/\s+/, $b))[1] } @tmp2;
    @tmp=();
    push(@tmp,"$seq\n");
    push(@tmp,@sorted);

    return @tmp;
}
    
1;
