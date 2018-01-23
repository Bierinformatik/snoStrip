package searchForTargetSites;


use strict;
use vars qw(@ISA @EXPORT);
use Bio::AlignIO;
use Bio::SimpleAlign;
use Logger;

require Exporter;


@ISA = qw(Exporter);
@EXPORT = qw(searchTargetSites);


######################################################################################################################
sub searchTargetSites
## This method tries to identify putative target sites in a new snoRNA candidate.
## Therefore, a PWM of length 9 is created starting from every position in the snoRNA alignment.
## Each PWM is then used to score a sliding window in the new sequence.
## If the score lies above a certain threshold, the sequence is assumed to contain a conserved target region.
######################################################################################################################
{

    my ( $aln_ref, $seq, $max_start_a, $min_start_b, $max_start_prime, $locations ) = @_;
    my $LENGTH = 13;
    my $BOX = 0;
    my %sequences;
    my $header;
    my $new_sequence;
    my $aln = $$aln_ref;

    my $max_start_prime_cutted = $max_start_prime;

    
    ## before the alignment will be cleaned we will created PWMs for box D and D'
    ## create PWM from the alignment
    my @pwm = @{createPWM( $aln )};

    ## create slice for box D
    my ( @Dbox, @Dprimebox );
    for ( my $start = 0; $start < 4; $start++ ){
	push @Dbox, [ @{ $pwm[$start] } [ $min_start_b..$min_start_b+3 ] ];
    }

    for ( my $start = 0; $start < 4; $start++ ){
	push @Dprimebox, [ @{ $pwm[$start] } [ $max_start_prime..$max_start_prime+3 ] ];
    }


    ## remove all columns containing more than 60% gaps
    my $aln2 = cleanAlignment( $aln_ref, $max_start_a, $min_start_b, \$max_start_prime_cutted );


    ## create PWM from the alignment
    @pwm = @{createPWM( $aln2 )};

    my @slice_d;
    my $i = $aln2->length() - ( $LENGTH + 4 );
    
    ## create a two dimensional slice of the pwm
    ## this slice will then be used to score a sliding window on the new snoRNA
    for ( my $start = 0; $start < 4; $start++ ){
	push @slice_d, [ @{ $pwm[$start] } [ $i..$i+($LENGTH+$BOX-1) ] ];
    }
    
    my @slice_prime;

    if( $max_start_prime > 0 ){

	my $i = 0;
	if ( $max_start_prime_cutted >= $LENGTH ){ $i = $max_start_prime_cutted - $LENGTH }
	else{ print STDERR "ERROR in searchTargetSite\nNot enough nucleotides in front of box D prime!\t$max_start_prime_cutted\n"; exit;}
	
	## create a two dimensional slice of the pwm
	## this slice will then be used to score a sliding window on the new snoRNA
	for ( my $start = 0; $start < 4; $start++ ){
	    push @slice_prime, [ @{ $pwm[$start] } [ $i..$i+($LENGTH+$BOX-1) ] ];
	}
	
    }


    my $box_ref = \@slice_d;
    my $prime_ref = \@slice_prime;


    my ( $D_start, $Dprime_start, $D_target, $Dprime_target, $D_score, $Dprime_score, $D_best, $Dprime_best, $D_seq, $Dprime_seq );
    my @line; 


    ## score the region in fron of both boxes.
    ( $D_best, $Dprime_best ) = ( 0, 0 );
    ( $D_seq, $Dprime_seq ) = ( "", "" );
    ( $D_start, $Dprime_start ) = ( 0, 0 );

    ## extract region in front of box D
    for my $i ( -1..4 ){
	if ( $D_start + $i <= length( $seq ) ){
	    $D_target = substr( $seq, $D_start - $LENGTH - $i, $LENGTH + $BOX);
	    $D_score = scorePattern( $D_target, $box_ref );
	    $D_best = $D_score if $D_score >= $D_best;
	    $D_seq = $D_target if $D_score >= $D_best;
	}
    }
    

    if ( $D_best < 0.7 ){
	( $D_best, $D_seq ) = scoreWholeSeq( $seq, $box_ref, $LENGTH, $BOX );
    }
    

    ## extract region in front of box Dprime
    if ( $Dprime_start > 0 ){
	
	for my $i ( -1..4 ){
	    $Dprime_target = substr( $seq, $Dprime_start - $LENGTH - $i, $LENGTH + $BOX);
	    if ( length ( $Dprime_target ) == $LENGTH + $BOX && $prime_ref->[0][0] ){
		$Dprime_score = scorePattern( $Dprime_target, $prime_ref );
		$Dprime_best = $Dprime_score if $Dprime_score >= $Dprime_best;
		$Dprime_seq = $Dprime_target if $Dprime_score >= $Dprime_best;
	    }
	}
	
	## check if the prime box score is greater than 0.6
	## otherwise, compute a scan over the whole sequence
	## to ensure to get the best score,
	## since sometimes the box annotation is wrong
	if( $Dprime_best < 0.7 ){
	    
	    ( $Dprime_best, $Dprime_seq ) = scoreWholeSeq( $seq, $prime_ref, $LENGTH, $BOX );
	    
	}
	
    }	
    else{
	
	( $Dprime_best, $Dprime_seq ) = scoreWholeSeq( $seq, $prime_ref, $LENGTH, $BOX );
	
    }
    
    
    if ( $D_best >= 0.7 or $Dprime_best >= 0.7 ){

	return ( "1", $D_best, $D_seq, $D_best, $Dprime_best ) if $D_best > $Dprime_best;
	return ( "1", $Dprime_best, $Dprime_seq, $D_best, $Dprime_best ) if $D_best <= $Dprime_best;

    }
    else{

	return ( "0", $D_best."\t".$Dprime_best, $D_seq, $D_best, $Dprime_best );

    }

}


##########################################################################################################################################
sub createPWM
#this method calculates the logarithm to base 2
##########################################################################################################################################
{

    my $aln = shift;

    my ( @pwm, @sequence );

    my $alignment_length = ( $aln->length() ) - 1;
    
    ## initialize array with '0'
    for my $j (0..3){
	for my $i (0..$alignment_length){
	    $pwm[$j][$i] = 0;
	}
    }


    my @sequences =  $aln->each_seq() ;

    foreach my $seq ( @sequences ){
	
	@sequence =  split( //, $seq->seq() );

	#calculate the amount of each nucleotide at each position of the box
	for(my $i = 0; $i < scalar @sequence ; $i++){
	    if( $sequence[$i] eq "A" ){ $pwm[0][$i]++; next; }
	    if( $sequence[$i] eq "G" ){ $pwm[1][$i]++; next; }
	    if( $sequence[$i] eq "T" || $sequence[$i] eq "U" ){ $pwm[2][$i]++; next; }
	    if( $sequence[$i] eq "C" ){ $pwm[3][$i]++; next; }
	}
	
    }


    ## calculate log-odd ratios out of absolute numbers
    ## use a pseudo-count of 0.0001 to avoid logarithms of zero
    ## a background frequency of 0.25 for each nucleotide is used (multiplying by 4)
    ## would it make a difference if we use a other background frequency? --> no
    for my $j (0..3){
	for my $i (0..$alignment_length){
	    $pwm[$j][$i] = 0.0001 if $pwm[$j][$i] == 0;
#	    $pwm[$j][$i] = log2(($pwm[$j][$i]*4) / scalar @sequences );
	    $pwm[$j][$i] = $pwm[$j][$i] / scalar @sequences;
	}
    }


    return \@pwm;

    

}


##########################################################################################################################################
sub log2
#this method calculates the logarithm to base 2
##########################################################################################################################################
{
    my $arg = shift;
    return log($arg)/log(2);
}


##########################################################################################################################################
sub cleanAlignment
#this method calculates the logarithm to base 2
##########################################################################################################################################
{

    my ( $aln, $max_start_a, $min_start_b, $start_prime_ref ) = @_;

    my ( @remove, @sequences, @sequence );

    my $start_prime_original = ${$start_prime_ref};

    ## note columns at the beginning that can be removed
    for my $i (0..$max_start_a){
	push @remove, $i;
    }
    
    

    @sequences = $$aln->each_seq();
    
    ## go through alignment and search for columns with more than 60% gap characters
    for my $i ( $max_start_a+6..$min_start_b-1 ){

	my $gaps = 0;

	## count all gaps
	foreach my $seq ( @sequences ){
	    if ( substr( $seq->seq(), $i, 1) eq "-" ){ $gaps++ }
	}

	## check whether this column contains too mamy gaps
	if( $gaps / scalar @sequences > 0.6 ){ push @remove, $i }

    }

    ## note all columns at the end of the alignment
    ## after box D
    for my $i ( $min_start_b+4..$$aln->length()-1 ){
	push @remove, $i;
    }

#    print join ("\t", @remove), "\n";

    ## create new alignment
    my $new_aln = new Bio::SimpleAlign();
    my $new_seq;

    ## adjust the $start_prime_b position in the alignment
    ## therefore go through the remove array
    ## and decrease the $start_prime_b position by one
    ## if the remove-position is smaller as the original $start_prime_b
    if ( $start_prime_original > 0 ){
	foreach my $int ( reverse @remove ){
	    if ( $int <= $start_prime_original ){
		$start_prime_original--;
		${$start_prime_ref}--;
#		print STDERR "remove: $int\t\tnew start: $start_prime_original\n";
		
	    }
	}
    }

    ## delete the 'gapped' columns from each sequence
    foreach my $seq ( @sequences ){
	@sequence = split( //, $seq->seq() );
	foreach my $index ( reverse @remove ){
	    splice( @sequence, $index, 1 );
	}

	## add shortened seq to the new alignment
	$new_seq = Bio::LocatableSeq->new( -seq => join( "", @sequence ),
	                          -id  => $seq->id(), 
	                         );

	## add new sequence to the novel alignment object
	$new_aln->add_seq( $new_seq );

    }


    return $new_aln;

}


##########################################################################################################################################
sub scorePattern
#this method scores a given pattern with a given matrix and returns the achieved score or probability
##########################################################################################################################################
{
    my $sequence = shift;
    my $boxRef = shift;
    my $score;
    my @seq = split(//, $sequence);


    if($seq[0] eq "A"){$score = $boxRef->[0][0]}
    if($seq[0] eq "G"){$score = $boxRef->[1][0]}
    if($seq[0] eq "T" || $seq[0] eq "U"){$score = $boxRef->[2][0]}
    if($seq[0] eq "C"){$score = $boxRef->[3][0]}

    
    for (my $i = 1; $i < scalar @seq; $i++){
	if($seq[$i] eq "A"){$score += $boxRef->[0][$i]}
	if($seq[$i] eq "G"){$score += $boxRef->[1][$i]}
	if($seq[$i] eq "T" || $seq[$i] eq "U" ){$score += $boxRef->[2][$i]}
	if($seq[$i] eq "C"){$score += $boxRef->[3][$i]}
    }


    ## normalize the score
    $score /= scalar @seq;


    return $score;
}


##########################################################################################################################################
sub scoreWholeSeq
#this method scores a given pattern with a given matrix and returns the achieved score or probability
##########################################################################################################################################
{

    my ( $seq, $target, $best_target, $score, $ref, $best_score, $LENGTH, $BOX );
    $seq = shift;
    $ref = shift;
    $LENGTH = shift;
    $BOX = shift;
    $best_score = 0;


    for my $i ( 0..length( $seq ) - $LENGTH + $BOX ){
	$target = substr( $seq, $i, $LENGTH + $BOX );
	if ( length ( $target ) == $LENGTH + $BOX && $ref->[0][0] ){
	    $score = scorePattern( $target, $ref );
	    $best_target = $target if $score >= $best_score;
	    $best_score = $score if $score >= $best_score;
	}
    }

    $best_score = "0.0000000" if $best_score == 0;

    return ( $best_score, $best_target );

}


1;
