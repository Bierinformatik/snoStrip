package primeBoxes;


use strict;
use vars qw(@ISA @EXPORT);
use Bio::AlignIO;
use Bio::SimpleAlign;
use boxPositions;
use makeMuscleAlignment;
use sequences;

require Exporter;

 
@ISA = qw( Exporter );
@EXPORT = qw( primeBoxSearch primeBoxSearch_pairwise );



##########################################################################################################################################
sub createBoxProfiles
#creates profiles of boxes for a given set of fastafiles (they need to provide the new header-format)
##########################################################################################################################################
{

    my $Ref = shift;

    #boxes is a hash containing two keys (boxes), for C/D type it is C-box => array-ref of PWM and D-box => array-ref of PWM
    #and for H/ACA type it is H-box => array-ref of PWM and ACA-box => array-ref of PWM
    my %boxes;

    my ( $boxA, $boxB );
    my ( @box1, @box2 );
    my ( @A, @B );

    my $runs = scalar keys %{$Ref};
    my ( $length1, $length2 );

    #length describes the boxlength which has to be the same for box1 and box2 in the whole multifasta file
    #which means that it is not allowed that a c-box consists of 6nt in one sequence and of 7nt in another sequence in one multifasta file
    my @keys = keys %{$Ref};
    

    ###JAN Changed because it can become bad infinity-loop
    ###    Function of the while and $k is unclear???
    ###OLD:
    my $k = 0;

    #print "createBoxProfiles ",scalar(@keys),"\n"; #jan
    
    while( $k<scalar(@keys)){
    	$boxA = $$Ref{$keys[$k]}->{box1Prime};
    	$boxB = $$Ref{$keys[$k]}->{box2Prime};
    	last if $boxA;
    	$k++;
    }
    
    unless(defined($boxA)){
	print STDERR "Problem with prime boxes\n";
	exit(0);
    }
    
    ( $length1, $length2 ) = ( length($boxA) - 1, length($boxB) - 1 );
  


    ## initialize arrays with '0'
    for my $j (0..3){
	for my $i (0..$length1){
	    $box1[$j][$i] = 0;
	}
	for my $i (0..$length2){
	    $box2[$j][$i] = 0;
	}
    }

        
    ## filter boxes from fasta-header
    ## and fill in the absolute numbers into both PWMs (box1- and box2-PWM)
    foreach ( @keys ){

	if( $$Ref{$_}->{box1Prime} ){
	    ( $boxA, $boxB ) = ( $$Ref{$_}->{box1Prime}, $$Ref{$_}->{box2Prime} );
	    
	    ## check if the current box motif does contain N-character
	    if( $boxA !~ /N/ ){
	    
		@A = split(//, $boxA);
		
		#calculate the amount of each nucleotide at each position of the box
		for(my $i = 0; $i < scalar @A; $i++){
		    if($A[$i] eq "A"){$box1[0][$i]++; next;}
		    if($A[$i] eq "G"){$box1[1][$i]++; next;}
		    if($A[$i] eq "T"){$box1[2][$i]++; next;}
		    if($A[$i] eq "C"){$box1[3][$i]++; next;}
		}
	    }

	    ## check if the current box motif does contain N-character
	    if ( $boxB !~ /N/ ){

		@B = split(//, $boxB);
		for(my $i = 0; $i < scalar @B; $i++){
	    if($B[$i] eq "A"){$box2[0][$i]++; next;}
	    if($B[$i] eq "G"){$box2[1][$i]++; next;}
	    if($B[$i] eq "T"){$box2[2][$i]++; next;}
	    if($B[$i] eq "C"){$box2[3][$i]++; next;}
		}
	    }
    
	}

    }


    ## calculate log-odd ratios out of absolute numbers
    ## use a pseudo-count of 0.0001 to avoid logarithms of zero
    ## a background frequency of 0.25 for each nucleotide is used (multiplying by 4)
    ## would it make a difference if we use a other background frequency? --> no
    for my $j (0..3){
	for my $i (0..$length1){
	    $box1[$j][$i] = 0.0001 if $box1[$j][$i] == 0;
	    $box1[$j][$i] = log2(($box1[$j][$i]*4) / $runs);
	    
	}
   
	for my $i (0..$length2){
	    $box2[$j][$i] = 0.0001 if $box2[$j][$i] == 0;
	    $box2[$j][$i] = log2(($box2[$j][$i]*4) / $runs);
	}
    }

    ## put everything in a return-hash
    %boxes = (
	'box1' => \@box1,
	'box2' => \@box2,
	);

    return \%boxes;

}


##########################################################################################################################################
sub primeBoxSearch_pairwise
## this method searches for possible box-motifs while using a pairwise muscle alignment of
## the new snoRNA-candidate with the best-scoring query sequence
##########################################################################################################################################
{

    my ( $fasta, $Ref, $res, $offset, $verbose, $query ) = @_;

    my %tmp_hash;
    $tmp_hash{$query} = $Ref->{$query};

    my $seqID = substr( $res->{bla}, 0, 30 );
    my $typeA = substr( $res->{type}, 0, 1 );
    my $typeB = substr( $res->{type}, 1 );
    $seqID =~ s/[\(\),]/_/g;

    ( my $aln_file = $fasta ) =~ s/fa$/aln/;
    ( my $dnd_file = $fasta ) =~ s/fa$/dnd/;

    my ( $bool, $message ) = &makeAlignment( $fasta, $aln_file );
    if ( $bool == 0 ){
	print STDERR "\nSome ERROR occured during alignment computation!\n\n";
	return ( "", "", "", "", 0, "alignment computation error" );
    }

    #create a new alignment object using bio::AlignIO
    my $str = Bio::AlignIO->new(-file => "$aln_file", -format => "clustalw");
    my $aln = $str->next_aln();

    my ($start_min_B, $start_max_B, $mean_L, $length_box_B );
    my ($start_min_A, $start_max_A, $length_box_A );
    my ( $boxA, $boxB, $penaltyA, $penaltyB, $seq ); 
    ( $start_min_B, $start_max_B, $mean_L, $length_box_B ) = getBoxStart( \%tmp_hash , \$aln, "prime2" );
    ( $start_min_A, $start_max_A, $mean_L, $length_box_A ) = getBoxStart( \%tmp_hash , \$aln, "prime1" );

#    print STDERR "boxA: $start_min_A\tboxB: $start_min_B\n";

    foreach my $seq_obj ( $aln->each_seq_with_id( $seqID ) ){

	$seq = $seq_obj->seq();

	## in case the query sequence contains the C'box 'NNNNNNN'
	## change the boxA to 'NNNNNNN' as well and update the start position accordingly
	if ( $tmp_hash{$query}->{box1PrimeStart} != 0 ){
	    $boxA = &extractSubsequence($seq, $start_min_A, $length_box_A);
	    $penaltyA  = isCorrectBox($boxA, $typeA);
	}
	else{
	    $boxA = "";
	    $penaltyA = 2;
	}
	
	$boxB = &extractSubsequence( $seq, $start_min_B, $length_box_B );
	$penaltyB = isCorrectBox($boxB, $typeB);


    }

    my $boxRef = createBoxProfiles( $Ref );
    if ( $penaltyB <= 1 ){

	## get additional information of box A
	my $scoreA = scorePattern( $boxA, $boxRef->{'box1'} ) if $boxA;
	my $startA = startPos($seq, $start_min_A) + 1 if $boxA;

	## get additional information of box B
	my $scoreB = scorePattern( $boxB, $boxRef->{'box2'} );
	my $startB = startPos($seq, $start_max_B) + 1 ;

	if( $penaltyA > 1 ){
	    return ("NNNNNNN", 0, $boxB, $startB, 2, $scoreA, $scoreB, $penaltyA, $penaltyB );
	}

	return ( $boxA, $startA, $boxB, $startB, 1, $scoreA, $scoreB, $penaltyA, $penaltyB );
    }

    else {

	return ( $boxA, 0, $boxB, 0, 0, 0, 0, $penaltyA, $penaltyB );

    }

    

}


##########################################################################################################################################
sub primeBoxSearch
#this method searches for possible box-motifs while using clustalw alignments of this new snoRNA-candidate with all known snoRNAs of this
#family
#the box-corresponding region in this clustalw alignment is cut out and the box suggested by this alignment is analysed with respect to
#our box motif definitions
##########################################################################################################################################
{

    my ( $fasta, $Ref, $res, $offset ) = @_;
    ( my $aln_file = $fasta ) =~ s/fa$/aln/;
 
    my $seqID = substr( $res->{bla}, 0, 30 );
    my $typeA = substr( $res->{type}, 0, 1 );
    my $typeB = substr( $res->{type}, 1 );
    my ( $minStartA, $maxStartA, $minStartB, $maxStartB ) = ( 1000, 0, 1000, 0 );
    my ( $boxA, $boxB );
    my ( $boxLengthA, $boxLengthB );
    my $tmpSeqLength = 0;
    my ( $seq, $tmpSeq );
    my $meanL = 0;
    my $maxL = 0;
    my $minL = 1000;


    ## create muscle alignment
    ## -clwstrict --> clustalw format, including CLUSTALW header
    ## -quiet -> suppress output to STDERR
    my ( $bool, $message ) = &makeAlignment( $fasta, $aln_file );
    if ( $bool == 0 ){
	print STDERR "\nSome ERROR occured during alignment computation!\n\n";
	return( "", 0, "", 0, 0 );
    }
    

    ## create a new alignment object using bio::AlignIO
    my $str = Bio::AlignIO->new(-file => "$aln_file", -format => "clustalw");
    my $aln = $str->next_aln();

    my ( $extendRangeA, $extendRangeB ) = ( 0, 0 );
    my ( $maxExtendRangeA, $maxExtendRangeB ) = ( 0, 0 );

	my ($start_min_B, $start_max_B, $mean_L, $length_box_B );
    ( $start_min_B, $start_max_B, $mean_L, $length_box_B ) = getBoxStart( $Ref, \$aln, "prime2" );


    ## find the alignment position of both box start positions in the muscle alignment
    ## and create a window out of these positions extracted from all query sequences
    HASH: foreach ( keys %{$Ref} ){

	if( $$Ref{$_}->{box1Prime} ){

	    my $seq_name = $_;
	    $seq_name =~ s/[\(\),]/_/g;

	    ## retrieve sequence from the alignment
	    $tmpSeq  = ($aln->each_seq_with_id($seq_name))[0]->seq();

        ## create a gapless sequence in order to extract putative box motifs later on
	    ( $seq = $tmpSeq ) =~ s/-//g;

	    $tmpSeqLength = length $seq;
	    $meanL += $tmpSeqLength;
	    $maxL = $tmpSeqLength if $tmpSeqLength > $maxL;
	    $minL = $tmpSeqLength if $tmpSeqLength < $minL;


	    ## start to analyze box D' since this is the box were really interested in
	    ## sometimes box C' is not really consereved and hence we'll annotate such boxes with NNNNNN
	    ## box2 or boxB is normally box D'
	    ## box1 or boxA is normally box C'
	    my $startB = $$Ref{$_}->{box2PrimeStart} - 1;
	    my $B = $$Ref{$_}->{box2Prime};
	    $boxLengthB = length $$Ref{$_}->{box2Prime};
	    
	    ## retrieve the actual alignment position of start position of box B
	    my $posB = $aln->column_from_residue_number( "$seq_name", $startB );

	    ## check if alignment position is correct in order to avoid errors like "----GTGA" were the first - is choosen instead of G
	    my $tmpB = substr( $tmpSeq, $posB, $boxLengthB );
	    if( $tmpB ne $B ){
		if( $tmpB =~ /^-/ ){
		    my $correctPos = 0;
		    
		    ## increase the positionA counter while walking over the gaps towards the real box start
		    while ( !$correctPos ){
			$posB++;
			$tmpB = substr( $tmpSeq, $posB, $boxLengthB );
			if ( !(defined $tmpB ) ) {
			    print STDERR "tmpSeq: $tmpSeq\tboxB: $boxB\tposB: $posB\n";
			    print STDERR $$Ref{$_}->{box2PrimeStart} - 1, "\n";
			    print STDERR $seq_name, "\n";			    
			    exit(0);
			}
			if( $tmpB =~ /^\w/ ){ $correctPos = 1; }
		    }

		}
		## check if extracted seq contains dashes in order to enlarge the scoring-window
		#elsif( $tmpB =~ /-/ ){
		#    my $letters = length( $tmpB =~ s/-//g ) + 1;
		#    $extendRangeB = $posB + $boxLengthB - 1;
		#    while( $letters != $boxLengthB ){
		#	$extendRangeB++;
		#	if( substr( $tmpSeq, $extendRangeB, 1 ) ne "-" ){
		#	    $letters++;
		#	}
		#    }
		#    $extendRangeB = $extendRangeB - $posB - $boxLengthB + 1;
		#}
	    }

	    if( $extendRangeB > $maxExtendRangeB ){ $maxExtendRangeB = $extendRangeB; }
	    if( $posB > $maxStartB ){ $maxStartB = $posB; }
	    if( $posB < $minStartB ){ $minStartB = $posB; }


	    ## continue with the C' box 
	    ## but first of all check whether there is a prime box anyway
	    ## otherwise go on with the next sequence
	    if ( $$Ref{$_}->{box1Prime} =~ /N/ ){
		next HASH;
	    }

	    my $startA = $$Ref{$_}->{box1PrimeStart} -1;
	    my $A = $$Ref{$_}->{box1Prime};
	    $boxLengthA = length $$Ref{$_}->{box1Prime};

	    ## retrieve the actual alignment position of start position of box A
	    my $posA = $aln->column_from_residue_number( "$seq_name", $startA );
	    
	    ## check if alignment position is correct in order to avoid errors like "----GTGA" were the first - is choosen instead of G
	    my $tmpA = substr( $tmpSeq, $posA, $boxLengthA );
#	    print STDERR "tmpA: $tmpA\n";
#	    print STDERR "tmpSeq: $tmpSeq\n";

	    if( $tmpA ne $A ){
		if( $tmpA =~ /^-/ ){
		    my $correctPos = 0;
		    
		    ## increase the positionA counter while walking over the gaps towards the real box start
		    while ( !$correctPos ){
#			print STDERR "correctPos: $correctPos\n";
			$posA++;
			$tmpA = substr( $tmpSeq, $posA, $boxLengthA );
#			print STDERR "nwe tmpA: $tmpA\n";
			if( $tmpA =~ /^\w/ ){ $correctPos = 1 }
		    }
		}
		## check if extracted seq contains gaps in order to enlarge the scoring-window
		elsif( $tmpA =~ /-/ ){
		    my $letters = length( $tmpA =~ s/-//g ) + 1;
		    $extendRangeA = $posA + $boxLengthA - 1;
		    while( $letters != $boxLengthA ){
			$extendRangeA++;
			if( substr( $tmpSeq, $extendRangeA, 1 ) ne "-" ){
			    $letters++;
			}
		    }
		    $extendRangeA = $extendRangeA - $posA - $boxLengthA + 1;
		}
	    }
	    

	    ## enlarge the window boundaries if necessary
	    if( $posA > $maxStartA ){ $maxStartA=$posA; }
	    if( $posA < $minStartA ){ $minStartA=$posA; }

	}
    }

    $meanL /= ( scalar keys %{$Ref} );

#    print STDERR "prime boxen: \n";
#    print STDERR "minA: ", $minStartA, "\tmaxA: ", $maxStartA, "\tminB: ", $minStartB, "\tmaxB: ", $maxStartB, "\n";
#    print STDERR "extendB: $extendRangeB, maxExtendRangeB: $maxExtendRangeB\n";

    ## look if all boxes in our alignment start at the same position ($minStart == $maxStart)
    ## is they do so we can go on and cut this block from the sequence
	$seqID =~ s/[\(\),]/_/g;

    my ( $returnA, $returnStartA, $returnB, $returnStartB );
    my ( $tmpreturnA, $tmpreturnStartA, $tmppenaltyA, $tmpscoreA );
    my ( $tmpreturnB, $tmpreturnStartB, $tmppenaltyB, $tmpscoreB );
    my ( $penaltyA, $penaltyB ) = ( 10, 10 );
    my ( $scoreA, $scoreB ) = ( 0, 0 );

    foreach my $seq_obj ( $aln->each_seq_with_id( $seqID ) ){

	my $seq = $seq_obj->seq();
	(my $tmpSeq = $seq) =~ s/-//g;
	my $s;
	my $boxRef = createBoxProfiles( $Ref );

	## check if there is at least one C' box
	## that means check if the start positions have changed since their initialization
	if( $minStartA != 1000 && $maxStartA != 0 ){
	    if($minStartA == $maxStartA){
		
		$boxA = &extractSubsequence($seq, $minStartA, $boxLengthA);
		$penaltyA = isCorrectBox($boxA, $typeA );
		if( $penaltyA <= 1 ){
		    $returnA = $boxA;
		    $returnStartA = startPos($seq, $minStartA);
		    $scoreA = scorePattern( $returnA, $boxRef->{'box1'});
		}
		else{
		    ($returnA, $returnStartA, $penaltyA, $scoreA) = explorer($seq, $boxRef->{'box1'}, $typeA, $maxStartA, $minStartA);
		}
	    }
	    
	    elsif( $maxExtendRangeA != 0 && $minStartA == $maxStartA ){
		$maxStartA = $minStartA + $maxExtendRangeA;
		$minStartA += 3;
		$maxStartA -= 3;
		($returnA, $returnStartA, $penaltyA, $scoreA) = explorer($seq, $boxRef->{'box1'}, $typeA, $maxStartA, $minStartA);
	    }
	    
	    else{
		#cut out a larger sequence and score all possible box-sequences
		($returnA, $returnStartA, $penaltyA, $scoreA) = explorer($seq, $boxRef->{'box1'}, $typeA, $maxStartA, $minStartA);	
	    }
	}


	#do the same comparison for second box
	if( $minStartB == $maxStartB && $maxExtendRangeB == 0 ){

#	    print STDERR "start_ ist gleich\n";

	    $boxB = &extractSubsequence($seq, $minStartB, $boxLengthB);
	    $penaltyB = isCorrectBox( $boxB, $typeB );
	    if( $penaltyB <= 1 ){
		$returnB = $boxB;
		$returnStartB = startPos($seq, $minStartB);
		$scoreB = scorePattern( $returnB, $boxRef->{'box2'});
	    }
	    else{
		($returnB, $returnStartB, $penaltyB, $scoreB) = explorer($seq, $boxRef->{'box2'}, $typeB, $maxStartB, $minStartB);
	    }
	}

	elsif( $maxExtendRangeB != 0 && $minStartB == $maxStartB ){
	    $maxStartB = $minStartB + $maxExtendRangeB;
#	    print "extende Range um $maxExtendRangeB!\n";
	    $minStartB += 3;
	    $maxStartB -= 3;
	    ($returnB, $returnStartB, $penaltyB, $scoreB) = explorer($seq, $boxRef->{'box2'}, $typeB, $maxStartB, $minStartB);
	}

	else{
#	    print "there are some position probs at boxB\n";
	    ($returnB, $returnStartB, $penaltyB, $scoreB) = explorer($seq, $boxRef->{'box2'}, $typeB, $maxStartB, $minStartB);
	}
	
    }
    $returnStartB++;
    $returnStartA++;

    ## rules for accepting a sequence:
    ## reject this sequence if 1) one or both boxes provides more than one mutation
    ##                         2) if it is an outlier and the summarized score is larger than 1.5
    ##                         3) or  
#    print "INFOS:\n";
#    print "A: $scoreA\nB: $scoreB\nPA: $penaltyA\nPB: $penaltyB\n";
#    print "boxA: $returnA\nboxB: $returnB\n\n";
#    print "startA: $returnStartA\nstartB: $returnStartB\n";
    $returnStartB += $offset;
    $returnStartA += $offset;

    ## in case only box D' is suitable, return this box and rewrite box C' as 'NNNNNNN'
    if( $penaltyB <= 1 && $penaltyA > 1 ) {
	return ( "NNNNNNN", "0", $returnB, $returnStartB, 1, $scoreA, $scoreB, $penaltyA, $penaltyB );
    }

    if( $penaltyA <= 1 && $penaltyB <= 1 ){
	if( ( $scoreA >= -20 && $scoreB >= -20) || ( ( $penaltyA + $penaltyB<=1.5 ) && $scoreA >= -20 && $scoreB >= -20 ) ){
	    return ( $returnA, $returnStartA, $returnB, $returnStartB, 1, $scoreA, $scoreB, $penaltyA, $penaltyB );
	}
	elsif( $scoreB >= -20 && $penaltyB <= 1 ){
	    return ( "NNNNNNN", "0", $returnB, $returnStartB, 1, $scoreA, $scoreB, $penaltyA, $penaltyB );
	}
	else{
	    return ( $returnA, $returnStartA, $returnB, $returnStartB, 0, $scoreA, $scoreB, $penaltyA, $penaltyB );
	}
    }

    else{
	return ( $returnA, $returnStartA, $returnB, $returnStartB, 0, $scoreA, $scoreB, $penaltyA, $penaltyB );
    }
}




##########################################################################################################################################
sub explorer
#
##########################################################################################################################################
{

    my ($seq, $boxRef, $type, $maxStart, $minStart) = @_;

#    print STDERR "\n\nexplorer: Seq: $seq\n";

#    print STDERR "max $maxStart min $minStart\tscalar", scalar(@{${$boxRef}[0]}), "\n";

    my $range = $maxStart - $minStart + 6 + scalar(@{${$boxRef}[0]});
    my $sequence = substr($seq, $minStart -3, $range);
    my $substitute = 0;
    my $overall_start = $minStart -3;
    my $seq_start;
    my $returnStartBox = 0;

#    print STDERR "EXPLORER:\nmaxStart: $maxStart\tminStart $minStart\tBoxlength:", scalar(@{${$boxRef}[0]}), "\t";
#    print STDERR "range: $range\nseq: $sequence\n";

    #check if sequence contains '-'
    while ( $sequence =~ /-/ ) {
#	print STDERR "seq enthaelt '-'\n$sequence\n";
	if ( $sequence =~ /^-/ ) {
	    $substitute = ($sequence =~ s/-//g);
	    $sequence = substr($seq, $minStart - 3 - $substitute, $range + $substitute);
	    $overall_start = $minStart - 3 - $substitute;
#	    print STDERR "novel seq: $sequence";
	    my $tmpSub = ($sequence =~ s/-/-/g);
	    while ( $tmpSub > $substitute ) {
		$substitute++;
		$sequence = substr($seq, $minStart - 3 - $substitute, $range + $substitute);
		$overall_start = $minStart - 3 - $substitute;
#		print STDERR "tmpSub: $tmpSub\t sub: $substitute\n$sequence\n";
		$tmpSub = ($sequence =~ s/-/-/g);
	    }
#	    print STDERR "FINAL SEQ: $sequence\neliminating $substitute dashes\n\n";
	    while ( $substitute > 0 ) {
		$sequence =~ s/^([^-]*)-/$1/;
#		print STDERR "Seq zwischendrin: $sequence\n";
		$substitute--;
	    }
	}
	else {
	    $substitute = ($sequence =~ s/-//g);
#	    print STDERR "substitute: $substitute\n";
	    $sequence = substr($seq, $minStart - 3, $range + $substitute);
	    $overall_start = $minStart - 3;
	    my $tmpSub = ($sequence =~ s/-/-/g);
	    while ( $tmpSub > $substitute ) {
		$substitute++;
		$sequence = substr($seq, $minStart - 3, $range + $substitute);
		$overall_start = $minStart - 3;
#		print STDERR "tmpSub: $tmpSub\t sub: $substitute\n$sequence\n";
		$tmpSub = ($sequence =~ s/-/-/g);
	    } 
#	    print STDERR "Seq zwischendrin: $sequence\n";
	    while ( $substitute > 0 ) {
		$sequence =~ s/-([^-]*)$/$1/;
#		print STDERR "Seq zwischendrin: $sequence\n";
		$substitute--;
	    }
	}
	
#	print STDERR "Gefixte seq: $sequence\n";
	
    }
       
    my ($returnBox, $startInTmpSeq, $penaltyBox, $scoreBox) = searchBoxes($sequence, $boxRef, $type);

#    if ( $startInTmpSeq ){

#    print STDERR "overall_start: $overall_start\n";
	$seq_start = substr($seq, 0, $overall_start );
#    print STDERR "seq_start: $seq_start\n";
	$seq_start =~ s/-//g;
#    print STDERR "seq_start: $seq_start\n";
	$returnStartBox = length( $seq_start ) + $startInTmpSeq;
#    print STDERR "start_pos: $returnStartBox\n";

#    }

#    print STDERR "startT: $startInTmpSeq\n";
#    $startInTmpSeq += $minStart -3;
#    print STDERR "startT: $startInTmpSeq\n";
    
#    my $cutSeq = substr($seq, 0, $startInTmpSeq);
#    print STDERR "cutSeq: $cutSeq\n";
#    print STDERR "laenge: ", length $cutSeq, "\n";
#    $cutSeq =~ s/-//g;
#    print STDERR "cutSeq: $cutSeq\n";
#    print STDERR "laenge: ", length $cutSeq, "\n";
#    my $returnStartBox = length $cutSeq;
	    
#    print STDERR "returnBox: $returnBox\nreturnStartBox: $returnStartBox\npenaltyBox: $penaltyBox\n scoreBox: $scoreBox\n\n";
    return($returnBox, $returnStartBox, $penaltyBox, $scoreBox);
    
}



##########################################################################################################################################
sub searchBoxes
#this method searches for possible box-motifs in a given snoRNA-candidate using the profiles calculated previously
##########################################################################################################################################
{

    my ($sequence, $boxRef, $type) = @_;

    my $laenge = length $sequence;

#    print "laenge: ",length $sequence, "\n";
#    print "sequen: ", $sequence, "\n";

    #variables declaration
    #use a sliding window to score all possible box-motifs and store the current leader
    #my @seq = split(//, $sequence);
    my ($curScore, $bestScore) = (-100,-100);
    my ($curBox, $bestBox) = ("","");
    my $startBox;
    my ($start, $end) = (0, (length($sequence))-scalar(@{${$boxRef}[0]}));
    my ($result, $tmpSeq);
    my ($curPenalty, $penalty) = (2,2);
#    print STDERR "start: $start\tend: $end\n";    

    for my $i ($start..$end){
	
	if($sequence =~ /^.{$i,$i}-/){next;}
	
	#calculate current score for this certain window
	($tmpSeq = $sequence) =~ s/.{$i,$i}//;
	$tmpSeq =~ s/-//g;
	$curBox = substr($tmpSeq, 0, scalar(@{${$boxRef}[0]}));
	$curScore = scorePattern($curBox, $boxRef);
	
	$curPenalty = isCorrectBox($curBox, $type);
#	printf STDERR "%-2d: box: %-15s: % .2f\t%.2f\n", $i, $curBox, $curScore, $curPenalty;
	
	#compare current score with the best score achieved so far
	#look at both starting positions and compare them in order to decline a box when box starting points a too far away from each other
	if($curScore>=$bestScore && $curPenalty <= $penalty ){ #&& $start > $rangeRef->[0] && $start < $rangeRef->[1]){
	    #print "drin!\n";
	    $bestScore = $curScore;
	    $bestBox = $curBox;
	    $startBox = $i;
	    $penalty = $curPenalty
	}
    }
    
    
    #test if the best box suites our constraints, which means only one mutation per box is allowed
    
#    print STDERR "\n\n", $bestBox , "\t", $bestScore ,"\tstarting a tosition: ", $startBox ,"\n\n";
    
    return ($bestBox, $startBox, $penalty, $bestScore);
    #welche rueckgabe werte sollte man nehmen? array fuer die boxen? boxen als string und noch die positionen als skalare?
    #weitere ideen?
    
}



##########################################################################################################################################
sub isCorrectBox
#this method analyses the current box whether it fits our constraints or not
#we penalize different mutations in different ways
##########################################################################################################################################
{
    my $box = shift;
    my $type = shift;
    my $penalty = 0;

    #print substr("vergleiche box : $box....................",0,30);


    #mutations of one of the first three nucleotides are more seldom than on one of the last three nucleotides
    #resulting in penalties of 0.5 and 1
    #every box with more than one mutation is penalized with 2
    #boxes with mutated G or A nucleotides in position 2 and 3 will get a penalty of 2
    if($type eq "C"){
	if((length $box) == 7){$box =~ s/^\w//;}
	if($box =~ /TGATG./){$penalty  = 0;}
	elsif($box =~ /TGA[CGA]GA|TGAT[CTA]A|TGATG[CTG]|[GCA]GATGA|T[TCA]ATGA|TG[CGT]TGA/){$penalty = 0;}
	elsif($box =~ /..ATGA|.G.TGA|.GA.TGA|.GAT.A|.GATG.|T..TGA|T.A.GA|T.AT.A|T.ATG.|TG..GA|TG.T.A|TG.TG.|TGA..A|TGA.G.|TGAT../){$penalty = 0}

	else{$penalty = 2;}
    }


    #boxes with one mutation will get penalty of 1
    #more than one mutation will result in a penalty of 2
    #boxes with mutated G or A nucleotides will get a penalty of 2
    elsif($type eq "D"){
	if($box =~ /CTGA/){$penalty = 0;}
	elsif($box =~ /CTG.|CT.A|C.GA|.TGA/){$penalty = 0}
	elsif($box =~ /CT..|C.G.|C..A|.TG.|.T.A|..GA/){$penalty = 0}
	elsif($box =~ /...A|..G.|.T..|C.../){$penalty = 2;}
	else{$penalty = 2;}
    }

    if($box =~ /-/){$penalty = 2}

    #print "penalty: $penalty\n";
    return $penalty;
}



##########################################################################################################################################
sub log2
#this method calculates the logarithm to base 2
##########################################################################################################################################
{
    my $arg = shift;
    return log($arg)/log(2);
}


1;
