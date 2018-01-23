package boxing;

use strict;
use vars qw(@ISA @EXPORT);
use Bio::AlignIO;
use Bio::SimpleAlign;
use boxPositions;
use makeMuscleAlignment;
use sequences;

require Exporter;

 
@ISA = qw( Exporter );
@EXPORT = qw( scorePattern boxSearch badScore boxSearch_pairwise );


##########################################################################################################################################
sub createBoxProfiles
#creates profiles of boxes for a given set of fastafiles (they need to provide the new header-format)
##########################################################################################################################################
{

    my $Ref = shift;

    #boxes is a hash containing two keys (boxes), for C/D type it is C-box => array-ref of PWM and D-box => array-ref of PWM
    #and for H/ACA type it is H-box => array-ref of PWM and ACA-box => array-ref of PWM
    my %boxes;

    my ($boxA, $boxB);
    my (@box1, @box2);
    my ($startA, $startB);
    my ($minStartA, $maxStartA, $minStartB, $maxStartB) = (100,0,1000,0);
    my (@A, @B, @tmp);

    my $runs = scalar keys %{$Ref};
    my ($length1, $length2);

    #length describes the boxlength which has to be the same for box1 and box2 in the whole multifasta file
    #which means that it is not allowed that a c-box consists of 6nt in one sequence and of 7nt in another sequence in one multifasta file
    my @keys = keys %{$Ref};
    $boxA = $$Ref{$keys[0]}->{box1};
    $boxB = $$Ref{$keys[0]}->{box2};


    ($length1, $length2) = (length($boxA) - 1, length($boxB) - 1);
  

    #initialize arrays with '0'
    for my $j (0..3){
	for my $i (0..$length1){
	    $box1[$j][$i] = 0;
	}
	for my $i (0..$length2){
	    $box2[$j][$i] = 0;
	}
    }

        
    #filter boxes from fasta-header
    #and fill in the absolute numbers into both PWMs (box1- and box2-PWM)
    foreach (@keys){
	($boxA, $startA, $boxB, $startB) = ($$Ref{$_}->{box1},$$Ref{$_}->{box1Start},$$Ref{$_}->{box2},$$Ref{$_}->{box2Start});

	#print STDERR "array:", $boxA, "\t", $startA , "\t", $boxB,"\t", $startB, "\n";

	if($boxA =~ /N/ || $boxB =~ /N/){next;}

	@A = split(//, $boxA);
	@B = split(//, $boxB);


	#calculate the amount of each nucleotide at each position of the box
	for(my $i = 0; $i < scalar @A; $i++){
	    if($A[$i] eq "A"){$box1[0][$i]++; next;}
	    if($A[$i] eq "G"){$box1[1][$i]++; next;}
	    if($A[$i] eq "T"){$box1[2][$i]++; next;}
	    if($A[$i] eq "C"){$box1[3][$i]++; next;}
	}

	for(my $i = 0; $i < scalar @B; $i++){
	    if($B[$i] eq "A"){$box2[0][$i]++; next;}
	    if($B[$i] eq "G"){$box2[1][$i]++; next;}
	    if($B[$i] eq "T"){$box2[2][$i]++; next;}
	    if($B[$i] eq "C"){$box2[3][$i]++; next;}
	}
	
	#compare current starting points with previous minima and maxima
	if($startA > $maxStartA){$maxStartA = $startA}
	if($startA < $minStartA){$minStartA = $startA}
	if($startB > $maxStartB){$maxStartB = $startB}
	if($startB < $minStartB){$minStartB = $startB}
	#}
    }


    #calculate log-odd ratios out of absolute numbers
    #use a pseudo-count of 0.0001 to avoid logarithms of zero
    #a background frequency of 0.25 for each nucleotide is used (multiplying by 4)
    #would it make a difference if we use a other background frequency? --> no
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

    # print STDERR "\n\npwm:";

    for my $j (0..3){
	for my $i (0..$length2){
	    #print STDERR $box2[$j][$i] , "\t";
	}
	#print STDERR "\n";
    }
    #print STDERR "\n\n";


    #extend the range for a possible box by changing the upper and lower bounds
    my @rangeA = ($minStartA -= 6, $maxStartA += 6);
    my @rangeB = ($minStartB -= 6, $maxStartB += 6);


    #put everything in a hash
    %boxes = (
	'box1' => \@box1,
	'range1' => \@rangeA,
	'box2' => \@box2,
	'range2' => \@rangeB,
	);

    return \%boxes;

}


##########################################################################################################################################
sub boxSearch_pairwise
## this method searches for possible box-motifs while using a pairwise muscle alignment of
## the new snoRNA-candidate with the best-scoring query sequence
##########################################################################################################################################
{

    my ( $fasta, $Ref, $res, $verbose, $query ) = @_;

    my %tmp_hash;
    $tmp_hash{$query} = $Ref->{$query};

    my $seqID = substr($res->{bla},0,30);
    my $typeA = substr($res->{type},0,1);
    my $typeB = substr($res->{type},1);
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


    ## get the start position of both boxes
    my ( $minStartA, $maxStartA , $meanL, $boxLengthA, $seq );
    my ( $boxA, $boxB, $penaltyA, $penaltyB );
    my ( $minStartB, $maxStartB , $boxLengthB );
    ( $minStartA, $maxStartA, $meanL, $boxLengthA ) = getBoxStart( \%tmp_hash, \$aln, "1" );
    ( $minStartB, $maxStartB, $meanL, $boxLengthB ) = getBoxStart( \%tmp_hash, \$aln, "2" );

    foreach my $seq_obj ( $aln->each_seq_with_id( $seqID ) ){

	$seq = $seq_obj->seq();

	$boxA = &extractSubsequence($seq, $minStartA, $boxLengthA);
	$penaltyA  = isCorrectBox($boxA, $typeA);
	
	$boxB = &extractSubsequence( $seq, $minStartB, $boxLengthB );
	$penaltyB = isCorrectBox($boxB, $typeB);

    }

    my $boxRef = createBoxProfiles($Ref);

    if ( $penaltyA <= 1 && $penaltyB <= 1 ){

	## get additional information of box A
	my $scoreA = scorePattern( $boxA, $boxRef->{'box1'} );
	my $startA = startPos($seq, $minStartA) + 1;

	## get additional information of box B
	my $scoreB = scorePattern( $boxB, $boxRef->{'box2'} );
	my $startB = startPos($seq, $minStartB) + 1;

	if($typeA eq "C" && !&kinkTurn( $boxA, $boxB )){
	    return ($boxA, $startA, $boxB, $startB, 0, 0, "No K-turn found", $scoreA, $scoreB, $penaltyA, $penaltyB );
	}

	return ( $boxA, $startA, $boxB, $startB, 1, 0, "", $scoreA, $scoreB, $penaltyA, $penaltyB );
    }

    else {

	return ( $boxA, 0, $boxB, 0, 0, 0, "No box motif found", 0, 0, $penaltyA, $penaltyB );

    }

}


##########################################################################################################################################
sub boxSearch
#this method searches for possible box-motifs while using clustalw alignments of this new snoRNA-candidate with all known snoRNAs of this
#family
#the box-corresponding region in this clustalw alignment is cut out and the box suggested by this alignment is analysed with respect to
#our box motif definitions
##########################################################################################################################################
{

    my ($fasta, $Ref, $res, $muscle, $verbose) = @_;

    ( my $aln_file = $fasta ) =~ s/fa$/aln/;
    ( my $dnd_file = $fasta ) =~ s/fa$/dnd/;
    

    my $seqID = substr($res->{bla},0,30);
    my $typeA = substr($res->{type},0,1);
    my $typeB = substr($res->{type},1);
    my ($minStartA, $maxStartA, $minStartB, $maxStartB) = (1000,0,1000,0);
    my ($boxA, $boxB);
    my ($boxLengthA, $boxLengthB);
    my $tmpSeqLength = 0;
    my ($seq, $tmpSeq);
    my $meanL = 0;
    my $maxL = 0;
    my $minL = 1000;
    
    my @alignment;
    my @clustal;


    my ( $bool, $message ) = &makeAlignment( $fasta, $aln_file );
    if ( $bool == 0 ){
	print STDERR "\nSome ERROR occured during alignment computation!\n\n";
	return ( "", "", "", "", 0, "alignment computation error" );
    }

    #create a new alignment object using bio::AlignIO
    my $str = Bio::AlignIO->new(-file => "$aln_file", -format => "clustalw");
    my $aln = $str->next_aln();
    #print "\n\nALN:\n";

    my ($extendRangeA, $extendRangeB) = (0,0);

    #find the position of our query sequences in the clustalw alignment
    #that means we search for all alignment positions and create a window out of these positions
    ( $minStartA, $maxStartA, $meanL, $boxLengthA ) = getBoxStart( $Ref, \$aln, "1" );
    ( $minStartB, $maxStartB, $meanL, $boxLengthB ) = getBoxStart( $Ref, \$aln, "2" );


    $meanL /= (scalar keys %{$Ref});

    #print "meanL: $meanL\nmaxL:  $maxL\nminL:  $minL\n\n";

    #look if all boxes in our alignment start at the same position ($minStart == $maxStart)
    #is they do so we can go on and cut this block out
    #$seqID =~ s/(\d)\./$1\_/g;
    $seqID =~ s/[\(\),]/_/g;
#    print  "analyse\nid: ", $seqID, "\n";
#    print  "sc: ", scalar $aln->each_seq_with_id($seqID), "\n";

    #print "length: $seqLength\n";

    my ($returnA, $returnStartA, $returnB, $returnStartB);
    my ($tmpreturnA, $tmpreturnStartA, $tmppenaltyA, $tmpscoreA);
    my ($tmpreturnB, $tmpreturnStartB, $tmppenaltyB, $tmpscoreB);
    my ($penaltyA, $penaltyB);
    my ($scoreA, $scoreB) = (0,0);
    my $isOutlier = 0;


    if( !$muscle ){ 
	#if the option b is given
	#that means that we check our query sequences for paralogs -> we can be more restrictive
	if($res->{b}){
	    my $identity = checkSeqIdentity(\@clustal, scalar keys %{$Ref}, $res);
	    if($identity < 90){
		#`rm $aln_file`;
		return("", "", "", "", 0, "seq identity is too low");
	    }
	}
    }


    ######
    ## CHECK FOR OUTLIERs
    ######

    if ( !$muscle ) {
	my @outlier = badScore(\@clustal, scalar keys %{$Ref}) if scalar keys %{$Ref} > 2;
	if($outlier[0] && $outlier[0] == ((scalar keys %{$Ref})+1)){
	    $isOutlier = 1;
	}        
    }
    
    my $offset = 0;


    ######
    ## check if new sequences fits into the alignment 
    ######
    foreach my $seq_obj ( $aln->each_seq_with_id( $seqID ) ){

	my $seq = $seq_obj->seq();
#	print STDERR $seq_obj->id() , "\n$seq\n\n";
	(my $tmpSeq = $seq) =~ s/-//g;
	my $s;

	## if a sequence is shorter than 70% of the longest snoRNA we will reject it
	if(length $tmpSeq < $meanL * 0.7){
	    #`rm $aln_file`;
	    return("", "", "", "", "seq is too short");
	}

	my $boxRef = createBoxProfiles($Ref);

	$boxA = &extractSubsequence($seq, $minStartA, $boxLengthA);
	$penaltyA  = isCorrectBox($boxA, $typeA);
	if( $minStartA == $maxStartA && $penaltyA <= 1 ){
	    $returnA = $boxA;
	    $returnStartA = startPos($seq, $minStartA);
	    $scoreA =  scorePattern( $boxA, $boxRef->{'box1'} );
	}
	else{
#	    print STDERR "there are some position probs at boxA and/or boxA is not valid\n";
#	    print STDERR "explorer($seq, $boxRef->{'box1'}, $typeA, $maxStartA, $minStartA, \"6\")\n";
#	    print STDERR "EXPLORER $seq\n $maxStartA, $minStartA\n";
	    ($returnA, $returnStartA, $penaltyA, $scoreA) = explorer($seq, $boxRef->{'box1'}, $typeA, $maxStartA, $minStartA, "6");
	}

	
	## if there is no suitable box at the beginning of the analysed sequence
	## use fastacmd to enlarge the sequence (only if tha candidate is of type C/D)
	if( $penaltyA > 1 ){

#	    print STDERR "PENALTY A TOO BAD!:\nreturnA: $returnA\nreturnStartA: $returnStartA\npenaltyA: $penaltyA\n scoreA: $scoreA\n\n";
	    
	    if( $typeA eq "C" ){
		#enlarge range of extracted sequence by 5(?)
		my ($start, $end) = split(/,/, $res->{range});
		my $diff = 5;
		
		if($res->{strand} eq "+"){
		    if($start < 5){$start = 0; $diff = $start}
		    else{$start -= $diff}
		    ($s, $start, $end)  = extractSeq($tmpSeq, $res->{database}, $start, $end, "1", $res->{chr});;
		}
		else{
		    $end += $diff;
		    ($s, $start, $end) = extractSeq($tmpSeq, $res->{database}, $start, $end, "2", $res->{chr});
		}
		
		if($s ne $tmpSeq){
		    
		    $tmpSeq = $s;
		    $seq = substr($tmpSeq, 0, $diff).$seq;
		    ($tmpreturnA, $tmpreturnStartA, $tmppenaltyA, $tmpscoreA) = explorer($seq, $boxRef->{'box1'}, $typeA, $maxStartA, $minStartA, "6");
		    
		    #penaltyB equals 0 means that we found a perfect box at most 20n away of our alignment box position
		    #we now have to update our result hash
		    if($tmppenaltyA <= 0.5){
			#if a suitable box is found within the enlarged seq, adapt start positions for box B
			$minStartB += $diff;
			$maxStartB += $diff;
			$offset += $diff;
			($returnA, $returnStartA, $penaltyA, $scoreA) = ($tmpreturnA, $tmpreturnStartA, $tmppenaltyA, $tmpscoreA);
			$res->{seq} = $tmpSeq;
			$res->{range} = "$start,$end";
			$res->{bla} =~ s/\d+,\d+/$start,$end/;
		    }
		    else{
			#restore old seq if none suitable box was found
			$seq =~ s/^.{$diff}//;
		    }
		}
	    }
	    else{
		($returnA, $returnStartA, $penaltyA, $scoreA) = explorer($seq, $boxRef->{'box1'}, $typeA, $maxStartA, $minStartA, "10");
	    }
		
	}

	#do the same comparison for second box
	
	if($minStartB == $maxStartB){
	
	    $boxB = &extractSubsequence( $seq, $minStartB, $boxLengthB );
	    $penaltyB = isCorrectBox($boxB, $typeB);
	    if( $penaltyB <= 1){
		$returnB = $boxB;
		$returnStartB = startPos($seq, $minStartB);
		$scoreB = scorePattern( $boxB, $boxRef->{'box2'} );
	    }
	    else{
#		print STDERR "boxB is not valid: $boxB\n";
		($returnB, $returnStartB, $penaltyB, $scoreB) = explorer($seq, $boxRef->{'box2'}, $typeB, $maxStartB, $minStartB, "6");
	    }

	}
	else{

#	    print STDERR "there are some position probs at boxB\n";
#	    $offset = 0;
	    ($returnB, $returnStartB, $penaltyB, $scoreB) = explorer($seq, $boxRef->{'box2'}, $typeB, $maxStartB, $minStartB, "6");

	}

	if($penaltyB > 1){

#	    print STDERR "PENALTY B TOO BAD!:\nreturnB: $returnB\nreturnStartB: $returnStartB\npenaltyB: $penaltyB\n scoreB: $scoreB\n\n";
	    my ($start, $end) = split(/,/, $res->{range});
	    my $diff = 5;
	    
	    if($res->{strand} eq "+"){
		$end += $diff;
		($s, $start, $end) = extractSeq($tmpSeq, $res->{database}, $start, $end, "1", $res->{chr});
	    }
	    else{
		if($start < 5){$start = 0; $diff = $start;}
		else{$start -= $diff;}
		($s, $start, $end) = extractSeq($tmpSeq, $res->{database}, $start, $end, "2", $res->{chr});
	    }
	    if($s ne $tmpSeq){
		$tmpSeq = $s;
		$seq .= substr($tmpSeq, length($tmpSeq) - $diff);
		$maxStartB += $diff;
		$minStartB += $diff;
		
		my ($tmp_returnB, $tmp_returnStartB, $tmp_penaltyB, $tmp_scoreB) = explorer($seq, $boxRef->{'box2'}, $typeB, $maxStartB, $minStartB, "6");			
		
		#penaltyB equals 0 means that we found a perfect box at most 20n away of our alignment box position
		#we now have to update our result hash
		#accept a boxes if it has a penalty lower than 1 (or 1.5) and it is not further away than ?? from the range of all other sequences
		if($tmp_penaltyB <= 1){
		    ($returnB, $returnStartB, $penaltyB, $scoreB) = ($tmp_returnB, $tmp_returnStartB, $tmp_penaltyB, $tmp_scoreB);
		    $res->{seq} = $tmpSeq;
		    $res->{range} = "$start,$end";
		    $res->{bla} =~ s/\d+,\d+/$start,$end/;
		}
		else{$seq =~ s/.{$diff}$//;}
	    }
	}
	
    }
    $returnStartB++;
    $returnStartA++;

    #rules for accepting a sequence:
    #reject this sequence if 1) one or both boxes provides more than one mutation
    #                        2) if it is an outlier and the summarized score is larger than 1.5
    #                        3) or  
#    print STDERR "INFOS:\n";
#    print STDERR "A: $scoreA\nB: $scoreB\nO: $isOutlier\nPA: $penaltyA\nPB: $penaltyB\n";
#    print STDERR "boxA: $returnA\nboxB: $returnB\n\n";
#    print STDERR "startA: $returnStartA\nstartB: $returnStartB\noffset: $offset\n";

    ## some more cleaning
    #`rm $aln_file`;


    if($penaltyA <= 1 && $penaltyB <= 1){
	if(($scoreA >= -10 && $scoreB >= -10) || (!$isOutlier && ($penaltyA+$penaltyB<=1.5))){

	    #analyse extracted boxes whether they can form a kink-turn or not (CD-box snoRNAs)
	    if($typeA eq "C" && !&kinkTurn( $returnA, $returnB )){
		return ($returnA, $returnStartA, $returnB, $returnStartB, 0, $offset, "No K-turn found", $scoreA, $scoreB, $penaltyA, $penaltyB );
	    }

	    return ($returnA, $returnStartA, $returnB, $returnStartB, 1, $offset, "", $scoreA, $scoreB, $penaltyA, $penaltyB );
	}
	else{
	    return ($returnA, $returnStartA, $returnB, $returnStartB, 0, $offset, "no box motifs found", $scoreA, $scoreB, $penaltyA, $penaltyB );
	}
    }
    else{
	return ($returnA, $returnStartA, $returnB, $returnStartB, 0, $offset, "no box motifs found", $scoreA, $scoreB, $penaltyA, $penaltyB );
    }

}


##########################################################################################################################################
sub explorer
#
##########################################################################################################################################
{

    my ($seq, $boxRef, $type, $maxStart, $minStart, $width) = @_;
#    print STDERR "minStart: $minStart\tmaxStart: $maxStart\n";
    my ( $sequence, $cutoff_5prime );

    my $range = $maxStart - $minStart + $width + scalar(@{${$boxRef}[0]});

    if ( $minStart - 3 >= 0 ){ $cutoff_5prime = 3 }
    else{ $cutoff_5prime = $minStart }

    ## in order to compensate for badly aligned c-boxmotifs (especially for fungi)
    ## we elongate the search window by four bases on the 3'end
    if ( $type eq "C" ){
	$sequence = substr($seq, $minStart - $cutoff_5prime, $range + 4 );
    }
    else{
	$sequence = substr($seq, $minStart -3, $range );
    }
   
#    print  STDERR "Explorer: $sequence\t$type\n$cutoff_5prime\n\n";

    my ($returnBox, $startInTmpSeq, $penaltyBox, $scoreBox) = searchBoxes($sequence, $boxRef, $type);
#    print STDERR "startT: $startInTmpSeq\n";
    $startInTmpSeq += $minStart - $cutoff_5prime;
#    print "startT: $startInTmpSeq\n";
    
    my $cutSeq = substr($seq, 0, $startInTmpSeq);
#    print STDERR "cutSeq: $cutSeq\n";
#    print STDERR "laenge: ", length $cutSeq, "\n";
    $cutSeq =~ s/-//g;
#    print STDERR "cutSeq: $cutSeq\n";
#    print STDERR "laenge: ", length $cutSeq, "\n";
    my $returnStartBox = length $cutSeq;
	    
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

#   print STDERR "laenge: ",length $sequence, "\n";
#   print STDERR "sequen: ", $sequence, "\n";

    #variables declaration
    #use a sliding window to score all possible box-motifs and store the current leader
    #my @seq = split(//, $sequence);
    my ($curScore, $bestScore) = (-100,-100);
    my ($curBox, $bestBox) = ("","");
    my $startBox;
    my ($start, $end) = (0, (length($sequence))-scalar(@{${$boxRef}[0]}));
    my ($result, $tmpSeq);
    my ($curPenalty, $penalty) = (2,2);
#   print STDERR "start: $start\tend: $end\n";    

    for my $i ($start..$end){
	
	if($sequence =~ /^.{$i,$i}-/){next;}
	
	#calculate current score for this certain window
	($tmpSeq = $sequence) =~ s/.{$i,$i}//;
	$tmpSeq =~ s/-//g;
	$curBox = substr($tmpSeq, 0, scalar(@{${$boxRef}[0]}));
#	print STDERR $curBox ,"\t";
	$curScore = scorePattern($curBox, $boxRef) if $curBox !~ /[Nn]/;
#	print STDERR $curScore,"\n";
	$curScore = -100 if $curBox =~ /[Nn]/;
	
	$curPenalty = isCorrectBox($curBox, $type);
#	printf STDERR "%-2d: box: %-15s: % .2f\t%.2f\n", $i, $curBox, $curScore, $curPenalty;

	#compare current score with the best score achieved so far
	#look at both starting positions and compare them in order to decline a box when box starting points a too far away from each other
#	if( $curScore >= $bestScore && $curPenalty <= $penalty && $curPenalty <= 1 ){ #&& $start > $rangeRef->[0] && $start < $rangeRef->[1]){
	if( $curScore >= $bestScore && $curPenalty <= $penalty ){ #&& $start > $rangeRef->[0] && $start < $rangeRef->[1]){
#	    print STDERR "new $bestScore\n";
	    $bestScore = $curScore;
	    $bestBox = $curBox;
	    $startBox = $i;
	    $penalty = $curPenalty
	}
    }
    
    
    #test if the best box suites our constraints, which means only one mutation per box is allowed
    
#    print STDERR "\n\n", $bestBox , "\t", $bestScore ,"\tstarting at position: ", $startBox ,"\n\n";
    
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
	if($box =~ /TGATGA/){$penalty  = 0;}
	elsif($box =~ /TGA[CGA]GA|TGAT[CTA]A|TGATG[CTG]/){$penalty = 0.5;}
	elsif($box =~ /[CGA]GATGA/){$penalty = 1;}
	else{$penalty = 2;}
    }


    #boxes with one mutation will get penalty of 1
    #more than one mutation will result in a penalty of 2
    #boxes with mutated G or A nucleotides will get a penalty of 2
    elsif($type eq "D"){
	if($box =~ /CTGA/){$penalty = 0;}
	elsif($box =~ /[CTGA]TGA|C[AGC]GA/){$penalty = 1;}
	else{$penalty = 2;}
    }


    #when the third 'A' is mutated the box will get a penalty of 0.5, because this mutation appears quite often
    #otherwise the box will be penalized with 1 (if another 'A' is changed) or 2 (more than one mutation)    
	elsif($type eq "H"){
	if($box =~ /A[CTGA]A[CTGA]{2}A/){$penalty = 0;}
	elsif($box =~ /A[CTGA]A[CTGA]{3}/){$penalty = 0.5;}
	elsif($box =~ /[CTGA]{2}A[CTGA]{2}A|A[CTGA]{4}A/){$penalty = 1;}
	else{$penalty = 2;}
    }


    #no penalty is set if 'C' is changed to an 'T', cause it occurs frequently
    #other mutation at position 2 are penalized with 0.5
    #every other single mutation gets a penalty of 1
    #penalty of 2 is given to every box with more than one mutation
    elsif($type eq "ACA"){
	if($box =~ /A[CT]A/){$penalty = 0;}
	elsif($box =~ /A[GA]A/){$penalty = 1;}
	elsif($box =~ /[CTG][CT]A|A[CT][CTG]/){$penalty = 1.5;}
	else{$penalty = 2;}
    }


    if($box =~ /-/){$penalty = 2}

    return $penalty;

}




##########################################################################################################################################
sub kinkTurn
#returns true if the given boxes are able to form a kink turn protein binding motif (only for CD-box snoRNAs)
#otherwise false
##########################################################################################################################################
{

    my ($boxC, $boxD ) = @_;
    $boxC =~ tr/U/T/;
    $boxD =~ tr/U/T/;
    my @boxC = split(//,$boxC);
    my @boxD = split(//,$boxD);
    
    my $pair = $boxC[5]."-".$boxD[0];

#    print STDERR $boxC ,"\t", $boxD, "\t", $pair,"\n";

    #test if 5th nt of c box forms a watson-crick bp with 1st nt of d box
    if($pair eq "G-C" || $pair eq "C-G" || $pair eq "A-T" || $pair eq "T-A" || $pair eq "T-G" || $pair eq "G-T" ){
	if($boxC[4] eq "T" || $boxD[1] eq "T"){
	    #print "watson-crick: ",$pair,"\n";
	    #print "u-u pair: ",$boxC[4]," ",  $boxD[1],"\n";
	    return 1;
	}
    }
    
#    print "ERROR: No K-turn motif found!\n\n";
    return 0;


}


##########################################################################################################################################
sub checkSeqIdentity
#returns true if the the clustal score of this new sequence is larger than 60 according to our starting sequence
#otherwise false
##########################################################################################################################################
{
	my ($scoresRef, $nr, $res) = @_;
	$nr++;
#	my @scores = split(/\n/, $clustal);
	my $score = 0;
	#go through clustal output and 

#	print STDERR "paralogs: ", $res->{paralogs}, "\n";

	SEQ: for(my $i = 0; $i < scalar @{$scoresRef} - 1 ; $i++){ 
		PARA: for(my $h = 1; $h <= $res->{paralogs}; $h ++){
        	if($scoresRef->[$i] =~ /Sequences\s*\($h:$nr\)\s*Aligned.\s*Score:\s*(\d+)/){  #die interessanten Zeilen des Clustaloutputs
#				print STDERR $scores[$i] , "\n";
		    if($1 > $score){$score = $1;}
		    last PARA;
        	}
    	}
	}	
#	print STDERR "\n";

	return $score;

}


##########################################################################################################################################
sub badScore
#returns an array with positions bad aligned targetSequences in Alignment
#bad aligned is defined by clustalw-score
##########################################################################################################################################
{
    #arguments passed to this subroutine are:
    #clustal output
    #number of used sequences
    #boolean which suggests to analyse only the last sequence
    my ($clustalRef, $nr) = @_;    #der output den clustal beim neu alignen zurueckgibt, da stehen die scores fuer die Sequenzen drin
    $nr++;
#    my @scores = split(/\n/,$clustal);
    my %seqScore;
    my ($sum, $mean, $variance, $deviation);
    my $outlier;
    my @outlier;

     #initialisierung alles Sequenzen mit 0, start bei 1, da BioPerl erst bei 1 anfaengt zu zaehlen
    for(my $j = 1; $j < $nr; $j++){
        $seqScore{$j}=0;
    }
    for(my $i = 0; $i < scalar @{$clustalRef} - 1; $i++){ 
        if($clustalRef->[$i] =~ /Sequences\s*\((\d+):(\d+)\)\s*Aligned.\s*Score:\s*(\d+)/){  #die interessanten Zeilen des Clustaloutputs
	    $seqScore{$1}+=$3;   #zu jeder Sequenz werden alle pairwise scores addiert
            $seqScore{$2}+=$3;
        }
    }


    #calculate final score values and standard mean
    for my $i (1..$nr){
	$seqScore{$i} /= $nr;
    }
    zScore(\%seqScore);


    #use routine getOutlier() to check a set of values for possible outliers
    #do this as long as an outlier is found
    while($outlier = getOutlier(\%seqScore)){
	#print "OUTLIER found! $outlier!!\n";
	push @outlier, $outlier;
	delete $seqScore{$outlier};
	last if !(scalar keys %seqScore > 2)
    }

    return @outlier;

}


##########################################################################################################################################
sub zScore
#this method normalizes all values in a given hash to z-Scores
##########################################################################################################################################
{

    my $hashRef = shift;
    my $nr = scalar keys %{$hashRef};
    my @key = sort {$a<=>$b} keys %{$hashRef};
    my ($sum, $variance, $deviation, $mean);


    #calculate final score values and standard mean
    foreach(@key){
	$sum += $hashRef->{$_};
    }
    $mean =  $sum/$nr;


    #calculate variance
    $sum = 0.0;
    foreach (@key){
	$sum += ($hashRef->{$_}-$mean)*($hashRef->{$_}-$mean);
    }


    $variance = $sum/($nr-1);
    $deviation = sqrt $variance;
    $deviation = 0.00001 if $deviation == 0;


    #calculate z-scores as usual
    foreach (@key){
	$hashRef->{$_} = ($hashRef->{$_}-$mean)/$deviation;
    }

}


##########################################################################################################################################
sub getOutlier
#this method tries to identify outliers in a given normal distrubtion
##########################################################################################################################################
{

    #hashRef is a reference of a hash which contains key/value pairs as follows:
    #'key' is the sequence number used in the subjacent clustalw alignment
    #'value' is the normalized z-score achieved with this alignment
    my $hashRef = shift;
    my $nr = scalar keys %{$hashRef};
    my @key = sort {$a<=>$b} keys %{$hashRef};


    #calculate mean, variance and standard deviation of all z-scores
    my ($sum, $variance, $deviation, $mean) = (0.0,0.0,0.0,0.0);
    foreach (@key){
	#print "$_: ",$hashRef->{$_} , "\n";
	$sum += $hashRef->{$_};
    }
    $mean = $sum/$nr;
    $sum = 0.0;


    foreach (@key){
	$sum += ($hashRef->{$_}-$mean)*($hashRef->{$_}-$mean);
    }
    $variance = $sum/($nr-1);
    $deviation = sqrt $variance;

    if($deviation == 0){$deviation = 0.00001}

    #search for the minimal value and use this one for a outlier test
    my $min = $hashRef->{(keys %{$hashRef})[-1]};# $hashRef->{$nr};
    #print "nr: $nr\nscalar: ", scalar keys %{$hashRef}, "\n";
    my $minNr = (keys %{$hashRef})[-1];
    foreach (@key){
	if($hashRef->{$_}<$min){$min = $hashRef->{$_};$minNr = $_;}
    }

    
    #print "sum: $sum\t var: $variance\tdev: $deviation\tmean: $mean\tmin: $min\n\n";

    #use Grubbs test for outliers to check if the new sequence fits into our alignment
    #use a significance value of 1% for rejecteing our hypothesis that there is no outlier
    #print STDERR "qt(0.01/($nr),($nr-2))\n";
    my $r = qx(echo "qt(0.01/($nr),($nr-2))" | R --slave);
    $r =~ s/\[1\]\s(-?\d+\.?\d*)\D*/$1/;

    my $g = ($mean - $min)/$deviation;
    my $za = ($nr-1)/(sqrt $nr);
    $za *= sqrt(($r*$r)/($nr-2+($r*$r)));
    

    #print STDERR"\ngrubbs: $g\nCvalue: $za\n\n\n";
    if($g > $za){
	return $minNr;
    }
    else{return "";}

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
