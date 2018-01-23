package sequences;

use strict;
use CONFIG;
use vars qw(@ISA @EXPORT);
require Exporter;


@ISA = qw(Exporter);
@EXPORT = qw(startPos extractSubsequence extractSeq scorePattern);

my $FASTACMD = $CONFIG::FASTACMD;


##########################################################################################################################################
sub extractSubsequence
## given a sequence, extract the putative boxmotif 
## begining at a certain startposition
## be careful with gap characters
## return: sequence motif
##########################################################################################################################################
{

    ## VARIABLES
    my ( $seq, $start, $motif_length ) = @_;
    
    my $motif = "";
    my $remaining_length = length($seq) - $start - $motif_length;
    my $extend_range = 0;

    while ( $remaining_length > 0 ){

	## extract sequence motif
	$motif = substr( $seq, $start, $motif_length + $extend_range );

	## remove gap characters from motif
	$motif =~ s/\-//g;

	## check if motif has the required length and exit
	last if length($motif) == $motif_length;

	## addapt variables
	$extend_range++;
	$remaining_length--;

    } 

    return $motif;

}



##########################################################################################################################################
sub startPos
## calculate the position of a base in a sequence according to its position in an alignment
## return: sequence specific position
##########################################################################################################################################
{

    ## VARIABLES
    my ( $seq, $position ) = @_;

    my $part = substr($seq, 0, $position);
    $part =~ s/\-//g;

    return (length $part);

}



##########################################################################################################################################
sub scorePattern
## score a given pattern with a given matrix
## return: score or probability
##########################################################################################################################################
{

    ## VARIABLES
    my $sequence = shift;
    my $boxRef = shift;

    my $score;
    my @seq = split( //, $sequence );


    if( $seq[0] eq "A" ){ $score = $boxRef->[0][0] }
    if( $seq[0] eq "G" ){ $score = $boxRef->[1][0] }
    if( $seq[0] eq "T" ){ $score = $boxRef->[2][0] }
    if( $seq[0] eq "C" ){ $score = $boxRef->[3][0] }

    
    for (my $i = 1; $i < scalar @seq; $i++){
	if( $seq[$i] eq "A" ){ $score += $boxRef->[0][$i] }
	if( $seq[$i] eq "G" ){ $score += $boxRef->[1][$i] }
	if( $seq[$i] eq "T" ){ $score += $boxRef->[2][$i] }
	if( $seq[$i] eq "C" ){ $score += $boxRef->[3][$i] }
    }


    return $score;
}




##########################################################################################################################################
sub extractSeq
## extract a sequence from genome file
## if end point is larger than contig extract as much nucleotides as possible
## return: genomic sequence, start, end positions
##########################################################################################################################################
{

    ## VARIABLES
    my ($seq, $db, $start, $end, $strand, $chr) = @_;

    my $tmpSeq;
    my $bool = 0;
    
    ## use fastacmd to extract the genomic sequence
    $tmpSeq = `$FASTACMD -d $db -s $chr -L "$start,$end" -S $strand -l 500 2>&1`;

    ## check for errors
    if($tmpSeq =~ /ERROR/){ $bool = 1; $tmpSeq =~ s/^\[fastacmd\][^\n]+\n+//;}
    $tmpSeq =~ s/^>[^\n]+\n//;
    $tmpSeq =~ s/\s//g;

    #test if fastacmd returned the whole chromosome, which would be the case if an error would have occured
    if( !$tmpSeq || $end-$start <= (length $tmpSeq) - 2 || $bool){
	$end = (length $tmpSeq) + 2;
	$seq = substr($tmpSeq, $start - 1) if $strand == 1;
	$seq = substr($tmpSeq, 0, $end - $start ) if $strand == 2;
	return ($seq, $start, $end);
    }

    return ($tmpSeq, $start, $end);

}


1;
