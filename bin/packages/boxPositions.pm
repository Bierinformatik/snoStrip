package boxPositions;

use strict;
use vars qw(@ISA @EXPORT);
require Exporter;
use Bio::AlignIO;
use Bio::SimpleAlign;


@ISA = qw(Exporter);
@EXPORT = qw(getBoxStart);

##########################################################################################################################################
sub getBoxStart
## this method 
##########################################################################################################################################
{

    ## initialize the arguments
    my ( $Ref, $aln_ref, $box_type, $consensus ) = @_;

    ## VARIABLES
    my ( $box_length, $start, $box, $tmp_seq, $seq, $seq_length, $pos, $tmp_box, $correct_pos, $letters,
         $meanL, $extend_range, $max_start, $min_start );

    my %consensus_position if $consensus;
    my $consensus_start if $consensus;


    ## initialize several variables
    $extend_range = 0;
    $meanL = 0;
    $max_start = 0;
    $min_start = 1000;

    ## defince the correct Box
    if ( $box_type eq "1"){ $box_type = "box1" }
    elsif ( $box_type eq "2" ){ $box_type = "box2" }
    elsif ( $box_type eq "prime1" ){ $box_type = "box1Prime" }
    else{ $box_type = "box2Prime" }

##    print "\n\n\nbox analyse:\n\nscalar: " , scalar keys %{$Ref}, "\n\n";

    ## go through the alignment and search for the start positions
  HASH: foreach (keys %{$Ref}){

      #print STDERR "key: " , $_, "\n";
      #print STDERR "box: " , $$Ref{$_}->{$box_type . "Start"}, "\n";

        ## check if this sequences contains a prime box
        if( $$Ref{$_}->{$box_type . "Start"} ){

	     $start = $$Ref{$_}->{$box_type . "Start"} - 1;
	     if( $start <= 0 ){ $start=1; }
	
	     $box_length = length $$Ref{$_}->{$box_type};
	     $box = $$Ref{$_}->{$box_type};

	     #print STDERR "Summary: start: $start, box: $box, len: $box_length\n";
	     
	     ## change header in case the alignment was computed with clustalw
	     $_ =~ s/[\(\),]/_/g;
	
	     ## retrieve the gap-free candidate sequence
	     #print STDERR "SUCHE: $_ \n";
	     $tmp_seq  = ( $$aln_ref->each_seq_with_id($_) )[0]->seq();
	     ($seq = $tmp_seq) =~ s/-//g;
	
	     #print STDERR "candidate sequence: $seq\n";

	     $seq_length = length $seq;
	     
	     ## get the corresponding alignment position for the start position of the box
	     $pos = $$aln_ref->column_from_residue_number("$_", $start);
	
	     #print STDERR "position der box im aln: $pos\n";
	
	     #check if alignment position is correct in order to avoid errors like "----GTGA" were the first - is chosen instead of G
	     $tmp_box = substr( $tmp_seq, $pos, $box_length );
	     if( $tmp_box ne $box ){
	
		  if( $tmp_box =~ /^-/ ){
		      $correct_pos = 0;
		      while ( !$correct_pos ){
			  $pos++;
			  $tmp_box = substr( $tmp_seq, $pos, $box_length );
			  if( $tmp_box =~ /^\w/ ){ $correct_pos = 1; }
			  last if $pos + $box_length >= length( $tmp_seq ); 
		      }
		  }
	
		  #check if extracted seq contains dashes in order to enlarge the scoring-window
		  elsif( $tmp_box =~ /-/ ){
		      $letters = length( $tmp_box =~ s/-//g ) + 1;
		      $extend_range = $pos + $box_length - 1;
		      while( $letters != $box_length ){
			  $extend_range++;
			  if( substr( $tmp_seq, $extend_range, 1) ne "-" ){
			      $letters++;
			  }
		      }
		      $extend_range = $extend_range - $pos - $box_length + 1;
		  }
	
	     }
	     
	     #print STDERR "len4: $box_length\n";
	     #print STDERR "summary : current: $pos\n";
	     #print STDERR "max: $max_start\tmin: $min_start\n\n";

	     if( $pos > $max_start ){ $max_start = $pos; }
	     if( $pos < $min_start ){ $min_start = $pos; }
	
	     ## try to find the start point where most of the sequences have the box start
	     if ( $consensus ){
		 #print STDERR "len5: $box_length\n";
		  if ( !$consensus_position{$pos} ){ $consensus_position{$pos} = 1 }
		  else{ $consensus_position{$pos}++ }
		 #print STDERR "len5-2: $box_length\n";
	     }
	
	     $meanL += $seq_length;
        }
    }

    ## retrieve the position with the most counts
    if ( $consensus ){
	
	#print STDERR "len6: $box_length\n";
	## identify the highest count
	$consensus_start = (sort {$b<=>$a} values %consensus_position )[0];

	## extract the hash-key that corresponds to the highest count.
	my @keys = grep { $consensus_position{$_} == $consensus_start } keys %consensus_position;
	$consensus_start = $keys[0];

    }
    

    #print STDERR "ENDE: $min_start, $max_start, $meanL, $box_length, $consensus_start, box: $box_type\n\n";

    return ( $min_start, $max_start, $meanL, $box_length ) if !$consensus;
    
    return ( $consensus_start ) if $consensus;
    
      

}

1;
