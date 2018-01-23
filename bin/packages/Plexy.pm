package Plexy;


use strict;
use Bio::SimpleAlign;
use Bio::AlignIO;
use find_alignment_position;
use CONFIG;
use vars qw(@ISA @EXPORT);

require Exporter;


@ISA = qw(Exporter);
@EXPORT = qw(computeTargetCD findAlignPos);


####################################################################################################
## VARIABLES
####################################################################################################

my $RNAPLEX = $CONFIG::RNAPLEX;



####################################################################################################
sub computeTargetCD
#finds regions which can bind to the InteractionRegion using RNAple
####################################################################################################
{
    my ($res, $org, $locations ,$Dbox, $kingdom ) = @_;
    my $target = {};
    my $message = "";

    my $ia = &cutInteractionRegion( $Dbox, $res->{seq} );
    
    ## get abbreviation for current organism
    (my $viech = $org) =~ s/\./_/;

    ## collects all possible targetRNAs
    my @snRNAs = glob($locations->{targetRNApath}.$viech."_*U*.fa");
    my @rRNAs = glob($locations->{targetRNApath}.$viech."_*[0-9]S*.fa");
   
    ## in case no target sequences are known, target prediction is useless
    if( scalar @snRNAs == 0 && scalar @rRNAs==0 ){
	my $targets = [];
	$message = "no targetRNAs available";
	return ( $targets, $message );
    }

    my @targetRNAs = (@snRNAs,@rRNAs);
    my @results;

    ## for each possible targetRNA
    foreach my $tarRNA ( @targetRNAs ){

	## read targetSequence
	open(TAR,"<".$tarRNA);
	my @tar = <TAR>;
	close(TAR);
	
	## only name of targetRNA, so throw away the path
	my $tarRNAFile = $tar[0];
	chomp($tar[1]);
	
	## print file for RNAplex
	open(PLEX, ">".$locations->{tmpPath}."$$\_plex.in");
	print PLEX $tarRNAFile,$tar[1],"\n>InteractionRegion\n",$ia,"\n";
	close PLEX;
	
	chomp $tarRNAFile ;
	$tarRNAFile =~ s/>//;

	## rRNA-name
	if($tarRNAFile =~ /\dS.\d/){
	    ( $tarRNA = $tarRNAFile ) =~ s/.*_([0-9\.]+S.[0-9]+).*/$1/;
	}
	elsif( $tarRNAFile =~ /\dS$/ ){
	    ( $tarRNA = $tarRNAFile ) =~ s/.*_([0-9\.]+S).*/$1/;
	}
	## snRNA-name
	else{
	    if($tarRNAFile =~ /\./){
		$tarRNA = ( split( /_/, ( split(/\./,$tarRNAFile) )[0] ) )[2] . "." . ( split( /\./,$tarRNAFile ) )[1];
	    }
	    else{
		$tarRNA = (split(/_/,$tarRNAFile))[2];
	    }
	}

	my @hits = `$RNAPLEX < $locations->{tmpPath}$$\_plex.in -l 20 -e -7.00 -f 2 ;`;
	
	foreach my $hit ( @hits ){

	    chomp $hit;
	    
	    if( $hit =~ /^>/ || $hit =~ /^$/ ){
		next;
	    }
	    #the duplexes
	    else{
		my $hash;
		my @parts = split( /\s+/, $hit );
		$hash->{ase} = $ia;


		## destiguish between D and D' box motifs
		if( $Dbox == $res->{startB} ){
		    $hash->{boxMotif} = $res->{boxB};
		}
		else{
		    $hash->{boxMotif} = $res->{boxBPrime};
		}


		$hash->{tarFile} = $tarRNAFile;
		$hash->{targetRNA} = $tarRNA;
		($hash->{energy} = $parts[4]) =~ s/[\(\)]//g;
		
		$hash->{targetPos} = $parts[1];
		$hash->{bindPos} = $parts[3]; 
		$hash->{raw_structure} = $parts[0];
		$hash->{structure} = $parts[0];

		#structure of Binding
		my @struc = split(/&/,$hash->{raw_structure});
		my $tarStart = (split(/,/,$hash->{targetPos}))[0];
		my ($duplexStart, $duplexEnd) = (( split(/,/,$hash->{bindPos} ))[0], (split( /,/,$hash->{bindPos} ))[1]);
		
                #RNAplex-bug?
		if( $struc[0] =~ /^\./ && $struc[1] =~ /[^\.]$/ ){
		    $struc[0] =~ s/^\.//;
		    $tarStart++;
		}
		$hash->{structure} = $struc[0]."&".$struc[1];

                #duplex has to start at most 2 nts before methylated site
		if( $duplexEnd <= 18 ){
		    #print STDERR "DUPLEX TOO FAR FROM BOX:\t ",$duplexEnd,"\n";
		    next;
		}
		$duplexStart-- if($duplexStart >0);
		$hash->{bindSeq} = substr( $ia, $duplexStart , $duplexEnd - $duplexStart );
		$hash->{tarSeq} = substr( $tar[1], $tarStart - 1, (split(/,/,$hash->{targetPos}))[1] - $tarStart + 1 );
		my $distToBox = 20 - $duplexEnd;

		$hash->{mod} = $tarStart + 4 - $distToBox;
		$hash->{seq} = $hash->{tarSeq}."&".$hash->{bindSeq};

				
		## duplex-length has to be at least 7 basepairs long
		if( length($hash->{tarSeq}) < 7 ){
		    next;
		}
		
                ## at beginning and end of the targetSeq a bulge is allowed to occur
		my $tarStruc = substr( $struc[0], 1, length($struc[0]) - 2 );
		my $bindStruc =  substr( $struc[1], 1, length($struc[1]) - 2 );

		## throw away duplexes with right or left bulges
		if( &hasNoBulge( $tarStruc, $bindStruc ) == -1 || &hasNoBulge( $struc[1], $struc[0] ) == -1 ){
		    next;
		}
		
		## check if 5th position upstream of D/D'-box, so methylated nucleotide forms a Watson-Crick basepair
		$hash->{modnt} = ( split( //,$hash->{tarSeq} ) )[4 - $distToBox];
		my $ntPair = ( split( //,$hash->{bindSeq} ) )[-5 + $distToBox];
		if( ( $hash->{modnt} eq "A" && $ntPair ne "U" &&  $ntPair ne "T" ) || (( $hash->{modnt} eq "U" || $hash->{modnt} eq "T" ) && $ntPair ne "A" ) || ($hash->{modnt} eq "G" && $ntPair ne "C") || ($hash->{modnt} eq "C" && $ntPair ne "G") ){
		    next;
		}

		$hash->{tarSeq} = substr($hash->{tarSeq},0, 5 - $distToBox)."m".substr($hash->{tarSeq},5 - $distToBox, length($hash->{tarSeq}));
		$hash->{seq} = $hash->{tarSeq}."&".$hash->{bindSeq};
		
		## 20 nts long fragment of the target RNA, so the maximal possible interacting region
		my $cutStart =  ($tarStart - $distToBox -1);
		$cutStart = 0 if ($cutStart < 0 );
		
		my $cutLen = 20;
		if( length( $tar[1] ) < 20 + ($cutStart) ){
		    $cutLen = length( $tar[1] ) - ($cutStart);
		}
		else{
		    $cutLen = 20;
		}


		$hash->{tarcomp} = substr( $tar[1], $cutStart, $cutLen );
		$hash->{tarcomp} = substr($hash->{tarcomp},0,5)."m".substr($hash->{tarcomp},5,length($hash->{tarcomp}));

		## 20 nts of target RNA & 20 nts antisense element and adjacent box motif
		$hash->{seq2} = $hash->{tarcomp}."&".$hash->{ase}.$hash->{boxMotif};

		## allow only duplexes with at most 1(maybe change to two) mismatches in duplex core-region between 3rd and 11th nt upsteam ot D/D'-box
		if(grep( /\./,split(//,substr($struc[1],2,9))) > 1){
		    next;
		}

	        ## filter hits with energy lower than treshold
		if($hash->{energy} > -7){
		    next;
		}

		push(@results,$hash);
	    }
	}
    }
	
    ## sort results after energy-values
    my @sorted = sort {$a->{energy} <=> $b->{energy}} @results;
    
    ## make modificated sites unique
    my @uniq;
    my %seen = ();

    foreach  (@sorted ) {
	my $tmp = $_;
	
	## there are some variants for the snRNAs, but only the best for each snRNA should be regarded
	( my $tmptar = $tmp->{targetRNA} ) =~ s/\.\d+$//;
	my $item = $tmptar."-".$tmp->{mod};
	
	push( @uniq,$tmp ) unless $seen{$item}++;
    }
    
    my $targets;
    my $targetRNAs_missing = 0;

    if( scalar(@uniq) > 0 ){
	foreach( @uniq ){
	    $target = $_;
	    ## find position in the alignment
	    $target->{alnPos} = &map_position( $target->{tarFile}, $target->{targetRNA}, $target->{mod}, $viech, $locations );
	    ## print STDERR $target->{alnPos},"\n";
	    if( $target->{alnPos} =~ /\-Na/ ){
		$targetRNAs_missing++;
	    }
	    push(@{$targets},$target);
	}
    }
    else{
	($target->{targetRNA},$target->{tarSeq},$target->{mod},$target->{energy},$target->{alnPos},$target->{bindSeq},$target->{structure},$target->{tarSeq},$target->{raw_Structure}, $target->{seq}, $target->{seq2}) = ("NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA");
	push(@{$targets},$target);
    }
    
    if( $targetRNAs_missing != 0 ){
	$message = "at least one targetRNA is missing or incomplete";
    }

    ## clean
    `rm $locations->{tmpPath}/$$\_plex.in`;
    
    return ( $targets, $message );

}



####################################################################################################
sub hasNoBulge
#checks weather the duplex has an bulge
#returns 1 if there is no bulge and -1 otherwise
####################################################################################################
{

    my($struc1,$struc2) = @_;
    my $bool = 1;

    my $offset = 0;
    my $mismatch = index($struc1, ".", $offset);
    #for every occurence of a mismatch(dot)
    while ($mismatch != -1) {
	#there has to be a mismatch in the other sequence,too			
	if(substr(reverse($struc2),$mismatch,1) ne "."){
	    last;
	}
	$offset = $mismatch + 1;
	$mismatch = index($struc1, ".", $offset);
    }
    #jumped out of while-loop by last (otherwise mismatch would be -1), so a bulge was found
    if($mismatch != -1){
	$bool = -1;
    }
    
    return $bool;

}



####################################################################################################
sub cutInteractionRegion
#cuts the region upstream of the D or D' box, which does potentialy bind
#to RNA region which will be methylated
####################################################################################################
{

    my ( $boxPos, $seq ) = @_;

    ## cuts 20 nucleotides before D|D'-box from sequence
    my $interaction;
    if( $boxPos >  20 ){
	$interaction = substr( $seq, $boxPos - 21, 20 );
    } else{
	$interaction = substr( $seq, 0, $boxPos ); 
    }

    return $interaction;

}


1;
