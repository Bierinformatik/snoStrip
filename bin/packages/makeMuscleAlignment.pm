package makeMuscleAlignment;

use Cwd 'abs_path';


BEGIN{
    push @INC, abs_path() . "/packages/";
}

use strict;
use vars qw(@ISA @EXPORT);
use Bio::AlignIO;
use Bio::SimpleAlign;
use CONFIG;

require Exporter;



@ISA = qw(Exporter);
@EXPORT = qw(makeAlignment makeSingleSeqAlignment);



###############################################################
## VARIABLES 
###############################################################

my ( $bool,$message, $fasta, $outfile, $MUSCLE, $cutting );
my ( @alignment );

$MUSCLE = $CONFIG::MUSCLE; 



###############################################################
sub makeAlignment
## this method computes a muscle alignment
## and reformats the output to fit the correct clustal format
## requires:  multifasta file, outputfile
## returns: 1 (successful alignment)
##          2 (some problems occured)
###############################################################
{

    ( $fasta, $outfile, $cutting ) = @_;
    $bool = 0;

    ## create muscle alignment
    ## -clwstrict --> clustalw format, including CLUSTALW header
    ## -quiet -> suppress output to STDERR
    $message = `$MUSCLE -in $fasta -out $outfile -clwstrict -quiet -seqtype rna 2>&1`;
    #print STDOUT "WAS: $message\n";

    if ( $message eq "" ){
	$bool = 1;
    }


    if( !$cutting){

	## open outputfile and reformat
	## therefore, escape special characters like '(', ')', ','
	## and cut the header to a total length of 30
	open ALIGNMENT, $outfile or die $!;
	@alignment = <ALIGNMENT>;
	close ALIGNMENT;
	
	`rm $outfile`;
	
	open ALIGNMENT, ">" . $outfile or die $!;
	foreach my $line ( @alignment ){
	    chomp $line;
	    $line =~ s/[\(\),]/_/g;
	    my @aln_line = split( //, $line );
	    if ( $line =~ /CD|HACA/ ){
		splice( @aln_line, 30, 2 );
	    } 
	    elsif ( scalar @aln_line > 1 && $line !~ /CLUSTAL/ ){
		splice( @aln_line, 0, 2 );
	    }
	    print ALIGNMENT join( "", @aln_line ), "\n";
	}
	close ALIGNMENT;
    }
    
    return ( $bool, $message ) if $bool == 0;
    return ( $bool, "" );

}


###############################################################
sub makeSingleSeqAlignment
## this method generates CLUSTAL formatted alignment file
## with solely one sequence
## returns: 1 (successful alignment)
##          2 (some problems occured)
###############################################################
{

    ( $fasta, $outfile ) = @_;

    my ( $header, $seq ) = split( /\n/, `cat $fasta` );
    $header = substr( $header, 1, 30 );
    $header =~ s/[\(\),]/_/g;
    chomp $seq;
    my @seq = split( //, $seq );
    
    open OUT, ">". $outfile or die $!;
    
    ## print the clustal header
    print OUT "CLUSTAL W (1.81) multiple sequence alignment\n\n";
    
    while ( @seq ){
	
	print OUT $header, " " x 5, splice( @seq, 0, 60 ), "\n\n";

    }
 
    close OUT;

    return ( "1", "" );

}
