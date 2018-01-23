package find_alignment_position;

use strict;
use CONFIG;
use Bio::AlignIO;
use vars qw(@ISA @EXPORT);

require Exporter;


@ISA = qw( Exporter );
@EXPORT = qw( map_position );


################################################################################
sub map_position
## maps the given single sequence position onto the
## corresponding position in the rRNA alignment
################################################################################
{

    my ( $target_RNA_file, $target_RNA, $position, $organism, $locations ) = @_;
 
    ## declare VARIABLES
    my ( $alignment_position, $alignment_path, $sequence_name, $alignment, $ALN_PATH, $TARGET_RNA_PATH );

    
    $TARGET_RNA_PATH = $locations->{targetRNApath};
    $ALN_PATH = $locations->{targetAlignmentPath};
    
    ## chomp trailing numbers like '.1' in U6.1
    $target_RNA =~ s/\.\d+$//;
    
    ## get alignment name and path
    ## and try to construct the presumed alignment name and path
	$alignment_path = $ALN_PATH . "ALL." . $target_RNA . ".aln";
	$sequence_name = $target_RNA_file;

	## OLD VERSION
    #if( $target_RNA =~ /S$/ ){
    #	$alignment_path = $ALN_PATH . "ALL." . $target_RNA . ".aln";
    #	( $sequence_name = `grep '^>' $TARGET_RNA_PATH$organism"_"$target_RNA".fa"` ) =~ s/^>//;
    #}
    #else{
	#	$alignment_path = $ALN_PATH . "ALL." . $target_RNA . ".aln";
	#	$sequence_name = $target_RNA_file;
	#	$sequence_name =~ s/[\._]/-/g;
    #}
    chomp( $sequence_name );
    
    ##  try to read the presumed alignment with BioPerl
    eval{
	$alignment = Bio::AlignIO -> new ( -file => $alignment_path, -format => 'clustalw' ) -> next_aln();
    };

    ## in case there is no alignment
    ## append 'NA' as alignment position
    if( $@ ne "" ){
	$alignment_position = "NA";
    }

    ## try to map the sequence position to an alignment position
    eval{
	$alignment_position = $alignment -> column_from_residue_number($sequence_name, $position);

    };
    
    ## in case the current sequence is not contained in the current alignment
    ## append 'NA' as alignment position
    if($@ ne ""){
	$alignment_position = "NA";
    }
    
    ## return alignment position
    return $target_RNA . "-" . $alignment_position;
    
}


1;
