package tools;

use strict;
use Bio::SeqIO;
use vars qw(@ISA @EXPORT);

require Exporter;

@ISA = qw( Exporter );
@EXPORT = qw( get_sequence_length );


############
# ROUTINES #
############



########################

sub get_sequence_length
# extracts the sequence length of fasta entries in a given fasta file
{

  ## get the given file
  my ( $file, $entry)= @_;

  ## check if the fasta file exists
  if ( !-e $file){
    print STDERR "ERROR: given fasta file does not exist!\n$file\n\n";
  }

  ## result_hash
  my %lengths;
  
  ## open file, go through all entries
  ## in case no name of a specific entry is given, extract the sequence length of all entries in the fasta file
  ## otherwise, get only the length of the specified entry
  my $seqin = new Bio::SeqIO( -format => 'fasta', -file => $file);
  while ( my $seq = $seqin->next_seq() ){
    if( !$entry || $entry eq $seq->display_id ){
      $lengths{ $seq->display_id } = $seq->length();
    }
  }

  return \%lengths;

}
