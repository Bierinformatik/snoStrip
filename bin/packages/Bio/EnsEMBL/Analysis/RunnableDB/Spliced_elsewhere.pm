# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::Spliced_elsewhere;

=head1 SYNOPSIS

my $runnabledb = Bio::EnsEMBL::Analysis::RunnableDB::Spliced_elsewhere->new
  (			
   -db => $dbadaptor,
   -input_id => flag ids,	
   -analysis => $analysis
  );

$runnabledb->fetch_input();
$runnabledb->run();
$runnabledb->write_output();

=head1 DESCRIPTION

Spliced elsewhere checks single exon transcripts and looks for copies 
elsewhere in the genome that contain introns ie: processed pseudogenes.
The module runs on chunk files generated by scripts in: 
ensembl-personal/sw4/Scripts/Pseudogenes
which partition the database into single and multiexon genes.The single exon 
genes are written to flat files, the multiexon genes are written into the target
databse and also formatted as a blast db.

The module calls a gene spliced elsewhere if it is > 80% identical to another EnsEMBL
transcript with > 80 % coverage and the target gene has a genomic span > 3x the 
span of the query. These numbers are configurable through:
Bio::EnsEMBL::Analysis::Config::Pseudogene

The databses used by the module are configured through:
Bio::EnsEMBL::Analysis::Config::Databases;

Runs as a part of the larger pseudogene analysis. Uses flagged single exon genes
identified by Bio::EnsEMBL::Analysis::RunnableDB::Pseudogene_DB.pm

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=cut


package Bio::EnsEMBL::Analysis::RunnableDB::Spliced_elsewhere;


use strict;
use Bio::EnsEMBL::Analysis::RunnableDB::Pseudogene_DB;
use Bio::EnsEMBL::Analysis::Config::Pseudogene; 
use Bio::EnsEMBL::Analysis::Runnable::Spliced_elsewhere;
use Bio::EnsEMBL::Pipeline::DBSQL::FlagAdaptor;
use Bio::EnsEMBL::Analysis::Runnable::BaseExonerate; 
use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;
use vars qw(@ISA);

@ISA = qw( Bio::EnsEMBL::Analysis::RunnableDB::Pseudogene_DB Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild);




=head2 fetch_input

Title   :   fetch_input
  Usage   :   $self->fetch_input
 Function:   Fetches input data for Spliced_elsewhere.pm from the database and flat files
  Returns :   none
  Args    :   none

=cut

sub fetch_input{
  my ($self)=@_;
  my ($start, $end);
  my @genes;
  my $count=0;
  my %parameters;
  if ($self->parameters_hash) {
    %parameters = %{$self->parameters_hash};
  }
  my $runname = "Bio::EnsEMBL::Analysis::Runnable::Spliced_elsewhere"; 

  my $dna_db = $self->get_dbadaptor("REFERENCE_DB") ;

  #genes come from final genebuild database
  my $genes_db = $self->get_dbadaptor($self->PS_INPUT_DATABASE);
  $self->gene_db($genes_db); 

  my $ga = $genes_db->get_GeneAdaptor;
  my $fa = Bio::EnsEMBL::Pipeline::DBSQL::FlagAdaptor->new($self->db);
  my $ids = $fa->fetch_by_analysis($self->analysis);
  $self->throw("No flags found for analysis " .
	       $self->SPLICED_ELSEWHERE_LOGIC_NAME ."\n")  unless (scalar(@{$ids}>0));
  if ($self->input_id =~ /(\d+):(\d+)/) {
    $start = $1;
    $end = $2;
  } else {
    $self->throw("Input id not recognised\n");
  }
  # get ids
  foreach my $flag (@{$ids}) {
    if ($flag->dbID >= $start && $flag->dbID <= $end) {
      $count++;
      my $gene = $ga->fetch_by_dbID($flag->ensembl_id);
      push @genes, $self->lazy_load($gene);
    }
  }
  print "$count genes retrieved\n";
  my $runnable = $runname->new 
    (
     '-genes'             => \@genes,
     '-analysis'          => $self->analysis,
     '-PS_MULTI_EXON_DIR' => $self->PS_MULTI_EXON_DIR,
    );

  $self->runnable($runnable);
  return 1;
}

=head2 run

Title   :   run
  Usage   :   $self->run
 Function:   Overrides run method in parent class
  Returns :   none
  Args    :   none

=cut

sub run {
  my($self) = @_;
  foreach my $runnable (@{$self->runnable}) {
    $self->throw("Runnable module not set") unless ($runnable);
    $runnable->run;
    $self->parse_results($runnable->output);
    $self->output($self->real_genes);
    # If you want to store retrotransposed genes for psilc
    return 1 unless ($self->retro_genes);
    if ($self->RETROTRANSPOSED ) {
      # Store gene dbIDS for PSILC 
      my @id_list;
      foreach my $gene (@{$self->retro_genes}) {
	push @id_list, $gene->dbID;
      }
      $self->store_ids(\@id_list,"psilc");
    } else {
      # store retrotransposed genes as pseudo or real depending on config
      foreach my $gene (@{$self->retro_genes}) {
	if ($self->RETRO_TYPE eq 'pseudogene') {
	  $self->pseudo_genes($gene);
	} else {	
	  $gene->type($self->RETRO_TYPE);
	  $self->output([$gene]);
	}
      }
    }
    $self->output($self->pseudo_genes);
  }
  return 1;
}

# Blast parsing done in runnable db as it needs acess to database to determine
# gene span

=head2 parse_results

Title   :   parse_results
  Usage   :   $self->parse_results
 Function:   Parses blast results to find spliced copies.
  Tests percent ID, coverage and span.
  Span is calculated as the distance between two genomic cordinates 
  corresponding to the HSP subject start and end positions
  Returns :   none
  Args    :   none

=cut

sub  parse_results{
  my ($self,$results)=@_;
  my $ta = $self->gene_db->get_TranscriptAdaptor;
  my $ga = $self->gene_db->get_GeneAdaptor;
 RESULT:foreach my $result_hash_ref (@{$results}) {
    my %result_hash = %{$result_hash_ref};
    my %trans_type ;
    next RESULT unless ( %result_hash );
  TRANS: foreach my $id (keys %result_hash){
      my $retro_trans = $ta->fetch_by_dbID($id);
      my $retro_span = $retro_trans->cdna_coding_end- $retro_trans->cdna_coding_start;
      # catch results where the blast job failed
      if ( $result_hash{$id} eq 'NONE' ){
	# gene is very unique we have to call it real really
	push @{$trans_type{'real'}},$retro_trans;
	next TRANS;
      }
      if ( $result_hash{$id} eq 'REPEATS' ) {
	# gene is covered with low complexity repeats probably - let it through
	push @{$trans_type{'real'}},$retro_trans;
	next TRANS;	
      }

      my @dafs =  @{$result_hash{$id}};
      @dafs = sort {$a->p_value <=> $b->p_value} @dafs;

    DAF: foreach my $daf (@dafs) {
	my $retro_coverage = 0;
	my $real_coverage = 0;
	# is the percent id above threshold?
	next DAF unless ($daf->percent_id > $self->PS_PERCENT_ID_CUTOFF);
	next DAF unless ($daf->p_value <  $self->PS_P_VALUE_CUTOFF);
	# dont want reverse matches
	next DAF unless ($daf->strand == $daf->hstrand ) ;
	# tighten up the result set
	next DAF unless ($daf->percent_id > 80 or ( $daf->percent_id > 70 and $daf->p_value < 1.0e-140 ));

	my @align_blocks = $daf->ungapped_features;
        # is the coverage above the threshold?
	# There are a few things to consider the coverage of the PROTEIN sitting in the genomic say > 50 %
	# There is how much GENOMIC dna has been aligned as well - thats a good sign it is a retro gene
	# lets work out the coverage exon at a time
	foreach my $e ( @{$retro_trans->get_all_Exons} ) {
	  # put the exons onto the same slice as the dafs
	  my $slice = $retro_trans->feature_Slice;
	  my $exon = $e->transfer($slice);

	  foreach my $block ( @align_blocks ) {
	    # compensate for the 1kb padding on the retro transcript
	    $block->start($block->start-1000);
	    $block->end($block->end-1000);
	    $real_coverage+= $block->length;

	    if ( $exon->overlaps($block) ) {
	      if ( $exon->start < $block->start && $exon->end < $block->end ) {
		#            Bs|---------------|------Be
		#         Es---|===============|Ee
		$retro_coverage+= $exon->end - $block->start ;
	      } elsif ( $exon->start >= $block->start && $exon->end >= $block->end ) {
		#    Bs-----|---------------|Be
		#         Es|===============|-----Ee		
		$retro_coverage+= $block->end - $exon->start;
	      } elsif ( $exon->start < $block->start && $exon->end >= $block->end ) {
		#         Bs|---------------------|Be
		#     Es----|=====================|-----Ee
		$retro_coverage+= $block->end - $block->start;
	      } elsif ( $exon->start >= $block->start && $exon->end < $block->end ) {
		#     Bs|---|---------------------|---|Be
		#         Es|=====================|Ee
		$retro_coverage+= $exon->end - $exon->start;
	      }
	    }
	  }
	}
	my $coverage = int(($retro_coverage / $retro_trans->length)*100);

	my $aligned_genomic = $real_coverage - $retro_coverage;
	next DAF unless ($coverage > $self->PS_RETOTRANSPOSED_COVERAGE);
	next DAF unless ($aligned_genomic > $self->PS_ALIGNED_GENOMIC);

	my $real_trans;
	# Throw if transcript cannot be found
	eval{
	  $real_trans =   $ta->fetch_by_stable_id($daf->hseqname);
	  unless ( $real_trans ) { 
             $real_trans =   $ta->fetch_by_dbID($daf->hseqname);
          } 
        };
	if ($@) {
	  $self->throw("Unable to find transcript $daf->hseqname \n$@\n");
	  next;
	}

	##################################################################
	# real span is the genomic span of the subject HSP
	# hstart+3 and $daf->hend-3 move inwards at both ends of the HSP by
	# 3 residues, this is so that if 1 residue is sitting on a new exon
	# it wont include that in the span, needs to overlap by 3 residues 
	# before it gets included
	
	my $real_span;
	my @genomic_coords = $real_trans->cdna2genomic($daf->hstart+9,$daf->hend-9);
	@genomic_coords = sort {$a->start <=> $b->start} @genomic_coords;
	$real_span = $genomic_coords[$#genomic_coords]->end - $genomic_coords[0]->start;
	
	# Is the span higher than the allowed ratio?
	if ($real_span / $retro_span > $self->PS_SPAN_RATIO ) {
	  print STDERR "---------------------------------------------------------------------------------------------------\n";
	print STDERR "DAF: " .
	  $daf->start . " " .
	    $daf->end . " " .
	      $daf->strand . " " .
		$daf->hstart . " " .
		  $daf->hend . " " .
		    $daf->hstrand . " " .
		      $daf->hseqname  . " " .
			$daf->cigar_string . "\n";
	  print STDERR "transcript id ". $retro_trans->dbID." matches transcriptid " . $daf->hseqname . " at "
	    . $daf->percent_id . " %ID and $coverage % coverage with $aligned_genomic flanking genomic bases aligned pseudogene span of $retro_span vs real gene span of $real_span\n";
	  print STDERR "P-value = " . $daf->p_value . "\n";
	  print STDERR  "Coverage $coverage% - $retro_coverage bp of gene and $aligned_genomic bp of genomic $real_coverage total\n";
	  print STDERR ">".$retro_trans->dbID;
	  print STDERR "\n".$retro_trans->feature_Slice->expand(1000,1000)->get_repeatmasked_seq->seq."\n";
	  print STDERR ">".$real_trans->dbID;
	  print STDERR "\n".$real_trans->seq->seq."\n";	 
	  # Gene is pseudo, store it internally
	  push @{$trans_type{'pseudo'}},$retro_trans;
	  next TRANS;
	}
      }
      # transcript is real
      push @{$trans_type{'real'}},$retro_trans;	
    }

    unless (defined($trans_type{'real'})) {
      # if all transcripts are pseudo get label the gene as a pseudogene
      my $gene = $ga->fetch_by_transcript_id($trans_type{'pseudo'}[0]->dbID);
      my @pseudo_trans = @{$gene->get_all_Transcripts};
      @pseudo_trans = sort {$a->length <=> $b->length} @pseudo_trans;
      my $only_transcript_to_keep = pop  @pseudo_trans;
      $self->transcript_to_keep($only_transcript_to_keep);
      foreach my $pseudo_transcript (@pseudo_trans) {
	my $blessed = $self->_remove_transcript_from_gene($gene,$pseudo_transcript);
	if ( $blessed ) {
	  print STDERR "Blessed transcript " . $pseudo_transcript->display_id . 
	    " is a retro transcript \n";
	}
      }
      $gene->biotype($self->RETRO_TYPE);
      $self->retro_genes($gene);
      next RESULT;
    }
    # gene has at least one real transcript so it is real
    $self->real_genes($ga->fetch_by_transcript_id($trans_type{'real'}[0]->dbID));
  }
  return 1; 
}


######################################
# CONTAINERS


=head2 retro_genes

Arg [1]    : Bio::EnsEMBL::Gene
 Description: get/set for retrogenes 
  Returntype : Bio::EnsEMBL::Gene
  Exceptions : none
  Caller     : general

=cut

sub retro_genes {
  my ($self, $retro_gene) = @_;
  if ($retro_gene) {
    unless ($retro_gene->isa("Bio::EnsEMBL::Gene")){
      $self->throw("retro gene is not a Bio::EnsEMBL::Gene, it is a $retro_gene");
    }
    push @{$self->{'_retro_gene'}},$retro_gene;
  }
  return $self->{'_retro_gene'};
}





1;