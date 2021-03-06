#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Compara::RunnableDB::CreateHomology_dNdSJobs

=cut

=head1 SYNOPSIS

my $aa = $sdba->get_AnalysisAdaptor;
my $analysis = $aa->fetch_by_logic_name('CreateHomology_dNdSJobs');
my $rdb = new Bio::EnsEMBL::Compara::RunnableDB::CreateHomology_dNdSJobs(
                         -input_id   => [[1,2,3,14],[4,13],[11,16]]
                         -analysis   => $analysis);

$rdb->fetch_input
$rdb->run;

=cut

=head1 DESCRIPTION

This is a homology compara specific runnableDB, that based on an input
of arrayrefs of genome_db_ids, creates Homology_dNdS jobs in the hive 
analysis_job table.

=cut

=head1 CONTACT

abel@ebi.ac.uk, jessica@ebi.ac.uk

=cut

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Compara::RunnableDB::CreateHomology_dNdSJobs;

use strict;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Hive;

use Bio::EnsEMBL::Hive::Process;
our @ISA = qw(Bio::EnsEMBL::Hive::Process);

sub fetch_input {
  my( $self) = @_;

  $self->{'species_sets_aref'} = undef;
  $self->throw("No input_id") unless defined($self->input_id);

  #create a Compara::DBAdaptor which shares the same DBI handle
  #with the Pipeline::DBAdaptor that is based into this runnable
  $self->{'comparaDBA'} = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(-DBCONN=>$self->db->dbc);

  $self->get_params($self->input_id);
  return 1;
}

sub get_params {
  my $self         = shift;
  my $param_string = shift;

  return unless($param_string);
  print("parsing parameter string : ",$param_string,"\n");
  
  my $params = eval($param_string);
  return unless($params);

  foreach my $key (keys %$params) {
    print("  $key : ", $params->{$key}, "\n");
  }

  if (defined $params->{'species_sets'}) {
    $self->{'species_sets_aref'} = $params->{'species_sets'};
  }
  if (defined $params->{'method_link_type'}) {
    $self->{'method_link_type'} = $params->{'method_link_type'};
  }
  
  return;
}

sub run
{
  my $self = shift;
  return 1 unless($self->{'species_sets_aref'});

  $self->create_analysis_jobs($self->{'species_sets_aref'});
  
  return 1;
}

sub write_output {
  my $self = shift;
  return 1;
}

##########################################
#
# internal methods
#
##########################################

sub create_analysis_jobs {
  my $self = shift;
  my $species_sets_aref = shift;

  my $aa = $self->db->get_AnalysisAdaptor;
  my $Homology_dNdS_analysis = $aa->fetch_by_logic_name('Homology_dNdS');

  my $sql = "insert ignore into analysis_job (analysis_id,input_id,status) select " .
    $Homology_dNdS_analysis->dbID .
      ",homology_id,'READY' from homology where method_link_species_set_id = ?";

  my $sth = $self->db->dbc->prepare($sql);

  my $mlssa = $self->{'comparaDBA'}->get_MethodLinkSpeciesSetAdaptor;
  foreach my $species_set (@{$species_sets_aref}) {
    while (my $genome_db_id1 = shift @{$species_set}) {
      foreach my $genome_db_id2 (@{$species_set}) {
        my $mlss = $mlssa->fetch_by_method_link_type_genome_db_ids($self->{'method_link_type'},[$genome_db_id1,$genome_db_id2]);
	next unless($mlss);
        $sth->execute($mlss->dbID);
      }
    }
  }

  $sth->finish;
}

1;
