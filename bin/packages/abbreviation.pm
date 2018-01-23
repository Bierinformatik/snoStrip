package abbreviation;


use strict;
use vars qw(@ISA @EXPORT);
require Exporter;


@ISA = qw(Exporter);
@EXPORT = qw( get_abbreviation_Gepithet get_scientific_name );


###############################################################################################
sub get_abbreviation_Gepithet
## this method return the 3letter-abbreviation of a given organism
## the abbreviation is taken from the information file given as snostrip argument 
###############################################################################################
{

    my ( $organism, $information_file ) = @_;
    my ( $grep, $abbreviation, $genus, $epithet );
    

    ## split the organism into the first letter of the genus
    ## and the whole epithet
    ( $genus, $epithet ) = split( /\./, $organism ); 
   
    ## try to grep the corresponding line to our current organism
    ## within the information file
    $grep = qx(grep -P "$genus\[a-z\]*_$epithet" $information_file | cut -f 3);

    chomp $grep;

    return $grep;
    
    
}



###############################################################################################
sub get_scientific_name
## this method return the 3letter-abbreviation of a given organism
## the abbreviation is taken from the information file given as snostrip argument 
###############################################################################################
{

    my ( $abbreviation, $information_file ) = @_;
    my ( $species_name );
    

    ## try to grep the corresponding line to our current organism
    ## within the information file
    $species_name = qx(grep -P "\\s$abbreviation\\s" $information_file | cut -f 2 );

    chomp $species_name;
    $species_name =~ s/[a-z]+\s/\./;

    return $species_name;
    
    
}


1;
