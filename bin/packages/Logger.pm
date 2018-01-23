package Logger;

use strict;
use vars qw( @ISA @EXPORT );

require Exporter;


@ISA = qw( Exporter );
@EXPORT = qw( writeLog );


########################################################################
## VARIABLES 
########################################################################

my ( $message, $locations );



########################################################################
sub writeLog
#this method tries to analyse whether there is already a snoRNA or not
########################################################################
{
    
    ( $message, $locations ) = @_;

    chomp $message;

    if( !-e $locations->{workingDirectory}."LOG" ) {
	`touch $locations->{workingDirectory}"LOG"`;
    }

    open ( OUT, ">>".$locations->{workingDirectory}."LOG" ) or print $!;

    printf OUT "%02d-%02d-%04d\t",(localtime)[3],((localtime)[4] +1),((localtime)[5]+1900);
    printf OUT "%02d:%02d\t\t", (localtime)[2], (localtime)[1];
    print OUT $message,"\n";

    close OUT;

    return 1;

}

1;
