package Snoopy;


use strict;
use Bio::SimpleAlign;
use Bio::AlignIO;
use find_alignment_position;
use CONFIG;
use vars qw(@ISA @EXPORT);

require Exporter;


@ISA = qw(Exporter);
@EXPORT = qw(splitHP singlesnoop makePostscript);



####################################################################################################
## VARIABLES
####################################################################################################
my $RNASNOOP = $CONFIG::RNASNOOP;
my $SVM = $CONFIG::SVM;


####################################################################################################
sub splitHP
#splits HACA according to box-information
#this is needed, as RNAsnoop needs single hairpin as input
####################################################################################################
{

    #declare some variables and paramaters
    my ($seq, $startA, $startB) = @_;
    my ($HP1, $HP2);
   
    $HP1 = substr($seq,0,$startA-1);
    $HP2 = substr($seq,$startA+5,($startB-$startA-6));

    return ($HP1,$HP2);
    
}


####################################################################################################
sub singlesnoop
#single sequence target search
####################################################################################################
{
    my ( $HP, $res, $locations, $kingdom, $e ) = @_;
    my ( $snoop, $target) = '';
    my $message_target_prediction = "";

    
    if( !$locations->{targetRNApath} ){
	return( "", "path to targetRNAs is unspecified" ); 
    }


    #print STDOUT "Predict Targets...\n";
    open(SNOOP,">".$locations->{tmpPath}."$res->{org}\_snoop.fa");
    #print STDOUT "\n>$dir/",$$,"snoop.fa\n";
    print SNOOP ">",$res->{org},"_HP\n",$HP,"\n";
    close(SNOOP);

    #get abbreviation for current organism
    ( my $viech = $res->{org} ) =~ s/\./_/; 


    #collects all possible targetRNAs
    my @snRNAs = glob( $locations->{targetRNApath}.$viech."_*U*.fa" );
    my @rRNAs = glob( $locations->{targetRNApath}.$viech."_*[0-9]S*.fa" );

    #no target sequences are known, so target prediction is useless
    my $targets;
    if(scalar @snRNAs == 0 && scalar @rRNAs==0){
	#print "ERROR:\n NO SEQUENCES OF PUTATIVE TARGET RNAS AVAILABLE IN ",$viech,"\n";  
	$message_target_prediction = "no targetRNAs available";
	$target->{seq} = "NA";
	$target->{targetRNA} = "NA";
	$target->{pseudo} = "NA";
	$target->{tarSeq} = "NA";
	$target->{lhs} = "NA";
	$target->{rhs} = "NA";
	$target->{energy} = "NA";
	$target->{struc} = "NA";
	$target->{loop} = "NA";
	$target->{ending} = "NA";
	$target->{alnPos} = "NA";
	push(@{$targets},$target);
	return ( $targets, $message_target_prediction );
    }

    my $targetRNAs = join(' ', @snRNAs)." ".join(' ', @rRNAs);
    `cat $targetRNAs >$locations->{tmpPath}$res->{org}\_targets.fa` if $targetRNAs ne " "; 

    if( $targetRNAs eq " " ){
	$message_target_prediction = "unknown error";
	return  ( "", $message_target_prediction );
    }
    print STDERR "\n\n";
    print STDERR $locations->{profilesPath},"\n";

    if( defined($e) ){
	#svm-approach of RNAsnoop
	$snoop = `$SVM -t $locations->{tmpPath}$res->{org}\_targets.fa  -s $locations->{tmpPath}$res->{org}\_snoop.fa -U $locations->{profilesPath} -S u1_to_30.out -r yeast.svm.range -m yeast.svm.scale.12.model -o $locations->{tmpPath}/pred`;
    } else {
	#unfiltered RNAsnoop
	$snoop = `$RNASNOOP -t $locations->{tmpPath}$res->{org}\_targets.fa -s $locations->{tmpPath}$res->{org}\_snoop.fa -U $locations->{profilesPath} -S u1_to_30.out -i 120 -j 11 -k 16 -e 40  -o -270 -p -170 -q -1090 -x -1370 -L 25 -f 1 | grep -v "no target found" `;
	#$snoop = `$RNASNOOP -t $locations->{tmpPath}$res->{org}\_targets.fa -s $locations->{tmpPath}$res->{org}\_snoop.fa -i 120 -j 11 -k 16 -e 40  -o -270 -p -170 -q -1090 -x -1370 -L 25 -f 1 | grep -v "no target found" `;

    }
    #print STDOUT "QUE:\n",$snoop,"\n";

    #my $targets;
    my $hlp = index($snoop,"&");
    my $targetRNA_count = 0;
    #print STDERR "INDEX: ",$hlp,"\n";
    if($hlp == -1){
	$target->{seq} = "NA";
	$target->{targetRNA} = "NA";
	$target->{pseudo} = "NA";
	$target->{tarSeq} = "NA";
	$target->{lhs} = "NA";
	$target->{rhs} = "NA";
	$target->{energy} = "NA";
	$target->{struc} = "NA";
	$target->{loop} = "NA";
	$target->{ending} = "NA";
	$target->{alnPos} = "NA";
	push(@{$targets},$target);
    } else {
	my $best = parseSnoop($snoop,$e);
	#when there are predictions which have no lower stem it might happen that nothing is returned
	if(!defined($best) || scalar(@{$best}) < 1  ){
	    $target->{seq} = "NA";
	    $target->{targetRNA} = "NA";
	    $target->{pseudo} = "NA";
	    $target->{tarSeq} = "NA";
	    $target->{lhs} = "NA";
	    $target->{rhs} = "NA";
	    $target->{energy} = "NA";
	    $target->{struc} = "NA";
	    $target->{loop} = "NA";
	    $target->{ending} = "NA";
	    $target->{alnPos} = "NA";
	    push(@{$targets},$target);
	}
	else{
	    #print STDOUT "BEST: ",scalar(@{$best}),"\n";
	    foreach(@{$best}){
		$target = $_;
		## find position in the alignment
		$target->{alnPos} = &map_position( $target->{targetRNAFile}, $target->{targetRNA}, $target->{pseudo}, $viech, $kingdom );
		if( $target->{alnPos} =~ /-NA/ ){ 
		   $targetRNA_count++; 
		}
		push(@{$targets},$target);
		# print $target->{targetRNAFile},"\n";
	    }
	}
    }

    $message_target_prediction = "at least one targetRNA is missing or incomplete" if $targetRNA_count != 0;

    ## cleaning
    `rm $locations->{tmpPath}$res->{org}\_targets.fa $locations->{tmpPath}$res->{org}\_snoop.fa`;

    if ( !$targets ) { return ( "", $message_target_prediction ); }
    return ( $targets, $message_target_prediction ) ;

}


####################################################################################################
sub makePostscript
#prints postscript-file for given snoop-result
####################################################################################################
{
    my ($snoop, $locations) = @_;

    #`echo '$snoop' | $RNASNOOP -I -U $locations->{profilesPath} -S u1_to_30.out -O $locations->{postscriptPath}`;

}


####################################################################################################
sub parseSnoop
#parses RNASnoop-output and returns best found target
####################################################################################################
{
    my ($snoop,$e) = @_;
    #die snoop-ausgabe zeilenweise im Array
    my @snoop = split(/\n/,$snoop);  
    my @cand;
    my ($tarRNAFile,$tarRNA);
    my ($svm_score,$svm_head);

    #durchlaufen der hits die snoop zurueckgibt  
    for(my $i = 0; $i < scalar(@snoop); $i += 2){		
	my $hit;
	if($snoop[$i] =~ /^>/ && $snoop[$i+1] =~ /^>/){
	    #header information
	    #svm-evaluation or snoRNA-header
	    ($svm_head = $snoop[$i]) =~ s/>//; 
	    if(defined($e)){
		$svm_score = (split(/\./,$snoop[$i]))[4];
		$svm_score =~ s/[\D]+//;
		$svm_score = "0.".$svm_score;
	    }
	    #targetRNA-header
	    ($tarRNAFile = $snoop[$i+1]) =~ s/>//;
            
            #rRNA-name
	    if($tarRNAFile =~ /\dS.\d/){
		( $tarRNA = $tarRNAFile ) =~ s/.*_([0-9\.]+S.[0-9]+).*/$1/;
	    }
	    elsif($tarRNAFile =~ /\dS$/){
		( $tarRNA = $tarRNAFile ) =~ s/.*_([0-9\.]+S).*/$1/;
	    }
	    #snRNA-name
	    else{
		if($tarRNAFile =~ /\./){
		    $tarRNA = (split(/_/,(split(/\./,$tarRNAFile))[0]))[2].".".(split(/\./,$tarRNAFile))[1];
		}
		else{
		    $tarRNA = (split(/_/,$tarRNAFile))[2];
		}
	    }
	}
	else {
	    if($snoop[$i] =~ /^no/){
		#unsuccsessful target-search
		$i--;
		next;
	    }
	    $hit->{seq} = $snoop[$i+1];
	    $hit->{res} = $snoop[$i];

	    #Ausgabe von snoop
	    $hit->{res} =~ /^([^&]+)\&([^\d]+)(\d+),(\d+)\s*;\s*(\d+)\s*:\s*(\d+)\s*,\s*(\d+)\s*\(([^=]+)\s*=/;  
	    $hit->{tarStruc} = $1;
	    $hit->{stemStruc} = $2;
	    $hit->{start} = $3;
	    $hit->{end} = $4;
	    $hit->{pseudo} = $5;
	    $hit->{leftStart} = $6;
	    $hit->{rightStart} = $7;
	    $hit->{energy} = $8;
	    
	    #filter hits where no lower stem exists
	    #if($hit->{stemStruc} =~ /^\.*>+/){
	    if($hit->{stemStruc} =~ /^\.*>+/ || $hit->{stemStruc} =~ /^\.*\([\.+>+]/){
		#print STDERR $hit->{stemStruc},"-->NO LOWER STEM\n";
		next;
	    }
	    #else{
	    #	print STDERR $hit->{stemStruc},"-->O.K.\n";
	    #   }
	    
	    $hit->{struc} = $hit->{tarStruc}."&".$hit->{stemStruc};
	    $hit->{targetRNAFile} = $tarRNAFile;
	    $hit->{targetRNA} = $tarRNA;
	    if(defined($e)){
		$hit->{svm_score} = $svm_score;
	    }
	    $hit->{snoop} = ">".$svm_head."\n>".$hit->{targetRNAFile}."\n".$hit->{res}."\n".$hit->{seq};
	    $hit->{ps_file} = "sno_".$svm_head."_".$hit->{targetRNAFile}."_u_".$hit->{pseudo}."_".($hit->{pseudo} - $hit->{start}).".ps";
	    
	    push(@cand,$hit);
	}
    } 

    my $best;
    if(scalar(@cand)>0){
	my @sorted;
	#if(!(defined($e))){
        #nach svm_score absteigend geordnet
	#@sorted = sort {$b->{svm_score} <=> $a->{svm_score}} @cand;
	#} else {
	#nach Energien aufsteigend geordnet 
	@sorted = sort {$a->{energy} <=> $b->{energy}} @cand;
	#} 
	
	#make modificated sites unique
	my @uniq;
	my %seen = ();
	foreach (@sorted) {
	    my $tmp = $_;
	    #there are some variants for the snRNAs, but only the best for each snRNA should be regarded
	    (my $tmptar = $tmp->{targetRNA}) =~ s/\.\d+$//;
	    my $item = $tmptar."-".$tmp->{pseudo};
	    
	    push( @uniq, $tmp ) unless $seen{$item}++;
	}

	#returns the best 10 targets from the svm 
	for(my $i = 0; ($i<scalar(@uniq) && $i <= 10 ); $i++){ 
	    $best->[$i] = &getBindingProperties($uniq[$i]);
	}
    }
    
    return $best;
}


####################################################################################################
sub getBindingProperties
#analyses snoRNA-targetRNA binding in more detail
#extract the properties listed in snoBoard.HACA-target
####################################################################################################
{
    my $target = shift;
    my ( @tmp, $relPseudo );

    #get bound target sequence
    @tmp = split(/&/,$target->{seq});
    $relPseudo = index($target->{tarStruc},"|");
    substr( $tmp[0], $relPseudo, 1, "Y" );

    $target->{tarSeq} = $tmp[0];
    $target->{rhs} = substr($tmp[1],$target->{rightStart},($target->{pseudo}-$target->{start}));
    $target->{lhs} = substr($tmp[1],$target->{leftStart},($target->{end}-$target->{pseudo}));
    $target->{loop} = $target->{rightStart} - $target->{leftStart} + length($target->{lhs}) ;
    $target->{ending} = length($target->{seq}) - $target->{leftStart} + length($target->{rhs}) ;

    return $target;
}



1;
