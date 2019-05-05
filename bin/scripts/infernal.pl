#!/usr/bin/env perl

use Cwd 'abs_path';


BEGIN{
    my $path = abs_path($0);
    $path =~ s/scripts.*/packages\//;
    push @INC, $path;
}

use strict;
use warnings;
use Getopt::Std;
use makeMuscleAlignment;
use CONFIG;
use tools;


my %opts;
getopts('e:d:i:m:r:w:g:t:k:n:shvb',\%opts);



####################################################################################################
## VARIABLES
####################################################################################################

my ( $master, $id, $copy, $dir, $total, $wild );
my ( $CMSEARCH, $CMCALIBRATE, $CMBUILD, $CMSTAT, $CMPATH, $STOCKHOLM, $FASTACMD );
my ( @data, @files );          #array to store the given multifastafiles (just in case wildcard expressions are used)


$copy = 0;
$total = 1;
$wild = "*.fa";


$CMSEARCH = $CONFIG::CMSEARCH;
$CMBUILD = $CONFIG::CMBUILD;
$CMCALIBRATE = $CONFIG::CMCALIBRATE;
$CMSTAT = $CONFIG::CMSTAT;

$STOCKHOLM = $CONFIG::STOCKHOLM;
$FASTACMD = $CONFIG::FASTACMD;

my $data_path = abs_path($0);
$data_path =~ s/bin\/.*/data\//;

if ( $opts{k} eq "deu" ){ $CMPATH = $data_path . $CONFIG::DEUTEROSTOMIA_CM_PATH; }
elsif( $opts{k} eq "pro" ){ $CMPATH = $data_path . $CONFIG::PROTOSTOMIA_CM_PATH; }
elsif( $opts{k} eq "fun" ){ $CMPATH = $data_path . $CONFIG::FUNGI_CM_PATH; }
elsif( $opts{k} eq "pla" ){ $CMPATH = $data_path . $CONFIG::PLANT_CM_PATH; }
else{ print STDERR "ERROR - No path to cm-models specified!"; exit(0); }



####################################################################################################
## START 
####################################################################################################

if( $opts{h} ){
    usage();
    exit(0);
}

if( defined $opts{g} && defined($opts{w}) && (defined($opts{d}) || (defined($opts{m})))){

    
    if( !defined($opts{e}) ){ $opts{e} = "1e-2"; }
    if( !defined($opts{i}) ){ $opts{i} = 1 };
    if( !$opts{t} ){ $opts{t} = "CD" }

    #open( LOG,">".$opts{w}."/log.run" );
    
    #path to multifasta file is defined by option m
    if( defined $opts{m} ){

	#check if multifasta file is located in current directory or not
	if( $opts{m} =~ /\// ){
	    my @temp = split( /\//, $opts{m} );
	    ( $dir = $opts{m} ) =~ s/$temp[-1]//g;
	}
	else{ $dir = "./"; }

	$total = `ls $dir/*.fa |wc -l`;
	push @files, $opts{m};
    }
    else{

	#directory and wildcard expression are defined by the options d and r, respectively
	$dir = $opts{d};
	$dir .= "/" if $dir !~ /\/$/;
	$wild = $opts{r} if $opts{r};
	@files = glob($dir."/".$wild);

    }

    if( !-e $dir ){
	# in case $dir doesn't exist, create it and change mod for other users
	`mkdir $dir`;
	`chmod 775 $dir`;
    }
    
    #read sequences from the given multifastafile(s)
    foreach( @files ){
	next if $_ =~ /rejected/;
    	$_ =~ /\/([^\/]+).fa$/;
    	push @data, ">".$1."_"; 
    }
    
    #open( STAT, ">".$opts{w}."/statistic.blast" );
    my $cnt = 0;
    my $hom = 0;

    for( my $i=$opts{i}; $i>0; $i-- ){
	my %queries;
	my $path = $opts{g}; #for the genomes in order of relativity
	$master = ();
	$cnt = 0;
	$hom = 0;
	$copy = 0;

	my $org = $opts{n};
	$org =~ /^(\w)\w* (\w+)/;
	my $org_out = uc($1).".".$2;

	my @genomes;
	$path =~ s/\n//g;
	my $test = `ls $path|wc -l`;
	( my $organism = $org ) =~ s/_/ /;
	chomp $organism;

	print STDERR "analyze $organism (infernal)";
	
	#more than one genome-file, specified with wildcard-expression
	if( $test > 1 || $opts{g} =~ /\*/){
	    @genomes =  glob( $opts{g} );
	}
	elsif( $test == 1 ){ #one genome-file is specified
	    push @genomes, $path;
	}
	else{
	    print STDERR " ... ... ERROR - genome file is unreadable!";
	    exit(0);	    
	}


        ## for all sequences in the fasta-file
	foreach my $dat ( @data ){ 

	    $dat =~ /^>*([^(_|\s)]+_[^(_|\s)]+)[_|\s]/;
	    $id = $1;
	    
	    
	    if( !-e  $dir."/".$id.".rejected.fa" ){
	    	`touch $dir"/"$id".rejected.fa"`;
	    }

	    if(!$opts{s}){
		if(`grep -i -c '$org_out' $dir$id".fa"` != 0 || `grep -i -c '$org_out' $dir$id".rejected.fa"` != 0){
		    print STDERR " ... ... ERROR - organism was already analyzed!";
		    next;
		}
	    }

	    
	    open(FILE, "<".$dir."/".$id.".fa") or die $!; #read fasta information
	    my @query = <FILE>;
	    close(FILE);
	
	    my $j = 0;
 	    my $jj = 0;
	    my $stop = "";
	    my @names;
	    my $gefunden = 0;
	    my ( $type, $species, $name, $query );
	    my $query_file = $opts{w} . "$$\_$id\_infernal.fa";
	    
	    open(OUT,">" . $query_file);

	    
	    ## write all sequences from the multifasta file in reversed order into a new file used for blast
	    ## rewrite the header
	    ## is this the main resaon for shifting the sequences into a new blast file?
	    ## store all provided information in a hash in order to eleminate blast hits which are already contained within the query file before returning blast results
	    for( my $i = (scalar @query) - 1; $i > 0; $i -= 2 ){
		$gefunden = 0;
		$dat =$query[$i-1];
		$dat =~ /^>([^_]+)_([^_]+)_([\w].[^_]+)_\(([^\)]+)\)_([^_]+)_([^_])/; #parse header of the file, not useful but let it as it came from the original script
		$type = $1;
		$name = $1 ."_". $2."_".$3;
		$species = $3;
		$stop = $species;
		$query = $query[$i];
		
		push @names, $name;
		if ( $org_out eq $3 ){
		    my %query;
		    $query{chr} = $4;
		    $query{position} = $5;
		    $query{strand} = $6;
		    
		    $queries{$1."_".$2} = \%query;
		}
		
		$query =~ s/\n//g;
		print OUT ">".$name."\n".$query."\n";
		
	    }
	    close OUT;
	    
	    if( (scalar @query) / 2 != scalar @names ){
		print STDERR " ... ... ERROR - not all queries could be used";
		exit(0);
	    }


	    ######
	    ## Search for an older calibrated cm-file that can be used
	    ## for an infernal search run.
	    ## This would save a huge amount of time, since calibrating
	    ## cm-file is awfully slow, but necessary.
	    ######

	    if( !&find_CM_file( $id, $query_file ) ){
		
		## in case only a single sequence is contained in the current snoRNA family.
		## this would cause no clustal-output, and hence no stockholm file and no cm model
		if ( `grep -c ">" $opts{w}"$$\_$id\_infernal.fa"` == 1 ){

		    ## snostrip will not search with cms built from single sequences
		    ## hence, print an error message and exit
		    print STDERR " ... ... ERROR - No CM file found or there is just one query sequence!";
		    exit(0);

      		}
		else{
		    print STDERR " ... ... cmbuild";
		    ## create clustal alignment from fastafile
		    &makeAlignment( $opts{w}.$$."\_".$id."\_infernal.fa", $opts{w}.$$."\_".$id."\_infernal.aln", 1 );
		    
		    ## generate a stockholm-formatted file that could be used as infernal input
		    `$STOCKHOLM $opts{w}"$$\_$id\_infernal.aln" $opts{w}"$$\_$id\_infernal.stk"`;
		}

		## build a covariance model from stk-file
		`$CMBUILD $opts{w}"$$\_$id\_infernal.cm" $opts{w}"$$\_$id\_infernal.stk"`;

		## calibrate cm-file
		## double the length of the random sequence that is used ( 1.5 --> 3.0 )
		## this will increase the runtime but refines the e-value calculation
#		print STDERR "... ... calibrating";
#		`$CMCALIBRATE -L 3.0 $opts{w}"$$\_$id\_infernal.cm"`;
	    }

	    my $blast = "";
	    my $tmpblast = "";
	    my ( $blasthit, $hit_length, $bitScore_threshold, $model_length, $target_length );
	    my ( @blasthit, @hits );


	    ## retrieve length of the model by the use of cmstat -m
	    ## this is necessary to distinguish between full-length hits and partial hits in the evaluation later on
	    $model_length = `$CMSTAT $opts{w}"$$\_$id\_infernal.cm" | grep -v "#" | awk '{print \$6}'`;
	    ## use query sequences to blast against all chromosomes specified in the @genomes-array, but skip random chromosomes
	    foreach my $genome ( @genomes ){
		    
		if( $genome =~ /random/ ){
		    next;
		}
		    
		#print STDERR "cmsearch --tabfile $opts{w}\"$$\_$id\_infernal.out\" $opts{w}\"$$\_$id\_infernal.cm\" $genome";

		`$CMSEARCH --tblout $opts{w}"$$\_$id\_infernal.out" $opts{w}"$$\_$id\_infernal.cm" $genome`;
		my @temp = `grep -v \"#\" $opts{w}"$$\_$id\_infernal.out" | grep "\!" | awk '{print \$3"\t"\$1"\t"\$8"\t"\$9"\t"\$6"\t"\$7"\t"\$15"\t"\$16"\t"\$13}' | sort -nrk 7`;
		
		`rm $opts{w}"$$\_$id\_infernal.out"`;

		foreach my $line ( @temp ){
	    
		    $line=~s/\n//;
		    next if $line eq "";

		    @blasthit = split( /\s+/, $line );


		    if ( $blasthit[7] > $opts{e} ){
			next;
		    }


		    ## the target length is described by the start and end positions of the target coords
		    ## the hit length is given by the query coords.
		    ## in case the hit length equals the model length, the whole model is mapped to the genome
		    ## the target length describes the length of the genomic region that was mapped.
		    ## a considerably longer or shorter targetregion indicates insertion and deletion events, respectively
		    $hit_length = $blasthit[5] - $blasthit[4] + 1;
		    if( $blasthit[3] > $blasthit[2] ){
			$target_length = $blasthit[3] - $blasthit[2] + 1;
		    }
		    else{
			$target_length = $blasthit[2] - $blasthit[3] + 1;
		    }

		    ## write a slightly different blastoutput, that contains the model_length instead of the percent identity
		    $blasthit = $blasthit[0] . "\t" . $blasthit[1] . "\t". $model_length. "\t" . $target_length ."\t0\t0\t" . $blasthit[4] . "\t" . $blasthit[5] . "\t" . $blasthit[2] . "\t" . $blasthit[3] . "\t" . $blasthit[7] . "\t" . $blasthit[6];

		    ## add formated infernal hit and the corresponding genome database to the infernal candidate set
		    $blast .= "$blasthit\t$genome\n";
		    push @hits, [split( /\s+/, $blasthit."\t".$genome )] 

		}

	    }


	    if( $blast eq "" ){
		
		## use two print outputs
		## print on STDERR is directly shown in the terminal
		print STDERR " ... ... ERROR - no significant hit found";
		

	    } 
	    else{

		## evaluate blast-like infernal output
		## do it seperatly for every sequence which retrieved a significant bitscore
	
		$j = 0; #anzahl der akkzeptierten Blasthits
		my ( $bestId, $bestLen, $bestScore );
	        HITS: while( scalar( @hits ) > $j ){

		    foreach my $tmp ( @{$hits[$j]} ){
			#print STDERR $tmp,"\n";
		    }

		    my $model_length = $hits[$j]->[2]; #baseidentity
		    my $perLen = $hits[$j]->[3] / ( $model_length / 100 ); #percentage of length of sequence
		    my $tmpScore = $hits[$j]->[11]; #score of this hit
		    my $evalue = $hits[$j]->[10];
		    my $ref;
		    
		    #save general blasthit-information
		    $ref->{database} = $hits[$j]->[-1];
		    $ref->{query_fam} = $id.".cm";
		    $ref->{found_with} = "cm_".$id;
		    $ref->{type} = $type;
		    $ref->{id} = $id;
		    $ref->{cnt} = $j;
		    $ref->{scores} = $model_length."-".sprintf ("%.2f", $perLen)."-".$tmpScore;
		    $ref->{chr} = $hits[$j]->[1];
			
		    ## save strand-specific blasthit-information
		    ## add a default of 10nt on both sides of the sequence
		    if( $hits[$j]->[8] > $hits[$j]->[9] ){
			$ref->{strand} = "-";
			$ref->{Sopt} = 2;
			$ref->{start} = $hits[$j]->[9] - 5 - ( $model_length - $hits[$j]->[7] );
			$ref->{start} = 0 if $ref->{start} < 0;
			$ref->{end} = $hits[$j]->[8] + 5 + $hits[$j]->[6];
		    }
		    else{
			$ref->{strand} = "+";
			$ref->{Sopt} = 1;
			$ref->{start} = $hits[$j]->[8] - 5 - $hits[$j]->[6];
			$ref->{end} = $hits[$j]->[9] + 5 + ( $model_length - $hits[$j]->[7] );
		    }
		    
		    if($j == 0){
			$bestLen = $perLen;
			$bestScore = $tmpScore;
		    }
		    
		    
		    #check if the current blasthit is already contained in the query file
		    my ($start, $end);
		    foreach my $blasthit ( keys %queries ){
			($start, $end) = split(/,/, $queries{$blasthit}->{position});
			
			if( $ref->{chr} eq $queries{$blasthit}->{chr} && 
			    $ref->{strand} eq $queries{$blasthit}->{strand} && 
			    ($ref->{start} <= $start && $ref->{end} >= $end || 
			     $ref->{start} >= $start && $ref->{end} <= $end ||
			     $ref->{start} <= $start && $ref->{end} >= $start ||
			     $ref->{start} >= $start && $ref->{start} <= $end) ){
#			    print STDERR "duplication found: $blasthit\n$start, $end\t\t$ref->{start}, $ref->{end}\n\n" if !$opts{v};
			    $j++;
			    next HITS;
			}
			
		    }
		    

		    #perfect blasthit
		    
		    if( $perLen >= 70 ){

			$ref->{range} = $ref->{start}.",".$ref->{end};
			( $ref->{seq}, $ref->{start}, $ref->{end} ) = cutSeq( $ref->{seq}, $ref->{chr},$ref->{range},$ref->{Sopt},$org,$hits[$j]->[-1] );
			my $fend = $ref->{end} + 200;
			$master->{$id}->{$jj} = $ref;
			$cnt++; #counter for all accepted blasthits
			$j++;  #counter for blasthit
			$jj++;
			$gefunden = 1;

		    }

		    else{
			$j++;
			print STDERR " ... ... ERROR - no significant hit found (too short)";
			next HITS;

		    }

		    
		}
		
		for(my $h=0; $h< $jj; $h++){
		    $hom++;
		    print ">".$id."_".$org_out."_(".$master->{$id}->{$h}->{chr}.")_".$master->{$id}->{$h}->{start}.",".$master->{$id}->{$h}->{end}."_".$master->{$id}->{$h}->{strand}."_".$master->{$id}->{$h}->{found_with}."_".$master->{$id}->{$h}->{scores}."\n".$master->{$id}->{$h}->{seq}."\n". $master->{$id}->{$h}->{database}."\n";
		}
	    }
	}

    }

    
    ## save cm files and the corresponding multifasta files in an extra subdirectory
    if( -e $opts{w}."$$\_$id\_infernal.stk" ){
	`rm $opts{w}"$$\_$id\_infernal.fa" $opts{w}"$$\_$id\_infernal.cm"`;
	`rm $opts{w}"$$\_$id\_infernal.stk"`;
	`rm $opts{w}"$$\_$id\_infernal.aln"` if -e $opts{w} . "$$\_$id\_infernal.aln";
	`rm $opts{w}"$$\_$id\_infernal.dnd" ` if -e $opts{w} . "$$\_$id\_infernal.dnd";
    }
    else{
	`rm $opts{w}"$$\_$id\_infernal.fa" $opts{w}"$$\_$id\_infernal.cm"` if -e $opts{w} . "$$\_$id\_infernal.fa";
    }

}
else{
    usage();
}


##########################################################################
sub usage
##########################################################################
{

    print STDERR "USAGE:\n";
    print STDERR "Given a file containing a list of ncRNA in fasta format or a directory with a fasta-files for each organism, this script\nwill search for homologs in user-specified organism\n";
    print STDERR "\t -h\tprints this help message\n";
    print STDERR "\t -d\tdirectory must exist and contain *.fa -files that are taken as queries\n";
    print STDERR "\t -r\tuse with -d and specify the multifasta to search for with wildcard expression, default is *.fa\n"; 
    print STDERR "\t -m\tone multifasta-file for which homolog sequences should be searched \n";
    print STDERR "\t -g\tpath that specifies the genome that should be used \n";
    print STDERR "\t -b\twhen searching for a potential cm file: each snoRNA of the current fasta file has to be included in the cm file and vice versa\n";
    print STDERR "\t -e\te-value for blast default 1e-3 \n";
    print STDERR "\t -i\thow many iterations of the blastsearch should be performed, default 1\n";
    print STDERR "\t -w\toutput directory\n";
    print STDERR "\t -v\tverbose-mode \n";

}



##########################################################################################
sub find_CM_file
## Search for an older but calibrated cm file that can be used for an infernal 
## search run. Therefore, the multifasta file that was used to construct the 
## covariance model hast to contain the same sequences which would be used now.
## If that's the case, the old model can be applied instead of generating a new one.
## This would save a huge amount of time, since the calibration of cm file is 
## awfully slow, but necessary for an e-value calculation.
##########################################################################################
{

    my ( $id, $found, $cmfile, $fastafile, $snoRNA, $query_file );
    my ( @cmfiles );

    $id = shift;
    $query_file = shift;

    ## Are there old cm files of the current snoRNA family?
    if ( `ls $CMPATH*$id*".cm" 2> /dev/null | wc -l` > 0 ){
	@cmfiles = split( /\s+/, `ls $CMPATH$id".cm" 2> /dev/null` );
    }

    ## check for each putative cm file if it is calibrated
    ## and if the corresponding multifasta file contains 
    ## the same snoRNA sequences
    ## Therefore, both directions have to be checked.
    FILES: while ( @cmfiles ){

	$cmfile = shift @cmfiles;
	chomp $cmfile;

#	next if !`grep calibrate $cmfile`;
	( $fastafile = $cmfile ) =~ s/\.cm/.fa/;
      ORGANISMS1: foreach $snoRNA ( split( /\s+/, `grep \">\" $CMPATH$id".fa"` ) ){
	  $snoRNA =~ s/_\(.*//;
	  next FILES if !`grep \"$snoRNA\" $query_file`;
	}

	if ( $opts{b}){
	  ORGANISMS2: foreach $snoRNA ( split( /\s+/, `grep \">\" $query_file` ) ){
	      next FILES if !`grep \"$snoRNA\" $CMPATH$id".fa"`;
	  }
	}
	
	$found = 1;
	last FILES;

    }

    ## if a suitable cm file is found, copy and rename this file in order to fit
    ## the current process id
    if ( $found ){
	`cp $cmfile $opts{w}"$$\_$id\_infernal.cm"`;
    }


    return $found;

}

##########################################################################################
sub createSingleSequenceStk
## This method creates an stockholm file conataining exactly one sequence. 
##########################################################################################
{
    my ( $id, $header, $sequence, $file );
    $id = shift;

    $sequence = `grep -v ">" $opts{w}"$$\_$id\_infernal.fa"`;
    chomp $sequence;
    $header = `grep ">" $opts{w}"$$\_$id\_infernal.fa" | sed 's/>//'`;
    chomp $header;
    $header =~ s/CD_[^_]+_|HACA_[^_]+_//;

    open( STK, ">" . $opts{w}. "$$\_$id\_infernal.stk") or die $!;
    print STK "# STOCKHOLM 1.0\n\n\n";
    printf STK "%-25s%s\n", $header, $sequence;
    printf STK "%-25s%s\n", "#=GC SS_cons", "." x length($sequence);
    print STK "//\n";
    close STK;

}


##########################################################################
sub cutSeq
#
##########################################################################
{

    my ($seq, $chr, $range, $Sopt,$genome,$bla) = @_;
    my $bool = 0;

    my $tmpSeq = `$FASTACMD -d $bla -s \"$chr\" -L $range -S $Sopt -l 500 2>&1`;


    if($tmpSeq =~ /ERROR/){ $bool = 1; $tmpSeq =~ s/^\[fastacmd\][^\n]+\n+//;}

    $tmpSeq =~ s/^>[^\n]+\n//;
    $tmpSeq =~ s/\s//g;
    
    my ($start, $end) = split(/,/ , $range);
    if(!$tmpSeq || $end - $start <= (length $tmpSeq) - 2 || $bool){
	$end = (length $tmpSeq) + 2;
	$seq = substr($tmpSeq, $start - 1) if $Sopt == 1;
	$seq = substr($tmpSeq, 0, $end - $start ) if $Sopt == 2;
	return ($seq, $start, $end);
    }

    return ($tmpSeq, $start, $end);

}




