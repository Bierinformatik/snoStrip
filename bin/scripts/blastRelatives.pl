#!/usr/bin/perl -w

use Cwd 'abs_path';

BEGIN{
    my $path = abs_path($0);
    $path =~ s/scripts.*/packages\//;
    push @INC, $path;
}

use strict;
use warnings;
use abbreviation;
use Getopt::Std;
use CONFIG;

my %opts;
getopts('f:e:d:i:m:r:w:g:t:k:shv',\%opts);


##################################################
## VARIABLES
##################################################

my ( $master, $message, $copy, $dir, $total, $wild, $BLAST, $FASTACMD, $RSCRIPT, $CLUSTER, $CD_CUTOFF, $HACA_CUTOFF, $CD_MEAN, $HACA_MEAN, $CD_SD, $HACA_SD );
my ( @data, @files );


$copy = 0;
$total = 1;                          #how many snoRNAs do I have in my dataset
$wild = "*.fa";                      #includes wildcard expression for multifasta-files in case directory is given
$BLAST = $CONFIG::BLASTALL;
$RSCRIPT = $CONFIG::RSCRIPT;
$CLUSTER = $CONFIG::CLUSTER;
$FASTACMD = $CONFIG::FASTACMD;

$CD_CUTOFF = $opts{k} eq "deu" ? $CONFIG::CD_CUTOFF_DEU : $CONFIG::CD_CUTOFF_FUNGI;
$HACA_CUTOFF = $opts{k} eq "deu" ? $CONFIG::HACA_CUTOFF_DEU : $CONFIG::HACA_CUTOFF_FUNGI;
$CD_MEAN = $CONFIG::CD_MEAN;
$CD_SD = $CONFIG::CD_SD;
$HACA_MEAN = $CONFIG::HACA_MEAN;
$HACA_SD = $CONFIG::HACA_SD;



##################################################
## START
##################################################

if( $opts{h} ){
    usage();
    exit(0);
}


if( defined $opts{g} && defined($opts{w}) && defined($opts{k}) && (defined($opts{d}) || (defined($opts{m})))){

    if( !defined($opts{e}) ){ $opts{e} = "1e-3"; }
    if( !defined($opts{i}) ){ $opts{i} = 1 };
    if( !$opts{t} ){ $opts{t} = "CD" }


    ## open( LOG,">".$opts{w}."/log.run" );
    
    ## path to multifasta file is defined by option m
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
    ## or by option d
    else{

	## directory and wildcard expression are defined by the options d and r, respectively
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
	my @temp = split( /\//, $path );
	

	my $tmppath = $path;
	$tmppath =~ s/\*/\\\*/;
	my $org = `grep "$tmppath" $opts{f} | awk '{print \$2"_"\$3}' `;
	chomp $org;


	$org =~ /^(\w)\w*_(\w+)/;
	my $org_out = uc($1).".".$2;
	my $end = $temp[-1];
	$end =~ s/\n//;
	my @genomes;
	$path =~ s/\n//g;
	my $test = `ls $path|wc -l`;
	( my $organism = $org ) =~ s/_/ /;
	chomp $organism;

	print STDERR "\nanalyze $organism (blast)";
	
	#more than one genome-file, specified with wildcard-expression
	if( $test > 1 ){
	    @genomes =  glob( $opts{g} );
	}
	elsif( $test == 1 ){ #one genome-file is specified
	    push @genomes,$path;
	}
	else{
	    print STDERR " ... ... ERROR - genome file is unreadable!";
	    exit(0);	    
	}

	foreach my $dat ( @data ){ #for all sequences in the fasta-file

	    $dat =~ /^>*([^(_|\s)]+_[^(_|\s)]+)[_|\s]/;
	    my $id = $1;
	    
	    if( !-e  $dir."/".$id.".rejected.fa" ){
	    	`touch $dir"/"$id".rejected.fa"`;
	    }

	    if(!$opts{s}){
		if(`grep -i -c '$org_out' $dir$id".fa"` != 0 || `grep -i -c '$org_out' $dir$id".rejected.fa"` != 0){
		    print STDERR " ... ... ERROR - organism was already analyzed!\n";
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
	    my %seq_length;
	    
	    open(OUT,">".$opts{w}."$$\_blast.fa");
	    
	    
	    ##================================================================================================
	    ## write all sequences from the multifasta file in reversed order into a new file used for blast
	    ## rewrite the header
	    ## 
	    ## store information of sequences that are alreadt contained in the current organism
	    ## in order to eleminate duplicated blast hits before returning the blast results to snostrip
	    ##================================================================================================
	    for( my $i = (scalar @query) - 1; $i > 0; $i -= 2 ){
		$gefunden = 0;
		$dat =$query[$i-1];
		$dat =~ /^>([^_]+)_([^_]+)_([\w].[^_]+)_\(([^\)]+)\)_([^_]+)_([^_])/; #parse header of the file, not useful but let it as it came from the original script
		$type = $1;
		$name = $1 ."_". $2."_".$3;
		$species = $3;
		$stop = $species;
		$query = $query[$i];
		$seq_length{$name} = length($query);

		push @names, $name;
		if($org_out eq $3){
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
		exit;
	    }
		

	    my $blast = "";
	    my $tmpblast = "";

	    #use query sequences to blast against all chromosomes specified in the @genomes-array, but skip random chromosomes
	    foreach my $genome ( @genomes ){
		    
		if( $genome =~ /random/ ){
		    next;
		}
		    
#		print STDERR  "\ncommand:\nblastall -p blastn -d $genome -i $opts{w}$$\_blast.fa -m 8 -e $opts{e} -b 100 -W 8 -r 1 -q -1 -G 2 -E 1\n";
		my @temp = `$BLAST -p blastn -d $genome -i $opts{w}$$\_blast.fa -m 8 -e 1e-1 -b 100 -W 8 -r 1 -q -1 -G 2 -E 1 -a 4`;

		foreach my $line ( @temp ){
		    $line=~s/\n//;
		    $line.="\t$genome\n";
		    $blast.=$line;
		}

	    }

	    if( $blast eq "" ){
		print STDERR " ... ... ERROR - no blasthits found";
		
	    } 
	    else{
		
		#check if there are overlaps in blastoutput
		#using christians cluster script
		#filter all blasthits with score lower 60
		
		my $tmpb;
		my @line;
		my @tmpb = split( /\n/, $blast );
		foreach my $line ( @tmpb ){
		    $tmpb .= $line."\n";
		}
		$blast = $tmpb;
		
		if(!$blast){
		    print STDERR " ... ... ERROR - no blasthits with score greater than 60 found"; 
		    next;
		}
		chomp $blast;
		
		open(TMP, ">".$opts{w}.$$."_blast.out") or die $!;
		print TMP $blast,"\n";
		close TMP;
		
		my $numberQueries = `cat $opts{w}$$\_blast.out | cut -f 1 | sort | uniq | wc -l`;
		if( !$opts{v} ){
#		    print STDERR "\n\nClustering ...\n";
#		    print STDERR "number of elements : ", scalar @names, "\n";			
#		    print STDERR "number of blasthits: ", scalar @{[split(/\n/, $blast)]},"\n";
#		    print STDERR "number of families : ", $numberQueries,"\n";
		}
		
		`$RSCRIPT $CLUSTER $opts{w}$$\_blast.out $opts{w}$$\_blast.chn 2> /dev/null`;

		if( -e $opts{w}.$$."_blast.chn" ){
		    $tmpblast = `cat $opts{w}$$\_blast.chn`;
		}
		    
		#evaluate blastoutput
		#do it seperatly for every sequence which retrieved a blasthit
		foreach ( @names ){

		    my @hits;

		    my $blast2 = `cat $opts{w}$$\_blast.chn | grep $_ | sort -rk 11,11g `;
		    $blast2 =~ s/^$//;
		
		    my @blasts = split( /\n/,$blast2 );
		    foreach my $bla ( @blasts ){
			
			next if( !($bla=~/\d/) );
			$bla=~s/\n+//g;
			chomp $bla;
			push @hits,[split(/\s+/, $bla)]

		    }
			
			
		    $j = 0; #anzahl der akkzeptierten Blasthits
		    my ( $bestId, $bestLen, $bestScore );
		    HITS: while( scalar(@hits) > $j ){

			my $tmp_query = $hits[$j]->[0];
			my $missing_upstream = $seq_length{$tmp_query} - $hits[$j]->[7];
			$tmp_query =~ s/_\(.*//;
			my $tmpId = $hits[$j]->[2]; #baseidentity
			my $length_ratio = $hits[$j]->[3] / $seq_length{$tmp_query};
			my $perLen = $hits[$j]->[3] / (length($query)/100); #percentage of length of sequence
			my $tmpScore = $hits[$j]->[11]; #score of this hit
			my $score_ratio = $hits[$j]->[11] / $seq_length{$tmp_query} ;
			my $ref;
			my @tmpFound = split( /_/, $hits[$j]->[0] );


			#save general blasthit-information
			$ref->{database} = $hits[$j]->[-1];
			$ref->{query_org} = $species;
			$ref->{found_with} = &get_abbreviation_Gepithet($tmpFound[-1], $opts{f} )."_".$tmpFound[0]."_".$tmpFound[1];
			if( $ref->{found_with} =~ /^_/){
			    print STDERR "WARNING: $tmpFound[-1] is not found in the information file!\n";
			}
			$ref->{type} = $type;
			$ref->{id} = $id;
			$ref->{cnt} = $j;
			$ref->{scores} = $tmpId . "-" . sprintf ("%.2f", $length_ratio) . "-" . sprintf( "%.2f", $score_ratio );
			$ref->{chr} = $hits[$j]->[1];
			
			#save strand-specific blasthit-information
			if( $hits[$j]->[8] > $hits[$j]->[9] ){
			    $ref->{strand} = "-";
			    $ref->{Sopt} = 2;
			    $ref->{start} = $hits[$j]->[9];
			    $ref->{end} = $hits[$j]->[8];
			}
			else{
			    $ref->{strand} = "+";
			    $ref->{Sopt} = 1;
			    $ref->{start} = $hits[$j]->[8];
			    $ref->{end} = $hits[$j]->[9];
			}
			
			if($j == 0){
			    $bestId = $tmpId;
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
				print STDERR " ... ... ERROR - no novel sequence found!";
				$message = "duplication found";
				$j++;
				next HITS;
			    }
			    
			}
			
			$ref->{score} = &is_above_cutoff( $tmpId, $length_ratio, $score_ratio );

			## blasthit with a sufficient cumulated score
			if( $ref->{score} ){
			    
			    ## enlarge BLASTHITs about 10nt on each side
			    ## beforehand: the query length was filled with missing nucleotides

#			    print STDERR "START: $ref->{start}\tEND: $ref->{end}\n";

			    if( $ref->{start} < $ref->{end} ){
 				$ref->{start} -= $hits[$j]->[6] ;
				$ref->{start} -= 10 if $opts{t} eq "CD";
				$ref->{start} = 0 if $ref->{start} < 0;
				$ref->{end} += $missing_upstream;
				$ref->{end} += 10 if $opts{t} eq "CD";
			    }
			    else{
				$ref->{start} += $missing_upstream;
				$ref->{start} += 10 if $opts{t} eq "CD";
				$ref->{end} -= $hits[$j]->[6];
				$ref->{end} -= 10 if $opts{t} eq "CD";				
				$ref->{end} = 0 if $ref->{end} < 0;
			    }

#			    print STDERR "START: $ref->{start}\tEND: $ref->{end}\n";
			    

			    $ref->{range} = $ref->{start}.",".$ref->{end};
			    ( $ref->{seq}, $ref->{start}, $ref->{end} ) = cutSeq( $ref->{seq}, $ref->{chr},$ref->{range},$ref->{Sopt},$org,$hits[$j]->[-1] );
			    $master->{$id}->{$jj} = $ref;
			    $cnt++; #counter for all accepted blasthits
			    $j++;  #counter for blasthit
			    $jj++;
			    $gefunden = 1;
			    $message = "found!";
			}
  

			else{
			    $message = "not found!";
			    $j++;
			}


		    }
		}

		if ($jj == 0){
		    print STDERR " ... ... ERROR - no blasthit with sufficient score found!\n";
		}

		

		for(my $h=0; $h< $jj; $h++){
		    $hom++;
		    print ">".$id."_".$org_out."_(".$master->{$id}->{$h}->{chr}.")_".$master->{$id}->{$h}->{start}.",".$master->{$id}->{$h}->{end}."_".$master->{$id}->{$h}->{strand}."_".$master->{$id}->{$h}->{found_with}."_".$master->{$id}->{$h}->{scores} . "-" . sprintf( "%.2f", $master->{$id}->{$h}->{score})."\n".$master->{$id}->{$h}->{seq}."\n". $master->{$id}->{$h}->{database}."\n";
		    
		}
		if ( $gefunden != 1 && $message eq "duplication found" ){
		    #print "skip due to a detected duplication";
		}
	    }
	}
	#print STAT $cnt." davon ".$hom." homologe und ".$copy." Kopien von insgesant ".$total." in: ".$org." gefunden\n";
    }
    #close(STAT);
    
    ## remove tmp files
    `rm $opts{w}$$\_blast.out $opts{w}$$\_blast.chn` if -e $opts{w}.$$."_blast.chn";
    `rm $opts{w}$$\_blast.fa` if -e $opts{w}.$$."_blast.fa";
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
    print STDERR "\t -e\te-value for blast default 1e-3 \n";
    print STDERR "\t -i\thow many iterations of the blastsearch should be performed, default 1\n";
    print STDERR "\t -w\toutput directory\n";
    print STDERR "\t -v\tverbose-mode \n";

}


##########################################################################
sub is_above_cutoff
## calculate score for blast hit
## and compare it to the CD or HACA cutoff-values
##########################################################################
{
    
    my $identity = shift;
    my $length = shift;
    my $bitscore = shift;


    my $score = ( $length*$length ) * $identity * $bitscore;
    $score = ( $score - $CD_MEAN ) / $CD_SD if $opts{t} eq "CD";
    $score = ( $score - $HACA_MEAN ) / $HACA_SD if $opts{t} eq "HACA";

    if ( $opts{t} eq "CD" ){
	return $score if $score > $CD_CUTOFF;
	return 0;
    }
    else{
	return $score if $score > $HACA_CUTOFF;
	return 0;
    }

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


