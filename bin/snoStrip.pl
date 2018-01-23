#!/usr/bin/perl -w


use Cwd 'abs_path';

BEGIN{
    my $path = abs_path($0);
    $path =~ s/snoStrip.pl/packages\//;
    push @INC, $path;
}

use warnings;
use strict;
use Getopt::Long;
use Fold;
use Snoopy;
use Plexy;
use primeBoxes;
use File::is;
use searchForTargetSites;
use makeMuscleAlignment;
use boxPositions;
use boxing;
use CONFIG;
use Logger;
use abbreviation;




####################################################################################################
## VARIABLES
####################################################################################################

my ( $opts_h, $opts_d, $opts_i, $opts_k, $opts_j, $opts_v, $opts_f, $opts_t, $opts_a, $opts_p, $opts_s, $opts_g, $opts_n ) = 
   ( "", "", "", "", "", "", "", "", "", "", "", "", "" ); 
my ( @opts_c, @opts_n);

my ( $tmp, $boxes, $message, $passed, $hasPrimeBox, $usage );
my ( $BLAST, $FASTACMD, $MUSCLE, $STOCKHOLM, $INFERNAL, $TARGET_PATH, $TARGET_ALN_PATH, $TARGET_PROFILE_PATH, $FASTA_PATH, $INFO_FILE );
my ( @multifasta, @tmp, @lst );
my ( %seqProperties, %locations, %genome_information, %bigHash, %opts );


$hasPrimeBox = 0;

&setOptions();

&setConstants();

( $passed, $message ) = &checkOptions();




####################################################################################################
## START
####################################################################################################

if( $passed ){ #search for homologs in this genome (based on snoBoard-data)


    ######
    ##
    ## CREATE OUTPUT FOLDER STRUCTURE 
    ##
    ## therein, separate folder for multifasta files, gff files, alignments, and target prediction files has to be created
    ## additionally, a tmp-directory is needed as well
    ######


    if( !-e $opts_d ){
	`mkdir $opts_d`;
    }
    if( !-e $opts_d. $$ . "_tmp" ){
	`mkdir $opts_d$$"_tmp"`;
    }
    if( !-e $opts_d."mfasta" ){
	`mkdir $opts_d"mfasta"`;
    }
    if( !-e $opts_d."alignments" ){
	`mkdir $opts_d"alignments"`;
    }
    if( !-e $opts_d."gff" ){
	`mkdir $opts_d"gff"`;
    }
    if( $TARGET_PATH && !-e $opts_d."targets" ){
	`mkdir $opts_d"targets"`;
    }
    if( $TARGET_PATH && !-e $opts_d."targets/psfiles/"){
	`mkdir $opts_d"targets/psfiles/"`;
    }
    if( !-e $opts_d."statistics/"){
	`mkdir $opts_d"statistics/"`;
    }

    
    $locations{workingDirectory} = $opts_d;
    $locations{tmpPath} = $opts_d . $$ . "_tmp/";
    $locations{mfastaPath} = $opts_d."mfasta/";
    $locations{gffPath} = $opts_d."gff/";
    $locations{alignmentPath} = $opts_d."alignments/";
    $locations{targetPath} = $opts_d."targets/";
    $locations{postscriptPath} = $opts_d."targets/psfiles/";
    $locations{statisticsPath} = $opts_d."statistics/";
    
    $locations{targetRNApath} = $TARGET_PATH;
    $locations{targetAlignmentPath} = $TARGET_ALN_PATH;
    $locations{profilesPath} = $TARGET_PROFILE_PATH;
    $locations{information_file} = $opts_d . $$ . "_tmp/" . $$ . "\_" . $opts_k . ".csv";




    ######
    ##
    ## EXTRACT GENOME AND SEQUENCE INFORMATION
    ##
    ######

    if ( $opts_i){
      open( INF, "<" . $opts_i ) or die $!;
      while ( <INF> ){
	chomp $_;
	my @line = split( /\t/, $_ );
	my %genome_name;
	$genome_information{$line[1]} = \%genome_name;
	$genome_name{genomefile} = $line[0];
	$genome_name{organism} = $line[1];
	$genome_name{abbreviation} = $line[2];
      }
      close INF;
    }
    else{

      my %genome_name;
      $genome_information{ $opts_n} = \%genome_name;
      $genome_name{genomefile} = $opts_g;
      $genome_name{organism} = $opts_n;
      my ($genus, $epithet) = split( / /, $opts_n);
      $genome_name{abbreviation} = lc( substr( $genus, 0, 1)) . substr( $epithet, 0, 2);

    }

    
    ## combine the original information file with the user given file of organism information that should be analyzed with snostrip
    `cp $INFO_FILE $opts_d$$"_tmp/"$$"_"$opts_k".csv"`;
    if ($opts_i){
      `cat $opts_i >> $opts_d$$"_tmp/"$$"_"$opts_k".csv"`;
    }
    else{
	`echo -e $genome_information{ $opts_n}->{genomefile}"\t"$opts_n"\t"$genome_information{ $opts_n}->{abbreviation} >> $opts_d$$"_tmp/"$$"\_"$opts_k".csv"`;
    }

 
    ######
    ##
    ## EXTRACT THE snoRNA FAMILIES WHICH SHOULD BE ANALYSED
    ##
    ## search through the given options for snoRNA families and extract the mfasta locations appropriately
    ## that means, look them up in snoBoard in case of option 'l', 'f', or no specified option
    ## or use the path that is specified with option 'c'
    ######

    my $lst;
    my $no_such_file = 0;
    ## with the use of option 'c', the user specifies snoRNA families that should be analysed
    ## this option requires the multi fasta file location as argument, hence there is no snoBoard interaction at this point
    if( @opts_c ){

	foreach my $opts ( @opts_c){
	    my @tmp =  glob( $FASTA_PATH . $opts. ".fa");
	    foreach my $fam ( @tmp ){
		if( ! -e $fam ){
		    print STDERR "\nERROR:\n" , $fam, ": No such file or directory\n\n";
		    $no_such_file = 1;
		}
		else{
		    push @lst, [$fam];
		}
	    }
	}
    }

    ## in case no family is specified with option -c
    ## use all available snoRNA families
    else{
	my @tmp =  glob( $FASTA_PATH . "*.fa");
	foreach my $fam ( @tmp ){
	    push @lst, [$fam];
	}
    }


    ## stop application in case at least one file could not be found
    if ( $no_such_file){
	print STDERR "snoStrip stopped due to missing fastafiles!\n\n";
	exit(0);
    }


    ## when no multifasta location could be retrieved, due to some unpredictable errors, or user inability
    if( scalar( @lst ) < 1 ){ print STDERR "ERROR:\nDid not find multifasta!"; exit(0) }


    ######
    ##
    ## ANALYZE snoRNA FAMILIES
    ##
    ## analyze each snoRNA family that is specified in the array '@lst'
    ######

    foreach ( @lst ){

	%bigHash = ();
	$locations{mf} = $_->[0];
	
	( my $fam = $_->[0] ) =~ s/^(.*\/)*([^\/\.]+)\.fa$/$2/;

	print STDERR "current snoRNA family: $fam\n\n";
	&writeLog("starting analysis of Family $fam", \%locations);

	## open statistics file
	open ( my $stat , ">" . $locations{statisticsPath} . $fam . ".stat" ) or die $!; 
	open ( my $stat_rej , ">" . $locations{statisticsPath} . $fam . ".rejected.stat" ) or die $!; 
	print $stat "#C score\tC penalty\tD score\tD penalty\tD target\tC'score\tC'penalty\tD' score\tD' penalty\tD' target\tdist(C-D)\tdist(C-D')\tdist(C'-D)\tKturn\tstem\tsnoRNA\n" if $fam =~ /CD/;
	print $stat "#H score\tH penalty\tACA score\tACA penalty\tmfe HP1\tlen HP1\tmfe HP2\tlen HP2\tsnoRNA\n" if $fam =~ /HACA/;
	print $stat_rej "#C score\tC penalty\tD score\tD penalty\tD target\tC'score\tC'penalty\tD' score\tD' penalty\tD'target\tdist(C-D)\tdist(C-D')\tdist(C'-D)\tKturn\tstem\tsnoRNA\n" if $fam =~ /CD/;
	print $stat_rej "#H score\tH penalty\tACA score\tACA penalty\tmfe HP1\tlen HP1\tmfe HP2\tlen HP2\tsnoRNA\n" if $fam =~ /HACA/;

	my $res;  

	if( $locations{mf} ne '' ){
	   
	    ## extract path and filename of the current multi fasta file
	    $locations{mf} =~ /^\.*([\w\/]*?)([\w\d\-]+.fa)$/;
	    $locations{path} = $1;
	    $locations{dat} = $2;

	    ## extract snoRNA family
	    ( $locations{id} = $locations{dat} ) =~ s/\.fa$//;

	    $locations{rej} = $locations{id} . ".rejected.fa";


	    $locations{mfasta_output} = $locations{mfastaPath} . $locations{dat};
	    $locations{mfasta_rejected_output} = $locations{mfastaPath} . $locations{rej};

	    ## copy multifasta to output directory and delete a former file containing rejected snoRNAs
	    `cp $locations{mf} $locations{mfasta_output}`;
	    if( -e $locations{mfasta_rejected_output} ){
		`rm $locations{mfasta_rejected_output}`;
	    }

	    ## delete old hairpin files
	    if(-e $locations{targetPath}.$locations{id}."_HP1.fa" ){
		`rm $locations{targetPath}$locations{id}"_HP1.fa"`;
		`rm $locations{targetPath}$locations{id}"_HP2.fa"`;
	    }


	    ## declare some more variables
	    my ($bla, $bool, $genome, $count, $found, $foundPrime, $opts_b_blast );
	    my (@hits, @gen);
	    my ($head, $seq);

	    
	    %seqProperties = %{ &readMultiFasta( \%locations ) };
	    

	    my @genomes;
	    if( $opts_i){
		open (INFO, $opts_i) or die print STDERR "ERROR:\n File defined by option -i does not exist!\n";
		foreach my $line ( <INFO>){
		    chomp $line;
		    push @genomes, (split /\t/, $line)[1];
		}
		close INFO;
	    }
	    else{
		push @genomes, $opts_n;
	    }

	    #this foreach loop goes through the given genome file in order to treat every new genome seperatly
	    #this means that each genome in this file is analysed on its own using the current multifasta file
	  BLASTME:
	    foreach $genome ( @genomes ){

		&writeLog( "analyze $fam in " . $genome, \%locations);
		
		if( $genome =~ /^#/ || $genome =~ /^\s*$/ ){ next; }
		chomp $genome;

		$count = "";
		$opts_b_blast = 1;
		$bla = '';


		## blast the current organism, although we already have organism specific snoRNAs.
		if( $opts_f ){

		    #count all paralogs of the current organism in order to make each header unique
		    ( my $tmpviech = $genome_information{$genome}->{organism} ) =~ s/^([A-Z])[a-z]*\s/$1./;
		    $res->{paralogs} = `grep -c $tmpviech $locations{mf}`;
		    chomp $res->{paralogs};

		    $count = $res->{paralogs} + 1 ;

		    $bla = `$BLAST -m $locations{mfasta_output} -g "$genome_information{$genome}->{genomefile}" -w $opts_d/$$\_tmp/ -t CD -f $locations{information_file} -s -v -k $opts_k`;
		    $bla = `$BLAST -m $locations{mfasta_output} -g "$genome_information{$genome}->{genomefile}" -w $opts_d/$$\_tmp/ -t HACA -f $locations{information_file} -s -v -k $opts_k` if $fam =~ /HACA/;
		}

		## in case the was no multifasta input for the current organism,
		## use different homology-based search methods
		## START WITH BLAST
		if ( !$bla ){

		    $bla = `$BLAST -m $locations{mfasta_output} -g "$genome_information{$genome}->{genomefile}" -w $opts_d/$$\_tmp/ -t CD -f $locations{information_file} -v -k $opts_k` if $fam =~ /CD/;
		    $bla = `$BLAST -m $locations{mfasta_output} -g "$genome_information{$genome}->{genomefile}" -w $opts_d/$$\_tmp/ -t HACA -f $locations{information_file} -v -k $opts_k` if $fam =~ /HACA/;

		}

		######
		## USE INFERNAL IF BOTH PREVIOUS METHODS FAILED TO DETECT HOMOLOGS
		## AND THE OPTION 'j' IS SPECIFIED
		######

		if ( $bla eq "" && $opts_j ){
		    print STDERR "\n";
		    $bla = `$INFERNAL -m $locations{mfasta_output} -g "$genome_information{$genome}->{genomefile}" -w $opts_d/$$\_tmp/ -n "$genome_information{$genome}->{organism}" -k $opts_k -s -v` if $opts_f;
		    $bla = `$INFERNAL -m $locations{mfasta_output} -g "$genome_information{$genome}->{genomefile}" -w $opts_d/$$\_tmp/ -n "$genome_information{$genome}->{organism}" -k $opts_k -v` if !$opts_f;

		}



#		print STDERR "\nBLASTOUT:\n", $bla, "\n\n";
		
		## if the genomefile wasn't found an error-message will be returned
		## in such cases go on with next given genome
		if( $bla =~ /^ERROR/ || $bla eq "" || $bla =~ /skip/ ){
		    print STDERR "$bla\n" if !$opts_v;

		    &writeLog( "No candidate sequences found for " . $fam . " in " . $genome, \%locations );

		    next BLASTME;
		}


		@hits = split( />/,$bla );
		if (!$hits[0] ){ shift @hits; }

		print STDERR  " ... ... ", scalar @hits, " sequence(s) found\n\n" if !$opts_v;
		&writeLog( scalar @hits . " candidate sequence(s) found for " . $fam . " in " . $genome, \%locations );

		$bool = 0;     #boolean for sequence accepted
		$count = 1;
		
		#this foreach loop goes through all blasthits retrieved before
		my $seq_count = 1;
		my $correct_snoRNAs = 1;
	      HITS: foreach( @hits ){
		    
			## save information about snoRNA specific propeties
			## that were calculated along the way
		    my %statistics_hash;

		    my $res;
		    my $hit = $_;
		    $boxes = "";
		    $found = 0;
		    $foundPrime = 0;
		    
		    if( $hit ne '' ){
			@tmp = split( /\n/,$hit );
			$res->{bla} = $tmp[0];
			$res->{genome} = $genome_information{$genome}->{genomefile};
			$res->{database} = $tmp[2];
			$res->{path} = $locations{path};
			$res->{seq} = $tmp[1];
			
			$res->{abbr} = $genome_information{$genome}->{abbreviation};
			$res->{organism} = $genome_information{$genome}->{organism};
			$res->{source} = $genome_information{$genome}->{source};
			$res->{db} = $genome_information{$genome}->{db};


			if( $res->{bla} =~ /__/){
			    print STDERR "\n\nERROR: One of the query organisms is not specified in the information file!\n";
			    exit(0)
			}
			
			if( $opts_b_blast ){
			    $res->{bla} =~ /^([^_]+)_([^_]+)_([\w].[^_]+)_\(([^\)]+)\)_([^_]+)_([^_])_([^_]+)_([^_]+)_([^_]+)_([\w].[^_]+)/;
			    $res->{name} = $1."_".$2;
			    ( $res->{family} = $1 . "_" . $2 ) =~ s/\-.*//;
			    $res->{org} = $3;
			    $res->{type} = $1;
			    $res->{chr} = $4;
			    $res->{range} = $5;
			    $res->{strand} = $6;
			    $res->{found_with} = $7."_".$8."_".$9;
			    $res->{found_with_seq} = $8 . "_" . $9 . "_" . &get_scientific_name( $7, $locations{information_file} );
			    $res->{score} = (split (/-/, $10))[-1];

			    $res->{count} = $count;
			    $res->{name} = $res->{name}."-".$res->{count};
			    $res->{bla} = $res->{name}."_".$res->{org}."_(".$res->{chr}.")_".$res->{range}."_".$res->{strand}."_".$res->{found_with}."_".$10;
			    
			}
			else{
			    $res->{bla} =~ /^([^_]+)_([^_]+)_([\w].[^_]+)_\(([^\)]+)\)_([^_]+)_([^_])_/;
			    $res->{name} = $1."_".$2;
			    $res->{org} = $3;
			    $res->{type} = $1;
			    $res->{chr} = $4;
			    $res->{range} = $5;
			    $res->{strand} = $6;
			    $res->{count} = $count;
			    if( $res->{name} !~ /\-\d+$/ ){
				$res->{name} = $res->{name}."-".$res->{count};
				$res->{bla} =~ s/$1_$2/$res->{name}/;
				print STDERR $res->{bla}
			    }

			    ## extract box motifs
			    my @tmp = split( /_/, $res->{bla} );
			    $res->{boxA} = $tmp[-4];
			    $res->{startA} = $tmp[-3];
			    $res->{boxB} = $tmp[-2];
			    $res->{startB} = $tmp[-1];
			    if( $tmp[-5] =~ /\d+$/ && $tmp[-6] =~ /[GTACU]{4}/){
				$res->{boxAPrime} = $tmp[-8];
				$res->{startAPrime} = $tmp[-7];
				$res->{boxBPrime} = $tmp[-6];
				$res->{startBPrime} = $tmp[-5];
			    }

			    if( $res->{genome} =~ /\*/ ){

				my @genome_file = split( /\n/, `grep -P \"$res->{chr}\[\\s\_\-\]|$res->{chr}\$\" $res->{genome}` );
				if( scalar @genome_file == 1 ){
				    $genome_file[0] =~ s/:.*$//;
				    $res->{database} = $genome_file[0];
				}
				else{
				    print STDERR "ERROR - more than one suitable genome files:\n";
				    foreach my $file ( @genome_file ){
					print STDERR $file,"\n";
				    }
				}

			    }
			    else{
				$res->{database} = $res->{genome};
			    }

			}

			print STDERR  "\tsequence: ", $seq_count++ , "\n" if !$opts_v;
			$correct_snoRNAs++;
			$statistics_hash{snoRNA} = $res->{name} . "_" . $res->{org};
			
			######
			##
			## CHECK FOR AT LEAST ONE CONSERVED TARGET SITE
			##  
			## run a PWM approach to search for for highly conserved regions in the novel snoRNA candidate
			## sequences with a a score below 12 will be rejected
			## 
			######

			## skip this analysis if there is only one sequence in the multi fasta file
			if ( $res->{type} eq "CD" && $opts_b_blast ){

			    print STDERR "\t\tcheck for conserved target sites ", "." x 6, " " if !$opts_v;
			    
			    ## copy the current mfasta file to the tmp-directory
			    `cp $locations{mfasta_output} $locations{tmpPath}`;

			    my ( $bool, $message ) = (1, "");
			    ## create the muscle Alignment
			    if ( `grep -c ">" $locations{mfasta_output} ` >  1 ){
				( $bool, $message ) = makeAlignment( $locations{tmpPath}.$locations{id}.".fa", $locations{tmpPath}.$locations{id}.".aln" );
			    }
			    else{
				( $bool, $message ) = makeSingleSeqAlignment( $locations{tmpPath}.$locations{id}.".fa", $locations{tmpPath}.$locations{id}.".aln" );
			    }

			    ## check for errors
			    if ( $bool == 0 ){
				print STDERR "\nSome ERROR occured during alignment computation!\n\n";
			    }
			    
			    ## create a new alignment object using bio::AlignIO
			    ## to access the alignment
			    my $str = Bio::AlignIO->new(-file => $locations{tmpPath}.$locations{id}.".aln", -format => "clustalw");
			    my $aln = $str->next_aln();

			    ## get the consensus start position of box C and D
			    my ( $start_a ) = getBoxStart(  \%seqProperties, \$aln, "1", "1" );
			    my ( $start_b ) = getBoxStart(  \%seqProperties, \$aln, "2", "1" );

			    my ( $min_start_prime_b, $max_start_prime_b ) = ( 0, 0 );
			    my ( $start_prime_b ) = 0;

			    if( $hasPrimeBox ) {
#				( $min_start_prime_b, $max_start_prime_b ) = getBoxStart(  \%seqProperties, \$aln, "prime2" );
				$start_prime_b = getBoxStart(  \%seqProperties, \$aln, "prime2", "1" );
			    } 

#			    print STDERR "start_primeB: ", $start_prime_b , "\n";

			    ## search for conserved site of length 13nt
			    my ( $has_target_site, $best_score, $target_site, $target_box );
			    ( $has_target_site, $best_score, $target_site, $statistics_hash{scoreDTarget}, $statistics_hash{scoreDprimeTarget} ) = searchTargetSites( \$aln, $res->{seq}, $start_a, $start_b, $start_prime_b, \%locations );

			    if ( !$has_target_site ){
				my @best_score = split( /\s/, $best_score );
				if ( $best_score[0] == $statistics_hash{scoreDTarget} ) {
				    $target_box = "Dbox";
				    $best_score = $best_score[0];
				}
				else{
				    $target_box = "D'box";
				    $best_score = $best_score[1];
				}
			    }
			    else{
				$target_box = $best_score == $statistics_hash{scoreDTarget} ? "Dbox" : "D'box";
			    }

			    ## remove temp files
			    `rm $locations{tmpPath}$locations{id}".fa" ` if ( -e $locations{tmpPath}.$locations{id}.".fa" ); 
			    `rm $locations{tmpPath}$locations{id}".aln"`  if ( -e $locations{tmpPath}.$locations{id}.".aln" );
			    
			    if ( !$has_target_site ){

				## Lower the threshold to distinguish between true and false positive target sites for 
				## snoRNA families with certain target switches, losses and duplications
				## Currently considered: CD_5, CD_19
				
				if ( $opts_k eq "deu" && ( $res->{name} =~ /CD_5\-/  || $res->{name} =~ /CD_19\-/ )){
				   
				    my ( $D_score, $prime_score ) = split ( /\s/, $best_score );
				    if ( $D_score > 0.6 || $prime_score > 0.6 ){
					
					printf STDERR "done. ( type: $target_box, score: %.3f)\n", $best_score if !$opts_v;
					
				    } 
				    else{
					
					printf STDERR "none found. ( best score: %.3f, threshold: 0.600). ", $best_score if !$opts_v;					
					&removeCandidate( "Candidate was rejected.", $locations{mfasta_rejected_output}, $res->{bla}, $res->{seq} );
					
					$res->{removed} = 1;
					$correct_snoRNAs--;
					#next HITS;
					
				    }
				    
				}
				else{

				    printf STDERR "none found. ( best score: %.3f, threshold: 0.700). ", $best_score if !$opts_v;
				    &removeCandidate( "Candidate was rejected.", $locations{mfasta_rejected_output}, $res->{bla}, $res->{seq} );
				    
				    $res->{removed} = 1;
				    $correct_snoRNAs--;
				    #next HITS;

				}
			    }
			    else{ 
				printf STDERR "done. ( type: $target_box, score: %.3f)\n", $best_score if !$opts_v;
			    }
			    
			}


			######
			## create a temporary count-number
			## in order to ensure correct computation of box motifs
			## by avoiding duplicated sequence names
			## this number will be checked again later on
			######
			if( $opts_b_blast ){
			
			    my $para = 1;
			    open ( FASTA,"<".$locations{mfasta_output} );
			    while ( <FASTA> ){
				if( $_ =~ />.*$res->{org}/ ){
				    $para++;
				}
			    }
			    close FASTA;
			    
			    if( $res->{name} !~ /-$para$/ ){
				my $oldname = $res->{name};
				$res->{name} =~ s/-\d+$/-$para/;
				$res->{bla} =~ s/$oldname/$res->{name}/;
			    }			    
			}


			######
			##
			## CHECK FOR CONSERVED BOX MOTIFS
			##
			## copy original multifasta in tmp directory, add novel snoRNA candidate
			## align all sequences by use of clustalw and search for conserved start positions of all boxtypes
			## in case any sequence of the original fastafile contains prime boxes, an additional prime box search will be done
			## in case no suitable box motif can be detected, the candidate will be rejected
			######


			my $offset;
			my $key_query;
			my $message_box_search = "";
			my ( $scoreA, $scoreB, $penaltyA, $penaltyB ) = ( 0, 0, 0, 0 );
			if( $opts_b_blast ){

			    print STDERR "\t\tbox motif search ", "." x 22, " " if !$opts_v;

			    ## try a pairwise alignment step solely for blast-generated snoRNA candidates
			    ## since there is no 'query' in the infernal-approach
			    if ( $res->{found_with} !~ /^cm_/ ){
			    
				## first, try a pairwise alignment step, 
				## including the novel candidate and the
				## the original best-scoring query
				`grep -A 1 $res->{found_with_seq} $locations{mfasta_output} > $locations{tmpPath}$locations{dat}`;
				my $query_header = `grep $res->{found_with_seq} $locations{mfasta_output} | sed 's/>//' `;
				chomp $query_header;
				$key_query = substr( $query_header, 0, 30 );
				
				open ( OUT, ">>$locations{tmpPath}$locations{dat}" ) or die $!;
				print OUT ">".$res->{bla}."\n".$res->{seq}."\n";
				close OUT;
				
				( $res->{boxA}, $res->{startA}, $res->{boxB}, $res->{startB}, $found, $offset, $message_box_search, $scoreA, $scoreB, $penaltyA, $penaltyB ) = boxSearch_pairwise( $locations{tmpPath}.$locations{dat}, \%seqProperties, $res, $opts_v, $key_query );

			    }
			    else{ $found = 0; }

			    ## in case neither appropriate box motifs nor a kink-turn could be found
			    ## try the old 'whole alignment' way
			    if ( !$found ){
			    
				`cp $locations{mfasta_output} $locations{tmpPath}$locations{dat}`;

				open ( OUT, ">>$locations{tmpPath}$locations{dat}" ) or die $!;
				print OUT ">".$res->{bla}."\n".$res->{seq}."\n";
				close OUT;

				## in this version each alignment is computed by 'muscle'
				## indicated by the boolean value "True" that is passed to boxSearch
				( $res->{boxA}, $res->{startA}, $res->{boxB}, $res->{startB}, $found, $offset, $message_box_search, $scoreA, $scoreB, $penaltyA, $penaltyB ) = boxSearch( $locations{tmpPath}.$locations{dat}, \%seqProperties, $res, "True", $opts_v );	

			    }

			    printf STDERR "done. ( C-box: $res->{boxA}, score: %.3f, penalty: $penaltyA and D-box: $res->{boxB}, score: %.3f, penalty: $penaltyB )\n", $scoreA, $scoreB if $found && !$opts_v;


			    $statistics_hash{scoreAbox} = $scoreA;
			    $statistics_hash{scoreBbox} = $scoreB;
			    $statistics_hash{penaltyAbox} = $penaltyA;
			    $statistics_hash{penaltyBbox} = $penaltyB;
			    $res->{penaltyA} = $penaltyA;
			    $res->{penaltyB} = $penaltyB;
			    $statistics_hash{Kturn} = 1;
			    $statistics_hash{Kturn} = 0 if $message_box_search =~ /K-turn/;
			    
			    
			    if( $found && $res->{type} eq "CD" && $hasPrimeBox ){

				print STDERR "\t\tprime box motif search ", "." x 16, " " if !$opts_v;
				my ( $tmp_boxB, $tmp_startB, $tmp_penaltyB, $tmp_scoreB );

				if ( $res->{found_with} !~ /^cm_/ ){ 
				    ## again try a pairwise approach at first
				    `grep -A 1 $res->{found_with_seq} $locations{mfasta_output} > $locations{tmpPath}$locations{dat}`;
				    
				    open ( OUT, ">>$locations{tmpPath}$locations{dat}" ) or die $!;
				    print OUT ">".$res->{bla}."\n".$res->{seq}."\n";
				    close OUT;
				    
				    ( $res->{boxAPrime}, $res->{startAPrime}, $res->{boxBPrime}, $res->{startBPrime}, $foundPrime, $statistics_hash{scoreAprimeBox}, $statistics_hash{scoreBprimeBox}, $statistics_hash{penaltyAprimeBox}, $statistics_hash{penaltyBprimeBox} ) = primeBoxSearch_pairwise( $locations{tmpPath}.$locations{dat}, \%seqProperties, $res, $offset, $opts_v, $key_query );
				    
				    ## in case a suitable D' box was found
				    ## save this one and apply alignment approach
				    ## to find a suitable C' box
				    if( $foundPrime == 2 ){
					$tmp_boxB = $res->{boxBPrime};
					$tmp_startB = $res->{startBPrime};
					$tmp_scoreB = $statistics_hash{scoreBprimeBox};
					$tmp_penaltyB = $statistics_hash{penaltyBprimeBox};
					$foundPrime = 0;
				    }
				}
				else{ $foundPrime = 0; }

				
				if ( !$foundPrime ){
				    
				    ## copy whole fasta file 
				    `cp $locations{mfasta_output} $locations{tmpPath}$locations{dat}`;
				    open ( OUT, ">>$locations{tmpPath}$locations{dat}" ) or die $!;
				    print OUT ">".$res->{bla}."\n".$res->{seq}."\n";
				    close OUT;
			
				    ( $res->{boxAPrime}, $res->{startAPrime}, $res->{boxBPrime}, $res->{startBPrime}, $foundPrime, $statistics_hash{scoreAprimeBox}, $statistics_hash{scoreBprimeBox}, $statistics_hash{penaltyAprimeBox}, $statistics_hash{penaltyBprimeBox} ) = primeBoxSearch( $locations{tmpPath}.$locations{dat}, \%seqProperties, $res, $offset, $opts_v );

				    ## keep the D' box that was found with the pairwise approach
				    if( $tmp_boxB ){
					$res->{boxBPrime} = $tmp_boxB ;
					$res->{startBPrime} = $tmp_startB;
					$statistics_hash{scoreBprimeBox} = $tmp_scoreB;
					$statistics_hash{penaltyBprimeBox} = $tmp_penaltyB;
					$foundPrime = 1;
				    }
				    
				}

				printf STDERR "done. ( C'-box: $res->{boxAPrime}, score: %.3f, penalty: $statistics_hash{penaltyAprimeBox} and D-box: $res->{boxBPrime}, score: %.3f, penalty: $statistics_hash{penaltyBprimeBox} )\n", $statistics_hash{scoreAprimeBox}, $statistics_hash{scoreBprimeBox}  if !$opts_v && $foundPrime && $res->{boxAPrime} ne "NNNNNNN";
				printf STDERR "done. ( D'-box: $res->{boxBPrime}, score: %.3f, penalty: $statistics_hash{penaltyBprimeBox} )\n", $statistics_hash{scoreBprimeBox} if !$opts_v && $foundPrime && $res->{boxAPrime} eq "NNNNNNN";
				print STDERR "none found.\n" if !$foundPrime;

				if ( $res->{boxAPrime} eq "NNNNNNN" ){
				    $statistics_hash{scoreAprimeBox} = "?";
				    $statistics_hash{penaltyAprimeBox} = "?";
				}


			    }

			    `rm $locations{tmpPath}$locations{id}.aln $locations{tmpPath}$locations{id}.fa `;
			    
			}
			else{
			    $found = 1;
			}			


			######
			##
			## REJECT snoRNA IN CASE NO BOXMOTIFS COULD BE FOUND
			##
			## write rejected snoRNA candidate into the `rejected` multifasta file
			## this file is located in the same directory as the original multifasta file
			######
			
			if ( !$found && !defined( $res->{removed} ) ){
			    $bool = 0;
			    if ( $message_box_search =~ /K-turn/ ){
				print STDERR "no K-Turn found. ";
				&removeCandidate( "Candidate was rejected.", $locations{mfasta_rejected_output}, $res->{bla}, $res->{seq} );
			    }
			    elsif( $message_box_search =~ /box/i ) {
				if( $res->{type} eq "CD"){
				    printf STDERR "none found. ( C-box: $res->{boxA}, score: %.3f, penalty: $statistics_hash{penaltyAbox} and D-box: $res->{boxB}, score: %.3f, penalty: $statistics_hash{penaltyBbox}). ", $statistics_hash{scoreAbox}, $statistics_hash{scoreBbox} if !$opts_v; 
				}
				else{
				    printf STDERR "none found. ( H-box: $res->{boxA}, score: %.3f, penalty: $statistics_hash{penaltyAbox} and ACA-box: $res->{boxB}, score: %.3f, penalty: $statistics_hash{penaltyBbox}). ", $statistics_hash{scoreAbox}, $statistics_hash{scoreBbox} if !$opts_v;
				}
				&removeCandidate( "Candidate was rejected.", $locations{mfasta_rejected_output}, $res->{bla}, $res->{seq} );
			    }

			    $correct_snoRNAs--;
			    $res->{removed} = 1;
			    
			}


			######
			##
			## EXTRACT snoRNA-SPECIFIC PROPERTIES
			##
			## this includes basics like renaming, folding, target prediction
			######

			if( $found && !defined( $res->{removed} )){
			    
			    $bool = 1;

			    ######
                            ## check if paraolog number is correct
			    ## therefore, count paralogous snoRNAs in the original multifasta file
			    ## rename current snoRNA-candidate if the presumed paralog-number is incorrect
			    ######

			    if( $opts_b_blast ){
			
				my $para = 1;
				open ( FASTA,"<".$locations{mfasta_output} );
				while ( <FASTA> ){
				    if( $_ =~ />.*$res->{org}/ ){
					$para++;
				    }
				}
				
				if( $res->{name} =~ /-\d+$/ || ( $res->{name} !~ /-/ && $para > 1 ) ){
				    if( $res->{name} !~ /-$para$/ ){
					my $oldname = $res->{name};
					$res->{name} =~ s/-\d+$/-$para/;
					$res->{bla} =~ s/$oldname/$res->{name}/;
				    }
				}
				close FASTA;
			    
				if($para == 2){
				    if( ( my $tmptmp = `grep $res->{org} $locations{mfasta_output}` ) !~ /-1_$res->{org}/ ){
					$tmptmp =~ s/^>//;
					`cat $locations{mfasta_output} | sed 's/_$res->{org}/-1_$res->{org}/' > $locations{mfasta_output}.tmp; mv $locations{mfasta_output}.tmp $locations{mfasta_output}`;
					$tmptmp = substr($tmptmp, 0, 30);
					( my $newtmp = $tmptmp ) =~ s/_$res->{org}/-1_$res->{org}/;
					$newtmp = substr( $newtmp, 0, 30 );
					$seqProperties{$newtmp} = $seqProperties{$tmptmp};
					delete $seqProperties{$tmptmp};
				    }
				}

			    }
			
    
			    ######
			    ## create anchor-constraints
			    ## all boxes are marked with an x and will be aligned underneath each other
			    ######

			    $boxes = "." x ( $res->{startA} - 1 ) . "x" x ( length $res->{boxA} );
			    $boxes .= "." x ( $res->{startB} - ( length $boxes ) - 1 ) . "x" x ( length $res->{boxB} );
			    $boxes .= "." x ( (length $res->{seq}) - (length $boxes) ) if $res->{type} eq "CD"; 
			    $boxes .= "x" x ( (length $res->{seq}) - (length $boxes) ) if $res->{type} eq "HACA";
			    $res->{box_cons}= $boxes;
			    $res->{length} = length($res->{seq});
			    
			    
			    ######
			    ##
			    ## CONSTRAINT FOLDING 
			    ##
			    ## typespecific folding of snoRNAs is computed by RNAsubopt of the Vienna RNA package
			    ## CD-snoRNAs: unpaired regions include both terminal boxes (box C and D) and the enclosed region
			    ## HACA-snoRNA: both boxes are unpaired, the regions ` 5'start - box H ` and ` box H - box ACA ` should form independent hairpins
			    ######

			    
			    print STDERR "\t\tstructure prediction ", "." x 18 ," " if !$opts_v;
			    
			    my ( $targets1, $targets2 );
			    my $message_target_prediction = "";
			    my $target_count = 0;
			    my $tmporg = $res->{org};
			    $tmporg =~ s/\./_/;

			    if( $res->{type} eq "CD" ){

				($res->{folding}, $res->{mfe}) = getFoldingCD($res->{seq}, $res->{startA}, ($res->{startB}+length($res->{boxB})) );
				
				$res->{length} = length($res->{seq});
				&cutFolding($res);

				## in case no C' box was found, which is indicated by "NNNNNN", 
				## reset the startposition of that box to '0'
				if ( $foundPrime && $res->{boxAPrime} =~ /NNN/ ){
				    $res->{startAPrime} = 0;
				}

				if( $res->{mfe} eq "atypical structure" ){
				    print STDERR "atypical structure found.\n" if !$opts_v;
				    $statistics_hash{mfe} = 0;
				}
				else{ 
				    print STDERR "done.\n" if !$opts_v; 
				    $statistics_hash{mfe} = 1;
				}

				
				######
				##
                                ## TARGET PREDICTION OF CD-snoRNAs
				##
				######
				if( $TARGET_PATH ){

				  print STDERR "\t\ttarget prediction ", "." x 21 , " " if !$opts_v;
				  
				  ( $message_target_prediction, $target_count ) = &CDtargetPrediction( $res, $tmporg, "B" ); 

				  print STDERR "none found" if !$opts_v && $target_count == 0;
				  print STDERR "$target_count putative targets found" if !$opts_v && $target_count > 0;
				  print STDERR " ($message_target_prediction).\n" if $message_target_prediction ne "" && !$opts_v;
				  print STDERR ".\n" if $message_target_prediction !~ /target/ && !$opts_v;


				  ######
				  ##
				  ## PRIME-TARGET PREDICTION OF CD-snoRNAs
				  ##
				  ######
				
				  if( defined( $res->{boxBPrime} ) ) {
				    
				    print STDERR "\t\tprime target prediction ", "." x 15, " " if !$opts_v;
				    
				    ( $message_target_prediction, $target_count ) = &CDtargetPrediction( $res, $tmporg, "BPrime" );
				    
				    print STDERR "none found" if !$opts_v && $target_count == 0;
				    print STDERR "$target_count putative targets found" if !$opts_v && $target_count > 0;
				    print STDERR " ($message_target_prediction).\n" if $message_target_prediction ne "" && !$opts_v;
				    print STDERR ".\n" if $message_target_prediction !~ /target/ && !$opts_v;
				    
				  }
				}

			    } else {

				&cutFoldingHACA( $res );

				( $res->{folding}, $res->{mfe} ) = getFoldingHACA( $res->{seq}, $res->{box_cons} ) ;
				( $res->{alternative_folding}, $res->{alternative_mfe}, $statistics_hash{mfe_HP1}, $statistics_hash{mfe_HP2}, $statistics_hash{length_HP1}, $statistics_hash{length_HP2}) = getIndependentFoldingHACA( $res->{seq}, $res->{box_cons}, $res->{startA}, $res->{startB} ) ;


				if( $res->{mfe} =~ /hairpin|strange/ ){ print STDERR "atypical structure found.\n" if !$opts_v }
				else{ print STDERR "done.\n" if !$opts_v; }

  
				######
				##
				## HACA TARGET PREDICTION FOR HAIRPIN 1
				##
				## split sequence into two distinct hairpins, use boxes as end or start positions
				## perform single hairpin target prediction by the use of RNAsnoop
				## write results in a hairpin-specific outputfile
				## visualize the targetinteraction with a postscript picture
				######
				if( $TARGET_PATH ){

				  print STDERR "\t\thp1 target prediction ", "." x 17," "if !$opts_v;
				  
				  ( $res->{HP1},$res->{HP2} ) = splitHP( $res->{seq}, $res->{startA}, $res->{startB} );
				  &writeSingleHairpins( $res, $tmporg );
				  ( $message_target_prediction, $target_count ) = &singleHairpinPrediction( $res, "1" );

				  print STDERR "none found" if !$opts_v && $target_count == 0;
				  print STDERR "$target_count putative targets found" if !$opts_v && $target_count > 0;
				  print STDERR " ($message_target_prediction).\n" if !$opts_v && $message_target_prediction ne "";
				  print STDERR ".\n" if !$opts_v && $message_target_prediction eq "";

				  print STDERR "\t\thp2 target prediction ", "." x 17, " " if !$opts_v;
				
				  ######
				  ##
				  ## HACA TARGET PREDICTION FOR HAIRPIN 2
				  ##
				  ## using the exact same procedure as for hairpin 1
				  ######
				
				  ( $message_target_prediction, $target_count ) = &singleHairpinPrediction( $res, "2" );
				    
				  print STDERR "none found" if !$opts_v && $target_count == 0;
				  print STDERR "$target_count putative targets found" if !$opts_v && $target_count > 0;
				  print STDERR " ($message_target_prediction).\n" if !$opts_v && $message_target_prediction ne "";
				  print STDERR ".\n" if !$opts_v && $message_target_prediction eq "";
				  
				}
			  }

			    ######
			    ##
			    ## UPDATE RESULT HASH AND FASTA-HEADER
			    ##
			    ######
			    
			    if( $opts_b_blast ){
				if ( $foundPrime ) { 
				    $res->{bla} .= "_" . $res->{boxAPrime} . "_" . $res->{startAPrime} . "_" . $res->{boxBPrime} . "_" . $res->{startBPrime};
				}
				else{
				    $res->{boxAPrime} = "";
				    $res->{startAPrime} = "";
				    $res->{boxBPrime} = "";
				    $res->{startBPrime} = "";
				}
				$res->{bla} .= "_" . $res->{boxA} . "_" . $res->{startA} . "_" . $res->{boxB} . "_" . $res->{startB};
			    }
			   

			    
			    ######
			    ##
			    ## WRITE SEQUENCES TO MULTIFASTA FILES
			    ##
			    ## write into mfasta file that is denoted in the db, and into the file given with option -c (unless they are not the same)
			    ######
			    
			    if( $opts_b_blast ){

				if( $locations{mfasta_output} ne $locations{mf} ){
				    open( FASTA, ">>".$locations{mfasta_output} ) or print $!;
				    print FASTA ">".$res->{bla}."\n".$res->{seq}."\n";
				    close FASTA;
				}

			    }
			    
			    ######
			    ## change Header in multifasta files if a paralog number was added
			    ######
			    else{

				#print "FEHLER: ".$res->{name}."\t";
				if ( $res->{name} =~ /-\d+$/ ){
				    ( my $tmp_name = $res->{name} )=~ s/-\d+$//;
				    #print $tmp_name."\t";
				    open( MF, "<" . $locations{mfasta_output} );
				    my @mf = <MF>;
				    close MF;

				    open( MF, ">" . $locations{mfasta_output} );
				    foreach my $line ( @mf ){
					if( $line =~ /$tmp_name\_[A-Z]/ ){
					    #print $line."\t";
					    $line =~ s/$tmp_name/$res->{name}/;
					    #print $line."?\n";
					    print MF $line;
					}
					else{
					    print MF $line;
					}
				    }
				    close MF;
				}
	
			
			    }


			    ######
			    ##
			    ## UPDATE THE HASH OF SEQUENCES AND THEIR PROPERTIES
			    ##
			    ## update hash of sequence properties
			    ## it is used for clustalw alignments and boxproperties
			    ######

			    my $key = substr( $res->{bla}, 0, 30 );
			    my %key;
			    if( !$seqProperties{$key} ){
				$seqProperties{$key} = \%key;
				$key{box1} = $res->{boxA};
				$key{box1Start} = $res->{startA};
				$key{box2} = $res->{boxB};
				$key{box2Start} = $res->{startB};
				$key{box1L} = substr( $res->{seq},$key{box1Start} - 3, ( length $key{box1} ) + 4 );
				$key{box2L} = substr( $res->{seq},$key{box2Start} - 3, ( length $key{box2} ) + 4 );

				if( $foundPrime ){
				    $key{box1Prime} = $res->{boxAPrime};
				    $key{box1PrimeStart} = $res->{startAPrime};
				    $key{box2Prime} = $res->{boxBPrime};
				    $key{box2PrimeStart}= $res->{startBPrime};
				    $key{box1LPrime} = substr( $res->{seq},$key{box1PrimeStart} - 3, ( length $key{box1Prime} ) + 4 );
				    $key{box2LPrime} = substr( $res->{seq},$key{box2PrimeStart} - 3, ( length $key{box2Prime} ) + 4 );
				}

			    }
			    
			    
			    $res = extractProperties( $res );

			    #######
			    ## PRINT SNORNA STATISTICS
			    ##
			    ## print the elements in the statistics hash in the family specific file
			    ##
			    #######
			    if ( $res->{type} ){
				$statistics_hash{CD} = $res->{startB} - $res->{startA} - 7;
				$statistics_hash{CDp} = $res->{startBPrime} - $res->{startA} - 7 if $res->{startBPrime};
				$statistics_hash{CpD} = $res->{startB} - $res->{startAPrime} - 4 if $res->{startAPrime} && $res->{boxAPrime} ne "NNNNNNN";
			    }
			    &printStatistics( \%statistics_hash, $stat, $res->{type} ) if $opts_b_blast;

			    
			    #######
			    ##
			    ## PRINT SUMMARY OF NOVEL snoRNA
			    ##
			    ## print short summary, not all retrieved properties
			    ## this includes header, sequence, folding, mfe, best targets
			    #######

			    print STDERR "\n\tsummary\n" if !$opts_v;
			    print STDERR "\n\n" if $opts_v;
			    printf STDERR "\t%-13s %s\n", "header:", $res->{bla};
			    printf STDERR "\t%-13s %s\n", "sequence:", $res->{seq};
			    printf STDERR "\t%-13s %s %s\n", "structure:", $res->{folding}, $res->{mfe} if $res->{mfe} !~ /atypical|hairpin|strange/;
			    printf STDERR "\t%-13s %s\n", "structure:", "atypical structure found" if $res->{mfe} =~ /atypical|hairpin|strange/;
			    printf STDERR "\t%-13s %s %s\n", "alt folding:", $res->{alternative_folding}, $res->{alternative_mfe} if $res->{alternative_mfe};

			    if( $TARGET_PATH ){
			      if( $res->{type} eq "CD" ){
				printf STDERR "\t%-13s %s %s %s\n", "target1:", $res->{alnPos1}, $res->{structure1}, $res->{energy1} if $res->{energy1} && $res->{structure1} && $res->{alnPos1};
				printf STDERR "\t%-13s %s\n", "target1:", "no putative target found" if !$res->{energy1} && !$res->{structure1} && !$res->{alnPos1};
				printf STDERR "\t%-13s %s %s %s\n", "target2:", $res->{alnPos2}, $res->{structure2}, $res->{energy2} if $res->{energy2} && $res->{structure2} && $res->{alnPos2};
				printf STDERR "\t%-13s %s\n", "target2:", "no putative target found" if !$res->{energy2} && !$res->{structure2} && !$res->{alnPos2};
			    
			      }
			      else{
			      
				printf STDERR "\t%-13s %s %s %s\n", "target1:", $res->{alnPos1}, $res->{tarStruc1}, $res->{energy1} if $res->{energy1} && $res->{tarStruc1} && $res->{alnPos1};
				printf STDERR "\t%-13s %s\n", "target1:", "no putative target found" if !$res->{energy1} && !$res->{tarStruc1} && !$res->{alnPos1};
				printf STDERR "\t%-13s %s %s %s\n", "target2:", $res->{alnPos2}, $res->{tarStruc2}, $res->{energy2} if $res->{energy2} && $res->{tarStruc2} && $res->{alnPos2};
				printf STDERR "\t%-13s %s\n", "target2:", "no putative target found" if !$res->{energy2} && !$res->{tarStruc2} && !$res->{alnPos2};
			      }
			    }
			    print STDERR "\n\n";



			    ######
			    ##
			    ## UPDATE SPECIES-SPECIFIC GFF-FILE
			    ##
			    ## first, check whether this snoRNAs is already contained in the gff-file
			    ## otherwise, write a new entry
			    ######

			    my ( $start_tmp, $end_tmp ) = split( /,/,$res->{range} );
			    `touch $locations{gffPath}$res->{org}".gff"` if !-e $locations{gffPath}.$res->{org}.".gff";
			    if ( `grep '$res->{abbr}\_$res->{name}' $locations{gffPath}/$res->{org}.gff | grep -cP '$start_tmp\t$end_tmp'` == 0 ){  
				open( GFF, ">>" . $locations{gffPath} . $res->{org} . ".gff" );
				print GFF $res->{chr}."\tSNOSTRIP\t".$res->{type}."-snoRNA\t".(split(/,/,$res->{range}))[0]."\t".(split(/,/,$res->{range}))[1]."\t".$res->{score}."\t".$res->{strand}."\t.\t".$res->{abbr}."_".$res->{name}."\n";	
				close GFF;
			    }

			    $count++ if $count;
				
			}
			else{

			    print STDERR "\n\t\tAnalyze remaining features.\n\n" if !$opts_v;

			    ## analze structure of snoRNAs that were removed previously
			    if( $res->{type} eq "CD" ){
				if ( $res->{penaltyA} <= 1 && $res->{penaltyB} <= 1 ){ 
				    ($res->{folding}, $res->{mfe}) = getFoldingCD($res->{seq}, $res->{startA}, ($res->{startB}+length($res->{boxB})) );
				    $res->{length} = length($res->{seq});
				    $res->{mfe} eq "atypical structure" ? $statistics_hash{mfe} = -1 : $statistics_hash{mfe} = $res->{mfe};
				}
				else{
				    $res->{length} = length($res->{seq});
				    $res->{mfe} = "atypical structure";
				    $res->{mfe} eq "atypical structure" ? $statistics_hash{mfe} = -1 : $statistics_hash{mfe} = $res->{mfe};
				}
			    }
			    else{

				&cutFoldingHACA( $res );
				( $res->{alternative_folding}, $res->{alternative_mfe}, $statistics_hash{mfe_HP1}, $statistics_hash{mfe_HP2}, $statistics_hash{length_HP1}, $statistics_hash{length_HP2}) = getIndependentFoldingHACA( $res->{seq}, $res->{box_cons}, $res->{startA}, $res->{startB} ) ;


			    }

			    &printStatistics( \%statistics_hash, $stat_rej, $res->{type} );
			    next HITS;

			}

		    }
		    
		}

		&writeLog( --$correct_snoRNAs . " candidate sequence(s) found for " . $fam . " on " . $genome . " fulfill each snoRNA criteria", \%locations );
		

	    }
	    

	    if ( `grep -c ">" $locations{mfasta_output}` > 1 ){
		
		## calculate multiple sequence alignment
		&writeLog( "realigning of family $locations{id}", \%locations );
		`$MUSCLE -in $locations{mfasta_output} -out $locations{alignmentPath}$locations{id}.muscle.aln -clwstrict -quiet`;
		`$STOCKHOLM $locations{alignmentPath}$locations{id}.muscle.aln $locations{alignmentPath}$locations{id}.muscle.stk`;
	    
	    }

	    
	}
	

	&writeLog( "finished analysis of Family $locations{id}\n", \%locations); 
	
	## close statistic file
	close $stat;
	close $stat_rej;
	
    }

    `chmod -R g+w $locations{tmpPath}`;


    &cleanup();
    
}
else{
    &usage();
    $message ? die "ERROR: $message\n" : die "$message\n";
    
}







####################################################################################################
####################################################################################################
##
## SUBROUTINES
##
####################################################################################################
####################################################################################################







####################################################################################################
sub extractProperties
## analyses found snoRNA
####################################################################################################
{
    my $res = shift;
    $res->{length} = length( $res->{seq} );
    if( $res->{type} eq "HACA" ){
	analyseHACABoxes( $res );
    }
    else{
	analyseCDBoxes( $res );
    }
    return $res;
}


####################################################################################################
sub analyseHACABoxes
## define HACA & ACA boxes and their distance
####################################################################################################
{
    my $res = shift;

    $res->{distAB} = $res->{startB} - $res->{startA} - (length ( $res->{boxA} ));

    return $res;
}



####################################################################################################
sub analyseCDBoxes
## define C & D & C' & D' boxes and their distances
####################################################################################################
{
    my $res = shift;

    if (!$res->{startB} or !$res->{startA}){
      print STDERR "boxA: ", $res->{boxA};
      print STDERR "startA: ", $res->{startA}, "\n";
      print STDERR "boxB: ", $res->{boxB}, "\n";
      print STDERR "startB: ", $res->{startB}, "\n";
      if($tmp[-5] =~ /^\d+$/ && $tmp[-6] =~ /[GTACU]{4}/){
        print STDERR "Aprime: ", $res->{boxAPrime}, "\n";
        print STDERR "AprimeStart: ", $res->{startAPrime}, "\n";
        print STDERR "Bprime: ", $res->{boxBPrime}, "\n";
        print STDERR "BprimeStart: ", $res->{startBPrime}, "\n";
      }
      print STDERR $res->{bla}, "\n";
      print STDERR $res->{name}, "\n";
      print STDERR $res->{org}, "\n";
      print STDERR $res->{type}, "\n";
      print STDERR $res->{chr}, "\n";
      print STDERR $res->{range}, "\n";
      print STDERR $res->{strand}, "\n";
      print STDERR $res->{found_with}, "\n";      
    }
    $res->{distAB} = $res->{startB} - $res->{startA} - ( length ( $res->{boxA} ) );
    if( $res->{boxAPrime} ){
	$res->{distAPrimeB} = $res->{startB} - $res->{startAPrime} - (length ( $res->{boxAPrime} ));
	$res->{distAPrimeBPrime} = $res->{startAPrime} - $res->{startBPrime} - (length ( $res->{boxBPrime} ));
	$res->{distABPrime} = $res->{startBPrime} - $res->{startA} - (length ( $res->{boxA} ));
    }

    return $res;
}


####################################################################################################
sub cutFolding
## analyze the given folding and cut the sequence appropriately
####################################################################################################
{
    
    my $res = shift;
    my ($tmpStart, $tmpEnd) = split(/,/, $res->{range});


    my $B = $res->{length} - 3 - $res->{startB};
    my $A = $res->{startA} - 6;
    
    if( $res->{startA} > 6 ){
	$res->{folding} =~ s/^.{$A}//;
	$res->{seq} =~ s/^\w{$A}//;
	$res->{startA} -= $A;
	$res->{startB} -= $A;
	$res->{startAPrime} -= $A;
	$res->{startBPrime} -= $A;
	$res->{length} -= $A;
	$res->{box_cons} =~ s/^\.{$A}//;
	
    }
    if( $B > 5 ){
	$B -= 5;
	$res->{folding} =~ s/.{$B}$//;
	$res->{seq} =~ s/\w{$B}$//;
	$res->{box_cons} =~ s/\.{$B}$//;
	$res->{length} -= $B;
	if($res->{strand} eq "+"){
	    $tmpStart += $A;
	    $tmpEnd -= $B;
	}
	else{
	    $tmpStart += $B;
	    $tmpEnd -= $A;
	}
    }
    
    $res->{range} = $tmpStart.",".$tmpEnd;

    #update fasta-header
    $res->{bla} =~ s/_\d+,\d+_/_$res->{range}\_/;
    
}


####################################################################################################
sub cutFoldingHACA
## this method cuts the sequence of HACA snoRNA exactly 3nt after ACA-box
####################################################################################################
{

    my $res = shift;
    my $cutoff = $res->{length} - 5 - $res->{startB};
    my ($start, $end) = split(/,/, $res->{range});
    my ($newStart, $newEnd);
    if( $cutoff >= 1 ){
	$res->{seq} =~ s/\w{$cutoff}$//;

	if( $res->{strand} eq "+" ){
	    $newStart = $start;
	    $newEnd = $end - $cutoff;
	}
	else{
	    $newEnd = $end;
	    $newStart = $start + $cutoff;
	}
	$res->{range} = "$newStart,$newEnd";
	$res->{bla} =~ s/$start,$end/$newStart,$newEnd/;
	$res->{box_cons} =~ s/x{$cutoff}$//;
	$res->{length} -= $cutoff;
	
    }
    elsif( $cutoff < 0 ){

	if( $res->{strand} eq "+" ){
	    $newStart = $start;
	    $newEnd = $end - $cutoff;
	}
	else{
	    $newEnd = $end;
	    $newStart = $start + $cutoff;
	}

	my ($tmpSeq, $tmpStart,$tmpEnd) = extractSeq($res->{seq}, $res->{database}, $newStart, $newEnd, getStrand($res->{strand}), $res->{chr});
	if($tmpSeq ne $res->{seq}){
	    if((length ($tmpSeq) - 5 - $res->{startB}) == 0){
		$res->{seq} = $tmpSeq;
		$res->{bla} =~ s/$start,$end/$newStart,$newEnd/;
		$res->{box_cons} .= "x" x abs($cutoff);
		$res->{length} = length $tmpSeq;
		$res->{folding} .= "." x abs($cutoff);
		$res->{range} = "$newStart,$newEnd";
	    }
	}
    }

}




####################################################################################################
sub getStrand
## returns the strand of a hit in numerical forma
####################################################################################################
{

    return 1 if $_[0] eq "+";
    return 2 if $_[0] eq "-";

}



####################################################################################################
sub extractSeq
## extract the seq from genome file
## if end point is larger than contig extract as much nucleotides as possible
####################################################################################################
{

    my ($seq, $db, $start, $end, $strand, $chr) = @_;
    my $tmpSeq;
    my $bool = 0;
    
    $tmpSeq = `$FASTACMD -d $db -s \"$chr\" -L "$start,$end" -S $strand -l 500 2>&1`;

    if($tmpSeq =~ /ERROR/){ $bool = 1; $tmpSeq =~ s/^\[fastacmd\][^\n]+\n+//;}
    $tmpSeq =~ s/^>[^\n]+\n//;
    $tmpSeq =~ s/\s//g;

    #test if fastacmd returned the whole chromosome, which would be the case if an error would have occured
    if(!$tmpSeq || $end-$start <= (length $tmpSeq) - 2 || $bool){
	$end = (length $tmpSeq) + 2;
	$seq = substr($tmpSeq, $start - 1) if $strand == 1;
	$seq = substr($tmpSeq, 0, $end - $start ) if $strand == 2;
	return ($seq, $start, $end);
    }

    return ($tmpSeq, $start, $end);

}




####################################################################################################
sub readMultiFasta
## read a given multifasta file and fill the hash seqProperties
####################################################################################################
{

    ######
    ## VARIABLES 
    ######
    my %locations = %{ $_[0] };
    my %seqProperties = ();

    my ( $keyRef, $tmpRef, $head, $seq, $boxes );
    my @tmp;


    open ( MF, $locations{mf} ) or die $!;
 
    ######
    ## analyse sequences already contained in the given multifasta-file
    ## go through the input fastafile and extract all necessary information about boxes and queryname
    ######
    while ( <MF> ) {

	chomp $_;
	if( $_ =~ /^>/ ) {

	    $_ =~ s/^>//;
	    my $res;
	    $head = $_;
	    @tmp = split(/_/, $_);

	    $res->{bla} = $_;
	    my %key;
	    my $key = substr($_, 0, 30);
	    $keyRef = \%key;
	    if(!$seqProperties{$key}){
		$seqProperties{$key} = \%key;
		$key{box1} = $tmp[-4];
		$key{box1Start} = $tmp[-3];
		$key{box2} = $tmp[-2];
		$key{box2Start} = $tmp[-1];
		if($tmp[-5] =~ /^\d+$/ && $tmp[-6] =~ /[GTACU]{4}/){
		    $key{box1Prime} = $tmp[-8];
		    $key{box1PrimeStart} = $tmp[-7];
		    $key{box2Prime} = $tmp[-6];
		    $key{box2PrimeStart}= $tmp[-5];
		    $hasPrimeBox = 1;
		}
	    }

	    $res->{boxA} = $key{box1};
	    $res->{startA} = $key{box1Start};
	    $res->{boxB} = $key{box2};
	    $res->{startB} = $key{box2Start};
	    if($tmp[-5] =~ /^\d+$/ && $tmp[-6] =~ /[GTACU]{4}/){
		$res->{boxAPrime} = $key{box1Prime};
		$res->{startAPrime} = $key{box1PrimeStart};
		$res->{boxBPrime} = $key{box2Prime};
		$res->{startBPrime} = $key{box2PrimeStart};
	    }
	    $res->{bla} =~ /^([^_]+)_([^_]+)_([\w].[^_]+)_\(([^\)]+)\)_([^_]+)_([^_])/;
	    $res->{name} = $1."_".$2;
	    $res->{org} = $3;
	    $res->{type} = $1;
	    $res->{chr} = $4;
	    $res->{range} = $5;
	    $res->{strand} = $6;
	    $res->{found_with} = "";


	    ######
	    ## check whether this snoRNA was found by another query
	    ######
	    if ( $res->{bla} =~ /^([^_]+)_([^_]+)_([\w].[^_]+)_\(([^\)]+)\)_([^_]+)_([^_])_([a-z]{3})_([A-Z]{2,})_([^_]+)_([\w].[^_]+)/ ){
		$res->{found_with} = $7 . "_" . $8 . "_" . $9;
	    }

	    ( my $expression = $res->{org} ) =~ s/^(\w)\./$1.*\\s/;
#	    ( $res->{genome}, $res->{organism}, $res->{abbr}, $res->{source}, $res->{db} ) = split( /\t/, `grep -P '$expression' $locations{information_file}`);
	     
	    $bigHash{$res->{name}."_".$res->{org}} = $res;
	    $tmpRef = $res;
	    

	    ######
	    ## SAVE ORIGINAL snoRNA NAME 
	    ######
	    
	    ( $tmpRef->{alternative_name} = $tmpRef->{name} ) =~ s/^[^_]+_//;
	    
	    ######
	    ## CHECK IF SEQUENCE CONTAINS PARALOG NUMBER
	    ######
	    
	    if( $tmpRef->{name} !~ /-\d+$/){
		my $tmp_name = $tmpRef->{name};
		my $tmp_key = substr( $tmpRef->{bla}, 0, 30 );
		$tmpRef->{name} .= "-1";
		$tmpRef->{bla} =~ s/$tmp_name/$tmp_name-1/;
		$bigHash{$res->{name}."_".$res->{org}} = $tmpRef;
		delete($bigHash{$tmp_name."_".$res->{org}});
		$seqProperties{substr( $tmpRef->{bla}, 0, 30 )} = \%key;
		delete $seqProperties{$tmp_key};
	    }


	} 

	else {

	    $seq = $_;
	    $tmpRef->{seq} = $seq;
	    chomp $seq;
	    chomp $head;


	    ######
	    ## CALCULATE BOX DISTANCES
	    ######
	    &extractProperties( $tmpRef );


	    ######
	    ## create anchor-constraints
	    ## all boxes are marked with an x and will be aligned underneath each other
	    ######
	    
	    $boxes = "." x ( $tmpRef->{startA} - 1 ) . "x" x ( length $tmpRef->{boxA} );
	    $boxes .= "." x ( $tmpRef->{startB} - ( length $boxes ) - 1 ) . "x" x ( length $tmpRef->{boxB} );
	    $boxes .= "." x ( (length $tmpRef->{seq}) - (length $boxes) ) if $tmpRef->{type} eq "CD"; 
	    $boxes .= "x" x ( (length $tmpRef->{seq}) - (length $boxes) ) if $tmpRef->{type} eq "HACA";
	    $tmpRef->{box_cons}= $boxes;


	    #add wider boxmotifs to hash 'seqProperties'
	    $keyRef->{box1L} = substr($_, $keyRef->{box1Start} - 3, (length $keyRef->{box1}) + 4);
	    $keyRef->{box2L} = substr($_, $keyRef->{box2Start} - 3, (length $keyRef->{box2}) + 4);

	    if($keyRef->{box1Prime}){
		$keyRef->{box1LPrime} = substr($_, $keyRef->{box1PrimeStart} - 3, (length $keyRef->{box1Prime}) + 4);
		$keyRef->{box2LPrime} = substr($_, $keyRef->{box2PrimeStart} - 3, (length $keyRef->{box2Prime}) + 4);
	    }


	}
	
    }
    
    close MF;
    return \%seqProperties;
    
}



####################################################################################################
sub writeSingleHairpins
## writes both hairpins into separated files
## these files are needed for single hairpin target prediction
####################################################################################################
{

    my ( $res, $organism ) = @_;

    open( HP1, ">>" . $locations{targetPath} . $locations{id} . "_HP1.fa" );
    open( HP2, ">>" . $locations{targetPath} . $locations{id} . "_HP2.fa" );

    print HP1 ">", $organism, "_", $res->{name}, "_HP1\n", $res->{HP1}, "\n";
    print HP2 ">", $organism, "_", $res->{name}, "_HP2\n", $res->{HP2}, "\n";

    close HP1;
    close HP2;

}



####################################################################################################
sub singleHairpinPrediction
## run target prediction on single hairpins with RNAnoop
####################################################################################################
{

    my ( $res, $nr ) = @_;

    my ( $targets, $message ) = singlesnoop( $res->{"HP".$nr}, $res, \%locations, $opts_k );
    my $target_count = 0;
    
    if( $targets->[0] ) {
	($res->{"RNA".$nr}, $res->{"pseudo".$nr}, $res->{"tarSeq".$nr}, $res->{"lhs".$nr}, $res->{"rhs".$nr}, $res->{"energy".$nr}, $res->{"tarStruc".$nr}, $res->{"distLoop".$nr}, $res->{"distEnd".$nr}, $res->{"alnPos".$nr}, $res->{"rawSeq".$nr} ) = ( $targets->[0]->{targetRNA}, $targets->[0]->{pseudo}, $targets->[0]->{tarSeq}, $targets->[0]->{lhs}, $targets->[0]->{rhs}, $targets->[0]->{energy}, $targets->[0]->{struc}, $targets->[0]->{loop}, $targets->[0]->{ending}, $targets->[0]->{alnPos}, $targets->[0]->{seq});
	
	open( TAR, ">>" . $locations{targetPath} . $locations{id} . "_targets_" . $nr );
	foreach ( @{$targets} ){
	    if( $_->{targetRNA} ne "NA" && `grep -P '$res->{org}\\s+$res->{name}-HP$nr' $locations{targetPath}$locations{id}_targets_$nr | grep -cP '$_->{targetRNA}-$_->{pseudo}\\s'` == 0){
		printf TAR "%-20s\t%-15s\t%-10s\t%-10s\t%-6s\t%-100s\t%-100s\n", $res->{org}, $res->{name} . "-HP" . $nr, $_->{targetRNA} . "-" . $_->{pseudo}, $_->{alnPos}, $_->{energy}, $_->{struc}, $_->{seq};
		&makePostscript( $_->{snoop}, \%locations );
		$target_count++;
	    }
	}
	close TAR;
    }
 
    return ( $message, $target_count );
   
}



####################################################################################################
sub CDtargetPrediction
## run target prediction with RNAplex 
####################################################################################################
{

    my ( $res, $organism, $which_box ) = @_;

    my $target_count = 0;
    my ( $nr, $box, $b );

    if( $which_box eq "B" ){
	$nr = "1";
	$box = "D";
	$b = "D";
    }
    else{
	$nr = "2";
	$box = "D prime";
	$b = "Dprime";
    }

    my ( $targets, $message ) = &computeTargetCD( $res, $organism, \%locations, $res->{"start".$which_box}, $opts_k );

    if( $targets->[0] ) {
	
	( $res->{"tarRNA".$nr}, $res->{"meth".$nr}, $res->{"tarSeq".$nr}, $res->{"bindSeq".$nr}, $res->{"box".$nr}, $res->{"structure".$nr}, $res->{"energy".$nr}, $res->{"alnPos".$nr} ) = ( $targets->[0]->{targetRNA}, $targets->[0]->{mod}, $targets->[0]->{tarSeq}, $targets->[0]->{bindSeq}, $box, $targets->[0]->{raw_structure}, $targets->[0]->{energy}, $targets->[0]->{alnPos} ); 
	open( TAR, ">>". $locations{targetPath} . $locations{id} . "_targets_" . $nr );

	foreach( @{$targets} ){
	    if( $_->{targetRNA} ne "NA" ){
		printf TAR "%-20s\t%-15s\t%-10s\t%-10s\t%-6s\t%-45s%-45s\n", $res->{org}, $res->{name} . "-" . $b, $_->{targetRNA} . "-" . $_->{mod}, $_->{alnPos}, $_->{energy}, $_->{structure}, $_->{seq2};
		$target_count++;
	    }
	}
	close TAR;
	
    }
    
    return ( $message, $target_count );
    
}



####################################################################################################
sub printStatistics
## this methods prints the elements from the given hash
## into a tab-separated file
####################################################################################################
{

    my $hash_ref = shift;
    my $fh = shift;
    my $type = shift;


    if ( $type eq "CD" ){
	$hash_ref->{scoreAbox} /= 7;
	$hash_ref->{scoreBbox} /= 4;
	$hash_ref->{scoreAprimeBox} /= 7 if $hash_ref->{scoreAprimeBox} && $hash_ref->{scoreAprimeBox} ne "?";
	$hash_ref->{scoreBprimeBox} /= 4 if $hash_ref->{scoreBprimeBox};
    } 
    elsif( $type eq "HACA" ){
	$hash_ref->{scoreAbox} /= 6;
	$hash_ref->{scoreBbox} /= 3;
    }


    printf $fh "%.3f\t" , $hash_ref->{scoreAbox};
    print $fh $hash_ref->{penaltyAbox}, "\t";
    printf $fh "%.3f\t", $hash_ref->{scoreBbox};
    print $fh $hash_ref->{penaltyBbox}, "\t";
    
    if ( $type eq "HACA" ){

	$hash_ref->{mfe_HP1} ? print $fh $hash_ref->{mfe_HP1} : print $fh "?";
	print $fh "\t";
	print $fh $hash_ref->{length_HP1} if $hash_ref->{length_HP1};
	print $fh "\t";
	$hash_ref->{mfe_HP2} ? print $fh $hash_ref->{mfe_HP2} : print $fh "?";
	print $fh "\t";
	print $fh $hash_ref->{length_HP2} if $hash_ref->{length_HP2};


    }
    else{

	printf $fh "%.3f", $hash_ref->{scoreDTarget} if defined($hash_ref->{scoreDTarget});
	print $fh "\t";

	if ( "scoreAprimeBox" ~~ %{$hash_ref} && $hash_ref->{scoreAprimeBox} ne "?" ) {
	    printf $fh "%.3f\t", $hash_ref->{scoreAprimeBox};
	    print $fh $hash_ref->{penaltyAprimeBox}, "\t";
	}
	else{
	    print $fh "?\t?\t";
	    $hash_ref->{CpD} = "?";
	}

	if ( "scoreBprimeBox" ~~ %{$hash_ref} ){
	    printf $fh "%.3f\t", $hash_ref->{scoreBprimeBox};
	    print $fh $hash_ref->{penaltyBprimeBox}, "\t";
	}
	else{
	    print $fh "?\t?\t";
	    $hash_ref->{CDp} = "?";
	}
	printf $fh "%.3f", $hash_ref->{scoreDprimeTarget} if defined($hash_ref->{scoreDprimeTarget});
	print $fh "\t";
	
	if ( $hash_ref->{CD}){ print $fh $hash_ref->{CD} , "\t" }
	else{ print $fh "?\t" }
	if ( $hash_ref->{CDp}){ print $fh $hash_ref->{CDp}, "\t" }
	else{ print $fh "?\t" }
	if (  $hash_ref->{CpD} ){ print $fh $hash_ref->{CpD}, "\t" }
	else{ print $fh "?\t" }

	print $fh $hash_ref->{Kturn}, "\t" if $hash_ref->{Kturn} == 0 || $hash_ref->{Kturn} == 1;
	print $fh $hash_ref->{mfe};
    }

    print $fh "\t", $hash_ref->{snoRNA} , "\n";

}



####################################################################################################
sub correctInformationfile
## checks whether the information file is formatted correctly
####################################################################################################
{

    my ( @line );

    if ( !-e $opts_i ){
	return 0;
    }

    open INF, $opts_i or die $!;
    while ( <INF> ){
	@line = split ( /\t/, $_ );
	if( scalar @line != 3 ){
	    return 0;
	}
    }

    return 1;

}



####################################################################################################
sub removeCandidate
## writes the current sequence that will be deleted into the rejected multi fasta file
####################################################################################################
{

    my ( $message, $location, $header, $sequence ) = @_;
    
    print STDERR  "$message\n" if !$opts_v;
			    
    open (FASTA, ">>".$location );
    print FASTA ">", $header, "\n", $sequence, "\n";
    close FASTA;

}



####################################################################################################
sub setOptions
## checks whether all necessary arguments are specified or not
####################################################################################################
{

    GetOptions(

	"h|help"              =>  \$opts_h,
	"d|dir=s"             =>  \$opts_d,
	"i|info=s"            =>  \$opts_i,
	"g|genome=s"          =>  \$opts_g,
	"n|name=s{2}"         =>  \@opts_n,
	"k|kingdom=s"         =>  \$opts_k,
	"c|clan=s{,}"         =>  \@opts_c,
	"targets"             =>  \$opts_t,
	"s|sequences=s"       =>  \$opts_s,       
	"a|alignments=s"      =>  \$opts_a,
	"p|profiles=s"        =>  \$opts_p,
	"q|quiet"             =>  \$opts_v,
	"f|force"             =>  \$opts_f,
	"j|infernal"          =>  \$opts_j,

    );


}



####################################################################################################
sub setConstants
## checks whether all necessary arguments are specified or not
####################################################################################################
{

    my $abs_path = abs_path($0);
    $abs_path =~ s/snoStrip.pl//;
    my $data_path = $abs_path;
    $data_path =~ s/bin/data/;
  
    $STOCKHOLM = $CONFIG::STOCKHOLM;
    $MUSCLE = $CONFIG::MUSCLE;
    $INFERNAL = $CONFIG::INFERNAL;
    $BLAST = $CONFIG::BLAST;
    $FASTACMD = $CONFIG::FASTACMD;

    if( $opts_s ){ $TARGET_PATH = $opts_s }
    if( $opts_p ){ $TARGET_PROFILE_PATH = $opts_p }
    if( $opts_a ){ $TARGET_ALN_PATH = $opts_a }

    if( !$TARGET_PATH && $opts_t ){
	if ( $opts_k eq "deu" ){ $TARGET_PATH = $data_path . $CONFIG::DEUTEROSTOMIA_TARGET_RNA_PATH }
	elsif ( $opts_k eq "pro" ){ $TARGET_PATH = $data_path . $CONFIG::PROTOSTOMIA_TARGET_RNA_PATH }
	elsif ( $opts_k eq "pla" ){ $TARGET_PATH = $data_path . $CONFIG::PLANT_TARGET_RNA_PATH }
	else{ $TARGET_PATH = $data_path . $CONFIG::FUNGI_TARGET_RNA_PATH }
    }
    if( !$TARGET_ALN_PATH && $opts_t ){
	$TARGET_ALN_PATH = $TARGET_PATH . "alignments/";
    }
    if( !$TARGET_PROFILE_PATH && $opts_t ){
	$TARGET_PROFILE_PATH = $TARGET_PATH . "profiles/";
    }

    if( $opts_k eq "fun" ){ 
	$FASTA_PATH = $data_path . $CONFIG::FUNGI_FASTA_PATH;
	$INFO_FILE = $data_path . $CONFIG::FUNGI_INFORMATION_FILE;
	
    }
    elsif( $opts_k eq "pro" ){
	$FASTA_PATH = $data_path . $CONFIG::PROTOSTOMIA_FASTA_PATH;
	$INFO_FILE = $data_path . $CONFIG::PROTOSTOMIA_INFORMATION_FILE;
	
    }
    elsif( $opts_k eq "pla" ){
	$FASTA_PATH = $data_path . $CONFIG::PLANT_FASTA_PATH;
	$INFO_FILE = $data_path . $CONFIG::PLANT_INFORMATION_FILE;
	
    }
    elsif( $opts_k eq "deu" ){
	$FASTA_PATH = $data_path . $CONFIG::DEUTEROSTOMIA_FASTA_PATH;
	$INFO_FILE = $data_path . $CONFIG::DEUTEROSTOMIA_INFORMATION_FILE;
	
    }


}



####################################################################################################
sub checkOptions
## checks whether all necessary arguments are specified or not
####################################################################################################
{

    ## print help message in case option 'h' is set
    if ( $opts_h ) {
	return ( 0, "" );
    }

    ## check if all mandatory arguments are given, i.e., the o-option-file
    ## the information file, and output directory, kingdom
    if ( !$opts_d || !$opts_k ){
	return ( 0, "Seems like you forgot a mandatory argument!" );
      }

    
    ## check if the organism-specific arguments are given
    ## either a file providing genom source, name, and abbreviation
    ## or two command line options: genome source and name.
    if( !$opts_i && !( $opts_g && @opts_n )){
      	return ( 0, "Seems like you forgot to specify the organism(s)!" );
    }

    ## analyze the format of the given organism name
    if ( @opts_n){

	if( scalar( @opts_n) != 2 ){
	    return ( 0, "You have to provide exactly two arguments for option '-n|--name'!" );
	}
	if( $opts_n[0] !~ /^[A-Z][a-z]*$/ || $opts_n[1] !~ /^[a-z\-]+$/ ){
	    return ( 0, "The format of your input name has to be 'Genus epithet', e.g., Mus musculus!");
	}

	$opts_n = $opts_n[0] . " " . $opts_n[1];
    }
    
    ## check if the mandatory arguments have the correct format
    if ( $opts_i && !&correctInformationfile ){
	return ( 0, "Badly formatted information file" );
    }

    if ( !-w $opts_d ){
	return ( 0, "You're not permitted to write in the given directory!" );
    }
    
    

    ## sanity checks
    $opts_d .= "/" if $opts_d !~ /\w\/$/;

    return ( 1, "" );

}

####################################################################################################
sub cleanup
## delete all temporary files that have been created by snostrip
####################################################################################################
{

    `rm $locations{information_file}`;
    `rm -r $locations{tmpPath}`;

	

}


####################################################################################################
sub usage
## prints help message
####################################################################################################
{

  my $path = abs_path();
  $path =~ s/snoStrip$/data\/targets/;
  
    $usage = << "FIN";

    usage: perl $0 -d|dir DIRECTORY -k|kingdom STRING -i|info FILE [OPTIONS]
    or :   perl $0 -d|dir DIRECTORY -k|kingdom STRING -g|genome FILE -n|name Genus epithet [OPTIONS]

    purpose: 
	This script runs the snoStrip-pipeline on different kingdoms
	of Eukarya in order to detect novel snoRNA candidates. 
	

     
    options:

	-h|help            Prints this help message.
	                   [OPTIONAL]



      GENERAL REQUIRED OPTIONS:

	-d|dir             Output directory where all retrieved information will
                           be stored.
			   [REQUIRED]

	-k|kingdom         Specify the kingdom which shall be analyzed.
	                   \'pro\' ... Protostomia
                           \'deu\' ... Deuterostomia
                           \'fun\' ... Fungi
                           \'pla\' ... Plants
                           [REQUIRED]



      SPECIFY THE ORGANISMS TO BE SEARCHED [REQUIRED]:
        
        -i|info            File containing information of the organisms to be analyzed
			   file format: genome\tGenus epithet\t3-letter abbreviation (tab-separated)
                           Especially useful in case more than one organism shall be searched.

        provide information for a single organism:
        (in this case, option -i can be skipped)

        -g|genome          Genome source to search for potential snoRNAs.

        -n|name            Name of the organism: 'Genus epithet'.


 
      SNORNA FAMILY OPTIONS (OPTIONAL, default: analyse ALL snoRNA families):

        -c|clan            Single snoRNA families to search for, e.g.,
                           -c CD_12 HACA_100  [analyze CD_12.fa and HACA_100.fa]
                           -c CD_1* HACA*     [analyze CD_10.fa, CD_11.fa, etc. and ALL HACA families] 



      FURTHER OPTIONAL SETTINGS:

        -j|infernal        Incorporate infernal in the homology-based 
                           search procedure.

        -q|quiet           Suppress unnecessary output.

        -f|force           Force to search for novel snoRNA candidates in every
	                   organism, even though there are already species-specific
                           sequences in the current family.


			   
      OPTIONAL TARGET PREDICTION:

        --targets          Enable target prediction and make use of the target RNAs
                           that were shipped with snoStrip.
                           Location: $path

        Change these default directories to own settings:

        -s|sequences       Directory of targetRNAs.

	-a|alignments      Directory of targetRNA alignments.

	-p|profiles        Directory of targetRNA profiles.



FIN


print STDERR $usage;


}
