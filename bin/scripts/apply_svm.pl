#!/usr/bin/env/ perl

use Cwd 'abs_path';

BEGIN{
    my $path = abs_path($0);
    $path =~ s/scripts.*/packages\//;
    push @INC, $path;
}


use strict;
use Getopt::Std;
use FileHandle;
use warnings;
use Data::Dumper;
use CONFIG;



my %opts;
getopt('mrtsPUSo', \%opts);

my $RNASNOOP = $CONFIG::RNASNOOP;


#h the homolog file
#n the ncRNAdb file
#p path where the genome files are
#b 1, both mRNA are folded for RNAup, else only the longest one
if( defined($opts{t}) && defined($opts{s}) && defined($opts{o}) && defined($opts{r}) && defined($opts{m})){ 
    compute_target();
}




sub readfile
{
	my $name=shift;
	open(TMP, "$name") || die "cant open $name\n";
	my @array=<TMP>;
	close(TMP);
	return \@array;
}



sub parse_structure
{
  my $structure = shift;
  my $gap_left;
  my $gap_right;
  my $U_gap;
  my $U_near_gap_l;
  my $U_near_gap_r;
  my $b_i_gap;
  my $i_t_gap;
  my $t_i_gap;
  my $i_b_gap;
  my $stem_ass;
  my $stem_max_ass;
  my $stem_length;
  my $count=0;
  ###################
  $structure=~/^([^\|]+)\|/;
  my $temp=$1;
  while($temp=~/<\.</gi){
    $count++;
  }
  $gap_left=$count;
  $count=0;
   ####################
  if($structure=~/<\.\|\.</){
    $U_gap=3;
  }else{
    $U_gap=2;
  }
  #################
  $structure=~/\|\.([^\&]+)/;
  $temp=$1;
  while($temp=~/<\.</gi){
    $count++;
  }
  $gap_right=$count;
  #################
  if($structure=~/<\.<\.{0,1}\|/){
    $U_near_gap_l =1 ;
  }
  else{
    $U_near_gap_l =0 ;
  }
  #################
  if($structure=~/\|\.<\.</){
    $U_near_gap_r =1 ;
  }
  else{
    $U_near_gap_r =0 ;
  }
  ##################
  $structure=~/&([^>]+>)/;
  $temp=$1;
  if(defined($1)){
    if($temp=~/(\(\.*\>)/){
      $b_i_gap = length($1)-2;
    }
    else{
      #    print $structure," ",$temp,"\n";exit(0);
      $b_i_gap="100";
    }
  }
  else{
    $b_i_gap="100";
  }
  ######################
  $structure=~/>(\.*)\(/g;
  $temp=$1;
  $i_t_gap = length($temp);
  #########################
  $structure=~/>\.*(\(.+\))\.*>/;
  $temp=$1;
  $stem_length=length($temp);
  if($stem_length==0){
    print $structure," ",$temp,"\n";
    exit(0);
  }
  my @temp_array;
  my $temp_pos;
  my $position = -2;

  while(($position)){
    $position=index($temp,"(",$position);
    $position++;
    if($position>0){
      push @{$temp_pos},$position;
    }

  }
  $position=-2;
  while(($position)){
    $position=index($temp,")",$position);
    $position++;
    if($position>0){
      push @{$temp_pos},$position;
    }
  }
  $stem_ass=0;
  $stem_max_ass=0;
  for(my $i=0; $i<$#$temp_pos/2; $i++){
    $stem_max_ass = $stem_max_ass > abs( $temp_pos->[$i+1]-$temp_pos->[$i] - ($temp_pos->[$#$temp_pos-$i] - $temp_pos->[$#$temp_pos-$i-1])) ? $stem_max_ass : abs( $temp_pos->[$i+1]-$temp_pos->[$i] - ($temp_pos->[$#$temp_pos-$i] - $temp_pos->[$#$temp_pos-$i-1]));
    $stem_ass += $temp_pos->[$i+1]-$temp_pos->[$i] - ($temp_pos->[$#$temp_pos-$i] - $temp_pos->[$#$temp_pos-$i-1]);
  }
  ############################
  $structure =~/(\)\.*\>)/;
  $temp=$1;
  $t_i_gap = length($temp)-2;
  ############################
  $structure =~/(\>[^>]*)$/;
  $temp=$1;
  if($temp=~/\>([^\)]*)\)/){
    $temp=$1;
    $i_b_gap=length($temp);
  }
  else{
    $i_b_gap="-1";
  }

  return [$gap_left, $gap_right, $U_gap, $U_near_gap_l, $U_near_gap_r, $b_i_gap, $i_t_gap, $stem_length, $stem_max_ass, $stem_ass, $t_i_gap, $i_b_gap];

}
    


sub compute_target
{
  my $sno=parse_sequence($opts{s});
  #my $tar=parse_sequence($opts{t});
  my $RNAup;
  my $results;
  my $selected_features;
  my $flag_access;
  my @keyz=("LE",  #left interaction
	    "RE",  #right interaction
	    "DE",  #lower stem energy
	    "TE",  #upper stem energy
	    "SE",  #Total duplex energy
	    "XE",  #Stem and duplex energy
	    "YE",  #both stem and duplex energy
	    "ZE", #duplex and lower stem
	    "OE", #total opening energy for RNAsnoop -U
	    "dSE", #delta duplex
	    "dXE", #delta XE
	    "dYE",
            "dZE", #delta ZE
	    "energy", #total energy
	    "gap_left",
	    "gap_right",
	    "U_gap",
	    "U_near_gap_l",
	    "U_near_gap_r",
	    "b_i_gap",
	    "i_t_gap",
	    "stem_length",
	    "stem_max_ass",
	    "stem_ass",
	    "t_i_gap",
	    "i_b_gap",
	    "Ll",
	    "Rl",
	    "leng");
  my @inter;
  mkdir $opts{o} unless -d $opts{o}; #create directory output
  ###produce output and save it in snoop.out
#  print "RNAsnoop...\n";
    $flag_access=1;
  if(defined($opts{U})) {
    if(defined($opts{S})){
	#if(defined($opts{l})){
	    #my @targetRNAs = split(/,/,$opts{l});
	    #foreach my $tarRNA(@targetRNAs){
		#print STDERR "TargetRNA: ",$tarRNA,"\n";
		#`./RNAsnoop -t $tarRNA -s $opts{s} -U $opts{U}  -S $opts{S} -i 120 -j 11 -k 16 -e 40  -o -270 -p -170 -q -1090 -x -1370 -L 25 -f 1 | grep -v "no target found" >> $opts{o}/snoop.out; echo "RNAsnoop -t \$i -s $opts{s} -U $opts{U}  -S $opts{S} -i 120 -j 11 -k 16 -e 40  -o -270 -p -170 -q -1090 -x -1370 -L 25 -f 1"  >> $opts{o}/log.snoop`;
		#print STDERR "./RNAsnoop -t $tarRNA -s $opts{s} -U $opts{U}  -S $opts{S} -i 120 -j 11 -k 16 -e 40  -o -270 -p -170 -q -1090 -x -1370 -L 25 -f 1\n";
	    #}
	#}else{
#	`./RNAsnoop -t $opts{t} -s $opts{s} -U $opts{U}  -S $opts{S} -o -1000 -p -1000 -q -1000 -j 1 -k 100 -l 0 -x -1000 -y -1000  -e 40 -f 1  | grep -v "no target found" > $opts{o}/snoop.out`;
    `$RNASNOOP -t $opts{t} -s $opts{s} -U $opts{U} -S $opts{S} -i 120 -j 11 -k 16 -e 40  -o -270 -p -170 -q -1090 -x -1370 -L 25 -f 1 | grep -v "no target found" > $opts{o}/snoop.out`;
	#}
    }
    else{
      printf("If you are using RNAup accessibility profiles, then you should give the suffix of the files -S and the directory where they are located -U \n"); exit(0);
    }
  }
  elsif(defined($opts{P})) {
    `$RNASNOOP -t $opts{t} -s $opts{s} -P $opts{P} -i 120 -j 11 -k 16 -e 40  -o -270 -p -170 -q -1090 -x -1370 -L 25 -f 1 | grep -v "no target found" > $opts{o}/snoop.out`;
  }
  else{
    `$RNASNOOP -t $opts{t} -s $opts{s}  -q -1090 -j 11 -k 16 -o -340 -p -220 -l -280 -e 30 -b -320 -v 3 -f 1 -w 5 -L 25 -x -1370 -i 120 | grep -v "no target found" > $opts{o}/snoop.out`;
    $flag_access=0;
  }
  ###
  ###read_in
#  print "svm\n";
  #open(FH,"<$opts{o}/snoop.out") or;
  my $result = readfile("<$opts{o}/snoop.out");
  #print "FILE: ",$opts{o},"/snoop.out\n";
  open(SNOOP,">$opts{o}/snoop.res");
  open(SNOOP_CSV,">$opts{o}/snoop.csv");
  open(SNOOP_SVM,">$opts{o}/snoop.svm");
  print SNOOP_CSV "\n";
  print SNOOP "\n";
  my $previous=();
  #really not nice and sureley not working in all cases, just a fast hack
  my $hack = `grep '>' $opts{s}`;
  #$hack =~ /^>([^_]+_[^_]+)/;
  $hack =~ /^>([^\b]+)$/;
  my $query = $1;
  chomp($query);

  my $target;
  for(my $i =0; $i<scalar(@{$result});$i += 2){
      my $previous = $result->[$i];
  #while($previous=<FH>){
    #if(!$results){last;} 
    #print $previous,$results,"\n";
    if($previous=~/^>/){
      #chomp($results);
      chomp($previous);
      #$query=$previous;
      $target=$previous;
      my @info=split(/_/,$previous);
      $target=~s/(>|\n|\s+)//g;
      $i--;
      #$query=~s/(>|\n|\s+)//g;
    }
    else{
	$results=$result->[$i+1];
	chomp($previous);
	my $hash;
	$previous=~/^([^\s]+)/;
	chomp($results);
	$hash->{struct} = $1;
	#print "REGEXP ",$previous,"\n";
	$previous=~/\s+([^,]+),([^\s^\;]+)\s*;([^:]+):([^,]+)\,([^\s]+)\s+\(([^=]+)=\s+([^\+]+)\+([^\+]+)\+([^\+]+)\+([^\+]+)\+([^\+]+)/;
	# print $previous,"\n";exit(0);
	$hash->{real_target}=$target;
	$hash->{query}=$query;
	$hash->{u} = $3;
	$hash->{energy} = $6;
	$hash->{LE} = $7;
	$hash->{RE}=$8;
	$hash->{DE}=$9;
	$hash->{TE}=$10;
	$hash->{OE}=0;
	if($flag_access){
	$hash->{OE}=$11;
	}
	$hash->{t_end} = $2;
	$hash->{t_begin}=$1;
	chomp($hash->{t_begin});
	chomp($hash->{u});
	$hash->{t_begin}=~s/\s//g;
	$hash->{t_end}=~s/\s//g;
	$hash->{q_begin}=$4;
	$hash->{q_end}=$5;
	$hash->{q_begin}=~s/\s//g;
	$hash->{q_end}=~s/\s//g;
	$hash->{leng} = $hash->{t_end} - $hash->{t_begin};
	$hash->{seq} = $results;
	$hash->{SE} = $hash->{RE} + $hash->{LE};
	$hash->{XE} = $hash->{SE} + $hash->{DE};
	$hash->{YE} = $hash->{XE} + $hash->{TE};
	$hash->{ZE} = $hash->{SE} + $hash->{TE};
	$hash->{u}=~s/\s//g;
	$hash->{Ll} = $hash->{u}-$hash->{t_begin}+1;
	$hash->{Rl} = $hash->{t_end}-$hash->{u}+1;
	$hash->{dSE} = $hash->{RE} + $hash->{LE} + $hash->{OE};
	$hash->{dXE} = $hash->{dSE} + $hash->{DE};
	$hash->{dYE} = $hash->{dXE} + $hash->{TE};
	$hash->{dZE} = $hash->{dSE} + $hash->{TE};
	my $temp=parse_structure($hash->{struct});
	($hash->{gap_left},$hash->{gap_right},$hash->{U_gap},
	 $hash->{U_near_gap_l},$hash->{U_near_gap_r},$hash->{b_i_gap},
	 $hash->{i_t_gap},$hash->{stem_length},$hash->{stem_max_ass},
	 $hash->{stem_ass},$hash->{t_i_gap},$hash->{i_b_gap})=@{$temp};
	print SNOOP $hash->{query},"#",$hash->{real_target},"#",$hash->{u},"#",$previous,"#",$results,"\n";
	print SNOOP_CSV $hash->{query},",",$hash->{real_target},",",$hash->{u};
	print SNOOP_SVM "1";
	my $count_key=1;
	foreach my $keys (@keyz){
	    $hash->{$keys}=~s/\s//g;
	    print SNOOP_CSV ",",$hash->{$keys};
	    print SNOOP_SVM " $count_key:$hash->{$keys}";
	    $count_key++;
	}
	print SNOOP_CSV "\n";
	print SNOOP_SVM "\n";
    }
  }
#  close(FH);
  close(SNOOP);
  close(SNOOP_CSV);
  close(SNOOP_SVM);
  ##################now apply the svm on it.
  ######1) scale
  my @scaled = `./svm-scale -r $opts{r} $opts{o}/snoop.svm`; # > $opts{o}/snoop.svm.scale`;

  open(SNOOP_SCALE,">$opts{o}/snoop.svm.scale");
  my $dummy_count=2;
NE:  foreach my $line (@scaled){
    $dummy_count++;
    my @splited = split(/\s+/,$line);
    if(@splited < 29){
      print SNOOP_SCALE "0 1:0\n";
      next;
    }
    if($flag_access){
      print SNOOP_SCALE $splited[0];
      foreach my $features ((3,6,7,12,14,16,17,21,22,24,25,26)){
	  my $split = $splited[$features];
	  $split =~ s/^\d+://g;
	if(($split>=-1 && $split<=1)){
	  print SNOOP_SCALE " ",$splited[$features];
	}
	else{print SNOOP_SCALE " ",$features,":",10000,"\n";next NE;}
	if(!defined($splited[$features])){
	  print $features,"\n";
	  print $dummy_count;
	  print $line,"\n";
	  exit(0);
	}
      }
      print SNOOP_SCALE "\n";
    }
    else{
      print SNOOP_SCALE $splited[0];
      foreach my $features ((16,17,21,22,25)){
	print SNOOP_SCALE " ",$splited[$features];
      }
      print SNOOP_SCALE "\n";
    }
  }
  ######2) predict
#  print "svm-predict -b 1 $opts{o}/snoop.svm.scale $opts{m} $opts{o}/snoop.pred\n";
  `./svm-predict  -b 1 $opts{o}/snoop.svm.scale $opts{m} $opts{o}/snoop.pred `;
  ######3) sort
  open(SNOOP_SELECTED, ">$opts{o}/snoop.selected");
  my @sorted = `paste  $opts{o}/snoop.pred  $opts{o}/snoop.res | grep -E "^1" | sort -nrk 3`;
  if(scalar @sorted == 0){
      print "SVM did evaluate everything as unprobable\n";
      exit(0);
  } else{
      foreach my $line (@sorted){
	  #print STDERR "LINE: ",$line,"\n";
	  my @splited = split(/[\t#]/, $line);
	  $splited[0]=~s/ /./g;
	  my $input = ">".$splited[0].$splited[1]."\n>".$splited[2]."\n".$splited[4]."\n".$splited[5];
	  print SNOOP_SELECTED $input;
      }
      close(SNOOP_SELECTED);
      if($flag_access){
	  #`cat $opts{o}/snoop.selected | ./RNAsnoop -I -U $opts{U} -S $opts{S} -O $opts{o}`;
	  print `cat $opts{o}/snoop.selected`;
      }
      else{
	  #`cat $opts{o}/snoop.selected | ./RNAsnoop -I -O $opts{o}`;
     print `cat $opts{o}/snoop.selected`;
      }
  }
    
  
 
  ######4) output final result
}



sub read_RNAup
{
	my $RNAup_data;
	my $directory=shift;
	my @files = <$directory/*.dat>;
	my $name;
	foreach my $file (@files){
		$file=~/^$directory(\/)+([^_]+)/;
		$name=$2;
		my $data=parse_RNAup($file);
		$RNAup_data->{$name}=$data;
	}
	return $RNAup_data;
}



sub parse_RNAup
{
    my $name=shift;
    my $data;
    my $file=readfile($name);#file containing the informationc
    my @temp=split(/\s+/, $file->[4]);#line containing the header of accessibility information
    my @sizes=@temp[3..$#temp];
	for(my $i=5; $i<=$#$file; $i++){
	    $file->[$i]=~s/\s+/\t/g;
	    my @access=split(/\s+/, $file->[$i]);
		$access[1]=~s/\s+//;
	    for(my $j=1; $j<=$#sizes; $j++){
		$sizes[$j-1]=~s/(u|S)//g;
		if($access[$j+1] eq 'NA'){
		    $access[$j+1] = 22;
		}
		$data->{$sizes[$j-1]}->{$access[1]}=$access[$j+1];
		#print $sizes[$j-1]," ",$access[1]," ",$access[$j+1],"\n";
		}
	}
    
    return $data;
}


sub parse_RNAsnoop
 {
  my $results=shift;
  my $master_array;
  my $rna;
  my $sno;
  my $target;
  my $query;
  my $u;
  my $accessibility=0;
  for(my $i=0; $i<$#$results; $i=$i+2){
    if($results->[$i]=~/^>/ && $results->[$i+1]=~/^>/)
      {
	chomp($results->[$i]);
	chomp($results->[$i+1]);
	my @info=split(/_/,$results->[$i+1]);
	$query=$results->[$i];
	$target=$results->[$i+1];
	$target=~s/(>|\n|\s+)//g;
	$query=~s/(>|\n|\s+)//g;
      }
    else{
      chomp($results->[$i]);
      if($results->[$i]=~/^no/){
	$i--;
	next;
      }
      my $hash;
      $results->[$i]=~/^([^\s]+)/;
      $hash->{struct} = $1;
      if($accessibility){
	$results->[$i]=~/\s+([^,]+),([^\s^\;]+)\s*;([^:]+):([^,]+)\,([^\s]+)\s+\(([^=]+)=\s+([^\+]+)\+([^\+]+)\+([^\+]+)\+([^\+]+)\+([^\+]+)/;
      }
      else{
	$results->[$i]=~/\s+([^,]+),([^\s^\;]+)\s*;([^:]+):([^,]+)\,([^\s]+)\s+\(([^=]+)=\s+([^\+]+)\+([^\+]+)\+([^\+]+)\+([^\+]+)/;
      }
      $hash->{real_target}=$target;
      $hash->{query}=$query;
      $hash->{real_u}=$u;	
      $hash->{u} = $3;
      $hash->{energy} = $6;
      $hash->{LE} = $7;
      $hash->{RE}=$8;
      $hash->{DE}=$9;
      $hash->{TE}=$10;
      if($accessibility){
	$hash->{OE}=$11;
      }
      else{
	$hash->{OE}=0;
      }
      $hash->{t_end} = $2;
      $hash->{t_begin}=$1;
      chomp($hash->{t_begin});
      chomp($hash->{u});
      $hash->{t_begin}=~s/\s//g;
      $hash->{t_end}=~s/\s//g;
      $hash->{q_begin}=$4;
      $hash->{q_end}=$5;
      $hash->{q_begin}=~s/\s//g;
      $hash->{q_end}=~s/\s//g;
      $hash->{leng} = $hash->{t_end} - $hash->{t_begin};
      $hash->{seq} = $results->[$i+1];
      push @{$master_array}, $hash;
    }
  }
  return $master_array;
}



sub parse_pseudo
{
  my $results;
  my $file=readfile($opts{p});
  foreach my $line (@{$file}){
    next if ($line=~/#/);
    my $hash;
    ($hash->{sno},$hash->{stem}, $hash->{target}, $hash->{u},$hash->{type})=split(/\s+/,$line);
    $results->{$hash->{target}}->{$hash->{u}}=1;
  }
  return $results;
}



sub parse_interaction
{
  my $file = readfile($opts{i});
  my $results;
  my $temp=();
  foreach my $line (@{$file}){
    next if ($line=~/#/);
    my $hash;
    ($hash->{sno},$hash->{stem}, $hash->{target}, $hash->{u},$hash->{type})=split(/\s+/,$line);
    push @{$results->{target}}, $hash;
  }
  return $results;
}



sub parse_sequence
{
	my $name=shift;
	my $file=readfile($name);
	my $hash;
	for(my $i=0; $i <=$#$file; $i=$i+2){
		chomp($file->[$i]);
		chomp($file->[$i+1]);
		$file->[$i]=~s/(>|\n|\s+)//g;
		$hash->{$file->[$i]}=$file->[$i+1];
	}
	return $hash;
}
		
	
