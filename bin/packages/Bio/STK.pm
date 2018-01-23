#!/usr/bin/perl -w
#####################################################################
#
# Copyright (c) 2009, Jon Hill, jon.hill@imperial.ac.uk
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#     * Redistributions of source code must retain the above copyright notice,
#       this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation
#       and/or other materials provided with the distribution.
#     * Neither the name of the University of Edinburgh or University of Oxford
#       the names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#####################################################################
#
# These functions and variables are used in a number of the Supertree
# Tool Kit executable scripts
#
#####################################################################

=head1 STK - Supertree Tool Kit

=head2 Description

B<STK> is a collection of tools to prepare data ready for supertree 
contruction with PAUP* or similar. It relies on a file-based database
based on tree and XML files. Functions include finding all tree files,
xml files and summarising data. STK also comes with a number of scripts
ready for use.

=head2 Requires

Perl 5.004, XML::Simple, Bio::NEXUS, File::Find, Carp, File::Copy;

=head2 Feedback

All feedback (bugs, feature enhancements, etc.) are all greatly appreciated. 

=head2 Authors

Jon Hill (jon.hill@imperial.ac.uk) and Katie Davis (katie.davis@oum.ox.ac.uk)

=head2 Methods

=cut

package Bio::STK;

use strict;
use XML::Simple;
use File::Find;
use File::Spec::Functions;
use File::Copy;
use Bio::NEXUS;
use Carp;
# Useful for debugging
#$SIG{__DIE__} =  \&confess;
#$SIG{__WARN__} = \&confess;

# ExtUtils::MakeMaker reads package global $VERSION
use vars qw($VERSION $AUTOLOAD);
$VERSION = "0.1.2";

# Character types set up a an enumeration (kind of)
# can use STK->MOLECULAR for example.
use constant {
    MOLECULAR     => 'Molecular',
    MORPHOLOGICAL => 'Morphological',
    BEHAVIOURAL   => 'Behavioural',
    OTHER         => 'Other',
};

# A global array so you can loop over the character types
# This should be added to if a new type is added
# This should be replaced with a XML schema in the near future
our @character_types;
$character_types[0] = "Molecular";
$character_types[1] = "Morphological";
$character_types[2] = "Behavioural";
$character_types[3] = "Other";

our @_files;    # needed for File::Find routines. Must be undef'd before use.

######################################################################
#                                                                    #
#                          FILE FUNCTIONS                            #
#                                                                    #
#  Functions to find, load and store data to and from files          #
#                                                                    #
#  Names should start with find_, read_, save_                       #
#                                                                    #
######################################################################

=head3 find_tree_files

 Title   : find_tree_files
 Usage   : my @files = Bio::STK::find_tree_files($dir);
 Function: Finds all tree files in specified directory and all sub-directories
 Returns : Array of complete filenames
 Args    : $dir, else croaks.
 
=cut

sub find_tree_files {

    my ($dir) = @_;
    undef @_files;

    # check if directory exists
    if ( -d $dir ) {
        find( \&Bio::STK::_find_trees, $dir );
    }
    else {
        croak("Error - specified directory does not exist: $dir\n");
    }

    @_files = sort @_files;
    return @_files;

}

=head3 find_xml_files

 Title   : find_xml_files
 Usage   : my @files = Bio::STK::find_xml_files($dir);
 Function: Finds all XML files in specified directory and all sub-directories
 Returns : Array of complete filenames
 Args    : $dir, else croaks.

=cut

sub find_xml_files {

    my ($dir) = @_;
    undef @_files;

    # check if directory exists
    if ( -d $dir ) {
        find( \&Bio::STK::_find_xml, $dir );
    }
    else {
        croak("Error - specified directory does not exist: $dir\n");
    }

    @_files = sort @_files;
    return @_files;

}

=head3 read_tree_file

 Title   : read_tree_file
 Usage   : my $tree = read_tree_file($file);
 Function: Reads in a correctly formatted NEXUS tree file
 Returns : An array of Newick tree formatted strings;
 Args    : $file, else croaks.

=cut

sub read_tree_file {

    my ($file) = @_;

    # _load_nexus check if the file exists
    my $nexus = _load_nexus($file);

    my @trees = _get_tree_array($nexus);

    return @trees;

}

=head3 read_xml_file

 Title   : read_xml_file
 Usage   : my $xml = read_xml_file($file);
 Function: Reads in an XML file
 Returns : The data structure from XML::Simple
 Args    : $file, else croaks. 

A sample data structure is:

 'Notes' => [
           '2,430bp of the mitochondrial genes ND2, ND3 and cyt b.'
         ],
 'Analysis' => [
              {
                'Type' => [
                            'MP'
                          ]
              }
            ],
 'Characters' => [
                {
                  'Other' => [
                               {}
                             ],
                  'Behavioural' => [
                                     {}
                                   ],
                  'Morphological' => [
                                       {}
                                     ],
                  'Molecular' => [
                                   {
                                     'Type' => [
                                                 'cytb',
                                                 'ND2',
                                                 'ND3'
                                               ],
                                     'number' => '3'
                                   }
                                 ]
                }
              ],
 'Taxa' => [
          {
            'List' => [
                        'Campylorhamphus falcularius',
                        'Campylorhamphus procurvoides',
                        'Campylorhamphus trochilirostris',
                        'Dendrexetastes rufigula',
                        'Dendrocolaptes certhia',
                        'Glyphorynchus spirurus',
                        'Hylexetastes perrotii',
                        'Lepidocolaptes albolineatus',
                        'Lepidocolaptes angustirostris',
                        'Lepidocolaptes fuscus',
                        'Nasica longirostris',
                        'Sittasomus griseicapillus',
                        'Xiphocolaptes promeropirhynchus',
                        'Xiphorhynchus erythropygius',
                        'Xiphorhynchus flavigaster',
                        'Xiphorhynchus guttatus',
                        'Xiphorhynchus kienerii',
                        'Xiphorhynchus lachrymosus',
                        'Xiphorhynchus obsoletus',
                        'Xiphorhynchus ocellatus',
                        'Xiphorhynchus pardalotus',
                        'Xiphorhynchus picus',
                        'Xiphorhynchus spixii',
                        'Xiphorhynchus susurrans',
                        'Xiphorhynchus triangularis'
                      ],
            'number' => '25',
            'fossil' => 'none'
          }
        ],
 'Source' => [
            {
              'Year' => [
                          '2000'
                        ],
              'Publisher' => [
                               {}
                             ],
              'Volume' => [
                            '119'
                          ],
              'Editor' => [
                            {}
                          ],
              'Pages' => [
                           '621-640'
                         ],
              'Booktitle' => [
                               {}
                             ],
              'Journal' => [
                             'Auk'
                           ],
              'Title' => [
                           'Molecular systematics and the role of the "Varzea"-"Terra-firme" ecotone in
 the diversification of Xiphorhynchus woodcreepers (Aves: Dendrocolaptidae).'
                         ],
              'Author' => [
                            'Bar, A.'
                          ]
            }
          ],
 'TreeFile' => [
              'tree1.tre'
            ]

You can access by using structures like this: C<@{$xml_contents->{Taxa}->[0]->{List}}>.

=cut

sub read_xml_file {

    my ($file) = @_;

    # check if file exists
    unless ( -e $file ) {
        croak("Error - specified XML file does not exist: $file\n");
    }
    return XMLin( $file, ForceArray => 1 );

}

=head3 save_tree_file

 Title   : save_tree_file
 Usage   : save_tree_file($treefile, \@tree_string, \@tree_names);
 Function: Save an array of tree strings to file. Will overwrite any existing file without checking
 Returns : 1 if sucessful, 0 if not. CHECK THIS
 Args    : $treefile, @tree_string, @tree_names

This function does not check if the treefile already exists - this is the responsibility of
the calling program. 

=cut

sub save_tree_file {

    my $file  = shift;
    my $trees = shift;
    my $names = shift;

    my $nexus_obj = _create_nexus_object( $trees, $names );
    my $t = $nexus_obj->get_block('Trees')->get_trees;

    # remove node names
    foreach my $tree ( @{$t} ) {
        foreach my $node ( @{ $tree->get_nodes } ) {
            if ( !$node->is_otu() ) {
                $node->set_name('');
            
            }
        }
    }

    # Bio::NEXUS spits out a very long and complex error message
    # when a file is not availble to write to.
    # Hide this from the user by checking it can be written to!
    open my $fh, ">$file" or die "Error saving tree file to $file: $!";
    close $fh;
    $nexus_obj->write($file);
    return 1;

}


=head3 save_tree_newick

 Title   : save_tree_newick
 Usage   : save_tree_newick($treefile, \@tree_string, \@tree_names)
 Function: Save an array of tree strings to file. Will overwrite any existing file without checking
 Returns : 1 if sucessful, 0 if not. CHECK THIS
 Args    : $treefile, \@tree_string, \@tree_names

This function does not check if the treefile already exists - this is the responsibility of
the calling program. The purpose of this function is to save a tree file in "raw" Newick format, 
i.e. without any other BLOCKS apart from TREE. This is the simplist form of the treefile and hence
the most compatible with other software.

=cut

sub save_tree_newick {

  # whilst we would like to use Bio::NEXUS to do this, I cannot
  # figure out a way to save trees without TAXA blocks
  # hence this function just saves to an ASCII file.

  my $file  = shift;
  my $trees = shift;
  my $names = shift;

  # Bio::NEXUS spits out a very long and complex error message
  # when a file is not availble to write to.
  # Hide this from the user by checking it can be written to!
  open my $fh, ">$file" or die "Error saving tree file to $file: $!";
 
  print $fh "#NEXUS\n\nBEGIN TREES;\n\n";

  my $i = 0;
  my $name;
  foreach my $tree ( @{$trees} ) {
    if ( !defined $$names[$i] ) {
      $name = "tree_" . ( $i + 1 );
    } else {
      $name = $$names[$i];
    }
    
    print $fh "tree $name = $tree;\n";
    $i++;
  }
 
  print $fh "\nEND;\n";
  close $fh;
  return 1;

}


=head3 save_xml_file

 Title   : save_xml_file
 Usage   : save_xml_file($xmlfile, $xml_hash);
 Function: Save a XML hash to file. Will overwrite any existing file without checking
 Returns : 1 if sucessful, 0 if not. CHECK THIS
 Args    : $xmlfile, $xml_hash

This function does not check if the xml file already exists - this is the responsibility of
the calling program.

$xml_hash is the same object as returned by L<read_xml_file>.

=cut

sub save_xml_file {

    my ( $file, $xml_hash ) = @_;

    open my $fh, ">$file" or die "Error saving file to $file: $!";
    XMLout( $xml_hash, OutputFile => $fh, XMLDecl => 1, RootName => 'Source' );
    close $fh;
    return 1;
}

=head3 read_taxa_file

 Title   : read_taxa_file
 Usage   : $nSubs = read_taxa_file($taxafile, \@old_taxa, \@new_taxa);
 Function: Read in a taxa substitution file, populating @old_taxa and @new_taxa.
 Returns : The number of substitutions.
 Args    : $taxafile, \@old_taxa, \@new_taxa

old_taxa contains each old taxa, with a corresponding entry in new_taxa, i.e.
 
 old_taxa[0] => new_taxa[0]
 old_taxa[1] => new_taxa[1]
 etc

new_taxa will contain a string which is a comma-seperated list of all taxa that should replace
the old_taxa

=cut

sub read_taxa_file {

    my $file = shift;
    my $old_taxa = shift;
    my $new_taxa = shift;

    # read in substitutions file and create arrays
    open( IN, "$file" ) or die "Read_Taxa_File: Can't read file\n";
    my $line;
    my $first_one = 0;
    my $subs = '';
    my $nSubs = 0;
    while ( defined( $line = <IN> ) ) {
        if ( $line =~ m/\s{1,}=\s{1,}/ ) {
            $nSubs++;
            # let's add the last one if there is one
            if ($first_one) {
                push( @{$new_taxa}, $subs );
            }

            #line has an equals sign so start a new higher taxa
            my @temp_line = split( /\s=\s/ , $line );
            # deal with old taxa
            # strip spaces
            $temp_line[0] =~ s/^\s+//;
            $temp_line[0] =~ s/\s+$//;
            # strip off any quotes - we must add them, but it's easier to strip then add than
            # check, then add
            $temp_line[0] =~ s/'//g;
            # might contain non-standard characters, quote these names
            $temp_line[0] =~ s/(\w*[=\+]\w*)/'$1'/g;

            # deal with the replacement taxa
            $temp_line[1] =~ s/^\s+//;
            $temp_line[1] =~ s/\s+$//;
            # might have spaces around the commas - remove them
            $temp_line[1] =~ s/\s*,\s*/,/g;
            # might contain non-standard characters, quote these names
            # strip off any quotes first - we must add them, but it's easier to strip then add than
            # check, then add
            $temp_line[1] =~ s/'//g;
            chomp $temp_line[1];
            $temp_line[1] =~ s/^\s+//;
            $temp_line[1] =~ s/\s+$//;
            $temp_line[1] =~ s/(,|^)(\w*[=\+]\w*)/$1'$2'/g;
            push( @{$old_taxa}, $temp_line[0] );
            $subs = $temp_line[1];
            $first_one = 1;
        }
        # a new line containing replacements, but no " = "
        else {
            # strip whitespace and add to current subs list
            my $temp_line = $line;
            chomp $temp_line;
            $temp_line =~ s/^\s+//;
            $temp_line =~ s/\s+$//;
            # empty line, keep going
            if (not $temp_line =~ m/\w/) {
                next;
            }

            # might have spaces around the commas - remove them
            $temp_line =~ s/\s*,\s*/,/g;
            # might contain non-standard characters, quote these names
            $temp_line =~ s/(,|^)(\w*[=\+]\w*)/$1'$2'/g;
            # does $subs end in a comma?, if so just add, if not
            # does this line start with one if not add one
            if ($subs =~ m/,$/) {
                $subs .= $temp_line;
            } else {
                if ($temp_line =~ m/^,/) {
                    $subs .= $temp_line;
                } else {
                    $subs .= "," . $temp_line;
                }
            }
        }
    }
    close IN;
    push( @{$new_taxa}, $subs );

    return $nSubs;
}


######################################################################
#                                                                    #
#                          DATA FUNCTIONS                            #
#                                                                    #
#  Functions that get data from a data structure or file             #
#                                                                    #
#  Names should start with X_from_Y, where X is the data required    #
#  and Y is the source, or get_                                      #
#                                                                    #
######################################################################

=head3 taxa_from_tree

 Title   : taxa_from_tree
 Usage   : my @taxa = taxa_from_tree($tree);
 Function: Obtains unique taxa from a tree file or tree string
 Returns : An array of taxa names or undef if error
 Args    : $tree, else croaks.

$tree can be a tree file or a NEXUS formatted tree string

=cut

sub taxa_from_tree {

    my ($input) = @_;
    my @taxlabels;
    my @_trees;

    # check if input is a file or tree string
    my $ret = _file_or_tree($input);
    if ( $ret == 1 ) {
        $_trees[0] = $input;
    }
    elsif ( $ret == 2 ) {
        @_trees = Bio::STK::read_tree_file($input);
    }
    else {

        # neither!
        return undef;
    }

    # Cannot figure out a way using Bio::NEXUS to get at
    # all the taxa names in a file with no TAXA block, but with multiple
    # trees that contain some, but not all, overlapping taxa
    # Therefore doing this the hard way...

    foreach my $tree (@_trees) {

        # we have a Newick string - extract name using reg exp magic
        # then use split to create array on ,
        # first let's remove ( & ), labels 
        # Note: Seems quicker to do these seperately, rather than in 1 regexp
        $tree =~ s/\(//g;
        $tree =~ s/\)//g;
        $tree =~ s/\:\w+(,|\))/$1/g;
        # check for any whitespace around commas
        $tree =~ s/\s*,\s*/,/g;
        #$tree =~ s/\'//g;

        push( @taxlabels, split( ',', $tree ) );
    }

    for ( my $i = 0; $i < @taxlabels; $i++ ) {
        $taxlabels[$i] =~ s/_/ /g;
    }

    # now get distinct names
    # The same taxa can be in multiple trees
    my %saw;
    undef %saw;
    @taxlabels = grep( !$saw{$_}++, @taxlabels );

    return @taxlabels;

}

=head3 get_taxa_list

 Title   : get_taxa_list
 Usage   : get_taxa_list($dir);
 Function: Gets a list of unique taxa from a directory (and subdirectories) of tree files.
 Returns : An array of unique names
 Args    : $dir, else croaks.

=cut

sub get_taxa_list {

    my ($dir) = @_;
    my @taxa;

    unless ( -d $dir ) {
        croak("Directory $dir does not exist\n");
    }

    #find all tree files first
    my @treefiles = Bio::STK::find_tree_files($dir);

    foreach my $file (@treefiles) {
        push( @taxa, Bio::STK::taxa_from_tree($file) );
    }

    # remove duplicates
    my %saw;
    undef %saw;
    @taxa = grep( !$saw{$_}++, @taxa );

    @taxa = sort(@taxa);

    return @taxa;

}

=head3 taxa_from_xml

 Title   : taxa_from_xml
 Usage   : my @taxa = taxa_from_xml($input)
 Function: Gets the taxa containedin XML file or XML hash
 Returns : An array of strings containing the taxa
 Args    : $input, else croaks. 

Note that input can be a filename or XML hash.

=cut

sub taxa_from_xml {

    my ($input) = @_;
    my @taxa;

    # check if input is a file or xml hash
    my $ret = Bio::STK::_file_or_xml($input);
    if ( $ret == 1 ) {
        @taxa = @{ $input->{Taxa}->[0]->{List} };
    }
    elsif ( $ret == 2 ) {

        # read_xml_file checks that the file exists
        my $xml_contents = read_xml_file($input);
        @taxa = @{ $xml_contents->{Taxa}->[0]->{List} };
    }
    else {

        # neither!
        return undef;
    }

    return @taxa;

}

=head3 get_analysis

 Title   : get_analysis
 Usage   : $analysis = get_analysis($xmlfile);
 Function: Gets the type of analysis used from an XML file
 Returns : A string containing the type of analysis used
 Args    : $xmlfile, else croaks.

=cut

sub get_analysis {

    my ($file) = @_;
    my $xml;

    # check if it's a file or data structure
    my $ret = Bio::STK::_file_or_xml($file);
    if ( $ret == 1 ) {
        $xml = $file;
    }
    elsif ( $ret == 2 ) {

        # read_xml_file checks that the file exists
        $xml = read_xml_file($file);
    }
    else {
        return undef;
    }

    return $xml->{Analysis}->[0]->{Type}->[0];

}

=head3 get_source_data

 Title   : get_source_data
 Usage   : %source = get_source_data($xmlfile);
 Function: Gets the tsource information from an XML file
 Returns : A hash containing the source information (journal, title, etc)
 Args    : $xmlfile, else croaks. 

The output from get_source_data returns a Hash with the following structure:

 my %expected = (year => '2000',
                volume => '119',
                pages => '621-640',
                journal => 'Auk',
                title => 'Molecular systematics and the role of the "Varzea"-"Terra-firme" ecotone in the diversification of Xiphorhynchus woodcreepers (Aves: Dendrocolaptidae).',
                author => 'Bar, A.',
                booktitle => 'Some title',
                editor => 'Some editors',
                publisher => 'Some publisher',
                );

Access as a normal Hash.

=cut

sub get_source_data {

    my ($input) = @_;
    my %source_data;

    my $xml_contents;

    # Check if it's a file or data structure
    my $ret = Bio::STK::_file_or_xml($input);
    if ( $ret == 1 ) {
        $xml_contents = $input;
    }
    elsif ( $ret == 2 ) {

        # read_xml_file checks that the file exists
        $xml_contents = read_xml_file($input);
    }
    else {
        return undef;
    }

    # for some reason XML::Simple gives an empty hash if the XML tags are
    # empty and a string if there is content. This was figured out by
    # trial and error - don't ask me to explain...
    if ( ref( $xml_contents->{Source}->[0]->{Title}->[0] ) ne "HASH" ) {
        $source_data{'title'} = $xml_contents->{Source}->[0]->{Title}->[0];
    }
    if ( ref( $xml_contents->{Source}->[0]->{Year}->[0] ) ne "HASH" ) {
        $source_data{'year'} = $xml_contents->{Source}->[0]->{Year}->[0];
    }
    if ( ref( $xml_contents->{Source}->[0]->{Publisher}->[0] ) ne "HASH" ) {
        $source_data{'publisher'} =
            $xml_contents->{Source}->[0]->{Publisher}->[0];
    }
    if ( ref( $xml_contents->{Source}->[0]->{Volume}->[0] ) ne "HASH" ) {
        $source_data{'volume'} = $xml_contents->{Source}->[0]->{Volume}->[0];
    }
    if ( ref( $xml_contents->{Source}->[0]->{Editor}->[0] ) ne "HASH" ) {
        $source_data{'editor'} = $xml_contents->{Source}->[0]->{Editor}->[0];
    }
    if ( ref( $xml_contents->{Source}->[0]->{Pages}->[0] ) ne "HASH" ) {
        $source_data{'pages'} = $xml_contents->{Source}->[0]->{Pages}->[0];
    }
    if ( ref( $xml_contents->{Source}->[0]->{Booktitle}->[0] ) ne "HASH" ) {
        $source_data{'booktitle'} =
            $xml_contents->{Source}->[0]->{Booktitle}->[0];
    }
    if ( ref( $xml_contents->{Source}->[0]->{Journal}->[0] ) ne "HASH" ) {
        $source_data{'journal'} = $xml_contents->{Source}->[0]->{Journal}->[0];
    }
    if ( ref( $xml_contents->{Source}->[0]->{Author}->[0] ) ne "HASH" ) {
        $source_data{'author'} = $xml_contents->{Source}->[0]->{Author}->[0];
    }

    return %source_data;

}

=head3 get_treefile

 Title   : get_treefile
 Usage   : get_treefile($xmlfile);
 Function: Get the tree files that this XML file is for
 Returns : A string containing the B<relative> path of the tree file
 Args    : $xmlfile, else croaks. 

=cut

sub get_treefile {

    my ($input) = @_;
    my $treefile;

    my $xml_contents;

    # Check if it's a file or data structure
    my $ret = Bio::STK::_file_or_xml($input);
    if ( $ret == 1 ) {
        $xml_contents = $input;
    }
    elsif ( $ret == 2 ) {

        # read_xml_file checks that the file exists
        $xml_contents = read_xml_file($input);
    }
    else {
        return undef;
    }

    $treefile = $xml_contents->{TreeFile}->[0];

    return $treefile;
}

=head3 get_characters

 Title   : get_characters
 Usage   : get_characters($xmlfile);
 Function: Get the chracters used in this analysis
 Returns : A Hash containing the character types
 Args    : $xmlfile, else croaks.

The output from get_characters returns a Hash with the following structure:

 my %expected = (Molecular => [ "cytb", "ND2", "ND3" ],
                Morphological => [ "feathers" ],
                Other => [ "stuff" ],
                Behavioural => [ "stuff", "more stuff" ],
                );

Access as a normal Hash and then array notation, i.e. C<$expected{molecular}[2]>;

=cut

sub get_characters {

    my ($file) = @_;
    my $xml_contents;

    my %characters;

    my $ret = Bio::STK::_file_or_xml($file);
    if ( $ret == 1 ) {
        $xml_contents = $file;
    }
    elsif ( $ret == 2 ) {

        # read_xml_file checks that the file exists
        $xml_contents = read_xml_file($file);
    }
    else {
        return undef;
    }

    # check if some of the characters exist
    for my $type (@character_types) {
        if (defined(
                $xml_contents->{Characters}->[0]->{$type}->[0]->{'number'}
            )
            )
        {
            $characters{$type} =
                $xml_contents->{Characters}->[0]->{$type}->[0]->{'Type'};
        }
    }

    return %characters;

}

######################################################################
#                                                                    #
#                       VALIDITY FUNCTIONS                           #
#                                                                    #
#  Functions to check that a data structure or file is valid         #
#                                                                    #
#  Names should start with check_                                    #
#                                                                    #
######################################################################

=head3 check_tree_file

 Title   : check_tree_file
 Usage   : check_tree_file($file);
 Function: Test that a tree file is a valid NEXUS file
 Returns : 1 on success, 0 on failure - REMEMBER TO CHECK THIS!
 Args    : $file, else croaks.

=cut

sub check_tree_file {

    my ($file) = @_;
    my @tree_loaded;

    unless ( -e $file ) {
        croak("File $file not found\n");
    }

    eval { @tree_loaded = Bio::STK::read_tree_file($file); };

    if (@tree_loaded) {
        return 1;
    }

    return 0;

}

######################################################################
#                                                                    #
#                   DATA MANIPULATION FUNCTIONS                      #
#                                                                    #
#  Functions that edit data structures                               #
#                                                                    #
######################################################################

=head3 replace_taxon_tree

 Title   : replace_taxon_tree
 Usage   : replace_taxon_tree(\@treeString, $old_taxon, \@new_taxa);
 Function: Swap a taxon in a tree with a new one within an array of tree strings. Leave new_taxon blank to remove old_taxon
 Returns : 1 on success, 0 on failure. \@treeString contains new array of trees
 Args    : \@treeString, $old_taxon, \@new_taxa. new_taxa may be null, in which case old_taxon is removed.
           @new_taxon is an array of names - if this is longer than 1 element, then the old_taxon is
           replaced with a polytomy
 
Note that this function uses "pass by reference" - remember to add the slashes

=cut

sub replace_taxon_tree {

    my $treeStrings = shift;
    my $old_taxon   = shift;
    my $new_taxa    = shift;

    # check for spaces
    $old_taxon =~ s/ /_/g;
    
    if ( defined($new_taxa) ) {
        for ( my $i = 0; $i < @$new_taxa; $i++ ) {
            $$new_taxa[$i] =~ s/ /_/g;
        }
    }


    if ( !defined $$new_taxa[0] ) {

        # remove old_taxon
        # if this case insensitive? Doubt it...
        my $nexus = _create_nexus_object($treeStrings);
        my $t_old_taxon = $old_taxon;
        $t_old_taxon =~ s/'//g;
        {
            # temp hide warnings as exclude_otus can sometimes give an
            # uninit var warning
            $SIG{__WARN__} = sub {};
            $nexus = $nexus->exclude_otus( [$t_old_taxon] );
            $SIG{__WARN__} = 'DEFAULT';
        }
        @{$treeStrings} = _get_tree_array($nexus);


    } else {

        # swap old for new
        
        # uniqu-ify incoming names
        my %seen = ();
        @$new_taxa = grep { ! $seen{$_} ++ } @$new_taxa;

        # we need to check if we want to replace with polytomy here
        my $length = @$new_taxa;

        my $i = 0;
        # we're not using $tree below as $treeStrings is a ref
        # and we pass this back, hence we need to alter treeStrings directly
        foreach my $tree ( @{$treeStrings} ) {
            # check for extra whitespace in NEXUS tree
            $$treeStrings[$i] =~ s/\s*,\s*/,/g;

            # construct string to replace taxon with polytomy
            my $new_string = "";
            foreach my $taxon (@$new_taxa) {
                unless ( Bio::STK::tree_contains( $taxon, $$treeStrings[$i] ) )
                {
                    $new_string = $new_string . "," . $taxon;
                }
            }
	        
            # strip first comma I've just added
            if ( length($new_string) > 1 ) {
                $new_string = substr( $new_string, 1 );
            }
            if ( length($new_string) == 0 ) {
                my $nexus = _create_nexus_object([$$treeStrings[$i]]);
                # Bio::NEXUS doesn't like quoted taxa so remove them
                my $t_old_taxon = $old_taxon;
                $t_old_taxon =~ s/'//g;
                {
                    # temp hide warnings as exclude_otus can sometimes give an
                    # uninit var warning
                    $SIG{__WARN__} = sub {};
                    $nexus = $nexus->exclude_otus( [$t_old_taxon] );
                    $SIG{__WARN__} = 'DEFAULT';
                }
                my @output = _get_tree_array($nexus);
                $$treeStrings[$i] = $output[0];
            } else {
                # do replacement on this string
                # note - take care with delimiters
                my $replace = quotemeta $old_taxon;
                $$treeStrings[$i] =~ s/(\(|,)$replace(\)|,)/$1$new_string$2/;
            }
            $i++;
        }
    }

    return 1;

}

=head3 replace_taxon_xml

 Title   : replace_taxon_xml
 Usage   : replace_taxon_xml($xml_data, $old_taxon, @new_taxon);
 Function: Swap a taxon in a XML data structure with a new one. Leave new_taxon blank to remove old_taxon
 Returns : 1 on success, 0 on failure. $xml_data is altered.
 Args    : $xml_data, $old_taxon, @new_taxon. new_taxon may be null, in which case old_taxon is removed.
           @new_taxon is an array of names - if this is longer than 1 element, then the old_taxon is
           replaced with a polytomy

=cut

# ToDo:
#
# - If we remove a fossil taxa, do we check the "fossil" part of the xml? Nope! This a job to
#   add when we link up to taxanomic databases
#

sub replace_taxon_xml {

    my ( $xml, $old_taxon, @new_taxon ) = @_;

    my @taxa = @{ $xml->{Taxa}->[0]->{List} };
        
    my $i = 0;

    # _ instead of spaces
    $old_taxon =~ s/_/ /g;
    # remove any ' - not needed in XML
    $old_taxon =~ s/'//g;
    $old_taxon = quotemeta $old_taxon;
    foreach my $taxon (@taxa) {
        # remove any ' - not needed in XML
        $taxon =~ s/'//g;  
        if ( $taxon =~ m/^$old_taxon$/ ) {

            # remove old taxon
            splice @taxa, $i, 1;
            $i--;
        }
        $i++;
    }

    # add new taxa if required
    if ( defined $new_taxon[0] ) {

        # uniqu-ify incoming names
        my %seen = ();
        @new_taxon = grep { ! $seen{$_} ++ } @new_taxon;

        # replace _ with spaces
        foreach my $taxon (@new_taxon) {
            $taxon =~ s/_/ /g; 
            $taxon =~ s/'//g; 
            my @found = grep(/^$taxon$/,@taxa);
            # check if in file already
            if ( @found == 0) {
                push @taxa, $taxon; 
            }
        }
    }

    @{ $xml->{Taxa}->[0]->{List} } = @taxa;

    # fix number
    my $nTaxa = @{ $xml->{Taxa}->[0]->{List} };
    $xml->{Taxa}->[0]->{number} = $nTaxa;

    return 1;
}

=head3 tree_equals

 Title   : tree_equals
 Usage   : tree_equals($tree1, $tree2);
 Function: Compare tree1 to tree2 and test for equality
 Returns : 1 if equal, 0 if not. CHECK THIS
 Args    : $tree1, $tree2 - both NEWICK formatted strings.

=cut

sub tree_equals {

    my ( $tree1, $tree2 ) = @_;

    my @t1;
    my @t2;

    $t1[0] = $tree1;
    $t2[0] = $tree2;

    my $nexus1 = _create_nexus_object( \@t1 );
    my $nexus2 = _create_nexus_object( \@t2 );

    my $treesblock_1 = $nexus1->get_block('trees');
    my $treesblock_2 = $nexus2->get_block('trees');

    if ( $treesblock_1->equals($treesblock_2) ) {
        return 1;
    }
    else {
        return 0;
    }

}

######################################################################
#                                                                    #
#                        SEARCH FUNCTIONS                            #
#                                                                    #
#  Functions to search for data within a data structure or file      #
#                                                                    #
#  Names should include _contains_                                   #
#                                                                    #
######################################################################

=head3 tree_contains

 Title   : tree_contains
 Usage   : tree_contains($taxon,$tree_data);
 Function: Check if $tree_data contains $taxon. Case insensitive. Accounts for _ or spaces in name
 Returns : 1 if it does, 0 if not. CHECK THIS
 Args    : $taxon, $tree_data

Note that $tree_data can be a tree file or a tree string

=cut

sub tree_contains {

    my ( $taxon, $input ) = @_;

    $taxon =~ s/_/ /g;
    $taxon = quotemeta $taxon;
    my @taxa = Bio::STK::taxa_from_tree($input);

    foreach my $t (@taxa) {
        if ( $t =~ m/^$taxon$/ ) {
            return 1;
        }

    }

    return 0;
}

=head3 xml_contains_taxon

 Title   : xml_contains_taxon
 Usage   : xml_contains_taxon($taxon,$xml_file);
 Function: Check if $xml_file contains $taxon. Case insensitive. Accounts for _ or spaces in name
 Returns : 1 if it does, 0 if not. CHECK THIS
 Args    : $taxon, $xml_file

$xml_file can be a file or XML data structure

=cut

sub xml_contains_taxon {

    my ( $taxon, $input, $partial ) = @_;

    # partial is optional - defaults to off!
    unless (defined $partial) {
        $partial = 0;
    }
    $taxon =~ s/ /_/g;
    $taxon =~ s/'//g;
    my @taxa = taxa_from_xml($input);
    $taxon = quotemeta $taxon;

    if ($partial) {
        foreach my $t (@taxa) {
            $t =~ s/'//g;
            $t =~ s/ /_/g;
            if ( $t =~ m/(_|^)$taxon(_|$)/) {
                return 1;
            }
        }
    } else {
        foreach my $t (@taxa) {
            $t =~ s/'//g;
            $t =~ s/ /_/g;
            if ( $t =~ m/^$taxon$/) {
                return 1;
            }
        }
    }

    return 0;
}

=head3 xml_contains_analysis

 Title   : xml_contains_analysis
 Usage   : xml_contains_analysis($analysis,$xml);
 Function: Check if $xml contains $analysis
 Returns : 1 if it does, 0 if not. CHECK THIS
 Args    : $analysis, $xml

$xml can be a file or XML data structure

=cut

sub xml_contains_analysis {

    my ( $analysis, $input ) = @_;

    my @analyses = Bio::STK::get_analysis($input);

    foreach my $a (@analyses) {
        if ( $a =~ m/$analysis/i ) {
            return 1;
        }
    }

    return 0;
}

=head3 xml_contains_character

 Title   : xml_contains_character
 Usage   : xml_contains_character($character,$xml,$only);
 Function: Check if $xml contains $character
 Returns : 1 if it does, 0 if not. CHECK THIS
 Args    : $character, $xml

$xml can be a file or XML data structure. C<$character> can be a specific charater type, e.g. 'cytb' or
can be one of "Molecular", "Morphological", "Behavioural" or "Other".

=cut

sub xml_contains_character {

    my ( $c, $input ,$only ) = @_;

    my $found = 0;
    my %characters = Bio::STK::get_characters($input);

    # Is it a main character type, or a sub-type
    if (grep {/$c/i} @character_types) {
        # look for type of character (Morphological, etc)
        for my $type (@character_types) {
            if ( $c =~ m/^$type$/i ) {
                if ( defined(${ $characters{$type} }[0]) ) {
                    $found++;
                }
            }
        }
        # OK, we've found something, no if $only, we need to check there aren't
        # any more lurking in the file...
        if ($only and $found > 0) {
            for my $type (@character_types) {
                if ( $c !~ m/^$type$/i ) {
                    if ( defined(${ $characters{$type} }[0]) ) {
                        $found++;
                    }
                }
            }
        } 
    } else {

        # now check for specific characters
        # feathers, bones, etc
        for my $type (@character_types) {
            for my $m ( @{ $characters{$type} } ) {
                if ( $c =~ m/$m/i ) {
                    $found++;
                }
            }
        }
        # as above
        if ($only and $found > 0) {
            for my $type (@character_types) {
                for my $m ( @{ $characters{$type} }) { 
                    if ( $c !~ m/$m/i ) {
                        $found++;
                    }
                }
            }
        }
    }

    if ($only) {
        if ($found == 1) {
            return 1;
        }
    } else {
        if ($found > 0) {
            return 1;
        }
    }

    return 0;
}

=head3 contains_data

 Title   : contains_data
 Usage   : contains_data(\@data,$xml)
 Function: Check if $xml contains $data. $data is case insensitive, but spaces in taxa names are not accounted for
 Returns : 1 if it does, 0 if not. CHECK THIS
 Args    : $data, $xml
 
C<$data> is a string or array of strings containing the data to be searched for within the C<$xml>. C<$data> can be one of:

=over 4

=item Characters
=item Analysis
=item Taxa

=back

If C<$data> is an array an OR search is carried out.

C<$xml> can be either a XML has (from read_xml_file) or a file name.

=cut

sub contains_data {

    my $needle = shift;
    my $input  = shift;
    my $xml;

    # Check if it's a file or data structure
    # Note that the individual functions below do
    # this, but let's read in the file once
    # if it is a file, rather than 3 times
    my $ret = Bio::STK::_file_or_xml($input);
    if ( $ret == 1 ) {
        $xml = $input;
    }
    elsif ( $ret == 2 ) {

        # read_xml_file checks that the file exists
        $xml = read_xml_file($input);
    }
    else {
        return undef;
    }

    # search through taxa
    for my $taxon ( @{$needle} ) {
        if ( Bio::STK::xml_contains_taxon( $taxon, $xml ) ) {
            return 1;
        }
    }

    # search though analyses
    for my $analysis ( @{$needle} ) {
        if ( Bio::STK::xml_contains_analysis( $analysis, $xml ) ) {
            return 1;
        }
    }

    # search through characters
    for my $char ( @{$needle} ) {
        if ( Bio::STK::xml_contains_character( $char, $xml ) ) {
            return 1;
        }
    }

    return 0;
}

=head3 contains_fossils

 Title   : contains_fossils
 Usage   : contains_fossils($xmlfile);
 Function: Tests if XML file contains fossil taxa
 Returns : 0 - no fossils, 1 all fossils, 2 mixed, -1 incorrect info
 Args    : $xmlfile, else croaks. 

=cut

sub contains_fossils {

    my ($input) = @_;

    my $xml_contents;

    # Check if it's a file or data structure
    my $ret = Bio::STK::_file_or_xml($input);
    if ( $ret == 1 ) {
        $xml_contents = $input;
    }
    elsif ( $ret == 2 ) {

        # read_xml_file checks that the file exists
        $xml_contents = read_xml_file($input);
    }
    else {
        return undef;
    }

    my $fossils = $xml_contents->{Taxa}->[0]->{fossil};

    if ( $fossils =~ m/none/i ) {
        return 0;
    }
    elsif ( $fossils =~ m/all/i ) {
        return 1;
    }
    elsif ( $fossils =~ m/some/i ) {
        return 2;
    }
    else {
        return -1;
    }

}

# give this the taxon and the tree string
sub in_polytomy {

    my ( $taxon, $tree ) = @_;

    $taxon =~ s/ /_/g;
    if ($tree =~ m/\w+,$taxon,\w+(,|\))/ or
        $tree =~ m/$taxon,\w+,\w+(,|\))/ or
        $tree =~ m/\w+,\w+,$taxon(,|\))/) {
        return 1
    } else {
        return 0
    }
}

######################################################################
#                                                                    #
#                          MISC FUNCTIONS                            #
#                                                                    #
#  Misc functions that don't fit into other categories               #
#                                                                    #
######################################################################

=head3 get_short_study_name

 Title   : get_short_study_name
 Usage   : get_short_study_name(@xmlFileList);
 Function: Return a unique (for the dataset defined by the input array) name for a source tree study
 Returns : A linked list containing files (key) and the short name (value)
 Args    : @xmlfile

Call this function before the main loop over XML files to generate the short name. You can then use
the resulting hash to grab the short name for whichever file you're currently dealing with.

=cut

sub get_short_study_name {

    my @xml_files = @_;
    my %short_names;
    my $last_author = '';
    my $last_title  = '';
    my $last_year   = '';
    my $year_a      = 'a';
    my $need_letter = 0;

    for my $file (@xml_files) {

        # load file
        my $xml = read_xml_file($file);

        # put together the author string
        # assume: "last, I.N., Last, I.N, ..." format
        # count commas

        # author_year
        my $author = $xml->{Source}->[0]->{Author}->[0];
        my $year   = $xml->{Source}->[0]->{Year}->[0];
        my $title  = $xml->{Source}->[0]->{Title}->[0];

        # temp vars for above
        my $y = $year;
        my $t = $title;
        my $a = $author;

        # put together the author string
        # assume: "last, I.N., Last, I.N, ..." format
        # count commas
        my $nCommas = $a =~ tr/,//;
        $a = substr( $a, 0, index( $a, ',' ) );
        if ( $nCommas > 2 ) {

            # add etal and year
            $a = $a . "_etal_" . $y;
        }
        elsif ( $author =~ m/and\s(\w+),/ ) {
            $a = $a . "_" . $1 . "_" . $y;
        }
        else {
            $a = $a . "_" . $y;
        }

        # just to be rid of any spaces...
        $a =~ s/ /_/g;
        $a = lc($a);

        # now figure out year, if the authors match and the title is
        # different, we need to add a letter
        if ( $title eq $last_title ) {
            if ( $year_a ne 'a' ) {
                $need_letter = 1;
            }
        }
        else {
            $need_letter = 0;
            $year_a      = 'a';
        }

        if (   $year eq $last_year
            && $a eq $last_author
            && $title ne $last_title )
        {
            $year_a++;
            $a .= $year_a;
        }
        elsif ($need_letter) {
            $a .= $year_a;
        }

        $last_title  = $title;
        $last_author = $a;
        $last_year   = $year;

        $short_names{$file} = $a;

    }

    return %short_names

}

######################################################################
#                                                                    #
#                       INTERNAL FUNCTIONS                           #
#                                                                    #
#  Functions used with the other functions above and not available   #
#  publically. Name must start with an _                             #
#                                                                    #
######################################################################

=begin comment

 Title   : _find_xml
 Usage   : find( \&Bio::STK::_find_xml, $dir );
 Function: The subroutine used by File::Find to locate xml
 Returns : nothing
 Args    : nothing
 
Uses .xml to search for XML files. Note that there are a few directories ignored. 

TODO: Make directories ignored some kind of configurable argument
 
=end comment 

=cut

sub _find_xml {

    if ( $File::Find::dir =~ m/original/ ) {
        return;
    }
    if ( $File::Find::dir =~ m/Adj for/i ) {
        return;
    }
    if ( $File::Find::dir =~ m/Adj_for/i ) {
        return;
    }

    
    if ( (/\.xml$/) ) {

        push( @_files, $File::Find::name );
    }
}

=begin comment

 Title   : _find_trees
 Usage   : find( \&Bio::STK::_find_trees, $dir );
 Function: The subroutine used by File::Find to locate tre
 Returns : nothing
 Args    : nothing
 
Uses .tre to search for Tree files. Note that there are a few directories ignored. 

TODO: Make directories ignored some kind of configurable argument
 
=end comment 

=cut

sub _find_trees {

    if ( $File::Find::dir =~ m/original/i ) {
        return;
    }
    if ( $File::Find::dir =~ m/Adj for/i ) {
        return;
    }
    if ( $File::Find::dir =~ m/Adj_for/i ) {
        return;
    }
    if ( $File::Find::name =~ m/temp.tre/i ) {
        return;
    }
    if ( $File::Find::name =~ m/combinedTrees/i ) {
        return;
    }
    if ( $File::Find::name =~ m/all_trees/i ) {
        return;
    }

    if ( (/\.tre$/) ) {

        push( @_files, $File::Find::name );
    }
}

=begin comment

 Title   : _load_nexus
 Usage   : my $nexus = _load_nexus($file);
 Function: Load a NEXUS file and return a nexus object for manipulation
 Returns : a Bio::NEXUS object
 Args    : filename to load. Croaks if non-existant.
 
=end comment 

=cut

sub _load_nexus {

    my ($file) = @_;

    # check if file exists
    unless ( -e $file ) {
        croak("Error - specified tree file does not exist: $file\n");
    }

    # Bio::NEXUS spits out a hideously complex error message
    # if it can't read the file. Trap this and spit out a friendlier
    # message instead.
    my $nexus;
    eval { 
        no warnings 'all';
    	$nexus = new Bio::NEXUS( $file, 0 );
    };
    if ($@ ne '') {
        croak("Error loading tree file: $file. It exists, but is not truely compliant with the NEXUS format");
    } 

    return $nexus;

}

=begin comment

 Title   : _get_tree_array
 Usage   : my @trees = _get_tree_array($nexus);
 Function: Grab the tree strings from a Bio::NEXUS object
 Returns : an array of Newick formatted tree strings
 Args    : a Bio::NEXUS object. Note no checks are done on this...
 
TODO: Check for Bio::NEXUS object
 
=end comment 
=cut

sub _get_tree_array {

    my ($nexus) = @_;

    my @trees;
    # generate new newick string
    # loop over all trees found in file and return
    foreach my $block ( @{ $nexus->get_blocks() } ) {
        my $type = $block->{'type'};
        if ( $type eq "trees" ) {
            foreach my $tree ( @{ $block->get_trees() } ) {
                # add to array after stripping off ;
                my $t =  substr( $tree->as_string_inodes_nameless(), 0, -1 );
                $t =~ s/\:\w*(,|\))/$1/g;
                push( @trees, $t );
            }
        }
    }
    

    return @trees;

}

=begin comment

 Title   : _create_nexus_object
 Usage   : $nexus_obj = _create_nexus_object(\@trees,\@names);
 Function: Creates a Bio::NEXUS object from an array of NEWICK formatted strings
 Returns : A Bio::NEXUS object
 Args    : An array of NEWICK formatted strings, array of names
 
=end comment 

=cut

sub _create_nexus_object {

    my $trees = shift;
    my $names = shift;
    my $name  = '';

    my $nexus_obj   = new Bio::NEXUS();
    my $trees_block = new Bio::NEXUS::TreesBlock('trees');

    my $i = 0;
    foreach my $tree ( @{$trees} ) {
        if ( !defined $$names[$i] ) {
            $name = "tree_" . ( $i + 1 );
        }
        else {
            $name = $$names[$i];
        }
        $trees_block->add_tree_from_newick( $tree, $name );
        $i++;
    }

    $nexus_obj->add_block($trees_block);

    # fixes bug with $nexus->get_otus, which needs a taxa block
    $nexus_obj->set_taxablock();

    return $nexus_obj;
}

=begin comment

 Title   : _file_or_tree
 Usage   : $ret = _file_or_tree($input);
 Function: Check a string to see if it's a file or a tree string
 Returns : 0: undef, 1: tree_string, 2: file
 Args    : A string which could be either file or tree_string
 
=end comment 

=cut

sub _file_or_tree {

    my $input = shift;

    # check if input is a file or tree string
    # a file cannot contain ( and ) and ,
    if ( defined($input) ) {
        if (-e $input) {

            # it's a file
            return 2;
        }
        else {
            if ( $input =~ m/\(.*,.*\)+/g ) {
              # it's a tree string (probably)
              return 1;
            } else {
              # random string
              return 0;
            }
        }
    }
    else {
        return 0;
    }

}

=begin comment

 Title   : _file_or_xml
 Usage   : $ret = _file_or_xml($input);
 Function: Check a string to see if it's a file or an xml hash
 Returns : 0: undef, 1: xml hash, 2: file
 Args    : A variable which could be either file or xml hash
 
=end comment 

=cut

sub _file_or_xml {

    my $input = shift;

    # check if input is a file or tree string
    # a file cannot contain ( and ) and ,
    if ( defined($input) ) {
        if ( ref($input) eq "HASH" ) {
            # our hash should contain these XML tags...
            if ( $input->{Taxa}  &&
                 $input->{Source} &&
                 $input->{TreeFile} ) {

                # it's a xml hash
                return 1;
            }
            else {

                # it's a hash, but not an xml one...
                return 0;
            }

         # note - could add check for the existance of a file, but
         # we asusme that is done elsewhere as this is menat to be a lightweight
         # function
        }
        else {
           if (-e $input) {
              # it's a file
              return 2;
            } else {
              return 0;
            }
        }
    }
    else {
        return 0;
    }

}

1;
