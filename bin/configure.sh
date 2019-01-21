#!/bin/bash
#
# File:   configure.sh    
# Author: sebastian@vodka.bioinf.uni-leipzig.de
# Date:   2013/09/03 08:41:11
#
# Changes: Jan Engelhardt, 21/1/19
#

# Globale Variablen
SCRIPTNAME=$(basename $0 .sh)
EXIT_SUCCESS=0
EXIT_FAILURE=1
EXIT_ERROR=2
EXIT_BUG=10

# --- your code BEGIN --------------------------------------------------------- #



## get the current working directory
WRK_DIR=`pwd`"/"
DATA_DIR=`pwd | sed 's/\/[^\/]*$/\/data\//'`
CONFIG=$WRK_DIR"packages/CONFIG.pm"

echo -n "write header ... "

####################################################################################################
## STEP 1  -  HEADER
## write header into the CONFIG package
####################################################################################################
echo -e "package CONFIG;\n\n" > $CONFIG
echo -e "use vars qw(@ISA @EXPORT);" >>  $CONFIG
echo -e "require Exporter;\n\n" >> $CONFIG
echo -e "@ISA = qw(Exporter);\n\n" >> $CONFIG

echo -e "done\n"
echo -n "collect all PATH information ... "


####################################################################################################
## STEP 2  -  PATHS
####################################################################################################

## write the working path into the CONFIG file
echo -e "##################################################\n###################################\n## PATHS\n" >> $CONFIG

## write TARGET paths
echo -e "## TARGET RNAs" >> $CONFIG
echo -e "\$FUNGI_TARGET_RNA_PATH = \"targetRNAs/fungi/\";" >> $CONFIG
echo -e "\$PROTOSTOMIA_TARGET_RNA_PATH = \"targetRNAs/protostomia/\";" >> $CONFIG
echo -e "\$DEUTEROSTOMIA_TARGET_RNA_PATH = \"targetRNAs/deuterostomia/\";" >> $CONFIG
echo -e "\$PLANT_TARGET_RNA_PATH = \"targetRNAs/plants/\";\n" >> $CONFIG

## write CM MODEL paths
echo -e "## MODELS" >> $CONFIG
echo -e "\$FUNGI_CM_PATH = \"models/fungi/\";" >> $CONFIG
echo -e "\$PROTOSTOMIA_CM_PATH = \"models/protostomia/\";" >> $CONFIG
echo -e "\$DEUTEROSTOMIA_CM_PATH = \"models/deuterostomia/\";" >> $CONFIG
echo -e "\$PLANT_CM_PATH = \"models/plants/\";\n" >> $CONFIG

## write MULTI FASTA paths
echo -e "## FASTA FILES" >> $CONFIG
echo -e "\$FUNGI_FASTA_PATH = \"fasta/fungi/\";" >> $CONFIG
echo -e "\$PROTOSTOMIA_FASTA_PATH = \"fasta/protostomia/\";" >> $CONFIG
echo -e "\$DEUTEROSTOMIA_FASTA_PATH = \"fasta/deuterostomia/\";" >> $CONFIG
echo -e "\$PLANT_FASTA_PATH = \"fasta/plants/\";\n" >> $CONFIG


## write INFORMATION paths
echo -e "## INFORMATION FILES" >> $CONFIG
echo -e "\$FUNGI_INFORMATION_FILE = \"information/fungi_information.csv\";" >> $CONFIG
echo -e "\$PROTOSTOMIA_INFORMATION_FILE = \"information/protostomia_information.csv\";" >> $CONFIG
echo -e "\$DEUTEROSTOMIA_INFORMATION_FILE = \"information/deuterostomia_information.csv\";" >> $CONFIG
echo -e "\$PLANT_INFORMATION_FILE = \"information/plant_information.csv\";\n" >> $CONFIG

echo -e "done\n"


####################################################################################################
## STEP 3  -  PROGRAMS
####################################################################################################

echo -e "##################################################\n###################################\n## PROGRAMS\n" >> $CONFIG


## search for muscle
echo -n "search for muscle ... "
if test `command -v muscle`; then
    echo "found at "`command -v muscle`
    MUSCLE=`command -v muscle`
else
    echo "not found"
    echo -n "Please enter path to muscle: "
    read MUSCLE
fi

echo "\$MUSCLE = \""$MUSCLE"\";" >> $CONFIG



## search for blastall
echo -n "search for blastall ... "
if test `command -v blastall`; then
    echo "found at "`command -v blastall`
    BLASTALL=`command -v blastall`
else
    echo "not found"
    echo -n "Please enter path to blastall: "
    read BLASTALL
fi

echo "\$BLASTALL = \""$BLASTALL"\";" >> $CONFIG



## search for RNAsubpot
echo -n "search for RNAsubopt ... "
if test `command -v RNAsubopt`; then
    echo "found at "`command -v RNAsubopt`
    RNASUBOPT=`command -v RNAsubopt`
else
    echo "not found"
    echo -n "Please enter path to RNAsubopt: "
    read RNASUBOPT
fi

echo "\$RNASUBOPT = \""$RNASUBOPT"\";" >> $CONFIG


## search for cmsearch
echo -n "search for cmsearch ... "
if test `command -v cmsearch`; then
    tmp=`command -v cmsearch`
    if test `$tmp -h | grep INFERNAL | grep -c 1.1` -eq 1; then
	echo "found at "`command -v cmsearch`
	CMSEARCH=`command -v cmsearch`
    else
	echo "cmsearch of version 1.1 not found"
	echo -n "Please enter path to cmsearch: "
	read CMSEARCH
    fi
else
    echo "not found"
    echo -n "Please enter path to cmsearch: "
    read CMSEARCH
fi

echo "\$CMSEARCH = \""$CMSEARCH"\";" >> $CONFIG



## search for cmstat
echo -n "search for cmstat ... "
if test `command -v cmstat`; then
    tmp=`command -v cmstat`
    if test `$tmp -h | grep INFERNAL | grep -c 1.1` -eq 1; then
	echo "found at "`command -v cmstat`
	CMSTAT=`command -v cmstat`
    else
	echo "cmstat of verion 1.1 not found"
	echo -n "Please enter path to cmstat: "
	read CMSTAT
    fi
else
    echo "not found"
    echo -n "Please enter path to cmstat: "
    read CMSTAT
fi

echo "\$CMSTAT = \""$CMSTAT"\";" >> $CONFIG



## search for cmbuild
echo -n "search for cmbuild ... "
if test `command -v cmbuild`; then
    tmp=`command -v cmbuild`
    if test `$tmp -h | grep INFERNAL | grep -c 1.1` -eq 1; then
	echo "found at "`command -v cmbuild`
	CMBUILD=`command -v cmbuild`
    else
	echo "cmbuild of verion 1.1 not found"
	echo -n "Please enter path to cmbuild: "
	read CMBUILD
    fi
else
    echo "not found"
    echo -n "Please enter path to cmbuild: "
    read CMBUILD
fi

echo "\$CMBUILD = \""$CMBUILD"\";" >> $CONFIG



## search for cmcalibrate
echo -n "search for cmcalibrate ... "
if test `command -v cmcalibrate`; then
    tmp=`command -v cmcalibrate`
    if test `$tmp -h | grep INFERNAL | grep -c 1.1` -eq 1; then
	echo "found at "`command -v cmcalibrate`
	CMCALIBRATE=`command -v cmcalibrate`
    else
	echo "cmcalibrate of verion 1.1 not found"
	echo -n "Please enter path to cmcalibrate: "
	read CMCALIBRATE 
    fi
else
    echo "not found"
    echo -n "Please enter path to cmcalibrate: "
    read CMCALIBRATE
fi

echo "\$CMCALIBRATE = \""$CMCALIBRATE"\";" >> $CONFIG



## search for fastacmd
echo -n "search for fastacmd ... "
if test `command -v fastacmd`; then
    echo "found at "`command -v fastacmd`
    FASTACMD=`command -v fastacmd`
else
    echo "not found"
    echo -n "Please enter path to fastacmd: "
    read FASTACMD
fi

echo -e "\$FASTACMD = \""$FASTACMD"\";\n" >> $CONFIG


## search for Rscript
echo -n "search for Rscript ... "
if test `command -v Rscript`; then
    echo "found at "`command -v Rscript`
    RSCRIPT=`command -v Rscript`
else
    echo "not found"
    echo -n "Please enter path to Rscript: "
    read RSCRIPT
fi

echo -e "\$RSCRIPT = \""$RSCRIPT"\";" >> $CONFIG


## search for necessary R packages
echo -n "search for necessary R package 'genomeIntervals' ... "
MESSAGE=`$RSCRIPT $WRK_DIR\scripts/check_R_packages.R genomeIntervals`
if test "$MESSAGE" == "package found"; then
    echo "found"
else
    echo "ERROR: PACKAGE 'genomeIntervals' not found"
    exit $EXIT_ERROR
fi



## search for RNAplex
RNAPLEX=$WRK_DIR"programs/RNAplex"
echo "\$RNAPLEX = \""$RNAPLEX"\";" >> $CONFIG


## search for RNAsnoop
RNASNOOP=$WRK_DIR"programs/RNAsnoop"
echo -e "\$RNASNOOP = \""$RNASNOOP"\";\n\n\n" >> $CONFIG


####################################################################################################
## STEP 4  -  SCRIPTS
####################################################################################################

echo -ne "\ncollecting scripts ... "

echo -e "##################################################\n###################################\n## SCRIPTS\n" >> $CONFIG

## search for STOCKHOLM
STK=$WRK_DIR"scripts/stockholm.pl"
echo "\$STOCKHOLM = \""$STK"\";" >> $CONFIG

## search for CLUSTER
CLUSTER=$WRK_DIR"scripts/clusterChains.R"
echo "\$CLUSTER = \""$CLUSTER"\";" >> $CONFIG

## search for SVM
SVM=$WRK_DIR"scripts/apply_svm.pl"
echo "\$SVM = \""$SVM"\";" >> $CONFIG

## search for INFERNAL
INF=$WRK_DIR"scripts/infernal.pl"
echo "\$INFERNAL = \""$INF"\";" >> $CONFIG

## search for BLAST
BLAST=$WRK_DIR"scripts/blastRelatives.pl"
echo -e "\$BLAST = \""$BLAST"\";" >> $CONFIG

## search for modification mapping
MAP=$WRK_DIR"scripts/map_modifications.pl"
echo -e "\$MAP_MODIFICATIONS = \""$MAP"\";\n\n" >> $CONFIG


echo "done"



####################################################################################################
## STEP 5 -  ADD INFORMATION ABOUT THE BLASTHIT PARAMETERS
####################################################################################################

echo -e "##################################################\n###################################\n## BLASTHIT PARAMETER\n" >> $CONFIG

echo -e "\$CD_CUTOFF_DEU = 3;" >> $CONFIG
echo -e "\$CD_CUTOFF_FUNGI = 2;" >> $CONFIG
echo -e "\$HACA_CUTOFF_DEU = 2;" >> $CONFIG
echo -e "\$HACA_CUTOFF_FUNGI = 1;" >> $CONFIG
echo -e "\$CD_MEAN = 10.55;" >> $CONFIG
echo -e "\$CD_SD = 16.05;" >> $CONFIG
echo -e "\$HACA_MEAN = 9.5;" >> $CONFIG
echo -e "\$HACA_SD = 17.57;\n\n1;" >> $CONFIG



exit $EXIT_SUCCESS


# End of file
