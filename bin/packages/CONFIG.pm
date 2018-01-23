package CONFIG;


use vars qw(@ISA @EXPORT);
require Exporter;


@ISA = qw(Exporter);


##################################################
###################################
## PATHS

## TARGET RNAs
$FUNGI_TARGET_RNA_PATH = "targetRNAs/fungi/";
$PROTOSTOMIA_TARGET_RNA_PATH = "targetRNAs/protostomia/";
$DEUTEROSTOMIA_TARGET_RNA_PATH = "targetRNAs/deuterostomia/";
$PLANT_TARGET_RNA_PATH = "targetRNAs/plants/";

## MODELS
$FUNGI_CM_PATH = "models/fungi/";
$PROTOSTOMIA_CM_PATH = "models/protostomia/";
$DEUTEROSTOMIA_CM_PATH = "models/deuterostomia/";
$PLANT_CM_PATH = "models/plants/";

## FASTA FILES
$FUNGI_FASTA_PATH = "fasta/fungi/";
$PROTOSTOMIA_FASTA_PATH = "fasta/protostomia/";
$DEUTEROSTOMIA_FASTA_PATH = "fasta/deuterostomia/";
$PLANT_FASTA_PATH = "fasta/plants/";

## INFORMATION FILES
$FUNGI_INFORMATION_FILE = "information/fungi_information.csv";
$PROTOSTOMIA_INFORMATION_FILE = "information/protostomia_information.csv";
$DEUTEROSTOMIA_INFORMATION_FILE = "information/deuterostomia_information.csv";
$PLANT_INFORMATION_FILE = "information/plant_information.csv";

##################################################
###################################
## PROGRAMS

$MUSCLE = "/usr/local/bin/muscle";
$BLASTALL = "/usr/local/bin/blastall";
$RNASUBOPT = "/usr/local/bin/RNAsubopt";
$CMSEARCH = "~/bin/cmsearch";
$CMSTAT = "~/bin/cmstat";
$CMBUILD = "~/bin/cmbuild";
$CMCALIBRATE = "~/bin/cmcalibrate";
$FASTALENGTH = "/usr/local/bin/fastalength";
$FASTACMD = "/usr/local/bin/fastacmd";

$RSCRIPT = "/usr/bin/Rscript";
$RNAPLEX = "/scr/pitu/sebastian/SNOSTRIP/snoStrip/bin/programs/RNAplex";
$RNASNOOP = "/scr/pitu/sebastian/SNOSTRIP/snoStrip/bin/programs/RNAsnoop";



##################################################
###################################
## SCRIPTS

$STOCKHOLM = "/scr/pitu/sebastian/SNOSTRIP/snoStrip/bin/scripts/stockholm.pl";
$CLUSTER = "/scr/pitu/sebastian/SNOSTRIP/snoStrip/bin/scripts/clusterChains.R";
$SVM = "/scr/pitu/sebastian/SNOSTRIP/snoStrip/bin/scripts/apply_svm.pl";
$INFERNAL = "/scr/pitu/sebastian/SNOSTRIP/snoStrip/bin/scripts/infernal.pl";
$BLAST = "/scr/pitu/sebastian/SNOSTRIP/snoStrip/bin/scripts/blastRelatives.pl";
$MAP_MODIFICATIONS = "/scr/pitu/sebastian/SNOSTRIP/snoStrip/bin/scripts/map_modifications.pl";


##################################################
###################################
## BLASTHIT PARAMETER

$CD_CUTOFF_DEU = 3;
$CD_CUTOFF_FUNGI = 2;
$HACA_CUTOFF_DEU = 2;
$HACA_CUTOFF_FUNGI = 1;
$CD_MEAN = 10.55;
$CD_SD = 16.05;
$HACA_MEAN = 9.5;
$HACA_SD = 17.57;

1;
