# snoStrip - an automatic snoRNA annotation pipeline
-----
Current release: v2.0


## Installation

Instruction for installation and configuration of "snoStrip" on Linux/Unix machines.


### REQUIRED PROGRAMS

The snoStrip pipeline utilizes several open source programs for detection and property extraction of putative snoRNAs.
Please ensure that the following programs are correctly installed:

1) snoRNA detection:
   
   * blastall		version 2.2.26 or newer
   * infernal		version 1.1 or newer

   * R			version 3.0.1 or newer
   * R packages:	genomeIntervals


2) property extraction:
   
   * MUSCLE		version 3.7 or newer
   * RNAsubopt		version 2.1.9 or newer

   * fastacmd		version 2.2.26 or newer
   * fastalength	version 2.0.0 or newer


### CONFIGURATION OF SNOSTRIP

The configuration file, the snoStrip pipeline, all necessary packages, and additional
scripts are located in the folder 'bin'.

Please move into the snoStrip directory and start the configuration script by

```
cd bin/
./configure.sh
```

The configure.sh script will search for all necessary programs and additional scripts. 
Several snoStrip scripts will be adapted to find the path required modules.

Once the configuration process is finished, the snoStrip pipeline should work properly.



## Run snoStrip

Plase read the help page of the snoStrip pipeline.
```
./snoStrip.pl -h
```

