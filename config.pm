#!/usr/bin/perl
package SNPiR::config;

use lib qw(../);

use base 'Exporter';


################################################################################
################################################################################
# $Revision: $
# Authors: Robert Piskol ( piskol@stanford.edu )
# Last modification $Author: piskol $
#
# configuration file

#our $BLATEXE = '/path/to/blat/executable'; #this is the path to the blat executable
#our $SAMTOOLSEXE = '/path/to/samtools/executable'; #this is the path to the samtools executable
#our $FASTAFROMBED = '/path/to/fastaFromBed'; #this is the path to the fastaFromBed executable in the bedtools package 

our $BLATEXE = '/srv/gs1/projects/li/shared-data/software/kent/src/blat/blat'; #this is the path to the blat executable
our $SAMTOOLSEXE = '/srv/gs1/projects/li/shared-data/software/samtools-0.1.16/samtools'; #this is the path to the samtools executable
our $FASTAFROMBED = '/srv/gs1/projects/li/shared-data/software/bedtools-2.17.0/bin/fastaFromBed'; #this is the path to the fastaFromBed executable in the bedtools package 

our @EXPORT = (
        '$BLATEXE',
        '$SAMTOOLSEXE',
        '$FASTAFROMBED',
);

1;