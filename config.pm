#!/usr/bin/perl
package snpir::config;

use lib qw(../);

use base 'Exporter';


################################################################################
################################################################################
# $Revision: $
# Authors: Robert Piskol ( piskol@stanford.edu )
# Last modification $Author: piskol $
#
# configuration file

our $BLATEXE = '/path/to/blat/executable'; #this is the path to the blat executable
our $SAMTOOLSEXE = '/path/to/samtools/executable'; #this is the path to the samtools executable
our $FASTAFROMBED = '/path/to/fastaFromBed'; #this is the path to the fastaFromBed executable in the bedtools package 

our @EXPORT = (
        '$BLATEXE',
        '$SAMTOOLSEXE',
        '$FASTAFROMBED',
);

1;