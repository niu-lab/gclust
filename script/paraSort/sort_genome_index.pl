#!/usr/bin/perl
# Author: Beifang Niu
# This script is used to sort genome sequences as sequence length from larger to small.
#

use strict;
use warnings;

use IO::File;
use Getopt::Long;

my ( $gf, $sgf, $help, $options, );
$options = GetOptions( 'index-file=s' => \$gf, '--sortedindex-file=s' => \$sgf, 'help' => \$help, );
unless( $options ) { die help_text(); } 
if ($help) { print STDOUT help_text(); exit 0; }
unless( $gf ) { warn "You must provide a genomes file to be sorted ! \n"; die help_text(); }
my $gfr = IO::File->new($gf) or die "can not open this file ! \n";
# reading input genomes file
my @genome;
my ($index, $value);
while(<$gfr>){ chomp; ($index, $value)=split/\t/; push(@genome, [$index, $value] ); }
undef $gfr;
my $genomenum = @genome;
print STDOUT "total genomes:\t".$genomenum."\n";
print STDOUT "sorting....\n";
@genome = sort{$b->[1] <=> $a->[1]} @genome;
print STDOUT "genomes sorting by size finished\n";
print STDOUT "writing output ...\n";
my $wsgf = IO::File->new("> $sgf") or die "can not create this file to write !\n";
foreach (@genome) {
    #chop($_->[0]);
    print $wsgf $_->[0]."\t".$_->[1]."\n";
}
undef $wsgf;

sub help_text{
		return <<HELP
Sortgenome - This script is used to sort genome sequences as sequence length from larger to small.

OPTIONS
--index-file          the input genome index file
--sortedindex-file    the output sorted genome index file
--help                  this message 

SUPPORT
For user support please mail neilniu.cn\@gmail.com

HELP
}

1;

