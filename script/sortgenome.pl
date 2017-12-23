#!/usr/bin/perl
# Author: Beifang Niu
# This script is used to sort genome sequences as sequence length from larger to small.
#

use strict;
use warnings;

use IO::File;
use Getopt::Long;

my ( $gf, $sgf, $help, $options, );
$options = GetOptions( 'genomes-file=s' => \$gf, 'sortedgenomes-file=s' => \$sgf, 'help' => \$help, );
unless( $options ) { die help_text(); } 
if ($help) { print STDOUT help_text(); exit 0; }
unless( $gf ) { warn "You must provide a genomes file to be sorted ! \n"; die help_text(); }
my $gfr = IO::File->new($gf) or die "can not open this file ! \n";
# reading input genomes file
my ( $i, $t, @genome, $name, $len, $str, $ll , );
$i = $t = 0;
while(<$gfr>){chomp; $_=$_."\n"; if (/^>/) { if ($i > 0) { push( @genome, [$name,$len,$str] ); $str = ""; $len = 0; } $name = $_; $i++; } else { $t = length($_) - 1; $str .= $_; $len += $t; } }
push( @genome,[$name,$len,$str] );
undef $gfr;
my $genomenum = @genome;
print STDOUT "total genomes:\t".$genomenum."\n";
print STDOUT "sorting....\n";
@genome = sort{$b->[1] <=> $a->[1]} @genome;
print STDOUT "genomes sorting by size finished\n";
print STDOUT "writing output ...\n";
my $wsgf = IO::File->new("> $sgf") or die "can not create this file to write !\n";
foreach (@genome) {
    chop($_->[0]);
    print $wsgf $_->[0]."\n".$_->[2];
}
undef $wsgf;

sub help_text{
		return <<HELP
Sortgenome - This script is used to sort genome sequences as sequence length from larger to small.

OPTIONS
--genomes-file          the input genome sequences file
--sortedgenomes-file    the output sorted genome sequences file
--help                  this message 

SUPPORT
For user support please mail neilniu.cn\@gmail.com

HELP
}

1;

