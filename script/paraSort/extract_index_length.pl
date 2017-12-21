#!/usr/bin/perl
# Author: Beifang Niu
# This script is used to sort genome sequences as sequence length from larger to small.
#

use strict;
use warnings;

use IO::File;
use Getopt::Long;

my ( $gf, $sgf, $help, $options, );
$options = GetOptions( 'genomes-file=s' => \$gf, 'index-file=s' => \$sgf, 'help' => \$help, );
unless( $options ) { die help_text(); } 
if ($help) { print STDOUT help_text(); exit 0; }
unless( $gf ) { warn "You must provide a genomes file to be sorted ! \n"; die help_text(); }
my $gfr = IO::File->new($gf) or die "can not open this file ! \n";
my $wsgf = IO::File->new("> $sgf") or die "can not create this file to write !\n";
# reading input genomes file
my ( $i, $t, $len, $ll ,$out );
$i = $t = 0;
print STDOUT "writing output ...\n";
#map{ if (/^>/) { if ($i > 0) { $out=join("\t",($i-1,$len));  print $wsgf $out."\n";  $len = 0; }  $i++; } else { $t = length($_) - 1;  $len += $t; } } <$gfr>;
while(<$gfr>){ if (/^>/) { if ($i > 0) { $out=join("\t",($i-1,$len));  print $wsgf $out."\n";  $len = 0; }  $i++; } else { $t = length($_) - 1;  $len += $t; } }

#push( @genome,join("\t",($i,$len)) );
$out=join("\t",($i-1,$len));  print $wsgf $out."\n";

undef $gfr;
undef $wsgf;

sub help_text{
		return <<HELP
Sortgenome - This script is used to sort genome sequences as sequence length from larger to small.

OPTIONS
--genomes-file          the input genome sequences file
--index-file            the output index and sequence length file
--help                  this message 

SUPPORT
For user support please mail neilniu.cn\@gmail.com

HELP
}

1;

