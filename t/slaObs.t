#!perl
# Test of slaObs

use strict;
use Test;

BEGIN { plan tests => 10 }

use Astro::SLA;

ok(1);

# First ask for telescope 1

my $i = 2;
my ($n, $name1, $w1, $p1, $h1);

slaObs($i, $n, $name1, $w1, $p1, $h1);
ok(1);
ok(length($n) > 0); # previous bug lost $n if number was specified

# Now ask for the parameters associated with the short telescope
# name associated with telescope 1
slaObs(-1, $n, my $name2, my $w2, my $p2, my $h2);

ok($name1, $name2);
ok($w1, $w2);
ok($p1, $p2);
ok($h1, $h2);

print "# $i $n $name2 $w2 $p2 $h2\n";

# Make sure we do not core dump if the second argument is undef
my $x;
slaObs(20,$x,my $a, my $b, my $c,my $d);
ok(1);

# Make sure we can specify a constant for the telescope name
slaObs(-1,'UKIRT',$a,$b,$c,$d);
ok(1);

# Same with undefined first arg
slaObs(undef,'JCMT',$a,$b,$c,$d);
ok(1);

