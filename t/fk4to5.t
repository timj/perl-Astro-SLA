#!/bin/perl

use Astro::SLA;

###########################################################################

# Only one test
print "1..1\n";

$ra = "6 10 23.9";
$dec = "-6 12 21.0";

print "Input (B1950): RA=$ra DEC=$dec\n";

($nra,$ndec) = &btoj($ra,$dec);

if ($nra eq "6 12 50.37" && $ndec eq "-6 13 11.76") {
  print "ok\n";
} else {
  print "not ok\n";
}

print "Output (J2000): RA=$ra DEC=$dec\n";


sub btoj {

  my ($ra,$dec) = (@_);
  my ($nra,$ndec,$bepoch,$dsign,$sign,$sign2,$status,$j);
  my (@idmsf,@ihmsf,$ra_rad,$dec_rad,$ra_j2000_rad,$dec_j2000_rad);


  my ($h,$m,$s) = split(/ /,$ra);
  my ($d,$dm,$ds) = split(/ /,$dec);

  slaDtf2r($h,$m,$s ,$ra_rad,$j);

  $dsign = ($d < 0 ? -1 : 1);
  $d = abs($d);

  slaDaf2r($d,$dm,$ds ,$dec_rad,$status);

  $dec_rad *= $dsign;

###########################################################################

$bepoch = 1950.0;
slaFk45z( $ra_rad, $dec_rad, $bepoch, $ra_j2000_rad, $dec_j2000_rad );

@idmsf = ();
&slaDr2af(2,$dec_j2000_rad,$sign,\@idmsf);

@ihmsf = ();
slaDr2tf(2,$ra_j2000_rad,$sign2,@ihmsf);


    $nra  = join(" ",@ihmsf[0..2]).".$ihmsf[3]";
    $ndec = $sign.join(" ",@idmsf[0..2]).".$idmsf[3]";


return $nra,$ndec;

}

