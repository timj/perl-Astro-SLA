package Astro::SLA;
 
require Exporter;
require DynaLoader;

use strict;
use Carp;
use vars qw(@ISA @EXPORT_OK $VERSION);

$VERSION = '0.10';

@ISA = qw(Exporter DynaLoader); 



@EXPORT_OK = qw( slaFk45z slaDtf2r slaDaf2r slaDr2af slaDr2tf slaPreces
	      slaEqeqx slaCldj slaDe2h slaDh2e slaDat slaGmst slaDjcl
	      slaDd2tf
	      slaAddet slaAfin slaAirmas slaAmp slaAmpqk slaAop
	      slaAoppa slaAoppat slaAopqk slaAtmdsp slaAv2m slaCaldj
	      slaCalyd slaClyd slaDafin slaDav2m slaDbear slaDbjin
	      slaDc62s slaDcc2s slaDcmpf slaDcs2c slaDeuler slaDfltin
	      slaDimxv slaDjcal slaDm2av slaDmoon slaDmxm slaDmxv
	      slaDpav slaDrange slaDranrm slaDs2c6 slaDs2tp
	      slaDsep slaDtf2d slaDtp2s slaDtp2v slaDtps2c
	      slaDtpv2c slaDtt slaDv2tp slaDvdv slaDvn slaDvn slaEarth
	      slaEcmat slaEcor slaEg50 slaEpb slaEpb2d slaEpco
	      slaEpj slaEpj2d slaEqecl slaEqgal slaEtrms slaEvp
	      slaFk425 sla524 slaFk54z  slaGaleq slaGalsup slaGe50
	      slaGeoc slaGmsta slaGresid slaImxv slaInvf
	      slaKbj slaMap slaMappa slaMapqk slaMapqkz slaMoon
	      slaNut slaNutc slaOap slaOapqk slaObs slaPa
	      slaPcd slaPda2h slaPdq2h slaPlanel slaPlanet slaPlante
	      slaPm slaPolmo slaPrebn slaPrec slaPrecl slaPrenut
	      slaPvobs slaRcc slaRdplan slaRefco slaRefcoq
	      slaRefv slaRefz slaRverot slaRvgalc slaRvlg slaRvlsrd
	      slaRvlsrk slaS2tp slaSubet slaSupgal slaUnpcd slaWait
	      slaXy2xy slaZd

		 DPI D2PI D1B2PI D4PI D1B4PI DPISQ DSQRPI DPIBY2 
		 DD2R DR2D DAS2R DR2AS DH2R DS2R D15B2P

		 lst_from_ut
	    ); 
 
 
bootstrap Astro::SLA;
 
=item Constants

Constants supplied by this module.

=cut

# Pi
use constant DPI => 3.1415926535897932384626433832795028841971693993751;
 
# 2pi 
use constant D2PI => 6.2831853071795864769252867665590057683943387987502;
 
# 1/(2pi)
use constant D1B2PI => 0.15915494309189533576888376337251436203445964574046;
 
# 4pi
use constant D4PI => 12.566370614359172953850573533118011536788677597500;
 
# 1/(4pi) 
use constant D1B4PI => 0.079577471545947667884441881686257181017229822870228;
 
# pi^2 
use constant DPISQ => 9.8696044010893586188344909998761511353136994072408;
 
# sqrt(pi) 
use constant DSQRPI => 1.7724538509055160272981674833411451827975494561224;
 
# pi/2:  90 degrees in radians 
use constant DPIBY2 => 1.5707963267948966192313216916397514420985846996876;
 
# pi/180:  degrees to radians 
use constant DD2R => 0.017453292519943295769236907684886127134428718885417;
 
# 180/pi:  radians to degrees 
use constant DR2D => 57.295779513082320876798154814105170332405472466564;
 
# pi/(180*3600):  arcseconds to radians 
use constant DAS2R => 4.8481368110953599358991410235794797595635330237270e-6;
 
# 180*3600/pi :  radians to arcseconds 
use constant DR2AS => 2.0626480624709635515647335733077861319665970087963e5;
 
# pi/12:  hours to radians 
use constant DH2R => 0.26179938779914943653855361527329190701643078328126;
 
# 12/pi:  radians to hours 
use constant DR2H => 3.8197186342054880584532103209403446888270314977709;
 
# pi/(12*3600):  seconds of time to radians 
use constant DS2R => 7.2722052166430399038487115353692196393452995355905e-5;
 
# 12*3600/pi:  radians to seconds of time 
use constant DR2S => 1.3750987083139757010431557155385240879777313391975e4;
 
# 15/(2pi):  hours to degrees x radians to turns 
use constant D15B2P => 2.3873241463784300365332564505877154305168946861068;


=item ($lst, $mjd) = lstnow($long);
 
Return current LST and MJD

=cut

sub lstnow {

   croak 'Usage: lstnow($long)' unless scalar(@_) == 1;

   my $long = shift;

   my ($sign, @ihmsf);

   # Get current UT time
   my ($sec, $min, $hour,$mday,$mon,$year,$wday,$yday,$isdst) = gmtime(time);
 
   # Calculate LST
   $year += 1900;
   $mon++;
   my ($lst, $mjd) = lst_from_ut($year, $mon, $mday, $hour, $min, $sec, $long);

   slaDr2tf(1, $lst, $sign, \@ihmsf);
   print "LST: $sign$ihmsf[0] $ihmsf[1] $ihmsf[2]\n";
     
   return ($lst, $mjd);

}


=item ($lst, $mjd) = lst_from_ut(yy, mn, dd, hh, mm, ss, long)

Given the UT time, calculate the Modified Julian date and the 
local sidereal time for the specified longitude.

Longitude should be negative if degrees west.

=cut

sub lst_from_ut {

  croak 'Usage: lst_from_ut(yy,mn,dd,hh,mm,ss,long)' 
    unless scalar(@_) == 7;

  my ($yy, $mn, $dd, $hh, $mm, $ss, $long) = @_;

  my ($rad, $j, $fd, $mjd, $slastatus, $gmst, $eqeqx, $lst);

  # Calculate fraction of day
  slaDtf2r($hh, $mm, $ss, $rad, $j);

  $fd = $rad / D2PI;

  # Calculate modified julian date of UT day
  slaCldj($yy, $mn, $dd, $mjd, $slastatus);
    
  if ($slastatus != 0) {
    croak "Error calculating modified Julian date with args: $yy $mn $dd\n";
  }

  # Calculate sidereal time of greenwich
  $gmst = slaGmsta($mjd, $fd);
  
  # Find MJD of current time (not just day)
  $mjd += $fd;
 
  # Equation of the equinoxes
  $eqeqx = slaEqeqx($mjd);

  # Local sidereal time = GMST + EQEQX + Longitude in radians

  $lst = $gmst + $eqeqx + ($long * DD2R);
  
  $lst += D2PI if $lst < 0.0;

  return ($lst, $mjd);

}



1;
