/*        -*- C -*-
  
  perl-SLA glue - 99% complete
                                        t.jenness@jach.hawaii.edu

  Copyright (C) 1998-2001 Tim Jenness.  All rights reserved.
  This program is free software; you can redistribute it and/or
  modify it under the same terms as Perl itself.

  Has been tested with the May 1998 release of SLALIB

 */
 
 
#include "EXTERN.h"   /* std perl include */
#include "perl.h"     /* std perl include */
#include "XSUB.h"     /* XSUB include */
 
#include "slalib.h"

#include "arrays.c"


MODULE = Astro::SLA   PACKAGE = Astro::SLA


# Add a few routines

void
slaAddet(rm, dm, eq, rc, dc)
  double rm
  double dm
  double eq
  double rc = NO_INIT
  double dc = NO_INIT
 PROTOTYPE: $$$$$
 CODE:
  slaAddet(rm, dm, eq, &rc, &dc);
 OUTPUT:
  rc
  dc

void
slaAfin(string, nstrt, reslt, jf)
  char * string
  int nstrt
  float reslt = NO_INIT
  int jf = NO_INIT
 PROTOTYPE: $$$$
 CODE:
  slaAfin(string, &nstrt, &reslt, &jf);
 OUTPUT:
  nstrt
  reslt
  jf



double
slaAirmas(zd)
  double zd
 PROTOTYPE: $
 CODE:
  RETVAL = slaAirmas(zd);
 OUTPUT:
  RETVAL

void
slaAltaz(ha, dec, phi, az, azd, azdd, el, eld, eldd, pa, pad, padd)
  double ha
  double dec
  double phi
  double az = NO_INIT
  double azd = NO_INIT
  double azdd = NO_INIT
  double el = NO_INIT
  double eld = NO_INIT
  double eldd = NO_INIT
  double pa = NO_INIT
  double pad = NO_INIT
  double padd = NO_INIT
 PROTOTYPE: $$$$$$$$$$$$
 CODE:
  slaAltaz(ha, dec, phi, &az, &azd, &azdd, &el, &eld, &eldd, &pa, &pad, &padd);
 OUTPUT:
  az
  azd
  azdd
  el
  eld
  eldd
  pa
  pad
  padd

void
slaAmp(ra, da, date, eq, rm, dm)
  double ra
  double da
  double date
  double eq
  double rm = NO_INIT
  double dm = NO_INIT
 PROTOTYPE: $$$$$$
 CODE:
  slaAmp(ra, da, date, eq, &rm, &dm);
 OUTPUT:
  rm
  dm

# FLAG: Need to add a check for number of components in amprms

void
slaAmpqk(ra, da, amprms, rm, dm)
  double ra
  double da
  double * amprms
  double rm = NO_INIT
  double dm = NO_INIT
 PROTOTYPE: $$\@$$
 CODE:
  slaAmpqk(ra, da, amprms, &rm, &dm);
 OUTPUT:
  rm
  dm

void
slaAop(rap,dap,date,dut,elongm,phim,hm,xp,yp,tdk,pmb,rh,wl,tlr,aob,zob,hob,dob,rob)
  double rap
  double dap
  double date
  double dut
  double elongm
  double phim
  double hm
  double xp
  double yp
  double tdk
  double pmb
  double rh
  double wl
  double tlr
  double aob = NO_INIT
  double zob = NO_INIT
  double hob = NO_INIT
  double dob = NO_INIT
  double rob = NO_INIT
 PROTOTYPE: $$$$$$$$$$$$$$$$$$$
 CODE:
  slaAop(rap,dap,date,dut,elongm,phim,hm,xp,yp,tdk,pmb,rh,wl,tlr,&aob,&zob,&hob,&dob,&rob);
 OUTPUT:
  aob
  zob
  hob
  dob
  rob

void
slaAoppa(date,dut,elongm,phim,hm,xp,yp,tdk,pmb,rh,wl,tlr,aoprms)
  double date
  double dut
  double elongm
  double phim
  double hm
  double xp
  double yp
  double tdk
  double pmb
  double rh
  double wl
  double tlr
  double * aoprms = NO_INIT
 PROTOTYPE: $$$$$$$$$$$\@
 CODE:
  aoprms = get_mortalspace(14,'d');
  slaAoppa(date,dut,elongm,phim,hm,xp,yp,tdk,pmb,rh,wl,tlr,aoprms);
  unpack1D( (SV*)ST(12), (void *)aoprms, 'd', 14);
 OUTPUT:
  aoprms

### FLAG: Can give 13 input arguments and receive 14 for slaAoppat
### Must make absolutely sure that we have 14 args going in.
### Too lazy for now

void
slaAoppat(date, aoprms)
  double date
  double * aoprms
 PROTOTYPE: $\@
 CODE:
  /* Should probably allocate [14] doubles here and copy the array
     myself */
  slaAoppat(date, aoprms);
  unpack1D( (SV*)ST(1), (void *)aoprms, 'd', 14);
 OUTPUT:
  aoprms

void
slaAopqk(rap, dap, aoprms, aob, zob, hob, dob, rob)
  double rap
  double dap
  double * aoprms
  double aob = NO_INIT
  double zob = NO_INIT
  double hob = NO_INIT
  double dob = NO_INIT
  double rob = NO_INIT
 PROTOTYPE: $$\@$$$$$
 CODE:
  slaAopqk(rap, dap, aoprms, &aob,&zob,&hob,&dob,&rob);
 OUTPUT:
  aob
  zob
  hob
  dob
  rob

void
slaAtmdsp(tdk, pmb, rh, wl1, a1, b1, wl2, a2, b2)
  double tdk
  double pmb
  double rh
  double wl1
  double a1
  double b1
  double wl2
  double a2 = NO_INIT
  double b2 = NO_INIT
 PROTOTYPE:
 CODE:
  slaAtmdsp(tdk, pmb, rh, wl1, a1, b1, wl2, &a2, &b2);
 OUTPUT:
  a2
  b2

#### FLAG: Need to check return array
#### Should really be a PDL function

void
slaAv2m(axvec, rmat)
  float * axvec
  float * rmat = NO_INIT
 PROTOTYPE: \@\@
 CODE:
  rmat = get_mortalspace(9,'f');
  slaAv2m(axvec, (void*)rmat);
  unpack1D( (SV*)ST(1), (void *)rmat, 'f', 9);
 OUTPUT:
  rmat


### SKIP: slaBear - use DOUBLE precisions version - slaDbear

### SKIP: slaCaf2r - use DOUBLE precisions version - slaDaf2r


void
slaCaldj(iy, im, id, djm, j)
  int iy
  int im
  int id
  double djm = NO_INIT
  int j = NO_INIT
 PROTOTYPE: $$$$$
 CODE:
  slaCaldj(iy, im, id, &djm, &j);
 OUTPUT:
  djm
  j

void
slaCalyd(iy, im, id, ny, nd, j)
  int iy
  int im
  int id
  int ny = NO_INIT
  int nd = NO_INIT
  int j = NO_INIT
 PROTOTYPE: $$$$$$
 CODE:
  slaCalyd(iy, im, id, &ny, &nd, &j);
 OUTPUT:
  ny
  nd
  j

### SKIP: slaCc2s - use Double precision version - slaDc2s
### SKIP: slaCc62s - use Double precision version - slaDc62s
### SKIP: slaCd2tf - use Double precision version - slaDd2tf

# Calendar to MJD

void
slaCldj(iy, im, id, djm, status)
  int iy
  int im
  int id
  double djm = NO_INIT
  int status = NO_INIT
 PROTOTYPE: $$$$$
 CODE:
  slaCldj(iy, im, id, &djm, &status);
 OUTPUT:
  djm
  status



void
slaClyd(iy, im, id, ny, nd, j)
  int iy
  int im
  int id
  int ny = NO_INIT
  int nd = NO_INIT
  int j = NO_INIT
 PROTOTYPE: $$$$$$
 CODE:
  slaClyd(iy, im, id, &ny, &nd, &j);
 OUTPUT:
  ny
  nd
  j


### SKIP: slaCr2af - use DOUBLE instead - slaDr2af
### SKIP: slaCr2tf - use DOUBLE instead - slaDr2tf
### SKIP: slaCs2c - use DOUBLE instead 
### SKIP: slaCs2c6 - use DOUBLE instead - slaDs2c6
### SKIP: slaCtf2d - use DOUBLE instead
### SKIP: slaCtf2r - use DOUBLE instead


## Up to slaDaf2r
#   Converts DMS to radians 

void
slaDaf2r(ideg, iamin, asec, rad, j)
  int ideg
  int iamin
  double asec 
  double  rad = NO_INIT
  int  j = NO_INIT
 ALIAS:
  slaCaf2r = 1
 PROTOTYPE: $$$$$
 CODE:
  slaDaf2r(ideg, iamin, asec, &rad, &j);
 OUTPUT:
 rad
 j


void
slaDafin(string, nstrt, dreslt, jf)
  char * string
  int nstrt
  double dreslt = NO_INIT
  int jf = NO_INIT
 PROTOTYPE: $$$$
 CODE:
  slaDafin(string, &nstrt, &dreslt, &jf);
 OUTPUT:
  nstrt
  dreslt
  jf


# Added 5/5/98

double
slaDat(utc)
  double utc
 PROTOTYPE: $
 CODE:
  RETVAL = slaDat(utc);
 OUTPUT:
  RETVAL




#### Should really be a PDL function

void
slaDav2m(axvec, rmat)
  double * axvec
  double * rmat = NO_INIT
 PROTOTYPE: \@\@
 CODE:
  rmat = get_mortalspace(9,'d');
  slaDav2m(axvec, (void*)rmat);
  unpack1D( (SV*)ST(1), (void *)rmat, 'd', 9);
 OUTPUT:
  rmat


double
slaDbear(a1, b1, a2, b2)
  double a1
  double b1
  double a2
  double b2
 ALIAS:
  slaBear = 1
 PROTOTYPE: $$$$
 CODE:
  RETVAL = slaDbear(a1, b1, a2, b2);
 OUTPUT:
  RETVAL

void
slaDbjin(string, nstrt, dreslt, j1, j2)
  char * string
  int nstrt
  double dreslt = NO_INIT
  int j1 = NO_INIT
  int j2 = NO_INIT
 PROTOTYPE: $$$$
 CODE:
  slaDbjin(string, &nstrt, &dreslt, &j1, &j2);
 OUTPUT:
  nstrt
  dreslt
  j1
  j2

void
slaDc62s(v, a, b, r, ad, bd, rd)
  double * v
  double a = NO_INIT
  double b = NO_INIT
  double r = NO_INIT
  double ad = NO_INIT
  double bd = NO_INIT
  double rd = NO_INIT
 ALIAS:
  slaCc62s = 1
 PROTOTYPE: \@$$$$$$
 CODE:
  slaDc62s(v, &a, &b, &r, &ad, &bd, &rd);
 OUTPUT:
  a
  b
  r
  ad
  bd
  rd


void
slaDcc2s(v,a,b)
  double * v
  double a = NO_INIT
  double b = NO_INIT
 ALIAS:
  slaCc2s = 1
 PROTOTYPE: \@$$
 CODE:
  slaDcc2s(v, &a, &b);
 OUTPUT:
  a
  b


void
slaDcmpf(coeffs, xy, yz, xs, ys, perp, orient)
  double * coeffs
  double xy = NO_INIT
  double yz = NO_INIT
  double xs = NO_INIT
  double ys = NO_INIT
  double perp = NO_INIT
  double orient = NO_INIT
 PROTOTYPE: \@$$$$$$
 CODE: 
  slaDcmpf(coeffs, &xy, &yz, &xs, &ys, &perp, &orient);
 OUTPUT:
  xy
  yz
  xs
  ys
  perp
  orient

void
slaDcs2c(a, b, v) 
  double a
  double b
  double * v = NO_INIT
 PROTOTYPE: $$\@
 CODE:
  v = get_mortalspace(3,'d');
  slaDcs2c(a, b, v);
  unpack1D( (SV*)ST(2), (void *)v, 'd', 3);
 OUTPUT:
  v

#   Converts decimal day to hours minutes and seconds

void
slaDd2tf(ndp, days, sign, ihmsf)
  int ndp
  double  days
  char sign = NO_INIT
  int * ihmsf = NO_INIT
 ALIAS:
  slaCd2tf = 1
 PROTOTYPE: $$$\@
 CODE:
  ihmsf = get_mortalspace(4,'i');
  slaDd2tf(ndp, days, &sign, ihmsf);
  unpack1D( (SV*)ST(3), (void *)ihmsf, 'i', 4);
 OUTPUT:
 sign
 ihmsf


# Equatorial to horizontal

void
slaDe2h(ha, dec, phi, az, el)
  double ha
  double dec
  double phi
  double az = NO_INIT
  double el = NO_INIT
 PROTOTYPE: $$$$$
 ALIAS:
  slaE2h = 1
 CODE:
  slaDe2h(ha, dec, phi, &az, &el);
 OUTPUT:
  az
  el


void
slaDeuler(order, phi, theta, psi, rmat)
  char * order
  double phi
  double theta
  double psi
  double * rmat = NO_INIT
 PROTOTYPE: $$$$\@
 ALIAS:
  slaEuler = 1
 CODE:
  rmat = get_mortalspace(9,'d');
  slaDeuler(order, phi, theta, psi, (void*)rmat);
  unpack1D( (SV*)ST(4), (void *)rmat, 'd', 9);
 OUTPUT:
  rmat


void
slaDfltin(string, nstrt, dreslt, jflag)
  char * string
  int nstrt
  double dreslt
  int jflag = NO_INIT
 PROTOTYPE: $$$$
 ALIAS:
  slaFloatin = 1
 CODE:
  slaDfltin(string, &nstrt, &dreslt, &jflag);
 OUTPUT:
  nstrt
  dreslt
  jflag

# Horizontal to equatorial

void
slaDh2e(az, el, phi, ha, dec)
  double az
  double el
  double phi
  double ha = NO_INIT
  double dec = NO_INIT
 PROTOTYPE: $$$$$
 ALIAS:
  slaH2e = 1
 CODE:
  slaDh2e(az, el, phi, &ha, &dec);
OUTPUT:
  ha
  dec

void
slaDimxv(dm, va, vb)
  double * dm
  double * va 
  double * vb = NO_INIT
 PROTOTYPE: \@\@\@
 CODE:
  vb = get_mortalspace(3,'d');
  slaDimxv((void*)dm, va, vb);
  unpack1D( (SV*)ST(2), (void *)vb, 'd', 3);
 OUTPUT: 
  vb


void
slaDjcal(ndp, djm, iymdf, j)
  int ndp
  double djm
  int * iymdf = NO_INIT
  int j
 PROTOTYPE: $$\@$
 CODE:
   iymdf =  get_mortalspace(4,'i');
   slaDjcal(ndp, djm, iymdf, &j);
   unpack1D( (SV*)ST(2), (void *)iymdf, 'i', 4);
 OUTPUT:
  iymdf
  j 

# MJD to UT

void
slaDjcl(mjd, iy, im, id, fd, j)
  double mjd
  int iy = NO_INIT
  int im = NO_INIT
  int id = NO_INIT
  double fd = NO_INIT
  int j = NO_INIT
 PROTOTYPE: $$$$$$
 CODE:
  slaDjcl(mjd, &iy, &im, &id, &fd, &j);
 OUTPUT:
  iy
  im
  id
  fd
  j

void
slaDm2av(rmat, axvec)
  double * rmat
  double * axvec = NO_INIT
 PROTOTYPE: \@\@
 ALIAS:
  slaM2av = 1
 CODE:
  axvec = get_mortalspace(3,'d');
  slaDm2av((void*)rmat, axvec);
  unpack1D( (SV*)ST(1), (void *)axvec, 'd', 3);
 OUTPUT:
  axvec


###### FLAG:   Do slaDmat at the end

void
slaDmoon(date, pv)
  double date
  double * pv = NO_INIT
 PROTOTYPE: $\@
 CODE:
   pv = get_mortalspace(6,'d');
   slaDmoon(date, pv);
   unpack1D( (SV*)ST(1), (void *)pv, 'd', 6);
 OUTPUT:
  pv


#### FLAG : Matrix manipulation should be using PDLs

void
slaDmxm(a, b, c)
  double * a
  double * b
  double * c = NO_INIT
 PROTOTYPE: \@\@\@
 ALIAS:
  slaMxm = 1
 CODE:
  c = get_mortalspace(9, 'd');
  slaDmxm((void*)a,(void*)b,(void*)c);
  unpack1D( (SV*)ST(2), (void *)c, 'd', 9);
 OUTPUT:
  c

void
slaDmxv(dm, va, vb)
  double * dm
  double * va
  double * vb = NO_INIT
 PROTOTYPE: \@\@\@
 ALIAS:
  slaMxv = 1
 CODE:
  vb = get_mortalspace(3, 'd');
  slaDmxv((void*)dm, va, vb);
  unpack1D( (SV*)ST(2), (void *)vb, 'd', 3);
 OUTPUT:
  vb


double
slaDpav(v1, v2)
  double * v1
  double * v2
 PROTOTYPE: \@\@
 ALIAS:
  slaPav = 1
 CODE:
  RETVAL = slaDpav(v1, v2);
 OUTPUT:
  RETVAL

#   Converts radians to HMS

void
slaDr2tf(ndp, angle, sign, ihmsf)
  int ndp
  double angle
  char sign = NO_INIT
  int * ihmsf = NO_INIT
 ALIAS:
  slaCr2tf = 1
 PROTOTYPE: $$$\@
 CODE:
  ihmsf = get_mortalspace(4,'i');
  slaDr2tf(ndp, angle, &sign, ihmsf);
  unpack1D( (SV*)ST(3), (void *)ihmsf, 'i', 4);
 OUTPUT:
 sign
 ihmsf

double
slaDrange(angle)
  double angle
 PROTOTYPE: $
 ALIAS:
  slaRange = 1
 CODE:
  RETVAL = slaDrange(angle);
 OUTPUT:
  RETVAL

double
slaDranrm(angle)
  double angle
 PROTOTYPE: $
 ALIAS:
  slaRanorm = 1
 CODE:
  RETVAL = slaDranrm(angle);
 OUTPUT:
  RETVAL


#   Converts radians to DMS 

void
slaDr2af(ndp, angle, sign, idmsf)
  int ndp
  double angle
  char sign = NO_INIT
  int * idmsf = NO_INIT
 ALIAS:
  slaCr2af = 1
 PROTOTYPE: $$$\@
 CODE:
  idmsf = get_mortalspace(4,'i');
  slaDr2af(ndp, angle, &sign, idmsf);
  unpack1D( (SV*)ST(3), (void *)idmsf, 'i', 4);
 OUTPUT:
 sign
 idmsf





void
slaDs2c6(a, b, r, ad, bd, rd, v)
  double a
  double b
  double r
  double ad
  double bd
  double rd
  double * v = NO_INIT
 ALIAS:
  slaCs2c6 = 1
 PROTOTYPE: $$$$$$\@
 CODE:
  v = get_mortalspace(6,'d');
  slaDs2c6(a, b, r, ad, bd, rd, v);
  unpack1D( (SV*)ST(6), (void *)v, 'd', 6);
 OUTPUT:
  v

void
slaDs2tp(ra, dec, raz, decz, xi, eta, j)
  double ra
  double dec
  double raz
  double decz
  double xi = NO_INIT
  double eta = NO_INIT
  int j = NO_INIT
 PROTOTYPE: $$$$$$$
 CODE:
  slaDs2tp(ra, dec, raz, decz, &xi, &eta, &j);
 OUTPUT:
  xi
  eta
  j

double
slaDsep(a1, b1, a2, b2)
  double a1
  double b1
  double a2
  double b2
 PROTOTYPE: $$$$
 ALIAS:
  slaSep = 1
 CODE:
  RETVAL = slaDsep(a1, b1, a2, b2);
 OUTPUT:
  RETVAL


double
slaDt(epoch)
  double epoch
 PROTOTYPE: $
 CODE:
  RETVAL = slaDt(epoch);
 OUTPUT:
  RETVAL



void
slaDtf2d(ihour, imin, sec, days, j)
  int ihour
  int imin
  double sec
  double  days = NO_INIT
  int  j = NO_INIT
 PROTOTYPE: $$$$$
 CODE:
  slaDtf2d(ihour, imin, sec, &days, &j);
 OUTPUT:
 days
 j


#  Converts HMS to radians 

void
slaDtf2r(ihour, imin, sec, rad, j)
  int ihour
  int imin
  double sec
  double  rad = NO_INIT
  int  j = NO_INIT
 PROTOTYPE: $$$$$
 CODE:
  slaDtf2r(ihour, imin, sec, &rad, &j);
 OUTPUT:
 rad
 j


void
slaDtp2s(xi, eta, raz, decz, ra, dec)
  double xi
  double eta
  double raz
  double decz
  double ra = NO_INIT
  double dec = NO_INIT
 PROTOTYPE: $$$$$$
 ALIAS:
  slaTp2s = 1
 CODE:
  slaDtp2s(xi, eta, raz, decz, &ra, &dec);
 OUTPUT:
  ra
  dec


void
slaDtp2v(xi, eta, v0, v)
  double xi
  double eta
  double * v0
  double * v = NO_INIT
 PROTOTYPE: $$\@\@
 ALIAS:
  slaTp2v = 1
 CODE:
  v = get_mortalspace(3, 'd');
  slaDtp2v(xi, eta, v0, v);
  unpack1D( (SV*)ST(3), (void *)v, 'd', 3);
 OUTPUT:
  v

void
slaDtps2c(xi, eta, ra, dec, raz1, decz1, raz2, decz2, n)
  double xi
  double eta
  double ra
  double dec
  double raz1 = NO_INIT
  double decz1 = NO_INIT
  double raz2 = NO_INIT
  double decz2 = NO_INIT
  int n = NO_INIT
 PROTOTYPE: $$$$$$$$$
 ALIAS:
  slaTps2c = 1
 CODE:
  slaDtps2c(xi, eta, ra, dec, &raz1, &decz1, &raz2, &decz2, &n);
 OUTPUT:
  raz1
  decz1
  raz2
  decz2
  n

void
slaDtpv2c(xi, eta, v, v01, v02, n)
  double xi
  double eta
  double * v
  double * v01 = NO_INIT
  double * v02 = NO_INIT
  int n = NO_INIT
 PROTOTYPE: $$\@\@\@
 ALIAS:
  slaTpv2c = 1
 CODE:
  v01 = get_mortalspace(3,'d');
  v02 = get_mortalspace(3,'d');
  slaDtpv2c(xi, eta, v, v01, v02, &n);
  unpack1D( (SV*)ST(3), (void *)v01, 'd', 3);
  unpack1D( (SV*)ST(4), (void *)v02, 'd', 3);
 OUTPUT:
  v01
  v02
  n


double
slaDtt(dju)
  double dju
 PROTOTYPE: $
 CODE:
  RETVAL = slaDtt(dju);
 OUTPUT:
  RETVAL

void
slaDv2tp(v, v0, xi, eta, j)
  double * v
  double * v0
  double xi = NO_INIT
  double eta = NO_INIT
  int j = NO_INIT
 PROTOTYPE: \@\@$$$
 ALIAS:
  slaV2tp = 1
 CODE:
  slaDv2tp(v, v0, &xi, &eta, &j);
 OUTPUT:
  xi
  eta
  j

double
slaDvdv(va, vb)
  double * va
  double * vb
 PROTOTYPE: \@\@
 ALIAS:
  slaVdv = 1
 CODE:
   RETVAL = slaDvdv(va, vb);
 OUTPUT:
  RETVAL

void
slaDvn(v, uv, vm)
  double * v
  double * uv = NO_INIT
  double vm
 PROTOTYPE: \@\@$
 ALIAS:
  slaVn = 1
 CODE:
  uv = get_mortalspace(3,'d');
  slaDvn(v, uv, &vm);
  unpack1D( (SV*)ST(1), (void *)uv, 'd', 3);
 OUTPUT:
  uv
  vm

void
slaDvxv(va, vb, vc)
  double * va
  double * vb
  double * vc = NO_INIT
 PROTOTYPE: \@\@\@
 ALIAS:
  slaVxv = 1
 CODE:
  vc = get_mortalspace(3,'d');
  slaDvxv(va,vb,vc);
  unpack1D( (SV*)ST(2), (void *)vc, 'd', 3);
 OUTPUT:
  vc

#### slaE2h - use Double precision


void
slaEarth(iy, id, fd, pv)
  int iy
  int id
  float fd
  float * pv = NO_INIT
 PROTOTYPE: $$$\@
 CODE:
   pv = get_mortalspace(6,'f');
   slaEarth(iy, id, fd, pv);
   unpack1D( (SV*)ST(3), (void *)pv, 'f', 6);
 OUTPUT:
  pv
  
void
slaEcleq(dl, db, date, dr, dd)
  double dl
  double db
  double date
  double dr = NO_INIT
  double dd = NO_INIT
 PROTOTYPE: $$$$$
 CODE:
  slaEcleq(dl, db, date, &dr, &dd);
 OUTPUT:
  dr
  dd

void
slaEcmat(date, rmat)
  double date
  double * rmat
 PROTOTYPE: $\@
 CODE:
  rmat = get_mortalspace(9,'d');
  slaEcmat(date, (void*)rmat);
  unpack1D( (SV*)ST(1), (void *)rmat, 'd', 9);
 OUTPUT:
  rmat

void
slaEcor(rm, dm, iy, id, fd, rv, tl)
  float rm
  float dm
  int iy
  int id
  float fd
  float rv = NO_INIT
  float tl = NO_INIT
 PROTOTYPE: $$$$$$$
 CODE:
  slaEcor(rm, dm, iy, id, fd, &rv, &tl);
 OUTPUT:
  rv
  tl

void
slaEg50(dr, dd, dl, db)
  double dr
  double dd
  double dl = NO_INIT
  double db = NO_INIT
 PROTOTYPE: $$$$
 CODE: 
   slaEg50(dr, dd, &dl, &db);
 OUTPUT:
  dl
  db


double
slaEpb(date)
  double date
 PROTOTYPE: $
 CODE:
  RETVAL = slaEpb(date);
 OUTPUT:
  RETVAL

double
slaEpb2d(epb)
  double epb
 PROTOTYPE: $
 CODE:
  RETVAL = slaEpb2d(epb);
 OUTPUT:
  RETVAL

double
slaEpco(k0, k, e)
  char  k0
  char  k
  double e
 PROTOTYPE: $$$
 CODE:
  RETVAL = slaEpco(k0, k, e);
 OUTPUT:
  RETVAL

double
slaEpj(date)
  double date
 PROTOTYPE: $
 CODE:
  RETVAL = slaEpj(date);
 OUTPUT:
  RETVAL


double
slaEpj2d(epj)
  double epj
 PROTOTYPE: $
 CODE:
  RETVAL = slaEpj2d(epj);
 OUTPUT:
  RETVAL

void
slaEqecl(dr, dd, date, dl, db)
  double dr
  double dd
  double date
  double dl = NO_INIT
  double db = NO_INIT
 PROTOTYPE: $$$$$
 CODE:
   slaEqecl(dr, dd, date, &dl, &db);
 OUTPUT:
  dl
  db

# Equation of the equinoxes

double
slaEqeqx(date)
  double date
 PROTOTYPE: $
 CODE:
  RETVAL = slaEqeqx(date);
 OUTPUT:
  RETVAL

void
slaEqgal(dr, dd, dl, db)
  double dr
  double dd
  double dl = NO_INIT
  double db = NO_INIT
 PROTOTYPE: $$$$
 CODE: 
   slaEqgal(dr, dd, &dl, &db);
 OUTPUT:
  dl
  db


void
slaEtrms(ep, ev)
  double ep
  double * ev
 PROTOTYPE: $\@
 CODE:
  ev = get_mortalspace(3, 'd');
  slaEtrms(ep, ev);
  unpack1D( (SV*)ST(1), (void *)ev, 'd', 3);
 OUTPUT:
  ev


#### FLAG:: slaEuler skipped in favcour of double prec version


void
slaEvp(date, deqx, dvb, dpb, dvh, dph)
  double date
  double deqx
  double * dvb = NO_INIT
  double * dpb = NO_INIT
  double * dvh = NO_INIT
  double * dph = NO_INIT
  PROTOTYPE: $$\@\@\@\@
  CODE:
   dvb = get_mortalspace(3,'d');
   dpb = get_mortalspace(3,'d');
   dvh = get_mortalspace(3,'d');
   dph = get_mortalspace(3,'d');
   slaEvp(date, deqx, dvb, dpb, dvh, dph);
   unpack1D( (SV*)ST(2), (void *)dvb, 'd', 3);
   unpack1D( (SV*)ST(3), (void *)dpb, 'd', 3);
   unpack1D( (SV*)ST(4), (void *)dvh, 'd', 3);
   unpack1D( (SV*)ST(5), (void *)dph, 'd', 3);
  OUTPUT:
  dvb
  dpb
  dvh
  dph

##### FLAG: Do slaFitxy some other time

void
slaFk425(r1950,d1950,dr1950,dd1950,p1950,v1950,r2000,d2000,dr2000,dd2000,p2000,v2000)
  double r1950
  double d1950
  double dr1950
  double dd1950
  double p1950
  double v1950
  double r2000 = NO_INIT
  double d2000 = NO_INIT
  double dr2000 = NO_INIT
  double dd2000 = NO_INIT
  double p2000 = NO_INIT
  double v2000 = NO_INIT
 PROTOTYPE: $$$$$$$$$$$$
 CODE:
  slaFk425(r1950,d1950,dr1950,dd1950,p1950,v1950,&r2000,&d2000,&dr2000,&dd2000,&p2000,&v2000);
 OUTPUT:
  r2000
  d2000 
  dr2000
  dd2000
  p2000
  v2000




#  B1950 to J2000

void
slaFk45z(r1950, d1950, bepoch, r2000, d2000)
  double r1950
  double d1950
  double bepoch
  double r2000 = NO_INIT
  double d2000 = NO_INIT
 PROTOTYPE: $$$$$
 CODE:
  slaFk45z(r1950, d1950, bepoch, &r2000, &d2000);
 OUTPUT:
 r2000
 d2000


void
slaFk524(r2000,d2000,dr2000,dd2000,p2000,v2000,r1950,d1950,dr1950,dd1950,p1950,v1950)
  double r2000
  double d2000
  double dr2000
  double dd2000
  double p2000
  double v2000
  double r1950 = NO_INIT
  double d1950 = NO_INIT
  double dr1950 = NO_INIT
  double dd1950 = NO_INIT
  double p1950 = NO_INIT
  double v1950 = NO_INIT
 PROTOTYPE: $$$$$$$$$$$$
 CODE:
  slaFk524(r2000,d2000,dr2000,dd2000,p2000,v2000,&r1950,&d1950,&dr1950,&dd1950,&p1950,&v1950);
 OUTPUT:
  r1950
  d1950 
  dr1950
  dd1950
  p1950
  v1950

void
slaFk54z(r2000, d2000, bepoch, r1950, d1950, dr1950, dd1950)
  double r2000
  double d2000
  double bepoch
  double r1950 = NO_INIT
  double d1950 = NO_INIT
  double dr1950 = NO_INIT
  double dd1950 = NO_INIT
 PROTOTYPE: $$$$$$$
 CODE:
  slaFk54z(r2000, d2000, bepoch, &r1950, &d1950, &dr1950, &dd1950);
 OUTPUT:
 r1950
 d1950
 dr1950
 dd1950


##### FLAG: SKIP slaFloatin - use slaDfltin instead

void
slaGaleq(dl, db, dr, dd)
  double dl
  double db
  double dr = NO_INIT
  double dd = NO_INIT
 PROTOTYPE: $$$$
 CODE: 
   slaGaleq(dl, db, &dr, &dd);
 OUTPUT:
  dr
  dd


void
slaGalsup(dl, db, dsl, dsb)
  double dl
  double db
  double dsl = NO_INIT
  double dsb = NO_INIT
 PROTOTYPE: $$$$
 CODE: 
   slaGalsup(dl, db, &dsl, &dsb);
 OUTPUT:
  dsl
  dsb

void
slaGe50(dl, db, dr, dd)
  double dl
  double db
  double dr = NO_INIT
  double dd = NO_INIT
 PROTOTYPE: $$$$
 CODE: 
   slaGe50(dl, db, &dr, &dd);
 OUTPUT:
  dr
  dd


void
slaGeoc(p, h, r, z)
  double p
  double h
  double r = NO_INIT
  double z = NO_INIT
 PROTOTYPE: $$$$
 CODE: 
   slaGeoc(p, h, &r, &z);
 OUTPUT:
  r
  z


# UT to GMST

double
slaGmst(ut1)
  double ut1
 PROTOTYPE: $
 CODE:
  RETVAL = slaGmst(ut1);
 OUTPUT:
  RETVAL


double
slaGmsta(date, ut1)
  double date
  double ut1
 PROTOTYPE: $$
 CODE:
  RETVAL = slaGmsta(date, ut1);
 OUTPUT:
  RETVAL


# slaGresid is not in C version

float
slaGresid(s)
  float s
 PROTOTYPE: $
 CODE:
  /* RETVAL = slaGresid(s); */
  croak("NOT implemented: slaGresid is not implemented\n");
 OUTPUT:
  RETVAL


##### SKIP::  slaH2e use slaDh2e instead

void
slaImxv(rm, va, vb)
  float * rm
  float * va 
  float * vb = NO_INIT
 PROTOTYPE: \@\@\@
 CODE:
  vb = get_mortalspace(3,'f');
  slaImxv((void*)rm, va, vb);
  unpack1D( (SV*)ST(2), (void *)vb, 'f', 3);
 OUTPUT: 
  vb


##### SKIP: slaIntin for now -- does perl need it?

void
slaInvf(fwds, bkwds, j)
  double * fwds
  double * bkwds = NO_INIT
  int j = NO_INIT
 PROTOTYPE: \@\@$
 CODE:
  bkwds = get_mortalspace(6,'d');
  slaInvf(fwds, bkwds, &j);
  unpack1D( (SV*)ST(1), (void *)bkwds, 'd', 6);
 OUTPUT:
  bkwds
  j


void
slaKbj(jb, e, k, j)
  int jb
  double e
  char * k = NO_INIT
  int j = NO_INIT
 PROTOTYPE: $$$$
 PREINIT:
  char string[256];
 CODE:
  k = string;
  slaKbj(jb, e, k, &j);
 OUTPUT:
  k
  j


#### SKIP:: slaM2av - use slaDm2av

void
slaMap(rm, dm, pr, pd, px, rv, eq, date, ra, da)
  double rm
  double dm
  double pr
  double pd
  double px
  double rv
  double eq
  double date
  double ra = NO_INIT
  double da = NO_INIT
 PROTOTYPE: $$$$$$$$$$
 CODE:
  slaMap(rm, dm, pr, pd, px, rv, eq, date, &ra, &da);
 OUTPUT: 
  ra
  da


void
slaMappa(eq, date, amprms)
  double eq
  double date
  double * amprms = NO_INIT
 PROTOTYPE: $$\@
 CODE:
  amprms = get_mortalspace(21, 'd');
  slaMappa(eq, date, amprms);
  unpack1D( (SV*)ST(2), (void *)amprms, 'd', 21); 
 OUTPUT:
  amprms

void
slaMapqk(rm, dm, pr, pd, px, rv, amprms, ra, da)
   double rm
  double dm
  double pr
  double pd
  double px
  double rv
  double * amprms
  double ra = NO_INIT
  double da = NO_INIT
 PROTOTYPE: $$$$$$\@$$
 CODE:
  slaMapqk(rm, dm, pr, pd, px, rv, amprms, &ra, &da);
 OUTPUT: 
  ra
  da

void
slaMapqkz(rm, dm, amprms, ra, da)
  double rm
  double dm
  double * amprms
  double ra = NO_INIT
  double da = NO_INIT
 PROTOTYPE: $$\@$$
 CODE:
   slaMapqkz(rm, dm, amprms, &ra, &da);
 OUTPUT:
  ra
  da
 

void
slaMoon(iy, id, fd, pv)
  int iy
  int id
  float fd
  float * pv = NO_INIT
 PROTOTYPE: $$$\@
 CODE:
   pv = get_mortalspace(6,'f');
   slaMoon(iy, id, fd, pv);
   unpack1D( (SV*)ST(3), (void *)pv, 'f', 6);
 OUTPUT:
  pv


#### FLAG: Miss slaMxm use slaDmxm instead

#### FLAG: Miss slaMxv use slaDmxv instead


void
slaNut(date, rmatn)
  double date
  double * rmatn = NO_INIT
 PROTOTYPE: $\@
 CODE:
  rmatn = get_mortalspace(9, 'd');
  slaNut(date, (void*)rmatn);
  unpack1D( (SV*)ST(1), (void *)rmatn, 'd', 9);
 OUTPUT:
  rmatn


void
slaNutc(date, dpsi, deps, eps0)
  double date
  double dpsi = NO_INIT
  double deps = NO_INIT
  double eps0 = NO_INIT
 PROTOTYPE: $$$$
 CODE:
  slaNutc(date, &dpsi, &deps, &eps0);
 OUTPUT:
  dpsi
  deps
  eps0

void
slaOap(type, ob1, ob2, date, dut, elongm, phim, hm, xp, yp, tdk, pmb, rh, wl, tlr, rap, dap)
  char * type
  double ob1
  double ob2
  double date
  double dut
  double elongm
  double phim
  double hm
  double xp
  double yp
  double tdk
  double pmb
  double rh
  double wl
  double tlr
  double rap = NO_INIT
  double dap = NO_INIT
 PROTOTYPE: $$$$$$$$$$$$$$$
 CODE:
   slaOap(type, ob1, ob2, date, dut, elongm, phim, hm, xp, yp, tdk, pmb, rh, wl, tlr, &rap, &dap);
 OUTPUT:
  rap
  dap

void
slaOapqk(type, ob1, ob2, aoprms, rap, dap)
  char * type
  double ob1
  double ob2
  double * aoprms
  double rap = NO_INIT
  double dap = NO_INIT
 PROTOTYPE: $$$\@$$
 CODE:
  slaOapqk(type, ob1, ob2, aoprms, &rap, &dap);
 OUTPUT:
  rap
  dap


# If c is undef we need to convert to a blank since slalib
# will generate segmentation violation if it receieves an undef
# value for the string (the strcpy fails for some reason).
# overcome this by providing a wrapper in the .pm file to check
# for this case. Have not got the time to work out a fix at the
# XS level. 'c' can be an input or output variable but must be
# guaranteed to contain a valid pointer to char.
#  slaObs is now defined in the .pm file

void
_slaObs(n, c, name, w, p, h)
  int n
  char * c
  char * name = NO_INIT
  double w = NO_INIT
  double p = NO_INIT
  double h = NO_INIT
 PROTOTYPE: $$$$$$
 PREINIT:
  char string[40];
 CODE:
  name = string;
  slaObs(n, c, name, &w, &p, &h);
 OUTPUT:
  c
  name
  w
  p
  h


double
slaPa(ha, dec, phi)
  double ha
  double dec
  double phi
 PROTOTYPE: $$$
 CODE:
  RETVAL = slaPa(ha, dec, phi);
 OUTPUT:
  RETVAL


#### SKIP: slaPav use slaDpav instead (is alias).

void
slaPcd(disco, x, y)
  double disco
  double x
  double y
 PROTOTYPE: $$$
 CODE:
  slaPcd(disco, &x, &y);
 OUTPUT:
  x
  y


void
slaPda2h(p, d, a, h1, j1, h2, j2)
  double p
  double d
  double a
  double h1  = NO_INIT
  int j1 = NO_INIT
  double h2 = NO_INIT
  int j2 = NO_INIT
 PROTOTYPE: $$$$$$$
 CODE:
  slaPda2h(p, d, a, &h1, &j1, &h2, &j2);
 OUTPUT:
  h1
  j1
  h2
  j2


void
slaPdq2h(p, d, q, h1, j1, h2, j2)
  double p
  double d
  double q
  double h1  = NO_INIT
  int j1 = NO_INIT
  double h2 = NO_INIT
  int j2 = NO_INIT
 PROTOTYPE: $$$$$$$
 CODE:
  slaPdq2h(p, d, q, &h1, &j1, &h2, &j2);
 OUTPUT:
  h1
  j1
  h2
  j2

void
slaPertel(jform,date0,date1,epoch0,orbi0,anode0,perih0,aorq0,e0,am0,epoch1,orbi1,anode1,perih1,aorq1,e1,am1,jstat)
  int jform
  double date0
  double date1
  double epoch0
  double orbi0
  double anode0
  double perih0
  double aorq0
  double e0
  double am0
  double epoch1 = NO_INIT
  double orbi1 = NO_INIT
  double anode1 = NO_INIT
  double perih1 = NO_INIT
  double aorq1 = NO_INIT
  double e1 = NO_INIT
  double am1 = NO_INIT
  int    jstat = NO_INIT
 PROTOTYPE: $$$$$$$$$$$$$$$$$$
 CODE:
  slaPertel(jform,date0,date1,epoch0,orbi0,anode0,perih0,aorq0,e0,am0,&epoch1,&orbi1,&anode1,&perih1,&aorq1,&e1,&am1,&jstat);
 OUTPUT:
  epoch1
  orbi1
  anode1
  perih1
  aorq1
  e1
  am1
  jstat

void
slaPertue(date,u,jstat)
  double date
  double * u
  int    jstat = NO_INIT
 PROTOTYPE: $\@$
 CODE:
  slaPertue(date,u,&jstat);
  unpack1D( (SV*)ST(1), (void *)u, 'd', 13);
 OUTPUT:
  u
  jstat


void
slaPlanel(date, jform, epoch, orbinc, anode, perih, aorq, e, aorl, dm, pv, jstat)
  double date
  int jform
  double epoch
  double orbinc
  double anode
  double perih
  double aorq
  double e
  double aorl
  double dm
  double * pv = NO_INIT
  int jstat = NO_INIT
 PROTOTYPE: $$$$$$$$$$\@$
 CODE:
  pv = get_mortalspace(6, 'd');
  slaPlanel(date, jform, epoch, orbinc, anode, perih, aorq, e, aorl, dm, pv, &jstat);
  unpack1D( (SV*)ST(10), (void *)pv, 'd', 6);
 OUTPUT:
  pv
  jstat

void
slaPlanet(date, np, pv, jstat)
  double date
  int np
  double * pv = NO_INIT
  int jstat = NO_INIT
 PROTOTYPE: $$\@$
 CODE:
   pv = get_mortalspace(6, 'd');
   slaPlanet(date, np, pv, &jstat);
   unpack1D( (SV*)ST(2), (void *)pv, 'd', 6);
 OUTPUT:
  pv
  jstat

void
slaPlante(date, elong, phi, jform, epoch, orbinc, anode, perih, aorq,e, aorl, dm, ra,dec, r, jstat)
  double date
  double elong
  double phi
  int jform
  double epoch
  double orbinc
  double anode
  double perih
  double aorq
  double e
  double aorl
  double dm
  double ra = NO_INIT
  double dec = NO_INIT
  double r = NO_INIT
  int jstat = NO_INIT
 PROTOTYPE: $$$$$$$$$$$$$$$$
 CODE:
  slaPlante(date, elong, phi, jform, epoch, orbinc, anode, perih, aorq,e, aorl, dm, &ra, &dec, &r, &jstat);
 OUTPUT:
  ra
  dec
  r
  jstat


void
slaPm(r0,d0,pr,pd,px,rv,ep0,ep1,r1,d1)
  double r0
  double d0
  double pr
  double pd
  double px
  double rv
  double ep0
  double ep1
  double r1 = NO_INIT
  double d1 = NO_INIT
 PROTOTYPE: $$$$$$$$$$
 CODE:
  slaPm(r0,d0,pr,pd,px,rv,ep0,ep1,&r1,&d1);
 OUTPUT:
  r1
  d1


void
slaPolmo(elongm, phim, xp, yp, elong, phi, daz)
  double elongm
  double phim
  double xp
  double yp
  double elong = NO_INIT
  double phi = NO_INIT
  double daz = NO_INIT
 PROTOTYPE: $$$$$$$
 CODE:
  slaPolmo(elongm, phim, xp, yp, &elong, &phi, &daz);
 OUTPUT:
  elong
  phi
  daz


##### Problem with slaPrebn - dont know the return args
# Think it is (3,3) - see slaPrec
void
slaPrebn(bep0, bep1, rmatp)
  double bep0
  double bep1
  double * rmatp
 PROTOTYPE: $$\@
 CODE:
  rmatp = get_mortalspace(9,'d');
  slaPrebn(bep0, bep1, (void*)rmatp);
  unpack1D( (SV*)ST(2), (void *)rmatp, 'd', 9);
 OUTPUT:
  rmatp


void
slaPrec(ep0, ep1, rmatp)
  double ep0
  double ep1
  double * rmatp
 PROTOTYPE: $$\@
 CODE:
  rmatp = get_mortalspace(9,'d');
  slaPrec(ep0, ep1, (void*)rmatp);
  unpack1D( (SV*)ST(2), (void *)rmatp, 'd', 9);
 OUTPUT:
  rmatp


# Precession
 
void
slaPreces(system, ep0, ep1, ra, dc)
  char *system
  double ep0
  double ep1
  double ra
  double dc
 PROTOTYPE: $$$$$
 CODE:
  slaPreces(system, ep0, ep1, &ra, &dc);
 OUTPUT:
 ra
 dc


void
slaPrecl(ep0, ep1, rmatp)
  double ep0
  double ep1
  double * rmatp
 PROTOTYPE: $$\@
 CODE:
  rmatp = get_mortalspace(9,'d');
  slaPrecl(ep0, ep1, (void*)rmatp);
  unpack1D( (SV*)ST(2), (void *)rmatp, 'd', 9);
 OUTPUT:
  rmatp

void
slaPrenut(epoch, date, rmatpn)
  double epoch
  double date
  double * rmatpn
 PROTOTYPE: $$\@
 CODE:
  rmatpn = get_mortalspace(9,'d');
  slaPrenut(epoch, date, (void*)rmatpn);
  unpack1D( (SV*)ST(2), (void *)rmatpn, 'd', 9);
 OUTPUT:
  rmatpn


void
slaPvobs(p, h, stl, pv)
  double p
  double h
  double stl
  double * pv = NO_INIT
 PROTOTYPE: $$$\@
 CODE:
   pv = get_mortalspace(6, 'd');
   slaPvobs(p, h, stl, pv);
   unpack1D( (SV*)ST(3), (void *)pv, 'd', 6);
 OUTPUT:
  pv


###### Skip slaPxy - do later


##### slaRandom is not implemented
float
slaRandom(seed)
  float seed
 PROTOTYPE: $
 CODE:
  /* RETVAL = slaRandom(&seed); */
  croak("NOT implemented: slaRandom is not implemented\n");
 OUTPUT:
  RETVAL
  seed

##### Skip: slaRange - use slaDrange
##### Skip: slaRanorm  use slaDranrm


double
slaRcc(tdb, ut1, wl, u, v)
  double tdb
  double ut1
  double wl
  double u
  double v
 PROTOTYPE: $$$$$
 CODE:
  RETVAL = slaRcc(tdb, ut1, wl, u, v);
 OUTPUT:
  RETVAL


void
slaRdplan(date, np, elong, phi, ra, dec, diam)
  double date
  int np
  double elong
  double phi
  double ra = NO_INIT
  double dec = NO_INIT
  double diam = NO_INIT
 PROTOTYPE: $$$$$$$
 CODE:
  slaRdplan(date, np, elong, phi, &ra, &dec, &diam);
 OUTPUT:
  ra
  dec
  diam


void
slaRefco(hm, tdk, pmb, rh, wl, phi, tlr, eps, refa, refb)
  double hm
  double tdk
  double pmb
  double rh
  double wl
  double phi
  double tlr
  double eps
  double refa = NO_INIT
  double refb = NO_INIT
 PROTOTYPE: $$$$$$$$$$
 CODE:
  slaRefco(hm, tdk, pmb, rh, wl, phi, tlr, eps, &refa, &refb);
 OUTPUT:
  refa
  refb

void
slaRefcoq(tdk, pmb, rh, wl, refa, refb)
  double tdk
  double pmb
  double rh
  double wl
  double refa = NO_INIT
  double refb = NO_INIT
 PROTOTYPE: $$$$$$
 CODE:
  slaRefcoq(tdk, pmb, rh, wl, &refa, &refb);
 OUTPUT:
  refa
  refb



void
slaRefro(zobs, hm, tdk, pmb, rh, wl, phi, tlr, eps, ref)
  double zobs
  double hm
  double tdk
  double pmb
  double rh
  double wl
  double phi
  double tlr
  double eps
  double ref = NO_INIT 
 PROTOTYPE: $$$$$$$$$$
 CODE:
  slaRefro(zobs, hm, tdk, pmb, rh, wl, phi, tlr, eps, &ref);
 OUTPUT:
  ref

void
slaRefv(vu, refa, refb, vr)
  double * vu
  double refa
  double refb
  double * vr = NO_INIT
 PROTOTYPE:  \@$$\@
 CODE:
  vr = get_mortalspace(3, 'd');
  slaRefv(vu, refa, refb, vr);
  unpack1D( (SV*)ST(3), (void *)vr, 'd', 3);
 OUTPUT:
  vr

void
slaRefz(zu, refa, refb, zr)
  double zu
  double refa
  double refb
  double zr = NO_INIT
 PROTOTYPE: $$$$
 CODE:
   slaRefz(zu, refa, refb, &zr);
 OUTPUT:
  zr

 
float
slaRverot(phi, ra, da, st)
  float phi
  float ra
  float da
  float st
 PROTOTYPE: $$$$
 CODE:
  RETVAL = slaRverot(phi, ra, da, st);
 OUTPUT:
  RETVAL


float
slaRvgalc(r2000, d2000)
  float r2000
  float d2000
 PROTOTYPE: $$
 CODE:
  RETVAL = slaRvgalc(r2000, d2000);
 OUTPUT:
  RETVAL

float
slaRvlg(r2000, d2000)
  float r2000
  float d2000
 PROTOTYPE: $$
 CODE:
  RETVAL = slaRvlg(r2000, d2000);
 OUTPUT:
  RETVAL


float
slaRvlsrd(r2000, d2000)
  float r2000
  float d2000
 PROTOTYPE: $$
 CODE:
  RETVAL = slaRvlsrd(r2000, d2000);
 OUTPUT:
  RETVAL

float
slaRvlsrk(r2000, d2000)
  float r2000
  float d2000
 PROTOTYPE: $$
 CODE:
  RETVAL = slaRvlsrk(r2000, d2000);
 OUTPUT:
  RETVAL

void
slaS2tp(ra, dec, raz, decz, xi, eta, j)
  float ra
  float dec
  float raz
  float decz
  float xi = NO_INIT
  float eta = NO_INIT
  int   j = NO_INIT
 PROTOTYPE: $$$$$$$
 CODE: 
  slaS2tp(ra, dec, raz, decz, &xi, &eta, &j);
 OUTPUT:
  xi
  eta
  j

###### SKIP slaSep - use SlaDsep instead

###### Skip slaSmat

void
slaSubet(rc, dc, eq, rm, dm)
  double rc
  double dc
  double eq
  double rm = NO_INIT
  double dm = NO_INIT
 PROTOTYPE: $$$$$
 CODE: 
  slaSubet(rc, dc, eq, &rm, &dm);
 OUTPUT: 
  rm
  dm

void
slaSupgal(dsl, dsb, dl, db)
  double dsl
  double dsb
  double dl = NO_INIT
  double db = NO_INIT
 PROTOTYPE: $$$$
 CODE: 
   slaSupgal(dsl, dsb, &dl, &db);
 OUTPUT:
  dl
  db


##### SDkip slaSVD

###### Skip slaSvdcov

###### Skip slaSvdsol

##### Skip slaTp2s - use slaDtp2s
##### Skip slaTp2v - use slaDtp2v
##### Skip slaTps2c - use slaDtps2c
##### Skip slaTpv2c - use slaDtpv2c

void
slaUnpcd(disco, x, y)
  double disco
  double x
  double y
 PROTOTYPE: $$$
 CODE:
  slaUnpcd(disco, &x, &y);
 OUTPUT:
  x
  y


##### Skip slaV2tp - use slaDv2tp
##### Skip slaVdv - use slaDvdv
##### Skip slaVn - use slaDvn
##### Skip slaVxv - use slaDvxv

# slaWait Not in C library -- implemented in perl via select()

##### Skip slaXy2xy for now

void
slaXy2xy(x1, y1, coeffs, x2, y2)
  double x1
  double y1
  double * coeffs
  double x2 = NO_INIT
  double y2 = NO_INIT
 PROTOTYPE: $$\@$$
 CODE: 
  slaXy2xy(x1, y1, coeffs, &x2, &y2);
 OUTPUT:
  x2
  y2


double
slaZd(ha, dec, phi)
  double ha
  double dec
  double phi
 PROTOTYPE: $$$
 CODE:
  RETVAL = slaZd(ha, dec, phi);
 OUTPUT:
  RETVAL
