use ExtUtils::MakeMaker;
use strict;
use warnings;

# A makefile can be generated for either the C or Fortran version
# of SLALIB. Change this flag to control the behaviour of the build.
# It does not yet determine the C/Fortran-ness of your library
# automatically.
# Defaults to using the C interface since that is the most complete
# implementation (albeit the C library is proprietary)
my $use_fortran = 1;

# Location of your slalib library (of whichever version)
# Please update as appropriate
my $sla_lib = '-L/home/timj/slalib/lib -lsla';

$sla_lib = " -L/home/timj/star/lib -lsla ";

# The #defines, libraries and any additional prerequisites
my ($defines, $libs, %prereqs);
if ($use_fortran) {
  print "Attempting to build against a Fortran SLALIB: $sla_lib\n";

  # Must switch modes
  $defines = " -DUSE_FORTRAN";

  # Must include fortran libraries but we do not want to force
  # ExtUtils::F77 on an unsuspecting C user
  $prereqs{'ExtUtils::F77'} = 0;

  # use eval here since we are listing ExtUtils as a prerequisite
  # anyway and so if this does not work we'll still get a message
  # later
  eval { require ExtUtils::F77; ExtUtils::F77->import };
  my $flibs = '';
  if ($@) {
    print "Error determining Fortran library requirements\n";
  } else {
    $flibs = ExtUtils::F77->runtime;
    $defines .= " -DHAS_UNDERSCORE"
      if ExtUtils::F77->trail_;
  }

  $libs = "$sla_lib $flibs";

} else {
  print "Attempting to build against a C SLALIB library: $sla_lib\n";

  # no special defines for the C version
  $defines = '';

  # Just need the math library
  $libs = "$sla_lib -lm";
}

WriteMakefile(
	      'NAME'      => 'Astro::SLA',
	      'VERSION_FROM'   => 'SLA.pm',
              'EXE_FILES' => [ 'stime' ],
	      'LIBS'      => [ $libs ],
	      'DEFINE'    => $defines,
              'PREREQ_PM' => {
			      'Test' => 0,
			      %prereqs,
			     },
	      'INC'       => '',
	      'dist'      => { COMPRESS => "gzip -9f"},
	      ($] >= 5.005 ?    ## Add these new keywords supported since 5.005
	       (ABSTRACT_FROM => 'SLA.pm',
		AUTHOR     => 'Tim Jenness <t.jenness@jach.hawaii.edu>') : ()),
	     );

