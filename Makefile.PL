use ExtUtils::MakeMaker;

WriteMakefile(
	      'NAME'      => 'Astro::SLA',
	      'VERSION_FROM'   => 'SLA.pm',
	      'LIBS'      => [
			      '-L/usr/local/lib -lcsla -lm',
			     ],
	      'INC'       => '',
	      'dist'      => { COMPRESS => "gzip -9f"}
	     );

