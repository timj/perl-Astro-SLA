use ExtUtils::MakeMaker;

WriteMakefile(
	      'NAME'      => 'Astro::SLA',
	      'VERSION_FROM'   => 'SLA.pm',
	      'LIBS'      => [
			      '-L/usr/local/lib -lcsla -lm',
			     ],
	      'INC'       => '',
	      'dist'      => { COMPRESS => "gzip -9f"},
	      ($] >= 5.005 ?    ## Add these new keywords supported since 5.005
	       (ABSTRACT_FROM => 'SLA.pm',
		AUTHOR     => 'Tim Jenness <t.jenness@jach.hawaii.edu>') : ()),
	     );

