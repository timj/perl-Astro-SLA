use ExtUtils::MakeMaker;
WriteMakefile(
    'NAME'      => 'Astro::SLA',
    'VERSION_FROM'   => 'SLA.pm',
    'LIBS'      => ['
-L/local/lib -lcsla -L/opt/SUNWspro/lib -lsunmath -lm
'],
    'DEFINE'    => '',
    'INC'       => ' ',
    );

