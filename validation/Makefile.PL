#!/usr/bin/env perl
use strict;
use warnings;
use ExtUtils::MakeMaker;

WriteMakefile (
       NAME         => 'ngs-validation-report',
       VERSION      => "1.0",
       EXE_FILES    => [ qw(ngs-validation-report ) ]
);