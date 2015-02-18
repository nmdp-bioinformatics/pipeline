#!/usr/bin/env perl
use strict;
use warnings;
use ExtUtils::MakeMaker;

WriteMakefile (
       NAME         => 'ngs-validation-report',
       VERSION      => "0.3",
       EXE_FILES    => [ qw(ngs-validation-report ) ]
);