#!/usr/bin/env perl
use strict;
use warnings;
use ExtUtils::MakeMaker;

WriteMakefile (
       NAME         => 'ngs-validation-report',
       VERSION      => "1.1",
       DISTNAME     => "backup",
       EXE_FILES    => [ qw(ngs-validation-report ) ]
);