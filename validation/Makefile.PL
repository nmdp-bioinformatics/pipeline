#!/usr/bin/env perl
=head1 LICENSE

    pipeline  Consensus assembly and allele interpretation pipeline.
    Copyright (c) 2015 National Marrow Donor Program (NMDP)

    This library is free software; you can redistribute it and/or modify it
    under the terms of the GNU Lesser General Public License as published
    by the Free Software Foundation; either version 3 of the License, or (at
    your option) any later version.

    This library is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; with out even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
    License for more details.
 
    You should have received a copy of the GNU Lesser General Public License
    along with this library;  if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA.

    > http://www.gnu.org/licenses/lgpl.html   
    
=cut
use strict;
use warnings;
use ExtUtils::MakeMaker;

WriteMakefile (

      NAME         => 'Validation',
      VERSION_FROM => 'ngs-validation-report',
      AUTHOR       => 'Mike Halagan mhalagan\@nmdp.org',

      DISTNAME     => 'validation-report',

      EXE_FILES    => [ qw(ngs-validation-report ngs-imgt-db) ],
      PMLIBDIRS    => [ qw(template lib) ],


      PREREQ_PM => {
        "LWP::UserAgent" => 0,
        "strict" => 0,
        "warnings" => 0,
        "Data::Dumper" => 0,
        "IO::Uncompress::AnyUncompress" => 0
      },

      TEST_REQUIRES => {
        "Test::More" => "0.96",
        "HTML::TreeBuilder" => 0,
        "strict" => 0,
        "warnings" => 0,
        "Data::Dumper" => 0
      },

      test => { TESTS => 't/*.t' },

      META_MERGE => {
          resources => {
            repository  => 'https://github.com/nmdp-bioinformatics/pipeline',
            bugtracker  => 'https://github.com/nmdp-bioinformatics/pipeline/issues',
          }
      }

);

