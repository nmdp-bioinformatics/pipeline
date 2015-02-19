#!/usr/bin/env perl
=head1 NAME

    report.t

=head1 SYNOPSIS


=head1 AUTHOR     Mike Halagan <mhalagan@nmdp.org>
    
    Associate Bioinformatics Scientist
    3001 Broadway Stree NE
    Minneapolis, MN 55413
    ext. 8225

=head1 DESCRIPTION

    This script takes in the output of ngs-validate-interp and the observed file and generates
    a static HTML website report.

 =head1 LICENSE

    pipeline  Consensus assembly and allele interpretation pipeline.
    Copyright (c) 2014-2015 National Marrow Donor Program (NMDP)

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
use Test::More;

my $number_of_tests_run = 0;

my $working  = `pwd`;chomp($working);
my $test_expected = $working."/t/ex0_expected.txt";
my $test_observed = $working."/t/ex0_observed.txt";
my $test_verified = $working."/t/ex0_verified.txt";

print `perl ngs-validation-report -d t -f 1 -t 1`;
my @a_directories = split(/,/,"report,report/css,report/ex0,report/img,report/js");
my @a_files = split(/,/,"report/css/bootstrap.min.css,report/experiment.html,report/help.html,report/log.html,report/ex0/errors.html,report/ex0/fails.html,report/ex0/index.html,report/ex0/results.html,report/ex0/subjects/subject1.html");

foreach my $s_dir (@a_directories){
    if(-d $s_dir){
        is(1,1);$number_of_tests_run++;
    }else{
        is(1,0,"$s_dir not created!\n");$number_of_tests_run++;
    }
}
foreach my $s_file (@a_files){
    if(-e $s_file){
        is(1,1);$number_of_tests_run++;
    }else{
        is(1,0,"$s_file not created!\n");$number_of_tests_run++;
    }
}

done_testing( $number_of_tests_run );



















