#!/usr/bin/env perl
=head1 NAME

    report.t

=head1 SYNOPSIS

    prove t/report.t

=head1 AUTHOR     Mike Halagan <mhalagan@nmdp.org>
    
    Bioinformatics Scientist
    3001 Broadway Stree NE
    Minneapolis, MN 55413
    ext. 8225

=head1 DESCRIPTION

    This is a test script tests that the expected html files are produced
    and the correct elements are within the html files. 

=head1 CAVEATS
    
    - In order for this to work you MUST have ngs-extract-expected-haploids
      installed on your machine.

=head1 OUTPUT

    bash-3.2$ make test
    t/report.t .. ok      
    All tests successful.
    Files=1, Tests=1844, 10 wallclock secs ( 0.20 usr  0.02 sys + 14.50 cusr  1.35 csys = 16.07 CPU)
    Result: PASS

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
use Test::More;
use Data::Dumper;
use HTML::TreeBuilder;
use vars qw($machine $user $working $b_npm $b_treeBuilder);
BEGIN{

    $working    = `pwd`;chomp($working);
    $user       = `whoami`;chomp($user);
    $machine    = `uname`;chomp($machine);
    my $npm     = `which npm`;chomp($npm);
    $b_npm      = defined $b_npm ? 1 : 0;
    if($working !~ /\/t/){
        my $lib = $working."/lib";
        push(@INC,$lib);
    }else{
        $working =~ s/\/t//;
        my $lib = $working."/lib";
        push(@INC,$lib);
    }

}

my $perl_v = $];             # Get the perl version
my $number_of_tests_run = 0; # Number of tests run

my $file_ending = $perl_v > 5.016002 ? "5.18" : "5.16";

my $experiment_file  = $working."/report/experiment.html";
my @a_directories    = split(/,/,"report,report/css,report/0cbdf011-79d7-4a36-897e-4fa4c17ec43c,report/img,report/js");

# Load the expected results
my %h_tests_report;
my $s_test_cfg = $working."/t/cfg/dom_results-".$file_ending."-uuid.cfg";
open(my $test_results,"<",$s_test_cfg) or die "CANT OPEN FILE $! $0";
while(<$test_results>){
    chomp;
    my($s_row,$s_report,$s_exp,$s_id,$s_locus,$loc_verdict,$allele_verdict,$s_expected) = split(/,/,$_);
    $h_tests_report{$s_exp}{$s_report}{$s_row}{$s_id}{$s_locus}{LOCUS}    = $loc_verdict;
    $h_tests_report{$s_exp}{$s_report}{$s_row}{$s_id}{$s_locus}{ALLELE}   = $allele_verdict;
    $h_tests_report{$s_exp}{$s_report}{$s_row}{$s_id}{$s_locus}{EXPECTED} = $s_expected;

}
close $test_results;

# Load the expected experiment results
my %h_tests_experiment;
my $s_test_experiment_cfg = $working."/t/cfg/dom_experiment-uuid.cfg";
open(my $test_experiment,"<",$s_test_experiment_cfg) or die "CANT OPEN FILE $! $0";
while(<$test_experiment>){
    chomp;
    my($s_row,$experiment,$exp,$num,$n_pass,$per_pass,$n_fail,$per_fail,$n_drop,$per_drop,$n_drbx,$per_drbx) = split(/,/,$_);
    
    #8,experiment,ex1,20,1,80.0%,2,10.0%,2,10.0%
    foreach($n_pass,$n_fail,$n_drop,$n_drbx){
        $_ = !defined $_ || $_ !~ /\S/ ? 0 : $_;
    }

    $h_tests_experiment{$s_row}{$exp}{TOTAL}   = $num;
    $h_tests_experiment{$s_row}{$exp}{PASS}    = $n_pass;
    $h_tests_experiment{$s_row}{$exp}{PERPAS}  = $per_pass;
    $h_tests_experiment{$s_row}{$exp}{FAIL}    = $n_fail;
    $h_tests_experiment{$s_row}{$exp}{PERFAIL} = $per_fail;
    $h_tests_experiment{$s_row}{$exp}{DROP}    = $n_drop; 
    $h_tests_experiment{$s_row}{$exp}{PERDROP} = $per_drop; 
    $h_tests_experiment{$s_row}{$exp}{DRBX}    = $n_drbx; 
    $h_tests_experiment{$s_row}{$exp}{PERDRBX} = $per_drbx;         

}
close $test_experiment;



my @a_files = split(/,/,"report/experiment.html,report/help.html,report/log.html,report/0cbdf011-79d7-4a36-897e-4fa4c17ec43c/drbx.html,report/0cbdf011-79d7-4a36-897e-4fa4c17ec43c/errors.html,report/0cbdf011-79d7-4a36-897e-4fa4c17ec43c/fails.html,report/0cbdf011-79d7-4a36-897e-4fa4c17ec43c/index.html,report/0cbdf011-79d7-4a36-897e-4fa4c17ec43c/results.html");

# Testing that it works with uuids
print `./ngs-validation-report -d t/uuid -f -t 1`;

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

&testDOM();
#print `rm -R report`;


done_testing( $number_of_tests_run );

sub testDOM{

    my $row=0;
    my $tree_experiment = HTML::TreeBuilder->new;
    $tree_experiment->parse_file($experiment_file);
    foreach my $tr ($tree_experiment->find('tr')) {
        my @td   = $tr->find('td');
        my @data = map($_->as_text, @td);

        if($row != 0 && $row != 3 && $row != 6){
            foreach(@data){
                if($_ =~ /\d{1,2}(.+)\(/){
                    $_ =~ s/$1//;
                }
                $_ =~ s/\&//g;$_ =~ s/\&nbsp;\&nbsp;//g;
                $_ =~ s/\s//g;$_ =~ s/ //g;
            }
            my($exp,$num,$passed,$failed,$dropout,$drbx) = @data;
            my($n_pass,$per_pass) = split(/\(/,$passed);$per_pass =~ s/\)//;
            my($n_fail,$per_fail) = split(/\(/,$failed);$per_fail =~ s/\)//;
            my($n_drop,$per_drop) = split(/\(/,$dropout);$per_drop =~ s/\)//;
            my($n_drbx,$per_drbx) = split(/\(/,$drbx);$per_drbx =~ s/\)//;
            
            foreach($n_pass,$n_fail,$n_drop,$n_drbx){
                $_ = !defined $_ || $_ !~ /\S/ ? 0 : $_;
            }

            if(defined $h_tests_experiment{$row}{$exp}{TOTAL} && $h_tests_experiment{$row}{$exp}{TOTAL} eq $num){
                is(1,1);$number_of_tests_run++;
            }else{
                is(0,1,"experiment.html - Doesnt match up! TOTAL $row! $exp | $h_tests_experiment{$row}{$exp}{TOTAL} eq $num ");$number_of_tests_run++;
            }
            if(defined $h_tests_experiment{$row}{$exp}{PASS} && $h_tests_experiment{$row}{$exp}{PASS} eq $n_pass){
                is(1,1);$number_of_tests_run++;
            }else{
                is(0,1,"experiment.html - Doesnt match up! PASS $row! $exp $n_pass | $h_tests_experiment{$row}{$exp}{PASS} eq $n_pass");$number_of_tests_run++;
            }
            if(defined $h_tests_experiment{$row}{$exp}{PERPAS} && $h_tests_experiment{$row}{$exp}{PERPAS} eq $per_pass){
                is(1,1);$number_of_tests_run++;
            }else{
                is(0,1,"experiment.html - Doesnt match up! PERPAS $row! $exp $per_pass | $h_tests_experiment{$row}{$exp}{PERPAS} eq $per_pass");$number_of_tests_run++;
            }
            if(defined $h_tests_experiment{$row}{$exp}{FAIL} && $h_tests_experiment{$row}{$exp}{FAIL} eq $n_fail){
                is(1,1);$number_of_tests_run++;
            }else{
                is(0,1,"experiment.html - Doesnt match up! FAIL $row! $exp | $h_tests_experiment{$row}{$exp}{FAIL} eq $n_fail");$number_of_tests_run++;
            }
            if(defined $h_tests_experiment{$row}{$exp}{PERFAIL} && $h_tests_experiment{$row}{$exp}{PERFAIL} eq $per_fail){
                is(1,1);$number_of_tests_run++;
            }else{
                is(0,1,"experiment.html - Doesnt match up! PERFAIL $row! $exp | $h_tests_experiment{$row}{$exp}{PERFAIL} eq $per_fail");$number_of_tests_run++;
            }        
            if(defined $h_tests_experiment{$row}{$exp}{DROP} && $h_tests_experiment{$row}{$exp}{DROP} eq $n_drop){
                is(1,1);$number_of_tests_run++;
            }else{
                is(0,1,"experiment.html - Doesnt match up! DROP $row! $exp | $h_tests_experiment{$row}{$exp}{DROP} eq $n_drop");$number_of_tests_run++;
            }
            if(defined $h_tests_experiment{$row}{$exp}{PERDROP} && $h_tests_experiment{$row}{$exp}{PERDROP} eq $per_drop){
                is(1,1);$number_of_tests_run++;
            }else{
                is(0,1,"experiment.html - Doesnt match up! PERDROP $row! $exp | $h_tests_experiment{$row}{$exp}{PERDROP} eq $per_drop");$number_of_tests_run++;
            } 
            if(defined $h_tests_experiment{$row}{$exp}{DRBX} && $h_tests_experiment{$row}{$exp}{DRBX} eq $n_drbx){
                is(1,1);$number_of_tests_run++;
            }else{
                is(0,1,"experiment.html - Doesnt match up! DRBX $row! $exp | $h_tests_experiment{$row}{$exp}{DRBX} ne $n_drbx");$number_of_tests_run++;
            }
            if(defined $h_tests_experiment{$row}{$exp}{PERDRBX} && $h_tests_experiment{$row}{$exp}{PERDRBX} eq $per_drbx){
                is(1,1);$number_of_tests_run++;
            }else{
                is(0,1,"experiment.html - Doesnt match up! PERDRBX $row! $exp | $h_tests_experiment{$row}{$exp}{PERDRBX} eq $per_drbx");$number_of_tests_run++;
            }  
        }

        $row++;
    }

    my @a_experiments    = ("0cbdf011-79d7-4a36-897e-4fa4c17ec43c","676cce41-1511-48e0-aa40-4a480d0a7198");
    foreach my $s_exp (@a_experiments){

        my $html_results     = $working."/report/".$s_exp."/results.html";
        my $drbx_results     = $working."/report/".$s_exp."/drbx.html";
        my $dropout_results  = $working."/report/".$s_exp."/errors.html";
        my $failed_results   = $working."/report/".$s_exp."/fails.html";
        my $qc_results       = $working."/report/".$s_exp."/qc.html";

        my $row=0;
        my $tree_results = HTML::TreeBuilder->new; 
        $tree_results->parse_file($html_results);
        foreach my $tr ($tree_results->find('tr')) {
            my @td   = $tr->find('td');
            my @data = map($_->as_text, @td);
            my($x,$s_id,$s_loc,$s_loc_verdict,$s_allele_verdict,$n_observed,$s_expected,$s_observed1,$s_observed2) = @data;

            if($row != 0){      
                if(defined $h_tests_report{"results"}{$s_exp}{$row}{$s_id}{$s_loc}{LOCUS} && $h_tests_report{"results"}{$s_exp}{$row}{$s_id}{$s_loc}{LOCUS} eq $s_loc_verdict){
                    is(1,1);$number_of_tests_run++;
                }else{
                    is(0,1,"$html_results - Doesnt match up! $s_exp | results $row $s_id $s_loc LOCUS | $h_tests_report{\"results\"}{$s_exp}{$row}{$s_id}{$s_loc}{LOCUS} ne $s_loc_verdict");$number_of_tests_run++;
                }

                if(defined $h_tests_report{"results"}{$s_exp}{$row}{$s_id}{$s_loc}{ALLELE} && $h_tests_report{"results"}{$s_exp}{$row}{$s_id}{$s_loc}{ALLELE} eq $s_allele_verdict){
                    is(1,1);$number_of_tests_run++;
                }else{
                    is(0,1,"$html_results - Doesnt match up! $s_exp | results $row $s_id $s_loc ALLELE | $h_tests_report{\"results\"}{$s_exp}{$row}{$s_id}{$s_loc}{ALLELE} ne $s_allele_verdict");$number_of_tests_run++;
                }

                if(defined $h_tests_report{"results"}{$s_exp}{$row}{$s_id}{$s_loc}{EXPECTED} && $h_tests_report{"results"}{$s_exp}{$row}{$s_id}{$s_loc}{EXPECTED} eq $s_expected){
                    is(1,1);$number_of_tests_run++;
                }else{
                    is(0,1,"$html_results - Doesnt match up! $s_exp | results $row $s_id $s_loc EXPECTED |$h_tests_report{\"results\"}{$s_exp}{$row}{$s_id}{$s_loc}{EXPECTED} ne $s_expected");$number_of_tests_run++;
                }
            }
            $row++;
        }

        $row=0;
        my $tree_failed = HTML::TreeBuilder->new; 
        $tree_failed->parse_file($failed_results);
        foreach my $tr ($tree_failed->find('tr')) {
            my @td   = $tr->find('td');
            my @data = map($_->as_text, @td);
            my($x,$s_id,$s_loc,$s_loc_verdict,$s_allele_verdict,$n_observed,$s_expected,$s_observed1,$s_observed2) = @data;

            if($row != 0){      
                if(defined $h_tests_report{"fail"}{$s_exp}{$row}{$s_id}{$s_loc}{LOCUS} && $h_tests_report{"fail"}{$s_exp}{$row}{$s_id}{$s_loc}{LOCUS} eq $s_loc_verdict){
                    is(1,1);$number_of_tests_run++;
                }else{
                    is(0,1,"$failed_results - Doesnt match up! LOCUS $s_exp $s_id $row! $s_loc | $h_tests_report{\"fail\"}{$s_exp}{$row}{$s_id}{$s_loc}{LOCUS} ne $s_loc_verdict");$number_of_tests_run++;
                }

                if(defined $h_tests_report{"fail"}{$s_exp}{$row}{$s_id}{$s_loc}{ALLELE} && $h_tests_report{"fail"}{$s_exp}{$row}{$s_id}{$s_loc}{ALLELE} eq $s_allele_verdict){
                    is(1,1);$number_of_tests_run++;
                }else{
                    is(0,1,"$failed_results - Doesnt match up! ALLELE $s_exp $s_id $row! $s_loc | $h_tests_report{\"fail\"}{$s_exp}{$row}{$s_id}{$s_loc}{ALLELE} ne $s_allele_verdict");$number_of_tests_run++;
                }

                if(defined $h_tests_report{"fail"}{$s_exp}{$row}{$s_id}{$s_loc}{EXPECTED} && $h_tests_report{"fail"}{$s_exp}{$row}{$s_id}{$s_loc}{EXPECTED} eq $s_expected){
                    is(1,1);$number_of_tests_run++;
                }else{
                    is(0,1,"$failed_results - Doesnt match up! EXPECTED $s_exp $s_id $row! $s_loc | $h_tests_report{\"fail\"}{$s_exp}{$row}{$s_id}{$s_loc}{EXPECTED} ne $s_expected");$number_of_tests_run++;
                }
            }
            $row++;
        }

        $row=0;
        my $tree_drbx = HTML::TreeBuilder->new; 
        $tree_drbx->parse_file($drbx_results);
        foreach my $tr ($tree_drbx->find('tr')) {
            my @td   = $tr->find('td');
            my @data = map($_->as_text, @td);
            my($x,$s_id,$s_loc,$s_loc_verdict,$s_allele_verdict,$n_observed,$s_expected,$s_observed1,$s_observed2) = @data;

            if($row != 0){      
                if(defined $h_tests_report{"drbx"}{$s_exp}{$row}{$s_id}{$s_loc}{LOCUS} && $h_tests_report{"drbx"}{$s_exp}{$row}{$s_id}{$s_loc}{LOCUS} eq $s_loc_verdict){
                    is(1,1);$number_of_tests_run++;
                }else{
                    is(0,1,"$drbx_results - Doesnt match up! LOCUS $s_exp $s_id $row! $s_loc | $h_tests_report{\"drbx\"}{$s_exp}{$row}{$s_id}{$s_loc}{LOCUS} ne $s_loc_verdict");$number_of_tests_run++;
                }

                if(defined $h_tests_report{"drbx"}{$s_exp}{$row}{$s_id}{$s_loc}{ALLELE} && $h_tests_report{"drbx"}{$s_exp}{$row}{$s_id}{$s_loc}{ALLELE} eq $s_allele_verdict){
                    is(1,1);$number_of_tests_run++;
                }else{
                    is(0,1,"$drbx_results - Doesnt match up! ALLELE $s_exp $s_id $row! $s_loc | $h_tests_report{\"drbx\"}{$s_exp}{$row}{$s_id}{$s_loc}{ALLELE} ne $s_allele_verdict");$number_of_tests_run++;
                }

                if(defined $h_tests_report{"drbx"}{$s_exp}{$row}{$s_id}{$s_loc}{EXPECTED} && $h_tests_report{"drbx"}{$s_exp}{$row}{$s_id}{$s_loc}{EXPECTED} eq $s_expected){
                    is(1,1);$number_of_tests_run++;
                }else{
                    is(0,1,"$drbx_results - Doesnt match up! EXPECTED $s_exp $s_id $row! $s_loc | $h_tests_report{\"drbx\"}{$s_exp}{$row}{$s_id}{$s_loc}{EXPECTED} ne $s_expected");$number_of_tests_run++;
                }
            }
            $row++;
        }

        $row=0;
        my $tree_drop = HTML::TreeBuilder->new; 
        $tree_drop->parse_file($dropout_results);
        foreach my $tr ($tree_drop->find('tr')) {
            my @td   = $tr->find('td');
            my @data = map($_->as_text, @td);
            my($x,$s_id,$s_loc,$s_loc_verdict,$s_allele_verdict,$n_observed,$s_expected,$s_observed1,$s_observed2) = @data;

            if($row != 0){      
                if(defined $h_tests_report{"error"}{$s_exp}{$row}{$s_id}{$s_loc}{LOCUS} && $h_tests_report{"error"}{$s_exp}{$row}{$s_id}{$s_loc}{LOCUS} eq $s_loc_verdict){
                    is(1,1);$number_of_tests_run++;
                }else{
                    is(0,1,"$dropout_results - Doesnt match up! LOCUS $s_exp $s_id $row! $s_loc | $h_tests_report{\"error\"}{$s_exp}{$row}{$s_id}{$s_loc}{LOCUS} ne $s_loc_verdict");$number_of_tests_run++;
                }

                if(defined $h_tests_report{"error"}{$s_exp}{$row}{$s_id}{$s_loc}{ALLELE} && $h_tests_report{"error"}{$s_exp}{$row}{$s_id}{$s_loc}{ALLELE} eq $s_allele_verdict){
                    is(1,1);$number_of_tests_run++;
                }else{
                    is(0,1,"$dropout_results - Doesnt match up! ALLELE $s_exp $s_id $row! $s_loc | $h_tests_report{\"error\"}{$s_exp}{$row}{$s_id}{$s_loc}{ALLELE} ne $s_allele_verdict");$number_of_tests_run++;
                }

                if(defined $h_tests_report{"error"}{$s_exp}{$row}{$s_id}{$s_loc}{EXPECTED} && $h_tests_report{"error"}{$s_exp}{$row}{$s_id}{$s_loc}{EXPECTED} eq $s_expected){
                    is(1,1);$number_of_tests_run++;
                }else{
                    is(0,1,"$dropout_results - Doesnt match up! EXPECTED $s_exp $s_id $row! $s_loc | $h_tests_report{\"error\"}{$s_exp}{$row}{$s_id}{$s_loc}{EXPECTED} ne $s_expected");$number_of_tests_run++;
                }
            }
            $row++;
        }

        


    }

}












