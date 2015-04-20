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
    Copyright (c) 2014 National Marrow Donor Program (NMDP)

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
use vars qw($machine $user $working $b_npm $b_treeBuilder);
BEGIN{
    $b_treeBuilder = 0;
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
    foreach(@INC){
        my @inc_files = glob "$_/*";
        foreach my $inc_file (@inc_files){
            $b_treeBuilder++ if $inc_file =~ /HTML/;
        }
    }
}

if($b_treeBuilder > 0){
    use HTML::TreeBuilder;
}else{
    print STDERR "The perl package HTML::TreeBuilder is not installed!!\n";
    print STDERR "No tests of the HTML content will be done!\n";
}

my $number_of_tests_run = 0;

my $test_expected    = $working."/t/ex0_expected.txt";
my $test_observed    = $working."/t/ex0_observed.txt";
my $test_verified    = $working."/t/ex0_verified.txt";
my $experiment_file  = $working."/report/experiment.html";
my $s_blast_file     = $working."/report/blast/ex0/Test1.A.0.html";
my @a_directories    = split(/,/,"report,report/css,report/ex0,report/img,report/js");

# Load the expected results
my %h_tests_report;
my $s_test_cfg = $working."/t/cfg/dom_results.cfg";
open(my $test_results,"<",$s_test_cfg) or die "CANT OPEN FILE $! $0";
while(<$test_results>){
    chomp;
    my($s_row,$s_report,$s_exp,$s_id,$s_locus,$loc_verdict,$allele_verdict,$s_expected) = split(/,/,$_);
    $h_tests_report{$s_exp}{$s_row}{$s_id}{$s_locus}{LOCUS}    = $loc_verdict;
    $h_tests_report{$s_exp}{$s_row}{$s_id}{$s_locus}{ALLELE}   = $allele_verdict;
    $h_tests_report{$s_exp}{$s_row}{$s_id}{$s_locus}{EXPECTED} = $s_expected;
}
close $test_results;

# Load the expected experiment results
my %h_tests_experiment;
my $s_test_experiment_cfg = $working."/t/cfg/dom_experiment.cfg";
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

# Load the expected experiment results
my %h_tests_qc;
my $s_test_qc_cfg = $working."/t/cfg/dom_qc.cfg";
open(my $test_qc,"<",$s_test_qc_cfg) or die "CANT OPEN FILE $! $0";
while(<$test_qc>){
    chomp;
    my($row,$exp,$s_id,$qc,$locus_verdict,$allele_verdict,$qc_hap,$s_e) = split(/,/,$_);

    $h_tests_qc{$exp}{$row}{$s_id}{QC}     = $qc;
    $h_tests_qc{$exp}{$row}{$s_id}{LOC}    = $locus_verdict;   
    $h_tests_qc{$exp}{$row}{$s_id}{ALLELE} = $allele_verdict;   
    $h_tests_qc{$exp}{$row}{$s_id}{QCHAP}  = $qc_hap;  
    $h_tests_qc{$exp}{$row}{$s_id}{EXP}    = $s_e;  

}
close $test_qc;



# Add the subject page back in ,report/ex0/subjects/subject3.html
print `./ngs-validation-report -d t/txt -f 1 -t 1`;
&testDOM() if($b_treeBuilder);&testJs() if $b_npm;

# Test to make sure all the directories exist
foreach my $s_dir (@a_directories){
    if(-d $s_dir){
        is(1,1);$number_of_tests_run++;
    }else{
        is(1,0,"$s_dir not created!\n");$number_of_tests_run++;
        die("Exiting do to critical failure!");
    }
}

my @a_files = split(/,/,"report/css/bootstrap.min.css,report/experiment.html,report/help.html,report/log.html,report/ex0/drbx.html,report/ex0/errors.html,report/ex0/fails.html,report/ex0/index.html,report/ex0/results.html");

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

print `./ngs-validation-report -d t/hml -f 1 -t 1`;
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

# Running blast test with txt input file
print `./ngs-validation-report -d t/txt -n t/blastn -f 1 -t 1`;
&testDOM() if($b_treeBuilder);&testJs() if $b_npm;
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
if( -e $s_blast_file){
    is(1,1);$number_of_tests_run++;
}else{
    is(1,0,"$s_blast_file BLAST file not created!\n");$number_of_tests_run++;
}

# Running blast test with hml input file
print `./ngs-validation-report -d t/hml -n t/blastn -f 1 -t 1`;
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
if( -e $s_blast_file){
    is(1,1);$number_of_tests_run++;
}else{
    is(1,0,"$s_blast_file BLAST file not created!\n");$number_of_tests_run++;
}
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
                if($_ =~ /\d(.+)\(/){
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
                is(0,1,"Doesnt match up! TOTAL $row! $exp ");$number_of_tests_run++;
            }
            if(defined $h_tests_experiment{$row}{$exp}{PASS} && $h_tests_experiment{$row}{$exp}{PASS} eq $n_pass){
                is(1,1);$number_of_tests_run++;
            }else{
                is(0,1,"Doesnt match up! PASS $row! $exp $n_pass ");$number_of_tests_run++;
            }
            if(defined $h_tests_experiment{$row}{$exp}{PERPAS} && $h_tests_experiment{$row}{$exp}{PERPAS} eq $per_pass){
                is(1,1);$number_of_tests_run++;
            }else{
                is(0,1,"Doesnt match up! PERPAS $row! $exp $per_pass");$number_of_tests_run++;
            }
            if(defined $h_tests_experiment{$row}{$exp}{FAIL} && $h_tests_experiment{$row}{$exp}{FAIL} eq $n_fail){
                is(1,1);$number_of_tests_run++;
            }else{
                is(0,1,"Doesnt match up! FAIL $row! $exp ");$number_of_tests_run++;
            }
            if(defined $h_tests_experiment{$row}{$exp}{PERFAIL} && $h_tests_experiment{$row}{$exp}{PERFAIL} eq $per_fail){
                is(1,1);$number_of_tests_run++;
            }else{
                is(0,1,"Doesnt match up! PERFAIL $row! $exp ");$number_of_tests_run++;
            }        
            if(defined $h_tests_experiment{$row}{$exp}{DROP} && $h_tests_experiment{$row}{$exp}{DROP} eq $n_drop){
                is(1,1);$number_of_tests_run++;
            }else{
                is(0,1,"Doesnt match up! DROP $row! $exp ");$number_of_tests_run++;
            }
            if(defined $h_tests_experiment{$row}{$exp}{PERDROP} && $h_tests_experiment{$row}{$exp}{PERDROP} eq $per_drop){
                is(1,1);$number_of_tests_run++;
            }else{
                is(0,1,"Doesnt match up! PERDROP $row! $exp ");$number_of_tests_run++;
            } 
            if(defined $h_tests_experiment{$row}{$exp}{DRBX} && $h_tests_experiment{$row}{$exp}{DRBX} eq $n_drbx){
                is(1,1);$number_of_tests_run++;
            }else{
                is(0,1,"Doesnt match up! DRBX $row! $exp $h_tests_experiment{$row}{$exp}{DRBX} ne $n_drbx");$number_of_tests_run++;
            }
            if(defined $h_tests_experiment{$row}{$exp}{PERDRBX} && $h_tests_experiment{$row}{$exp}{PERDRBX} eq $per_drbx){
                is(1,1);$number_of_tests_run++;
            }else{
                is(0,1,"Doesnt match up! PERDRBX $row! $exp ");$number_of_tests_run++;
            }  
        }

        $row++;
    }

    my @a_experiments    = ("ex0","ex1");
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
                if(defined $h_tests_report{"results"}{$row}{$s_id}{$s_loc}{LOCUS} && $h_tests_report{"results"}{$row}{$s_id}{$s_loc}{LOCUS} eq $s_loc_verdict){
                    is(1,1);$number_of_tests_run++;
                }else{
                    is(0,1,"Doesnt match up! LOCUS $row! $s_loc");$number_of_tests_run++;
                }

                if(defined $h_tests_report{"results"}{$row}{$s_id}{$s_loc}{ALLELE} && $h_tests_report{"results"}{$row}{$s_id}{$s_loc}{ALLELE} eq $s_allele_verdict){
                    is(1,1);$number_of_tests_run++;
                }else{
                    is(0,1,"Doesnt match up! ALLELE $row! $s_loc");$number_of_tests_run++;
                }

                if(defined $h_tests_report{"results"}{$row}{$s_id}{$s_loc}{EXPECTED} && $h_tests_report{"results"}{$row}{$s_id}{$s_loc}{EXPECTED} eq $s_expected){
                    is(1,1);$number_of_tests_run++;
                }else{
                    is(0,1,"Doesnt match up! EXPECTED $row! $s_loc");$number_of_tests_run++;
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
                if(defined $h_tests_report{"fail"}{$row}{$s_id}{$s_loc}{LOCUS} && $h_tests_report{"fail"}{$row}{$s_id}{$s_loc}{LOCUS} eq $s_loc_verdict){
                    is(1,1);$number_of_tests_run++;
                }else{
                    is(0,1,"Doesnt match up! LOCUS $row! $s_loc");$number_of_tests_run++;
                }

                if(defined $h_tests_report{"fail"}{$row}{$s_id}{$s_loc}{ALLELE} && $h_tests_report{"fail"}{$row}{$s_id}{$s_loc}{ALLELE} eq $s_allele_verdict){
                    is(1,1);$number_of_tests_run++;
                }else{
                    is(0,1,"Doesnt match up! ALLELE $row! $s_loc");$number_of_tests_run++;
                }

                if(defined $h_tests_report{"fail"}{$row}{$s_id}{$s_loc}{EXPECTED} && $h_tests_report{"fail"}{$row}{$s_id}{$s_loc}{EXPECTED} eq $s_expected){
                    is(1,1);$number_of_tests_run++;
                }else{
                    is(0,1,"Doesnt match up! EXPECTED $row! $s_loc");$number_of_tests_run++;
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
                if(defined $h_tests_report{"drbx"}{$row}{$s_id}{$s_loc}{LOCUS} && $h_tests_report{"drbx"}{$row}{$s_id}{$s_loc}{LOCUS} eq $s_loc_verdict){
                    is(1,1);$number_of_tests_run++;
                }else{
                    is(0,1,"Doesnt match up! LOCUS $row! $s_loc");$number_of_tests_run++;
                }

                if(defined $h_tests_report{"drbx"}{$row}{$s_id}{$s_loc}{ALLELE} && $h_tests_report{"drbx"}{$row}{$s_id}{$s_loc}{ALLELE} eq $s_allele_verdict){
                    is(1,1);$number_of_tests_run++;
                }else{
                    is(0,1,"Doesnt match up! ALLELE $row! $s_loc");$number_of_tests_run++;
                }

                if(defined $h_tests_report{"drbx"}{$row}{$s_id}{$s_loc}{EXPECTED} && $h_tests_report{"drbx"}{$row}{$s_id}{$s_loc}{EXPECTED} eq $s_expected){
                    is(1,1);$number_of_tests_run++;
                }else{
                    is(0,1,"Doesnt match up! EXPECTED $row! $s_loc");$number_of_tests_run++;
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
                if(defined $h_tests_report{"error"}{$row}{$s_id}{$s_loc}{LOCUS} && $h_tests_report{"error"}{$row}{$s_id}{$s_loc}{LOCUS} eq $s_loc_verdict){
                    is(1,1);$number_of_tests_run++;
                }else{
                    is(0,1,"Doesnt match up! LOCUS $row! $s_loc");$number_of_tests_run++;
                }

                if(defined $h_tests_report{"error"}{$row}{$s_id}{$s_loc}{ALLELE} && $h_tests_report{"error"}{$row}{$s_id}{$s_loc}{ALLELE} eq $s_allele_verdict){
                    is(1,1);$number_of_tests_run++;
                }else{
                    is(0,1,"Doesnt match up! ALLELE $row! $s_loc");$number_of_tests_run++;
                }

                if(defined $h_tests_report{"error"}{$row}{$s_id}{$s_loc}{EXPECTED} && $h_tests_report{"error"}{$row}{$s_id}{$s_loc}{EXPECTED} eq $s_expected){
                    is(1,1);$number_of_tests_run++;
                }else{
                    is(0,1,"Doesnt match up! EXPECTED $row! $s_loc");$number_of_tests_run++;
                }
            }
            $row++;
        }


        $row=0;
        my $tree_qc = HTML::TreeBuilder->new; 
        $tree_qc->parse_file($qc_results);
        foreach my $tr ($tree_qc->find('tr')) {
            my @td   = $tr->find('td');
            my @data = map($_->as_text, @td);
            my($x,$s_id,$qc,$locus_verdict,$allele_verdict,$qc_hap,$s_e)= @data;

            if($row != 0){      

                if(defined $h_tests_qc{$s_exp}{$row}{$s_id}{QC} && $h_tests_qc{$s_exp}{$row}{$s_id}{QC} eq $qc){
                    is(1,1);$number_of_tests_run++;
                }else{
                    is(0,1,"Doesnt match up! QC $row! ");$number_of_tests_run++;
                }
                if(defined $h_tests_qc{$s_exp}{$row}{$s_id}{LOC} && $h_tests_qc{$s_exp}{$row}{$s_id}{LOC} eq $locus_verdict){
                    is(1,1);$number_of_tests_run++;
                }else{
                    is(0,1,"Doesnt match up! LOC $row! ");$number_of_tests_run++;
                }
                if(defined $h_tests_qc{$s_exp}{$row}{$s_id}{ALLELE} && $h_tests_qc{$s_exp}{$row}{$s_id}{ALLELE} eq $allele_verdict){
                    is(1,1);$number_of_tests_run++;
                }else{
                    is(0,1,"Doesnt match up! ALLELE $row! ");$number_of_tests_run++;
                }
                if(defined $h_tests_qc{$s_exp}{$row}{$s_id}{QCHAP} && $h_tests_qc{$s_exp}{$row}{$s_id}{QCHAP} eq $qc_hap){
                    is(1,1);$number_of_tests_run++;
                }else{
                    is(0,1,"Doesnt match up! QCHAP $row! ");$number_of_tests_run++;
                } 
                if(defined $h_tests_qc{$s_exp}{$row}{$s_id}{EXP} && $h_tests_qc{$s_exp}{$row}{$s_id}{EXP} eq $s_e){
                    is(1,1);$number_of_tests_run++;
                }else{
                    is(0,1,"Doesnt match up! EXP $row! ");$number_of_tests_run++;
                }                              
            }
            $row++;
        }
        


    }

}


sub testJs{

    print STDERR `npm test &> npm.stderr`;
    my $js_test = `cat npm.stderr`;chomp($js_test);
    print STDERR `echo;tail -3 npm.stderr | head -1;rm npm.stderr`;
    if($js_test =~ /PASS/){
         is(1,1);$number_of_tests_run++;
    }else{
         is(1,0,"Javascript tests failed!! See npm-debug.log for details!");$number_of_tests_run++;
    }

}












