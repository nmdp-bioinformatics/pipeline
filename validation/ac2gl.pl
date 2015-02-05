#!/usr/bin/env perl
=head1 NAME

ac2gl.pl

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 SVN
	
=head1 TODO
	
=head1 AUTHOR - Mike Halagan <mhalagan@nmdp.org>

=head1 VERSIONS
	
=head1 APENDIX

=cut
use strict;
use warnings;
use Data::Dumper;
use LWP::UserAgent;

my $s_infile = shift @ARGV or die "No infile provided!!";


open(IN,$s_infile) or die "CANT OPEN FILE $! $0";
while(<IN>){
	chomp;
	#/mnt/common/data/incoming/ucla/ex014/final/1367-7100-3_S26_L001_RX_001.fastq.contigs.bwa.sorted.bam	HLA-A	./tutorial/regions/grch38/hla-a/hla-a.ars.txt	2	HLA-A*01:01	HLA-A*02:ANGA
	my($s,$h,$path,$num,$s_typ1,$s_typ2) = split(/\t/,$_);
	$s_typ1 = isAC($s_typ1) ? toGl($s_typ1) : $s_typ1;
	$s_typ2 = isAC($s_typ2) ? toGl($s_typ2) : $s_typ2;

	$s_typ1 = $s_typ1 !~ /\S/ ? $s_typ2 : $s_typ1;
	$s_typ2 = $s_typ2 !~ /\S/ ? $s_typ1 : $s_typ2;

	print join("\t",$s,$h,$path,$num,$s_typ1,$s_typ2),"\n";
}


sub toGl{

	my $typ1 = shift;
	$typ1 =~ /HLA-(\D+\d{0,1})\*/;
	my $s_gl = isAC($typ1) ? ac2gl($1,$typ1) : $typ1;
	return $s_gl;

}

sub isAC{
	my $allele = shift;
	return 1 if($allele =~ /\*\d{2,3}:\D{2,5}/);
	return 0;
}

################################################################################################################
=head2 g2p

	Title:     g2p
	Usage:     my $s_allele = g2p($s_genomic_allele)
	Function:  This turns genomic level alleles into protein level alleles 
			   This will r
	Returns:   protein level alleles
	Args:      genomic level allele

=cut
sub g2p{
	my $allele    = shift;
	my @a_alleles = split(/:/,$allele);
	return join(":",$a_alleles[0],$a_alleles[1]);
}


sub ac2gl{

	my($loc,$typ) = @_;
	$typ =~ s/HLA-//;
	$typ =~ s/$loc\*//;

	my $ua = new LWP::UserAgent;
    $ua->agent("AlleleCodeClient/0.1");
    my @allele_list_list;my @ret;

    my $url = "http://emmes-dev.nmdp.org:8080/ac/api/decode?typing=$loc*$typ";
    my $response = $ua->request(new HTTP::Request("GET", $url));
    my $code = $response->code;
    my $content = $response->content;
    my $headers = $response->headers_as_string;
    if ($code == 200) {  # OK
        #print STDERR "Request url:  $url\nStatus $code Content: \n$content\n"; 
        my @allele_list_list = split ("/", $content);
        push @ret, @allele_list_list;
    } elsif ($code == 400) { # Bad Request
        # Request syntax was bad, or the typing was bad
        # print error and keep original typing.
        print STDERR "Bad request: $content\n\turl:  $url\n";
        push @ret, ("$loc*$typ");
    } else {
        die "System error: code=$code $content\n";
    }
    foreach(@ret){ $_ = "HLA-".$_;}
    return join("/",@ret);

}




