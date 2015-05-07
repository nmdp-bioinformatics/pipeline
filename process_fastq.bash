#!/bin/bash

#
#    pipeline  Consensus assembly and allele interpretation pipeline.
#    Copyright (c) 2014-2015 National Marrow Donor Program (NMDP)
#
#    This library is free software; you can redistribute it and/or modify it
#    under the terms of the GNU Lesser General Public License as published
#    by the Free Software Foundation; either version 3 of the License, or (at
#    your option) any later version.
#
#    This library is distributed in the hope that it will be useful, but WITHOUT
#    ANY WARRANTY; with out even the implied warranty of MERCHANTABILITY or
#    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
#    License for more details.
#
#    You should have received a copy of the GNU Lesser General Public License
#    along with this library;  if not, write to the Free Software Foundation,
#    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA.
#
#    > http://www.gnu.org/licenses/lgpl.html
#
SOURCEFILE=$1
DEBUG=$2

echo "reading in file ${SOURCEFILE}"
if [[ ${DEBUG}"x" != "x" ]]; then
  echo "debug mode activated.  verbosity cranked."
  set -x
fi

TOOLDIR=/opt/ngs-tools/bin
REFDIR=/opt/data/grch38/seqs_for_alignment_pipelines

## note that REFDIR may require customization to your circumstances, as
## well as TOOLDIR
REFCHR=GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
## change above to just chromosome6 (chr6.fa) if you have less than 4GB or
## memory.  Frankly, we've seen issues even on machines _with_ 4GB, so
## for safeties sake, make sure you've got more memory.
TMPDIR=~/tmp

ABORT=0

if [[ ! -f ${REFDIR}/${REFCHR}.gz ]]; then
	echo "without ${REFDIR}/${REFCHR} I cannot proceed."
	let ABORT=ABORT+1
fi

###  validate that we have all the tools we need
if [[ ! -x ${TOOLDIR}/ngs-fastq-to-ssake ]]; then
	echo "did not find ${TOOLDIR}/ngs-fastq-to-ssake -- you must install NMDP ngs-tools."
	let ABORT=ABORT+1
fi

if [[ $(which SSAKE) == "" ]]; then
	echo "could not find SSAKE.  You must install the SSAKE package."
	let ABORT=ABORT+1
fi

#COMMENTS FROM RENE WARREN
#1. ssake.rb/tasr.rb/hlaminer.rb will be merged in homebrew today (26 sep 2014)
#2. TASR, which is a targeted assembler derived from SSAKE, reads a database file (eg. HLA allele seqs), deconstruct into kmers that are used to interrogate NGS reads for assembly.
#    I ran the tutorial data with TASR, using only HLA-I CDS sequences and assembled 3 large contigs, matching the A,B and C predicted alleles each with 100% sequence identity
#3. HLAminer, in targeted assembly mode uses TASR and produces an ouptut with following highest score/confidence
# HLA-A
#        Prediction #1 - A*24
#               A*24:02P,34168.1141,4.02e-13,124.0
#
#HLA-B
#        Prediction #1 - B*40
#                B*40:01P,55826.7173,3.57e-50,494.5
#
#HLA-C
#        Prediction #1 - C*03
#                C*03:04P,56610.9655,8.60e-40,390.7
#
#
#4. Additional assembly tools might be tried: Konnector (URL TBD, soon in homebrew) is a DBG assembler that will assemble only fragments and produce ambiguity scores when discordant allele data assembled, which is a nice feature.  It is also scalable
#5. If interested in micro assembly of overlapping reads, try abyss-mergepairs (linuxbrew/homebrew), will will output a long pseudoread with corrected bases if there's ground for merging

if [[ $(which bwa) == "" ]]; then
	echo "could not find bwa -- you can install bwa via homebrew.  Get at least version 0.7.8."
	let ABORT=ABORT+1
fi

if [[ $(which samtools) == "" ]]; then
	echo "could not find samtools -- you can install samtools via homebrew.  Get at least version 0.1.19."
	let ABORT=ABORT+1
fi

if [[ $(which bcftools) == "" ]]; then
	echo "could not find bcftools -- you can install bcftools via homebrew.  Get at least version 0.2.0-rc6."
	let ABORT=ABORT+1
fi

if [[ ${ABORT} -gt 0 ]]; then
	echo "exiting due to failures"
	exit 1
fi

if [[ ! -d $TMPDIR ]]; then
	echo "did not find $TMPDIR.  Creating"
	mkdir -p $TMPDIR || exit 1
fi

for entry in `cat ${SOURCEFILE}`; do
  path="${entry%/raw/*}"
  base="${entry##*/}"

  raw=${entry}
  intermediate="${path}/intermediate/${base}"
  final="${path}/final/${base}"
  SSAKE_FA="${intermediate}.ssake.reformatted.fa"
  SCAFFOLDFILE="${intermediate}.ssake.scaffolds"
  echo "path = ${path}"
  echo "base = ${base}"

## note that here, we get a filename with /path/to/file/PATTERN_"RX"_REMAININGPATTERN
## so we will create file1, file2 using sed, and replace RX with 
## R1 and R2 appropriately.

  FILE1=$(echo ${entry} | sed -e 's/RX/R1/g')
  FILE2=$(echo ${entry} | sed -e 's/RX/R2/g')
  PAIRED=1

  echo "Checking ${FILE1} and ${FILE2}"

  if [[ -f ${FILE2} && ${FILE1} != ${FILE2} ]];
  then
     echo "Found paired reads" 
     ${TOOLDIR}/ngs-fastq-to-ssake -1 ${FILE1} -2 ${FILE2} -o ${intermediate}.ssake.fa --insert-size 500
     sed -e 's/a/A/g' -e 's/c/C/g' -e 's/g/G/g' -e 's/t/T/g' -e 's/n/N/g' -e 's/^>.*$/>reformatted:500/' ${intermediate}.ssake.fa > $SSAKE_FA
  else
    PAIRED=0
    cp ${FILE1} ${intermediate}.ssake.reformatted.fa
  fi

  SSAKE -f ${SSAKE_FA} -b ${intermediate}.ssake -w 1 -h 1 -p ${PAIRED} -m 50 -o 30 -c 1 -e 0.90 -k 4 -a 0.1 -x 20

  echo "REFERENCE = ${REFDIR}/${REFCHR}"

BOILERPLATE_HEADER=$TMPDIR/boilerplate_header.bcl
if [[ ! -f ${BOILERPLATE_HEADER} ]]; then
# only recreate headerfile if it doesn't exist
cat > $BOILERPLATE_HEADER << EOF
##fileformat=VCFv4.2
##contig=<ID=gi|568336018|gb|CM000668.2|,length=171115067,assembly=hg19,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens">
#CHROM POS ID REF ALT QUAL FILTER INFO
EOF
fi
  bwa bwasw -b 1 ${REFDIR}/${REFCHR} ${intermediate}.ssake.contigs | samtools view -hub - | samtools sort -l 0 -O bam -T ${final}.contigs.bwa.tmp -o ${final}.contigs.bwa.sorted.bam
  bwa mem ${REFDIR}/${REFCHR} ${FILE1} ${FILE2} | samtools view -hub - | samtools sort -l 0 -O bam -T ${final}.reads.bwa.tmp -o ${final}.reads.bwa.sorted.bam
  ## Having major issues with samtools v1.1 -- will swap samtools mpileup for freebayes (unfortunately not in homebrew right now)
  # samtools index ${final}.reads.bwa.sorted.bam
  #samtools mpileup -RB -C 0 -Q 0 -f ${REFDIR}/${REFCHR}.gz ${final}.contigs.bwa.sorted.bam -v -u | cat $BOILERPLATE_HEADER -
  #samtools mpileup -f ${REFDIR}/${REFCHR}.gz ${final}.reads.bwa.sorted.bam -v -u > ${final}.reads.bwa.sorted.vcf
done

echo "completed processing file ${SOURCEFILE}"
