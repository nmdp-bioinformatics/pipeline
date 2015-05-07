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
OUTPUT_DIR=${1:-.}

ABORT=0
TOOLDIR=/opt/ngs-tools/bin
REFDIR=/mnt/common/data/reference/ipd-imgt-hla

###  validate that we have all the tools we need
if [[ ! -x ${TOOLDIR}/ngs-generate-paired-end-reads ]]; then
	echo "did not find ${TOOLDIR}/ngs-generate-paired-end-reads -- you must install NMDP ngs-tools."
	let ABORT=ABORT+1
fi
if [[ ! -x ${TOOLDIR}/ngs-split-fasta ]]; then
	echo "did not find ${TOOLDIR}/ngs-split-fasta -- you must install the NMDP ngs-tools."
	let ABORT=ABORT+1
fi

if [[ ${ABORT} -gt 0 ]]; then
	echo "exiting due to failures"
	exit 1
fi

MEAN_COVERAGE=1500
MEAN_LENGTH=250
MEAN_INSERT_SIZE=550

echo "generating reads for HLA-A..."
mkdir ${OUTPUT_DIR}/HLA-A
${TOOLDIR}/ngs-split-fasta -i ${REFDIR}/A_gen.fasta -d ${OUTPUT_DIR}/HLA-A -x fa.gz
for i in $( ls ${OUTPUT_DIR}/HLA-A/*.fa.gz ); do ${TOOLDIR}/ngs-generate-paired-end-reads -r $i -1 ${i/.fa.gz/_R1.fq.gz} -2 ${i/.fa.gz/_R2.fq.gz} --mean-length ${MEAN_LENGTH} --length-variation 0 --mean-insert-size ${MEAN_INSERT_SIZE} --insert-size-variation 0 --mean-coverage ${MEAN_COVERAGE} --quality illumina --mutation identity --seed 42; done

echo "generating reads for HLA-B..."
mkdir ${OUTPUT_DIR}/HLA-B
${TOOLDIR}/ngs-split-fasta -i ${REFDIR}/B_gen.fasta -d ${OUTPUT_DIR}/HLA-B -x fa.gz
for i in $( ls ${OUTPUT_DIR}/HLA-B/*.fa.gz ); do ${TOOLDIR}/ngs-generate-paired-end-reads -r $i -1 ${i/.fa.gz/_R1.fq.gz} -2 ${i/.fa.gz/_R2.fq.gz} --mean-length ${MEAN_LENGTH} --length-variation 0 --mean-insert-size ${MEAN_INSERT_SIZE} --insert-size-variation 0 --mean-coverage ${MEAN_COVERAGE} --quality illumina --mutation identity --seed 42; done

echo "generating reads for HLA-C..."
mkdir ${OUTPUT_DIR}/HLA-C
${TOOLDIR}/ngs-split-fasta -i ${REFDIR}/C_gen.fasta -d ${OUTPUT_DIR}/HLA-C -x fa.gz
for i in $( ls ${OUTPUT_DIR}/HLA-C/*.fa.gz ); do ${TOOLDIR}/ngs-generate-paired-end-reads -r $i -1 ${i/.fa.gz/_R1.fq.gz} -2 ${i/.fa.gz/_R2.fq.gz} --mean-length ${MEAN_LENGTH} --length-variation 0 --mean-insert-size ${MEAN_INSERT_SIZE} --insert-size-variation 0 --mean-coverage ${MEAN_COVERAGE} --quality illumina --mutation identity --seed 42; done

echo "generating reads for HLA-DRB1..."
mkdir ${OUTPUT_DIR}/HLA-DRB1
${TOOLDIR}/ngs-split-fasta -i ${REFDIR}/DRB1_gen.fasta -d ${OUTPUT_DIR}/HLA-DRB1 -x fa.gz
for i in $( ls ${OUTPUT_DIR}/HLA-DRB1/*.fa.gz ); do ${TOOLDIR}/ngs-generate-paired-end-reads -r $i -1 ${i/.fa.gz/_R1.fq.gz} -2 ${i/.fa.gz/_R2.fq.gz} --mean-length ${MEAN_LENGTH} --length-variation 0 --mean-insert-size ${MEAN_INSERT_SIZE} --insert-size-variation 0 --mean-coverage ${MEAN_COVERAGE} --quality illumina --mutation identity --seed 42; done

echo "generating reads for HLA-DPB1..."
mkdir ${OUTPUT_DIR}/HLA-DPB1
${TOOLDIR}/ngs-split-fasta -i ${REFDIR}/DPB_gen.fasta -d ${OUTPUT_DIR}/HLA-DPB1 -x fa.gz
for i in $( ls ${OUTPUT_DIR}/HLA-DPB1/*.fa.gz ); do ${TOOLDIR}/ngs-generate-paired-end-reads -r $i -1 ${i/.fa.gz/_R1.fq.gz} -2 ${i/.fa.gz/_R2.fq.gz} --mean-length ${MEAN_LENGTH} --length-variation 0 --mean-insert-size ${MEAN_INSERT_SIZE} --insert-size-variation 0 --mean-coverage ${MEAN_COVERAGE} --quality illumina --mutation identity --seed 42; done

echo "generating reads for HLA-DQB1..."
mkdir ${OUTPUT_DIR}/HLA-DQB1
${TOOLDIR}/ngs-split-fasta -i ${REFDIR}/DQB_gen.fasta -d ${OUTPUT_DIR}/HLA-DQB1 -x fa.gz
for i in $( ls ${OUTPUT_DIR}/HLA-DQB1/*.fa.gz ); do ${TOOLDIR}/ngs-generate-paired-end-reads -r $i -1 ${i/.fa.gz/_R1.fq.gz} -2 ${i/.fa.gz/_R2.fq.gz} --mean-length ${MEAN_LENGTH} --length-variation 0 --mean-insert-size ${MEAN_INSERT_SIZE} --insert-size-variation 0 --mean-coverage ${MEAN_COVERAGE} --quality illumina --mutation identity --seed 42; done

echo "done"
