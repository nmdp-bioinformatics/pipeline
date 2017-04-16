#!/bin/bash
#
# pipeline Consensus assembly and allele interpretation pipeline.
# Copyright (c) 2014-2015 National Marrow Donor Program (NMDP)
#
# This library is free software; you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation; either version 3 of the License, or (at
# your option) any later version.
#
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; with out even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
# License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library; if not, write to the Free Software Foundation,
# Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
#
# > http://www.gnu.org/licenses/lgpl.html
#

ABORT=0
TOOLDIR=/opt/ngs-tools/bin
FINAL=$1
TEMPFILE=temp.fasta
BREADTH_OF_COVERAGE=0.5
OUTFILE=$2
ZYGOSITY=1000
LOCI="HLA-A HLA-B HLA-C HLA-DRB1 HLA-DQB1"



### validate that we have all the tools we need
if [[ ! -x ${TOOLDIR}/ngs-filter-consensus ]]; then
  echo "did not find ${TOOLDIR}/ngs-filter-consensus -- you must install NMDP ngs-tools."
  let ABORT=ABORT+1
fi

if [[ ${ABORT} -gt 0 ]]; then
  echo "exiting due to failures"
  exit 1
fi

rm -f ${OUTFILE}
for BAMFILE in $(find $FINAL -name "*.contigs.bwa.sorted.bam" -type f -print); do
  	for LOCUS in ${LOCI}; do
        REGIONS_FILE=
        if [ "$LOCUS" == "HLA-A" ]; then
          REGIONS_FILE="./tutorial/regions/grch38/hla-a/hla-a.ars.txt"
        fi
        if [ "$LOCUS" == "HLA-B" ]; then
          REGIONS_FILE="./tutorial/regions/grch38/hla-b/hla-b.ars.txt"
        fi
        if [ "$LOCUS" == "HLA-C" ]; then
          REGIONS_FILE="./tutorial/regions/grch38/hla-c/hla-c.ars.txt"
        fi
        if [ "$LOCUS" == "HLA-DRB1" ]; then
          REGIONS_FILE="./tutorial/regions/grch38/hla-drb1/hla-drb1.ars.txt"
        fi
        if [ "$LOCUS" == "HLA-DQB1" ]; then
          REGIONS_FILE="./tutorial/regions/grch38/hla-dqb1/hla-dqb1.ars.txt"
        fi
        rm -f ${TEMPFILE}
        ${TOOLDIR}/ngs-filter-consensus -i ${BAMFILE} -x ${REGIONS_FILE} -g ${LOCUS} -b ${BREADTH_OF_COVERAGE} -c -r -p ${ZYGOSITY} -o ${TEMPFILE}
        for INTERPRETATION in $(curl -T ${TEMPFILE} http://interp.b12x.org/hla/api/rs/seqInterp); do
            printf '%s\t%s\n' ${BAMFILE} ${INTERPRETATION} >> ${OUTFILE}
        done
    done
done
