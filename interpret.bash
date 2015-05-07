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
SOURCEFILE=$1
TEMPFILE=temp.fasta
BREADTH_OF_COVERAGE=0.5
OUTFILE=$2

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

while read SAMPLE_FILE LOCUS REGIONS_FILE ZYGOSITY FIRST_ALLELE SECOND_ALLELE; do
  rm -f ${TEMPFILE}
  ${TOOLDIR}/ngs-filter-consensus -i ${SAMPLE_FILE} -x ${REGIONS_FILE} -g ${LOCUS} -b ${BREADTH_OF_COVERAGE} -c -r -p ${ZYGOSITY} -o ${TEMPFILE}
  for INTERPRETATION in $(curl -T ${TEMPFILE} http://interp.b12x.org/hla/api/rs/seqInterp); do
    printf '%s\t%s\n' ${SAMPLE_FILE} ${INTERPRETATION} >> ${OUTFILE}
  done
done < "${SOURCEFILE}"
