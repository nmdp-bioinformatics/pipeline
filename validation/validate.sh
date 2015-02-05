#!/bin/bash

#
#    pipeline  Consensus assembly and allele interpretation pipeline.
#    Copyright (c) 2014 National Marrow Donor Program (NMDP)
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

TOOLDIR=/opt/ngs-tools/bin
VALTOOL=/home/mhalagan/nmdp-bioinformatics/validation/pipeline/validation
ABORT=0
EXPECTED=$1
OBSERVED=$2
EXP=$3
EMAIL=$4

#reformat expected
perl ${VALTOOL}/ac2gl.pl ${EXPECTED} > ${EXP}_expected.txt
ngs-validate-interpretation -e ${EXP}_expected.txt -b ${OBSERVED} -r 2 > ${EXP}_validated.txt
perl ${VALTOOL}/ngs-validation-report -e ${EXP}_expected.txt -l ${EXP}_validated.txt -o ${OBSERVED} -v 1 -m ${EMAIL}
