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

function error_abort() {
	echo "$1"
	exit 1
}


TARGETDIR="$1"
echo "this script will setup the data directory and populate it with some sample fastq files"
echo "for experimenting with the pipeline."
echo "you have selected ${TARGETDIR} as the target."
echo "is this correct? (y|n)"
read goahead

if [[ ${goahead}"x" != "yx" ]]; then
  error_abort "exiting on user request."
fi

if [[ ${UID} == "0" ]]; then
  error_abort "do not run this script as root.  aborting"
fi

if [[ "${TARGETDIR}x" == "x" ]]; then
  error_abort "no target directory specified, refusing to proceed"
fi

if [[ -d ${TARGETDIR} ]]; then
  echo "${TARGETDIR} already exists."
  echo "Do you want to delete all data in it and recreate?(y|n)"
  echo "WARNING:  THIS WILL DELETE YOUR EXISTING ${TARGETDIR}"
  read goahead

  if [[ ${goahead}"x" != "yx" ]]; then
    error_abort "exiting on user request."
  else
    rm -rf ${TARGETDIR} || error_abort "unable to delete ${TARGETDIR}"
  fi 
fi
  

mkdir -p ${TARGETDIR} || error_abort "could not create directory"
for subdir in raw intermediate final; do
  mkdir ${TARGETDIR}/$subdir || error_abort "unable to create ${TARGETDIR}/$subdir"
done

TUTORIAL_DIR=/opt/data/tutorial/fastq

if [[ ! -d $TUTORIAL_DIR ]]; then
  echo "unable to find $TUTORIAL_DIR -- no sample files located."
  echo "you will have to manually copy your fastq files into ${TARGETDIR}/raw"
else
  cp -a $TUTORIAL_DIR/*.fastq ${TARGETDIR}/raw
  echo "copied sample files from ${TUTORIAL_DIR} into ${TARGETDIR}/raw"
fi

echo "$0 exiting."
