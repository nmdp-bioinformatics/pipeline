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

function bsd_linespersubfile {
  # this subroutine receives two arguments:  number of samples
  # and number of cores.
  # it then decides how to evenly divide the work into subfiles
  # checks to see if there are more subfiles than there are cores
  # and if there are, increases the line count in each subfile
  # and returns that value.

  # this is to avoid the case where we have more subfiles than cores

  let S=$1  #samples
  let C=$2  #cores
  let L=S/C  
  if [[ ${L} -eq 0 ]]; then
    let L=1  # we can't have less than one sample per file
  fi
  let E=S%C # see if we divided evenly
  if [[ ${E} -ne 0 ]]; then  # we did not divide evenly.
    let T=S/L
    if [[ ${T} -gt ${C} ]]; then  # check to see if files outnumber cores
      let L=L+1  # increment line-count-per-file by 1
    fi
  fi
  #finally, return how many lines should be in each file
  return $L
}




MYPID=$$
DEBUG=0  #change to 1 to increase verbosity and leave workfiles behind
MASTERWORKFILE=fullfile.$MYPID
SCRATCHFILE=scratch.$MYPID
TOOLPATH="$(pwd)"  #default to cwd if not specified
WORKDIR="$1"  #our data dir is argument one
DEBUG="$2"  #take debug mode as argument two 

#herein, we check to see if we're on osx or redhat or debian
# osx and debian require some changes as to how we split files
# into chunks, and how we detect number of CPUs
if [[ `uname` == "Darwin" ]]; then
  OS='osx'
else
  if [[ -f /etc/debian_version ]]; then
    OS='linuxDebian'
  elif [[ -f /etc/redhat-release ]]; then
    OS='linuxRedHat'
  else
    echo "this does not appear to be a version of Linux or OSX that I "
    echo "recognize.  Proceed?  (results may be unexpected)  (Y|N)"
    read plowahead
    if [[ $plowahead"x" != "Yx" ]]; then
      echo "exiting."
      exit 1
    fi
  fi 
fi

if [[ "${TOOLPATH}" == "$(pwd)" ]]; then 
  echo "TOOLPATH is the variable that must contain the path to your workerscript,"
  echo "which defaults to 'process_fastq.sh'.  This value can be altered by "
  echo "editing $0 and setting TOOLPATH correctly."
  echo "You have specified TOOLPATH as ${TOOLPATH}.  Is this REALLY correct? (y|n)"
  read goahead
  if [[ ${goahead}"x" != "yx" ]]; then
    echo "edit splitter.bash to have the desired TOOLPATH.  It should point to"
    echo "the directory containing process_fastq.bash"
    exit 1
  fi
fi
export goahead=""

if [[ ${WORKDIR}"x" == "x" ]]; then
  echo "location of work directory MUST be the first argument to this script."
  exit 1
fi

if [[ $DEBUG"x" != "x" ]]; then
  echo "debug mode enabled.  verbosity cranked."
  echo "WORKDIR = ${WORKDIR}"
fi

RAWDIR="${WORKDIR}/raw"
WORKERSCRIPT="${TOOLPATH}/process_fastq.bash"

if [[ ! -x ${WORKERSCRIPT} ]]; then
  echo "sorry, you specified a worker script (${WORKERSCRIPT}) which cannot be executed."
  exit 1
fi

INTERM="${WORKDIR}/intermediate"
FINAL="${WORKDIR}/final"
ABORT=0

##try to find all permission/directory existence issues
## in one pass, then spit out all the dir things that we found wrong.

for DIR in ${INTERM} ${FINAL} ${RAWDIR}; do
  if [[ ! -d ${DIR} ]]; then
    echo "directory ${DIR} must exist to proceed"
    let ABORT=ABORT+1
  fi
done

for DIR in ${INTERM} ${FINAL}; do
  if [[ ! -w ${DIR} ]]; then
    echo "directory ${DIR} must be writable (for you) to proceed"
    let ABORT=ABORT+1
  fi
  if [[ ! -x ${DIR} ]]; then
    echo "directory ${DIR} must be executable (for you) to proceed"
    let ABORT=ABORT+1
  fi
done

if [[ $ABORT -gt 0 ]]; then
  echo "exiting due to errors"
  exit 1
fi

# now, we accept fastq OR fq, but we want that as a suffix.
if [[ $(find ${RAWDIR} -type f | grep -cie '\(fq\|fastq\)$') -eq 0 ]]; then
  echo "no fastq files found in ${RAWDIR}.  Aborting."
  exit 1
fi

for DIR in ${INTERM} ${FINAL}; do
  if [[ $(ls ${DIR} | wc -l) -gt 0 ]]; then
    export goahead=""
    echo "Files already exist in ${DIR} -- enter 'y' to proceed and overwrite:"
    read goahead
    if [[ ${goahead}"x" != "yx" ]]; then
      exit 1
    fi
  fi
done

#  clean up before we start, in case we were interrupted last time
rm ${MASTERWORKFILE}  2>/dev/null
rm ${SCRATCHFILE}  2>/dev/null
if [[ ${DEBUG}"x" != "x" ]]; then
  echo "MASTERWORKFILE = ${MASTERWORKFILE}"
fi
for MYIDENTIFIER in $(find ${RAWDIR} \( -iname '*.fastq' -o -iname '*.fq' \) -print -type f | sed -e 's/R[12]/RX/g' | uniq); do if [[ ${DEBUG}"x" == "1x" ]]; then
    echo ${MYIDENTIFIER}
  fi
  echo ${MYIDENTIFIER} >> ${SCRATCHFILE}
  ## this should give us a unique list
  ## yes, this is a bit pathological, but duplicating
  ## entries could lead to two separate processing threads
  ## running on the same fastq file.  potentially ugly.
done

sort -u ${SCRATCHFILE} > ${MASTERWORKFILE}

# different unix variants have different ways of recording how many cores
# are in this cpu.
case $OS in 
  osx)
    let NUMCPU=$(sysctl -n machdep.cpu.core_count)  #OSX is a bit outdated.
  ;;
  *)
    #default is, so far, linux
    let NUMCPU=$(cat /proc/cpuinfo | grep -cw processor)
  ;;
esac

echo "preparing to divide work evenly among ${NUMCPU} cores"

echo "about to check $OS again"

## check number of cores.  if we're 1, we can skip ALL this nonsense
## likewise if we only have 1 sample.  TODO

SUBFILE_PATTERN=${MYPID}
## this will get overridden on osx, because osx hasn't been updated since '05.

case $OS in
  linuxDebian)
    #debian has an updated split command.  
    split ${MASTERWORKFILE} -n l/${NUMCPU} -d -e ${SUBFILE_PATTERN}
    ;;
  *)
    echo "ugly case"
    # this is for darwin, redhat, any other version where we're
    # not really sure that we get the -n l/ option (chunks)

    LINESINFILE=$(wc -l ${MASTERWORKFILE} | awk '{print $1}')
    ## above gets the length of our file, using the wc command and 
    # pulling just the linecount out.

    # honestly, the algorithm for computing lines per subfile got
    # complex enough that it's in the function at the top of the file
    # called "bsd_linespersubfile"

    bsd_linespersubfile ${LINESINFILE} ${NUMCPU}
    LINESPERSUBFILE=$?

    case $OS in
      linuxRedHat)
        split ${MASTERWORKFILE} -l ${LINESPERSUBFILE} -d ${SUBFILE_PATTERN}
        ;;
      osx)
        SUBFILE_PATTERN=x
        rm -f ${SUBFILE_PATTERN}[a-z][a-z]
        ## on OSX, we lose the protection of pid-named subfiles, so we're going to kill 
        ## any ^x files in an attempt to keep previous runs from munging up current runs
        split -l ${LINESPERSUBFILE} ${MASTERWORKFILE} 
        ;;
      default)
        split -l ${LINESPERSUBFILE} ${MASTERWORKFILE} 
        ;;
    esac
    ;;
esac

echo "about to call worker"

## now, time to call the worker script

if [ -f nohup.out ]; then
  rm nohup.out
fi

for file in $(find . -name "${SUBFILE_PATTERN}*" -type f -print); do
  BASEFILENAME=$(echo ${file} | awk -F\/ '{print $NF}')
  nohup time ${WORKERSCRIPT} ${file} ${DEBUG} >> timedata.${BASEFILENAME} &
done

if [[ ${DEBUG}"x" == "x" ]]; then
  #delete scratchfile and workfile unless we're in debug mode
  rm ${MASTERWORKFILE} ${SCRATCHFILE} 2>/dev/null
fi
