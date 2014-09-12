#!/bin/bash

MYPID=$$
DEBUG=0  #change to 1 to increase verbosity and leave workfiles behind
MASTERWORKFILE=fullfile.$MYPID
SCRATCHFILE=scratch.$MYPID
TOOLPATH="/mnt/scratch/janderson/parallel_genomic/"
WORKDIR="$1"
DEBUG="$2"  #take debug mode as argument one on command line

if [[ ${TOOLPATH}"x" == "x" ]]; then 
	echo "you need to reset TOOLPATH and then comment out the lines which have 'commentme' on them in ${0}"  #commentme
	exit #commentme
fi

if [[ $DEBUG"x" != "x" ]]; then
	echo "WORKDIR = ${WORKDIR}"
fi

RAWDIR="${WORKDIR}/raw"
WORKERSCRIPT="${TOOLPATH}/process_fastq.bash"

INTERM="${WORKDIR}/intermediate"

FINAL="${WORKDIR}/final"

for DIR in ${INTERM} ${FINAL}; do
	if [[ ! -d ${DIR} ]]; then
		echo "${DIR} must exist for this script to proceed"
		exit 1
	else
		if [[ ! -w ${DIR} ]]; then
			echo "${DIR} must be writable for this script to proceed"
			exit 1
		else
			if [[ ! -x ${DIR} ]]; then
				echo "${DIR} must be executable for this script to proceed"
				exit 1
			fi
		fi
	fi

done

if [[ ${WORKDIR}"x" == "x" ]]; then
	echo "location of work directory MUST be the first argument to this script."
	exit 1
fi

if [[ ! -d ${RAWDIR} ]]; then
	echo "Raw fastq files must exist in ${RAWDIR} for this script to proceed."
	exit 1
else
	if [[ $(ls $RAWDIR | grep -c fastq) -eq 0 ]]; then
		echo "no fastq files found in $RAWDIR.  Aborting."
		exit 1
	fi
fi

for DIR in ${INTERM} ${FINAL}; do
	if [[ $(ls ${DIR} | wc -l) -gt 0 ]]; then
		echo "Files already exist in ${DIR} -- enter 'y' to proceed and overwrite"
		read goahead
		if [[ ${goahead}"x" != "yx" ]]; then
			exit 1
		fi
	fi
done	



#  clean up before we start, in case we were interrupted last time
rm $MASTERWORKFILE  2>/dev/null
rm $SCRATCHFILE  2>/dev/null
if [[ $DEBUG"x" != "x" ]]; then
	echo "MASTERWORKFILE = ${MASTERWORKFILE}"
fi
for MYIDENTIFIER in $(find $RAWDIR -name '*.fastq' -print -type f | sed -e 's/R[12]/RX/g' | uniq); do if [[ $DEBUG"x" == "1x" ]]; then
		echo $MYIDENTIFIER
	fi
	echo $MYIDENTIFIER >> $SCRATCHFILE
	## this should give us a unique list
done

sort -u $SCRATCHFILE > $MASTERWORKFILE

let NUMCPU=$(cat /proc/cpuinfo | grep -cw processor)
echo "preparing to divide work evenly among $NUMCPU cores"


split $MASTERWORKFILE -n l/$NUMCPU -d -e $MYPID 

## now, time to call the script derived from caleb's script

### DO SOMETHING

if [ -f nohup.out ]; then
	rm nohup.out
fi
for file in $(find . -name "$MYPID*" -type f -print); do
	BASEFILENAME=$(echo $file | awk -F\/ '{print $NF}')
	nohup time ${WORKERSCRIPT} $file >> timedata.$BASEFILENAME &

done

if [[ $DEBUG"x" == "x" ]]; then
	#delete scratchfile and workfile unless we're in debug mode
	rm $MASTERWORKFILE $SCRATCHFILE 2>/dev/null
fi
