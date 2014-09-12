#!/bin/bash  -x

echo "reading in file $1"

TOOLDIR=/mnt/scratch/caleb/ngs-tools-1.0-SNAPSHOT/bin
REFDIR=/mnt/common/data/reference/grch38/Primary_Assembly/assembled_chromosomes/FASTA
REFCHR=all_chr.fa

if [[ ! -f ${REFDIR}/${REFCHR} ]]; then
	echo "without ${REFDIR}/${REFCHR} I cannot proceed.  aborting"
	exit 1
fi
## TODO:  quit pulling things from caleb's $HOME and get ngs-tools installed for real.

for entry in `cat $1`; do
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

echo "Checking ${FILE1} and ${FILE2}"

   if [[ -f ${FILE2} ]];
   then
     echo "Found paired reads" 
     ${TOOLDIR}/ngs-fastq-to-ssake -1 ${FILE1} -2 ${FILE2} -o ${intermediate}.ssake.fa --insert-size 500
     sed -e 's/a/A/g' -e 's/c/C/g' -e 's/g/G/g' -e 's/t/T/g' -e 's/n/N/g' -e 's/^>.*$/>reformatted:500/' ${intermediate}.ssake.fa > $SSAKE_FA
   else
     cp ${FILE1} ${intermediate}.ssake.reformatted.fa
   fi

	SSAKE -f ${SSAKE_FA} -b ${intermediate}.ssake -w 1 -h 1 -p 1 -m 50 -o 30 -c 1 -e 0.90 -k 4 -a 0.1 -x 20

   makeFastaFileFromScaffolds.pl ${SCAFFOLDFILE}

   echo "REFERENCE = ${REFDIR}/${REFCHR}"

BOILERPLATE_HEADER=/tmp/boilerplate_header.bcl
echo > $BOILERPLATE_HEADER << EOF
##fileformat=VCFv4.2
##contig=<ID=gi|568336018|gb|CM000668.2|,length=171115067,assembly=hg19,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens">
#CHROM POS ID REF ALT QUAL FILTER INFO
EOF

   bwa bwasw ${REFDIR}/${REFCHR} ${intermediate}.ssake.scaffolds.fa | samtools view -Sb - | samtools sort - ${final}.contigs.bwa.sorted
   bwa mem ${REFDIR}/${REFCHR} ${FILE1} ${FILE2} | samtools view -Sb - | samtools sort - ${final}.reads.bwa.sorted
   samtools index ${final}.reads.bwa.sorted
   samtools mpileup -RB -C 0 -Q 0 -f ${REFDIR}/${REFCHR} ${final}.contigs.bwa.sorted.bam | cat $BOILERPLATE_HEADER - | bcftools view -O z -o ${final}.vcf.gz
   samtools mpileup -f ${REFDIR}/${REFCHR}.gz ${final}.reads.bwa.sorted.bam
done

echo "completed processing file $1"
