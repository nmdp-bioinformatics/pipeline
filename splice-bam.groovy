/*
    pipeline  Consensus assembly and allele interpretation pipeline.
    Copyright (c) 2014 National Marrow Donor Program (NMDP)

    This library is free software; you can redistribute it and/or modify it
    under the terms of the GNU Lesser General Public License as published
    by the Free Software Foundation; either version 3 of the License, or (at
    your option) any later version.

    This library is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; with out even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
    License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this library;  if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA.

    > http://www.gnu.org/licenses/lgpl.html

*/

import java.io.File

import com.google.common.base.Joiner

import htsjdk.samtools.SAMFileReader

import htsjdk.samtools.util.SequenceUtil

import org.apache.commons.cli.Option

import org.biojava.bio.seq.DNATools

import org.biojava.bio.seq.impl.SimpleSequence

import org.biojava.bio.symbol.Edit
import org.biojava.bio.symbol.SymbolList
import org.biojava.bio.symbol.IllegalSymbolException

import org.nmdp.ngs.feature.Locus
import org.nmdp.ngs.feature.Allele

import org.nmdp.ngs.feature.parser.FeatureParser

CliBuilder cli = new CliBuilder()

cli.with {
  h longOpt: 'help', 'display usage and help message'
  x longOpt: 'intput-genomic-range-file', args: 1, 'input file of genomic ranges, one-per-line of form contig:xxxx-yyyy'
  i longOpt: 'input-file', args: 1, 'input file for unfiltered, aligned consensus sequences (.bam)'
  o longOpt: 'output-file', args: 1, 'output file for filtered consensus sequences'
  g longOpt: 'gene', args: 1, 'gene (to put in the FASTA header)'
  m longOpt: 'merge', 'strictly merge overlapping contigs'
  c longOpt: 'cdna', 'output cdna from the same contig (phased consensus sequence) in FASTA format (for interpretation)'
  b longOpt: 'minimum-breadth-of-coverage', args: 1, 'filter contigs less than minimum'
  v longOpt: 'verbose', 'produce verbose output'
}

def options = cli.parse(args)

assert options

if (options.h) {
 cli.usage()
}

def minimumBreadth = 0

if (options.b) {
  minimumBreadth = options.b as double
}

exons = [:]

new File(options.x).each { line ->
  def fields = line.split("\t")
  if(!(fields.size() >= 2)) {
    println "input file ${options.x} is not properly formatted"
    System.exit(-1)
  }
  def index = fields[0]
  def coordinate = fields[1]
  def locus = FeatureParser.parseLocus(coordinate)
  exons[index] = Allele.builder()
                  .withContig(locus.getContig())
                  .withStart(locus.getMin())
                  .withEnd(locus.getMax() + 1)
                  .build();
}

List<Allele> spliced = []

def regions = [:]

new SAMFileReader(new File(options.i)).each { record ->
  def alignment = new Locus("chr6", record.getAlignmentStart(), record.getAlignmentEnd());  
  def edits = cigarToEditList(record)
  def sequence = DNATools.createDNA(new String(record.getReadBases()))

  edits.each {
    sequence.edit(it)
  }

  def name = record.getReadName()  
  def contig = Allele.builder()
                 .withContig("chr6")
                 .withStart(record.getAlignmentStart())
                 .withEnd(record.getAlignmentEnd())
                 .withSequence(sequence)
                 .build();

  exons.each { index, exon ->
    def chr = exon.getContig()
    def min = exon.getMin()
    def max = exon.getMax()
    def range = "${chr}:${min}-${max}"
      
    if (alignment.overlaps(exon)) {
      def intersection = alignment.intersection(exon)
      def xover = exon.doubleCrossover(contig)

      def clipped = xover.leftHardClip("-").rightHardClip("-")
      clipped.setName(new String(">${name}|gene=${options.g}|exon=${index}|location=${range}|${max - min}"))

      if (!regions.containsKey(index)) {
        regions[index] = []
      }

      regions[index].add(clipped)

      def offset = 0
      edits.each { lesion ->
        if (lesion.replacement.equals(SymbolList.EMPTY_LIST) && exon.contains(lesion.pos)) {
          exon.sequence.edit(new Edit(lesion.pos + offset, 0, lesion.replacement))
          offset += lesion.length
        }
      }
    }
  }
}

def contigs = [:]

regions.each { region, list ->
  List<Allele> merged = list

  if (options.m) {
    merged = mergeAlleles(list)
  }

  merged.each { allele ->
    def sequenceLength = allele.sequence.seqString().length()
    def fields = allele.getName().tokenize('|')
    def locusLength = fields[-1] as int

    if (sequenceLength/locusLength >= minimumBreadth) {
      if (!contigs.containsKey(fields[0])) {
        contigs[fields[0]] = [:]
      }

      if (!contigs[fields[0]].containsKey(region)) {
        contigs[fields[0]][region] = []
      }

      contigs[fields[0]][region].add(allele)

      if (!options.c) {
        println "${allele.getName()}|${sequenceLength}"
        println allele.sequence.seqString().toUpperCase()
      }
    }
  }
}

if (options.c) {
  contigs.each { contig, exons ->    
    def cdna = ""
    exons.each { index, list ->
      best = list.sort { it.sequence.seqString().length() } .first()
      
      cdna += best.sequence.seqString()
    }

    println "${contig}"
    
    if(!(options.g.equals("HLA-A") || options.g.equals("HLA-DPB1"))) {
      cdna = DNATools.reverseComplement(DNATools.createDNA(cdna)).seqString() 
    }

    println "${cdna.toUpperCase()}"
  }
}

def cigarToEditList(record) {
  def edits = []
  def cigar = record.getCigar()
  def position = 0

  cigar.getCigarElements().each { element ->
    position += element.getLength();

    if (element.getOperator().toString().equals("I")) {
      edits += new Edit(position - element.getLength() + 1, element.getLength(), SymbolList.EMPTY_LIST)
      position -= element.getLength();
    }

    if (element.getOperator().toString().equals("D")) {
      def replace = DNATools.createDNA(Joiner.on("").join(Collections.nCopies(element.getLength(), "N")));
      edits += new Edit(position + element.getLength() - 1, 0, replace)
    }
  }
  return edits
}

def List<Allele> mergeAlleles(List<Allele> alleles) {
  def copy = alleles
  List<Allele> merged = new ArrayList<Allele>()

  while (copy.size() > 1) {
    Allele first = copy[0]

    def before = first
    copy.remove(0)
    def remove = []

    (0..copy.size() - 1).each { i ->
      Allele merge = first.merge(copy[i], 0)
      if (!merge.isEmpty()) {
        first = merge
        remove.add(i)
      }
    }

    remove.each { i->
      copy[i] = null
    }
    copy = copy.findAll { it != null }
    merged.add(first)
  }

  while (copy.size()) {
    merged.add(copy[0])
    copy.remove(0)
  }
  return merged
}
