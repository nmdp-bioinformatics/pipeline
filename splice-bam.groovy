import org.apache.commons.cli.Option

import java.io.File
import net.sf.samtools.SAMFileReader
import net.sf.samtools.util.SequenceUtil
import org.nmdp.ngs.feature.Locus
import org.nmdp.ngs.feature.Allele
import org.nmdp.ngs.feature.parser.FeatureParser
import org.biojava.bio.symbol.Edit
import org.biojava.bio.seq.DNATools
import org.biojava.bio.symbol.SymbolList
import org.biojava.bio.symbol.IllegalSymbolException
import org.biojava.bio.seq.impl.SimpleSequence
import com.google.common.base.Joiner
import static org.nmdp.ngs.feature.Allele.builder;

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

if(options.h) { cli.usage() }

def minimumBreadth = 0

if(options.b) {
  minimumBreadth = options.b as double
}

exons = [:]

new File(options.x).each { line ->
  def (index, coordinate) = line.split("\t")
  def locus = FeatureParser.parseLocus(coordinate)
  exons[index] = builder()
                .withContig(locus.getContig())
                .withMin(locus.getMin())
                .withMax(locus.getMax() + 1)
                .build();
}

ArrayList<Allele> spliced = []

def regions = [:]

new SAMFileReader(new File(options.i)).each { record ->
  def alignment = new Locus("chr6", record.getAlignmentStart(), record.getAlignmentEnd());
  // println "alignment = ${alignment}"
  

  
  def edits = CigarToEditList(record)
  def sequence = DNATools.createDNA(new String(record.getReadBases()))
  
  // println "${sequence.seqString()}"
  // println "EDITS:"
  
  edits.each {
    sequence.edit(it)
  }
  
  def name = record.getReadName()
  
  def contig = builder()
                        .withContig("chr6")
                        .withMin(record.getAlignmentStart())
                        .withMax(record.getAlignmentEnd())
                        .withSequence(sequence)
                        .build();
  
  // println "${record.getReadName()} ${record.getReferenceName()} ${record.getAlignmentStart()} ${record.getCigarString()} ${record.getAlignmentEnd()} ${new String(record.getReadBases())}"
  
  // maybe a bug in makeReferenceFromAlignment(record true) -- don't see 'N' for deleted reference bases!
  // println "${new String(SequenceUtil.makeReferenceFromAlignment(record, true))}"
  
  // println "EXONS:"
  
    
    def cdna = ""
    
    
    
    exons.each { index, exon ->
    // println "  ${exon.getContig()} ${exon.getMin()} ${exon.getMax()}"
      
    
      def chr = exon.getContig()
      def min = exon.getMin()
      def max = exon.getMax()
      def range = "${chr}:${min}-${max}"
      
    
    if(alignment.overlaps(exon)) {
      
      
      def intersection = alignment.intersection(exon);
      // println "INTERSECTION:"
      // println "  ${intersection.getMin()} ${intersection.getMax()}"
      //exonmap[exon] = new String(record.getReadBases()[exon.getMin()..exon.getMax()])
      
      
        

        def xover = exon.doubleCrossover(contig)



        //println "CLIP"

        def clipped = xover.leftHardClip("-").rightHardClip("-")
        clipped.setName(new String(">${name}|gene=${options.g}|exon=${index}|location=${range}|${max - min}"))

        //println clipped.getName()
        //println clipped.sequence.seqString()

        //println "END CLIP"

        spliced.add(clipped)

        if(!regions.containsKey(index)) {
          regions[index] = []
        }

        regions[index].add(clipped)

        def offset = 0
        edits.each { lesion ->
          if(lesion.replacement.equals(SymbolList.EMPTY_LIST)) {

            // println "INSERTION"

            if(exon.contains(lesion.pos)) {
              // println "IN EXON"
              // println "Edit(${lesion.pos + offset}, 0, ${lesion.replacement})"
              // this is untested; I think this should actually be xover.sequence.edit(...)
              exon.sequence.edit(new Edit(lesion.pos + offset, 0, lesion.replacement))
              offset += lesion.length
            }
          }
        }

        // println "  ${xover.sequence.seqString()}"


        def sequenceLength = xover.sequence.seqString().replaceAll("-", "").length()

        //println ">${name}|gene=${options.g}|exon=${index}|location=${range}|sequenceLength=${sequenceLength}"
        //println xover.sequence.seqString().toUpperCase()
        
 
    } else {
      
      //summarytable[range] += 0
      //summarytable[range] += max - min
      //header += "|${exon.getContig()}:${exon.getMin()}-${exon.getMax()}(0)"
      //cdna += exon.sequence.seqString()
    }
  }
 
}

/*
ArrayList<Allele> merged = AlleleUtils.mergeAlleles(spliced)

println "MERGED SIZE = ${merged.size()}"

merged.each { allele ->
  def sequenceLength = allele.sequence.seqString().length()
  
  def fields = allele.getName().tokenize('|')
  
  //if((fields[4] as int).equals(sequenceLength)) {
    println "${allele.getName()}|${sequenceLength}"
    println allele.sequence.seqString().toUpperCase()
  //}
}
*/

def contigs = [:]

regions.each { region, list ->
  

    
    ArrayList<Allele> merged = list
    
    
    if(options.m) {
      merged = AlleleUtils.mergeAlleles(list)
    }
    
    //println "EACH"
  
    merged.each { allele ->
      def sequenceLength = allele.sequence.seqString().length()
  
      def fields = allele.getName().tokenize('|')
      def locusLength = fields[-1] as int
      
      if(sequenceLength/locusLength >= minimumBreadth) {
        if(!contigs.containsKey(fields[0])) {
          contigs[fields[0]] = [:]
        }

        if(!contigs[fields[0]].containsKey(region)) {
          contigs[fields[0]][region] = []
        }

        contigs[fields[0]][region].add(allele)

        if(!options.c) {
          println "${allele.getName()}|${sequenceLength}"
          println allele.sequence.seqString().toUpperCase()
        }
      }
  }
}

if(options.c) {
  contigs.each { contig, exons ->
    
    def cdna = ""
    exons.each { index, list ->
      best = list.sort { it.sequence.seqString().length() } .first()
      cdna += best.sequence.seqString()
    }

    println "${contig}"
    println "${cdna.toUpperCase()}"
  }
}

def CigarToEditList(record) {
  def edits = []
  
  def cigar = record.getCigar()
  def position = 0
  
  cigar.getCigarElements().each { element ->
    position += element.getLength();
    
    if(element.getOperator().toString().equals("I")) {
      // println "  new Edit(${position - element.getLength() + 1}, ${element.getLength()}, SymbolList.EMPTY_LIST)"
      
      edits += new Edit(position - element.getLength() + 1, element.getLength(), SymbolList.EMPTY_LIST)
      position -= element.getLength();
    }
    
    if(element.getOperator().toString().equals("D")) {
      def replace = DNATools.createDNA(Joiner.on("").join(Collections.nCopies(element.getLength(), "N")));
      // println "  new Edit(${position + element.getLength() - 1}, 0, ${replace.seqString()})"
      edits += new Edit(position + element.getLength() - 1, 0, replace)
    }
  }
  return edits
}

class AlleleUtils {
  static ArrayList<Allele> mergeAlleles(ArrayList<Allele> alleles) {
    def copy = alleles
    
    //println "COPY = ${copy}"
    
    ArrayList<Allele> merged = new ArrayList<Allele>()

    while (copy.size() > 1) {
      Allele first = copy[0]
      
      //println "FIRST = ${first}"
      
      def before = first
      copy.remove(0)
      def remove = []

      (0..copy.size() - 1).each { i ->
        //println "i = ${i}"
        //println "${first}"
        //println ".merge(${copy[i]})"
        Allele merge = first.merge(copy[i], 0)
        if(!merge.isEmpty()) {
          first = merge

          //println "---MERGING---"
          //println "${before}"
          //println "${before.sequence.seqString()}"
          //println "${copy[i]}"
          //println "${copy[i].sequence.seqString()}"
          //println "---EQUALS---"
          //println "${first}"
          //println "---END---"

          remove.add(i)
        }
      }
      
      //println "REMOVING"

      remove.each { i->
        //println "remove({$i})"
        copy[i] = null
      }
      
      copy = copy.findAll { it != null }
      
      merged.add(first)
    }
    
    // println "END OF MERGE LIST SIZE IS " + copy.size()
    
    while(copy.size()) {
      merged.add(copy[0])
      copy.remove(0)
    }
    
    return merged
  }
}
