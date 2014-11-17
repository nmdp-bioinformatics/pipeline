import org.apache.commons.cli.Option

import groovy.io.FileType

CliBuilder cli = new CliBuilder()

cli.with {
  e longOpt: 'expected-file', args: 1, 'file of expected genotype calls'
  o longOpt: 'observed-file', args: 1, 'file of observed genotype calls'
  z longOpt: 'resolution', args: 1, 'minimum fields of resolution by which to match expected and observed alleles'
  s longOpt: 'summary-report', 'print report summary'
}

def options = cli.parse(args)

assert options

if (options.h) {
 cli.usage()
}

assert 1 == Match.byField("HLA-A\\*01", "HLA-A\\*01")
assert 1 == Match.byField("HLA-A\\*01:01", "HLA-A\\*01")
assert 1 == Match.byField("HLA-A\\*01:02:01", "HLA-A\\*01:01:01")
assert 2 == Match.byField("HLA-DQB1\\*05:02:01", "HLA-DQB1\\*05:02:04")
assert 3 == Match.byField("HLA-A\\*01:01:01", "HLA-A\\*01:01:01")
assert 4 == Match.byField("HLA-A\\*01:01:01:01", "HLA-A\\*01:01:01:01")
assert 0 == Match.byField("HLA-DRB1\\*13:01:01", "HLA-DRB3\\*01:01:02:01")
assert 0 == Match.byField("HLA-DRB1\\*13:01:01", "HLA-DRB3\\*01:01:02:02")

def minimumResolution = 2

if(options.z) {
  minimumResolution = options.z as int
}

def expected = [:]
new File(options.e).each { line ->
  (sample, locus, regionsFile, zygosity, firstAllele, secondAllele) = line.split("\\s+")
  if(!expected.containsKey(sample)) {
    expected[sample] = []
  }
    
  expected[sample].add(firstAllele)
  expected[sample].add(secondAllele)
}

def observed = [:]
new File(options.o).each { line ->
  (sample, interpretation) = line.split("\t")
  if(!observed.containsKey(sample)) {
    observed[sample] = []
  }
    
  observed[sample].add(interpretation)
}

/*
println "observed"
observed.each { sample, glstrings ->
  println "${sample} ${glstrings}"
}

println "expected"
expected.each { sample, alleles ->
  println "${sample} ${alleles[0]} ${alleles[1]}"
}
*/

def summary = ["PASS":0, "FAIL":0]
expected.each { sample, alleles ->
  
  
  /*
   * Add HLA- and * to allele names
   */

  /*
   * Split string on newline to give one phased interpreted allele set with ambiguity per consensus sequence
   */
  def interpretations = observed[sample]
  //println "contigGlStrings = ${contigGlStrings}"
  //
  def product = []
  
  alleles.each { expectedAllele ->
    def pass = "FAIL"
   interpretations.each { interpretation ->
      def interpretedAlleles = interpretation.split(/[\/|]/)

      //println "EXPECTED ALLELE = ${expectedAllele}"
      //println "INTERPRETED ALLELES = ${interpretedAlleles}"
     
      found = interpretedAlleles.findAll { Match.byField(expectedAllele, it) >= minimumResolution}

   
      //println "FOUND by match = ${found}"
      
      if(!found.isEmpty()) {
        pass = "PASS"
      }
    }
    
    summary[pass]++

    if(!options.s) {
      println "${pass}\t${sample}\t${expectedAllele}"
    }
  }
}

if(options.s) {
  summary.each { category, count ->
    println "${category}\t${count}"
  }
}

class Match {
  static int byField(first, second) {
    String[] firstAlleles = first.split(":")
    String[] secondAlleles = second.split(":")
    int smallest = firstAlleles.size() < secondAlleles.size() ? firstAlleles.size() : secondAlleles.size()

    for(int i = 0; i < smallest; i++) {
      if(!firstAlleles.getAt(i).equals(secondAlleles.getAt(i))) {
        return i
      }
    }
    
    return smallest
  }
}