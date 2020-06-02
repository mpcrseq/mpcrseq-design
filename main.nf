


Channel.fromPath(params.reference).into{ch_reference_for_maskreference; ch_reference_for_primer3; ch_reference_for_int_therm_matching}
ch_vcf = Channel.fromPath(params.vcf)
ch_targetbed = Channel.value(params.targets)

// First step is to use supplied VCF or BED file to mask the reference. Masked sites are lowercase, and primers should not overlap these sites.
process maskreference {
    tag "mask_${reference}_with_${vcf}"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    file(reference) from ch_reference_for_maskreference
    file(vcf) from ch_vcf

    output:
    file "*.masked.fasta.gz" into ch_fasta_for_primer3, ch_fasta_for_exact_matching, ch_fasta_for_thermo_matching

    script:
    """
    
    bgzip -d $reference
    bgzip -d $vcf
    
    bedtools maskfasta -fi ${reference.baseName} -bed ${vcf.baseName} -fo ${reference.simpleName}.masked.fasta
    bgzip -@${task.cpus} *.masked.fasta
    """
}


// This process uses primer3 to design primers for each supplied target. After designing the primers they are filtered to meet any further criteria specified in the pipeline parameters. Then a TSV file is output with the putative primer information.
// primer3_csv: target_name target_start target_end f_start f_end r_start r_end f_seq r_seq amplicon_start amplicon_end amplicon_seq f_gc r_gc f_tm r_tm pass_filter
process primer3 {
    tag "primerdesign"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    file fasta from ch_fasta_for_primer3
    file bed from ch_targetbed

    output:
    file "*.00_primer3.bed" into ch_primer3

    script:
    """
    
    mpcrutils primer3 --ref $fasta --target_bed $bed > ${params.run_prefix}.00_primer3.csv
    
    """
}


// This process filters primers based on some "rules" gathered from literature and the internet. A. Flag primers with GC at 3' end of primer (which stabilizes the binding site and may improve efficiency) B. Flag primers with too many GC at end of primer (too many can cause mispriming) C. Flag primers with single repeats D. Flag primers with di repeats
// heuristics_csv: target_name target_start target_end f_start f_end r_start r_end f_seq r_seq amplicon_start amplicon_end amplicon_seq f_gc r_gc f_tm r_tm pass_filter pass_GC3prime pass_GCclamp pass_SIrepeats pass_DIrepeats
process heuristics {
    tag "characteristics"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    file "primerpairs" from ch_primer3

    output:
    file "*.01_heuristics.pass.bed" into ch_heuristics

    script:
    """
    mpcrutils heuristics --in $primerpairs --pass ${run_prefix}.01_heuristics.pass.csv --fail ${run_prefix}.01_heuristics.fail.csv
    """
}
 
 
// Performs exact matching for off-target priming internally and in external genomes.


ch_internal_fasta = Channel.fromPath("$params.reference").map { it -> 
  ["genome", it]
}

ch_external_fastas = Channel.from(params.external_fastas).map{it ->
  name = it[0]
  fasta = it[1]
  [name, fasta]
}

ch_internal_fasta.combine(ch_external_fastas).into{ch_fastas_for_exact_matching; ch_fastas_for_thermal_matching}

process exact_matching {
    tag "extact_match"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    file "primerpairs" from ch_heuristics
    file fasta from ch_fastas_for_exact_matching

    output:
    file "*.02_int_exact_match.pass.tsv" into ch_int_exact_match

    script:
    """
    mpcrutils exact_match --fasta $reference --in $primerpairs --columnname "em_internal" --pass ${run_prefix}.02_exact_match.pass.csv --fail ${run_prefix}.02_exact_match.fail.csv
    """
}


// Performs thermodynamic matching for off-target priming internally, within the target species genome.
process int_therm_matching {
    tag "int_thm_match"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    file "primerpairs" from ch_int_exact_match
    file reference from ch_reference_for_int_therm_matching

    output:
    file "*.03_int_thermo_match.pass.tsv" into ch_thermo_exact_match

    script:
    """
    mpcrutils thermo_match --fasta $reference --in $primerpairs --columnname "thm_internal" --pass ${run_prefix}.03_thermo_match.pass.csv --fail ${run_prefix}.03_thermo_match.fail.csv
    """
}

  
// Checks for pairwise primer dimerization.
process dimerization {
    tag "dimer"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    file "primerpairs" from ch_thermo_exact_match

    output:
    file "*.04_dimerization.pass.tsv" into ch_dimerization

    script:
    """
    mpcrutils dimerization --fasta $reference --in $primerpairs --pass ${run_prefix}.04_dimerization.pass.csv --fail ${run_prefix}.04_dimerization.fail.csv
    """
}


// Performs thermodynamic matching for off-target priming externally, within other species genomes
if(params.ch_ext_references) {

  ch_ext_references = Channel.value(params.external_references)

  process ext_exact_matching {
      tag "ext_ext_match"
      publishDir "${params.outdir}", mode: 'copy'

      input:
      file "primerpairs" from ch_dimerization
      file references from ch_ext_references

      output:
      file "*.02_ext_exact_match.pass.tsv" into ch_ext_exact_match

      script:
      """
      mpcrutils exact_match --fasta $reference --in $primerpairs --columnname "em_external" --pass ${run_prefix}.05_exact_match.pass.csv --fail ${run_prefix}.05_exact_match.fail.csv
      """
  }

  // Performs thermodynamic matching for off-target priming externally, within other species genomes.
  process ext_therm_matching {
      tag "ext_thm_match"
      publishDir "${params.outdir}", mode: 'copy'

      input:
      file "primerpairs" from ch_int_exact_match
      file references from ch_ext_references

      output:
      file "*.03_ext_thermo_match.pass.tsv" into ch_ext_exact_match

      script:
      """
      mpcrutils thermo_match --fasta $references --in $primerpairs --columnname "thm_external" --pass ${run_prefix}.06_thermo_match.pass.csv --fail ${run_prefix}.06_thermo_match.fail.csv
      """
  }


}




