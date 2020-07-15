
Channel.fromPath(params.reference).into{ch_reference_for_maskreference }
ch_vcf = Channel.fromPath(params.vcf)
ch_targetbed = Channel.fromPath(params.targets)

// First step is to use supplied VCF or BED file to mask the reference. Masked sites are lowercase, and primers should not overlap these sites.
process maskreference {
    tag "mask_${reference}_with_${vcf}"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    file(reference) from ch_reference_for_maskreference
    file(vcf) from ch_vcf

    output:
    file "*.masked.fasta.gz" into ch_fasta_for_primer3, ch_fasta_for_exact_matching, ch_fasta_for_thermo_matching, ch_fasta_for_mfe_index, ch_fasta_for_mfeprimer

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

ch_targetbed.splitText( by: 10 ).set{ ch_target_chunks }

process primer3 {
    tag "primerdesign"
    publishDir "${params.outdir}", mode: 'copy'

    cache 'deep'

    input:
    file fasta from ch_fasta_for_primer3
    file bed from ch_target_chunks

    output:
    file "${params.run_prefix}.00_primer3.tsv" into ch_primer3_chunks

    script:
    """
    
    mpcrutils.py primer3 --product_size_range '0-200' $bed $fasta > ${params.run_prefix}.00_primer3.tsv
    
    """
}

ch_primer3_chunks.collectFile(name: "00_primer_3.csv", storeDir: params.outdir, sort: true ).set{ ch_primer3 }


// This process filters primers based on some "rules" gathered from literature and the internet. A. Flag primers with GC at 3' end of primer (which stabilizes the binding site and may improve efficiency) B. Flag primers with too many GC at end of primer (too many can cause mispriming) C. Flag primers with single repeats D. Flag primers with di repeats

process heuristics {
    tag "characteristics"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    file primerpairs from ch_primer3

    output:
    file "${params.run_prefix}.01_heuristics.tsv" into ch_heuristics

    script:
    """

    mpcrutils.py gcclamp $primerpairs > ${params.run_prefix}.01a_heuristics.tsv
    mpcrutils.py gc3prime ${params.run_prefix}.01a_heuristics.tsv > ${params.run_prefix}.01b_heuristics.tsv
    mpcrutils.py max_single_repeats ${params.run_prefix}.01b_heuristics.tsv > ${params.run_prefix}.01c_heuristics.tsv
    mpcrutils.py max_dinucleotide_repeats ${params.run_prefix}.01c_heuristics.tsv > ${params.run_prefix}.01_heuristics.tsv
    
    """
}


process mfe_index {
    tag "mfeprimer index"

    input:
    file fasta from ch_fasta_for_mfe_index

    output:
    set file("foo.fasta"), file('foo.fasta.fai'), file('foo.fasta.json'), file('foo.fasta.primerqc'), file('foo.fasta.primerqc.fai') into ch_mfe_index

    script:
    """
    gzcat $fasta > foo.fasta
    mfeprimer index -i foo.fasta
    """
}


process mfeprimer {
    tag "mfeprimer"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    file index from ch_mfe_index
    file primerpairs from ch_heuristics

    output:
    file("${params.run_prefix}.02_mfeprimer") into ch_mfeprimer

    script:
    """
    mpcrutils.py convert_to_fasta $primerpairs > primers.fasta

    mfeprimer -i primers.fasta -d foo.fasta -j -o ${params.run_prefix}.02_mfeprimer
    """
}


/*
 
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
    mpcrutils exact_match --fasta $reference --in $primerpairs --columnname "em_internal" --pass ${params.run_prefix}.02_exact_match.pass.csv --fail ${params.run_prefix}.02_exact_match.fail.csv
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
    mpcrutils thermo_match --fasta $reference --in $primerpairs --columnname "thm_internal" --pass ${params.run_prefix}.03_thermo_match.pass.csv --fail ${params.run_prefix}.03_thermo_match.fail.csv
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
    mpcrutils dimerization --fasta $reference --in $primerpairs --pass ${params.run_prefix}.04_dimerization.pass.csv --fail ${params.run_prefix}.04_dimerization.fail.csv
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
      mpcrutils exact_match --fasta $reference --in $primerpairs --columnname "em_external" --pass ${params.run_prefix}.05_exact_match.pass.csv --fail ${params.run_prefix}.05_exact_match.fail.csv
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
      mpcrutils thermo_match --fasta $references --in $primerpairs --columnname "thm_external" --pass ${params.run_prefix}.06_thermo_match.pass.csv --fail ${params.run_prefix}.06_thermo_match.fail.csv
      """
  }


}

*/




