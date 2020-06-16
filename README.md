## mpcrseq-design
A [nextflow](https://nextflow.io) based pipeline to design PCR primer pools for hundreds of targets.

## Introduction
Brief description here.

## Quick-start

  1. Install [Nextflow](https://www.nextflow.io).

  2. Use docker (very slow on macs), or install the following software:  
    htslib or tabix (for bgzip)  
    python3 and packages: biopython primer3-py pyfaidx  
    bedtools  

  2. Test the pipeline  
  with docker:
    nextflow run -profile test,docker

  without docker:
    nextflow run -profile test


  3. Run with your own data  
    nextflow run -profile docker -reference reference.fasta -targets targets.bed

## Pipeline steps

  1. Repeat Masking (important step to be completed before running the pipeline.)
  2. Variant Masking
  3. Mappability Masking
  4. Target Homopolymer Filtering
  5. Primer design
  6. Primer characteristics checking
  7. Primer off-target exact matching
  8. Primer off-target inexact matching
  9. Primer pair dimerization
  10. Non-target genome matching
  11. Filter primer pairs for compatability

### Repeat Masking
The input reference genome must be repeat masked. This pipeline does not perform repeat masking. [RepeatMasker](http://www.repeatmasker.org) will do the job.


### Variant Masking
(optional) Masks known variants so that primers do not overlap variable loci.

**Inputs** Fasta genome file, bed file with locations to mask  
**Outputs** Masked fasta (Ns in place of all polymorphic loci)


### Mappability masking
 (optional) https://github.com/cpockrandt/genmap

**Inputs** Variant masked genome file.
**Outputs** Masked fasta (Ns in place of all regions with low mappability)


### Target Homopolymer Filtering
(optional) Check for homo-polymer regions overlapping or very near target polymorphism. Homo-polymers have higher error rates when sequenced on Oxford Nanopore sequencers.

**Inputs** Target bed file
**Outputs**
  Filtered target bed file
  Rejected targets with reason (homopolymer overlap)

### Primer design
Uses primer3 to find primers for each target locus.

**Inputs** Masked fasta genome, bed file with targets, config variables  
  config variables: Tm range, primer size range, amplicon size range  
**Outputs**
  primer_pairs_tsv  
    primer_pairs_tsv: primer_pair_name(f_start,f_end,r_start,r_end) target_start target_end f_start f_end r_start r_end f_seq r_seq amplicon_start amplicon_end amplicon_seq f_gc r_gc f_tm r_tm
  Rejected targets with reason (no primers)

### Primer characteristics checking
A. Flag primers with GC at 3' end of primer (which stabilizes the binding site and may improve efficiency)
B. Flag primers with too many GC at end of primer (too many can cause mispriming)
C. Flag primers with single repeats
D. Flag primers with di repeats

**Inputs** primer_pairs_tsv  
**Outputs** primer_pairs_tsv with A-D flags


### Primer off-target exact matching
Flag any primer pairs that match off-target sequences. Exact matching is much faster than inexact matching, so we first remove exact matches.

**Inputs** primer_pairs_tsv  
**Outputs** primer_pairs_tsv flagged with exact matches


### Primer off-target inexact matching
Flag primer pairs that match thermodynamically to off-target sequences.

**Inputs** primer_pairs_tsv  
**Outputs** primer_pairs_tsv annotated with off-target Tm


### Primer pair dimerization
Check for pairwise primer dimers, including adapter sequences (primer dimers across pairs of primers). Primer3 (step 2) checks for within pair dimerization.

**Inputs** primer_pairs_tsv  
**Outputs** primer_pairs_tsv, primer_dimer_matrix
primer_dimer_matrix: a matrix with primer_dimer tms. [Not sure this is the best way to store this info]


### Non-target genome matching
Repeat off-target matching steps for non-target genomes.

**Inputs** primer_pairs_tsv  
**Outputs** primer_pairs_tsv annotated with non-target Tm


### Filter primer pairs for compatability
Filter primers based on the previous steps.

**Inputs** primer_pairs_tsv, primer_dimer_matrix  
**Outputs** compatible set(s) of primers

