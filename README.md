## mpcrseq-design
A pipeline to design PCR primer pools for hundreds of targets.

## Introduction
Brief description here.

## Pipeline steps

1. Variant Masking
2. Primer design
3. Primer characteristics checking
4. Primer off-target exact matching
5. Primer off-target inexact matching
6. Primer pair dimerization
7. Non-target genome matching
8. Filter primer pairs for compatability

### Variant Masking
(optional) Masks known variants so that primers do not overlap variable loci.

**Inputs** Fasta genome file, bed file with locations to mask  
**Outputs** Masked fasta (Ns in place of all variable loci)


### Mappability masking
 (optional) https://github.com/cpockrandt/genmap


### Variant Filtering
Check for homo-polymer regions


### Primer design
Uses primer3 to find primers for each target locus.

**Inputs** Masked fasta genome, bed file with targets, config variables  
  config variables: Tm range, primer size range, amplicon size range  
**Outputs** primer_pairs_tsv  
primer_pairs_tsv: primer_pair_name(f_start,f_end,r_start,r_end) target_start target_end f_start f_end r_start r_end f_seq r_seq amplicon_start amplicon_end amplicon_seq f_gc r_gc f_tm r_tm


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


## Quick-start guide

## Requirements

[Primer3](https://github.com/primer3-org/primer3)

[Python 3](https://www.python.org/downloads/)

[Jellyfish](https://github.com/gmarcais/Jellyfish)

[vsearch](https://github.com/torognes/vsearch)

### Install requirements on a comptuer running Mac OS 10.x
The easiest way is to install the [homebrew package manager](https://brew.sh). Then install the external software requirements.
Here we will use [pipenv](https://docs.python-guide.org/dev/virtualenvs/) to ensure that compatible versions of the python module requirements are installed.


    # Install python3
    brew install python3
    pip3 install pipenv
    
    # Now install the required python modules
    
    cd PATH/TO/mpcrseq
    pipenv install
    
    # Get access to the homebrew biology packages
    brew tap brewsci/bio
    # Install external programs using homebrew
    brew install primer3
    brew install jellyfish
    brew install vsearch
    

### Install requirements on a comptuer running Ubuntu

    # Install python3 and other external requirements
    sudo apt-get install -y python3 primer3 jellyfish vsearch
    # Install pipenv for easier python package management
    sudo pip3 install pipenv
    git clone https://github.com/rwtaylor/mpcrseq.git
    cd mpcrseq
    pipenv install

## Example
Compute primers for a small genome and a few targets.

### First pre-compute the jellyfish kmers (for exact off-target hits)
Note that unless you modify the makefile this will use all available CPUs.
The range of kmer sizes should contain all the primer sizes that you wish to search for.

    cd mpcrseq/examples/small_example/jellyfish_kmers && make


###Try running mpcrseq

    pipenv shell
    cd examples/small_example
    python3 ../../mpcrseq --config=mpcrseq.conf
    

## Input Files
1. Reference genome in fasta format.
2. A VCF or BED file that includes all polymorphisms (indels and snps) in the genome.
3. A VCF or BED file that includes potential SNP targets.
4. Any genomes (in fasta format) that should be avoided as primer targets. For example, prey genomes if primers are being designed to amplify a predator's dna from feces.

The paths to these files will be input with the config file.


## Outputs

## Files
mpcrdesign.py - The main functions and control flow are all in here for now...