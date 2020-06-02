#!/usr/bin/env python3

import click
import csv
import primer3

def global_args():
  



def run_primer3(fasta, targets_bed):
  # Set global arguments

    tbed = open(targets_bed, 'r'):
    fasta = open(fasta, 'r')
    for line in tbed:
        # first we get the sequence template for this target from the input fasta
        # The template is the maximum primer product size range before and after the target
        # This is a little larger than necessary, but it doesn't matter.
        seq_start    = target_loc - max(product_size_range)
        seq_end      = target_loc + max(product_size_range)
        # Maybe there is a python package to do this.
        seq_template = `samtools faidx chrom:seq_start-seq_end`
        # 
    
    
    bindings.designPrimers(
        {
            'SEQUENCE_ID': target_id,
            'SEQUENCE_TEMPLATE': seq_template,
            'SEQUENCE_INCLUDED_REGION': [0,len(seq_template)],
            'SEQUENCE_TARGET': seq_target
        },
        global_args)
    
    
    

def gcclamp(sequence):  # Returns True if 2 or fewer G or C nucleotides are within the last 5 bases.
    #Oligonucleotides with more than three G or C nucleotides within the last  five bases were  filtered.  This should help minimize mispriming at GC-rich binding sites
    # http://www.bio.net/mm/methods/2000-August/084438.html
    sequence = sequence.upper()
    sequence = sequence[-5:].replace("C","G")
    return(sequence.count('G') <= 2)

def max_single_repeats(seq):
    seq_length = int(len(seq))
    max = 1
    for i in range(seq_length-1):
        j = 1
        while (seq[i] == seq[i+j] and ((i+j) < (seq_length -1))):
            j += 1
        if j > max: max = j
    return(max)

def max_dinucleotide_repeats(seq):
    seq_length = int(len(seq))
    max = 1
    for i in range(seq_length-2):
        j = 2
        while (seq[i:i+1] == seq[i+j:i+j+1]) and (i+j+2 < seq_length - 2):
            j += 2
        if (j/2) > max: max = int(j/2)
    return(max)

def GC_3prime(sequence):  # Returns True if 3' end of oligo is G or C
    # Rationale:
    # For primers, a single G or C nucleotide at the 3'-end helps to stabilize
    # binding near the site of extension, which can reduce the possibility of
    # "breathing" and improves priming efficiency.  terefore, primers ending
    # in an A or a T base could optionally be  altered. In this study, the 
    # filter was applied. [Francis et al. 2017 - ThermoAlign: a genome-aware
    # primer design tool for tiled amplicon resequencing]
    # http://www.bio.net/mm/methods/2000-August/084438.html
    sequence = sequence.upper()
    return(set(sequence[-1:]) <= set('GC'))



@click.group()
def cli():
    pass

@cli.command('primer3')
@click.argument('ref')
@click.argument('target_bed')
def primer3_routine(path):
    


@cli.command('gcclamp')
@click.argument('path')
def gcclamp_routine(path):
    # Expects file to be tsv with f_seq and r_seq columns containing the forward and reverse sequences for primer pairs. Output will be the input file appended with a gcclamp column for each primer.
    f = open(path,'r')
    reader = csv.reader(f, delimiter="\t")
    header = next(reader)

    f_col = header.index('f_seq')
    r_col = header.index('r_seq')

    header.extend(['f_gc_clamp', 'r_gc_clamp'])
    
    print(*header, sep = '\t')

    for row in reader:
        f_seq = row[f_col]
        r_seq = row[r_col]

        f_clamp = gcclamp(f_seq)
        r_clamp = gcclamp(r_seq)
        row.extend([f_clamp, r_clamp])
        print(*row, sep = '\t')


@cli.command('gc3prime')
@click.argument('path')
def gcclamp_routine(path):
    # Expects file to be tsv with f_seq and r_seq columns containing the forward and reverse sequences for primer pairs. Output will be the input file appended with a gcclamp column for each primer.
    f = open(path,'r')
    reader = csv.reader(f, delimiter="\t")
    header = next(reader)

    f_col = header.index('f_seq')
    r_col = header.index('r_seq')

    header.extend(['f_gc_3prime', 'r_gc_3prime'])
    
    print(*header, sep = '\t')

    for row in reader:
        f_seq = row[f_col]
        r_seq = row[r_col]

        f_out = GC_3prime(f_seq)
        r_out = GC_3prime(r_seq)
        row.extend([f_out, r_out])
        print(*row, sep = '\t')

if __name__ == '__main__':
    cli()








