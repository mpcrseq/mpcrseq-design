#!/usr/bin/env python3

import click
import csv
from pyfaidx import Fasta

# Primer3-py bindings
import primer3

# import methods
import p3
from heuristics import *
from utils import *



# See https://click.palletsprojects.com/en/7.x/
@click.group()
def cli():
    pass

@cli.command('primer3')
@click.argument('targets', type=click.Path(exists=True))
@click.argument('fasta', type=click.Path(exists=True))
@click.option('--product_size_range', default = '0-100')
def primer3_routine(targets, fasta, product_size_range):
    p3_results = p3.run_primer3(targets, fasta, product_size_range)
    for t in p3_results:
        print('\t'.join(str(s) for s in t))


@cli.command('gcclamp')
@click.argument('path')
def gcclamp_routine(path):
    # Expects file to be tsv with pl_seq and pr_seq columns containing the forward and reverse sequences for primer pairs. Output will be the input file appended with a gcclamp column for each primer.
    f = open(path,'r')
    reader = csv.reader(f, delimiter="\t")
    header = next(reader)

    pl_col = header.index('pl_seq')
    pr_col = header.index('pr_seq')

    header.extend(['pl_gc_clamp', 'pr_gc_clamp'])
    
    print(*header, sep = '\t')

    for row in reader:
        pl_seq = row[pl_col]
        pr_seq = row[pr_col]

        pl_clamp = gcclamp(pl_seq)
        pr_clamp = gcclamp(pr_seq)
        row.extend([pl_clamp, pr_clamp])
        print(*row, sep = '\t')


@cli.command('gc3prime')
@click.argument('path')
def gc3prime_routine(path):
    # Expects file to be tsv with pl_seq and pr_seq columns containing the forward and reverse sequences for primer pairs. Output will be the input file appended with a gcclamp column for each primer.
    f = open(path,'r')
    reader = csv.reader(f, delimiter="\t")
    header = next(reader)

    pl_col = header.index('pl_seq')
    pr_col = header.index('pr_seq')

    header.extend(['pl_gc_3prime', 'pr_gc_3prime'])
    
    print(*header, sep = '\t')

    for row in reader:
        pl_seq = row[pl_col]
        pr_seq = row[pr_col]

        pl_out = GC_3prime(pl_seq)
        pr_out = GC_3prime(pr_seq)
        row.extend([pl_out, pr_out])
        print(*row, sep = '\t')


@cli.command('max_single_repeats')
@click.argument('path')
def max_single_repeats_routine(path):
    # Expects file to be tsv with pl_seq and pr_seq columns containing the forward and reverse sequences for primer pairs. Output will be the input file appended with a max_single_repeats column for each primer.
    f = open(path,'r')
    reader = csv.reader(f, delimiter="\t")
    header = next(reader)

    pl_col = header.index('pl_seq')
    pr_col = header.index('pr_seq')

    header.extend(['pl_si_repeats', 'pr_si_repeats'])
    
    print(*header, sep = '\t')

    for row in reader:
        pl_seq = row[pl_col]
        pr_seq = row[pr_col]

        pl_out = max_single_repeats(pl_seq)
        pr_out = max_single_repeats(pr_seq)
        row.extend([pl_out, pr_out])
        print(*row, sep = '\t')


@cli.command('max_dinucleotide_repeats')
@click.argument('path')
def max_dinucleotide_repeats_routine(path):
    # Expects file to be tsv with pl_seq and pr_seq columns containing the forward and reverse sequences for primer pairs. Output will be the input file appended with a max_dinucleotide_repeats column for each primer.
    f = open(path,'r')
    reader = csv.reader(f, delimiter="\t")
    header = next(reader)

    pl_col = header.index('pl_seq')
    pr_col = header.index('pr_seq')

    header.extend(['pl_di_repeats', 'pr_di_repeats'])
    
    print(*header, sep = '\t')

    for row in reader:
        pl_seq = row[pl_col]
        pr_seq = row[pr_col]

        pl_out = max_dinucleotide_repeats(pl_seq)
        pr_out = max_dinucleotide_repeats(pr_seq)
        row.extend([pl_out, pr_out])
        print(*row, sep = '\t')


@cli.command('convert_to_fasta')
@click.argument('path')
def convert_to_fasta_routine(path):
    # Expects file to be tsv with pl_seq and pr_seq columns containing the forward and reverse sequences for primer pairs. Output will be the input file appended with a max_dinucleotide_repeats column for each primer.
    f = open(path,'r')
    reader = csv.reader(f, delimiter="\t")
    header = next(reader)

    any_primers_col = header.index('any_primers')

    pl_seq_col   = header.index('pl_seq')
    pl_start_col = header.index('pl_start')
    pl_end_col   = header.index('pl_end')

    pr_seq_col   = header.index('pr_seq')
    pr_start_col = header.index('pr_start')
    pr_end_col   = header.index('pr_end')

    name_col = header.index('name')

    for row in reader:
        if row[any_primers_col] == 'TRUE':
            print('>'+ row[name_col] + ':' + 'pl' + ':' + row[pl_start_col] + '-' + row[pl_end_col])
            print(row[pl_seq_col])

            print('>'+ row[name_col] + ':' + 'pr' + ':' + row[pr_start_col] + '-' + row[pr_end_col])
            print(row[pr_seq_col])


if __name__ == '__main__':
    cli()








