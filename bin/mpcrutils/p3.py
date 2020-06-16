import primer3 as p3
from pyfaidx import Fasta


def run_primer3(targets_bed, fasta_path, product_size_range):
    product_size_range = [int(i) for i in product_size_range.split('-')] 
    
    # Set global arguments for primer 3
    # Need to set this up to pass variables to setGlobals, or write arguments to a file for setGlobals to read.
    p3_setGlobals()
    
    # targets_bed is a tsv with targets: chrom start end name
    targets = open(targets_bed, 'r')
    ref     = Fasta(fasta_path)
    
    # Place successful outputs here as a list of tuples
    p3_results = [('chrom','start','end','name','any_primers','pl_start', 'pl_end', 'pl_len', 'pl_seq', 'pr_start', 'pr_end', 'pr_len', 'pr_seq', 'product_size', 'amplicon_size')]
    
    for line in targets:
        target = line.strip('\n').split('\t')
        # If the name field doesn't exist, make it an empty string
        target_name = '' if len(target) < 4 else target[3]

        # first we get the sequence template for this target from the input fasta
        # The template is the maximum primer product size range before and after the target
        # This is a little larger than necessary, but it doesn't matter.
        template_padding  = max(product_size_range)
        template_start    = int(target[1]) - template_padding
        template_end      = int(target[2]) + template_padding
        template_seq      = ref[target[0]][template_start:template_end]
        template_len      = len(template_seq)
        target_length     = int(target[2]) - int(target[1])
        
        # Now design primers for the target using primer3
        p3out = p3.bindings.designPrimers(
            {
                'SEQUENCE_TEMPLATE': str(template_seq),
                'SEQUENCE_INCLUDED_REGION': [0, len(template_seq)],
                'SEQUENCE_TARGET': [template_padding, target_length]
            })
        
        # Primer3 may return multiple pairs of primers for the target.
        nreturned = p3out['PRIMER_PAIR_NUM_RETURNED']
        if nreturned > 0:
            for i in range(nreturned):
                i = str(i)
                # Primer3 returns the 0-based start of the primer, and the length. By adding the template start we get back the reference coordinates.
                pl_start = template_start + int(p3out['PRIMER_LEFT_' + i][0])
                pl_len   = p3out['PRIMER_LEFT_' + i][1]
                pl_end   = pl_start + int(pl_len)
                pl_seq   = p3out['PRIMER_LEFT_'+ i +'_SEQUENCE']
                
                pr_start = template_start + int(p3out['PRIMER_RIGHT_' + i][0])
                pr_len   = p3out['PRIMER_RIGHT_' + i][1]
                pr_end   = pr_start + int(pr_len)
                pr_seq   = p3out['PRIMER_RIGHT_'+ i +'_SEQUENCE']
                
                product_size = p3out['PRIMER_PAIR_'+ i +'_PRODUCT_SIZE']
                amplicon_size = int(product_size) - int(pl_len) - int(pr_len)
                
                # Append tuple with target information and found primer information to 'found_primers' list
                p3_results.append((target[0], target[1], target[2], target_name, 'TRUE', pl_start, pl_end, pl_len, pl_seq, pr_start, pr_end, pr_len, pr_seq, product_size, amplicon_size))
        else:
            p3_results.append((target[0], target[1], target[2], target_name, 'FALSE', '', '', '', '', '', '', '', '', '' ))

    return(p3_results)




def p3_setGlobals():
    global_args =  {
        'PRIMER_OPT_SIZE' : 20,
        'PRIMER_PICK_INTERNAL_OLIGO' : 1,
        'PRIMER_INTERNAL_MAX_SELF_END' : 8,
        'PRIMER_MIN_SIZE' : 18,
        'PRIMER_MAX_SIZE' : 25,
        'PRIMER_OPT_TM' : 60.0,
        'PRIMER_MIN_TM' : 57.0,
        'PRIMER_MAX_TM' : 63.0,
        'PRIMER_MIN_GC' : 20.0,
        'PRIMER_MAX_GC' : 80.0,
        'PRIMER_MAX_POLY_X' : 100,
        'PRIMER_INTERNAL_MAX_POLY_X' : 100,
        'PRIMER_SALT_MONOVALENT' : 50.0,
        'PRIMER_DNA_CONC' : 50.0,
        'PRIMER_MAX_NS_ACCEPTED' : 0,
        'PRIMER_MAX_SELF_ANY' : 12,
        'PRIMER_MAX_SELF_END' : 8,
        'PRIMER_PAIR_MAX_COMPL_ANY' : 12,
        'PRIMER_PAIR_MAX_COMPL_END' : 8,
        'PRIMER_PRODUCT_SIZE_RANGE' : [37,100]
    }
    p3.bindings.setP3Globals(global_args)

