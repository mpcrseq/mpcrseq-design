# Heuristics utilities

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


