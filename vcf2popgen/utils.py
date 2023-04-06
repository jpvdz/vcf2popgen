def to_nucleotides(variants, ref, alt):
    """Encode bi-allelic variants stored as 0 (reference allele) or 1 
    (alternate allele) as nucleotides.
    
    Parameters
    ----------
    variants
        NumPy array of bi-allelic variants encoded as 0's and 1's.
    ref
        NumPy array of reference allele nucleotides encoded as strings.
    alt
        NumPy array of alternate allele nucleotides encoded as strings.
        
    Returns
    -------
    nucs
        NumPy array of nucleotides encoded as strings.
    """
    is_ref = variants == 0
    is_alt = variants > 0

    nucs = (ref * is_ref) + (alt * is_alt)

    return nucs


def recode_nucleotides(nucs, miss_encode = -9):
    """Recodes nucleotides as integers.
    
    Parameters
    ----------
    nucs
        NumPy array of nucleotides encoded as strings.
    miss_encode
        Integer used to encode missing data. Defaults to -9.
    Returns
    -------
    recoded_nucs
        NumPy array of nucleotides encoded as integers.
    """
    nuc_codes = {'A': 1, 'T': 2, 'C': 3, 'G': 4, '' : miss_encode}

    is_A = nucs == 'A'
    is_T = nucs == 'T'
    is_C = nucs == 'C'
    is_G = nucs == 'G'
    is_N = nucs == ''

    recoded_nucs = (
        (is_A * nuc_codes['A']) +
        (is_T * nuc_codes['T']) +
        (is_C * nuc_codes['C']) +
        (is_G * nuc_codes['G']) +
        (is_N * nuc_codes[''])
    )

    return recoded_nucs


def _create_pop_dict(sample_map):
    """Create a population dictionary with unique population labels
    mapped to integers.
    
    Parameters
    ----------
    sample_map
        Pandas dataframe or series that contains a column labeled 'pop_id'
        containing population labels encoded as strings.
        
    Returns
    -------
    pop_dict
        Dictionary of population labels (keys) and population IDs (values)
        encoded as integers.
    """
    pop_dict = {}

    for pop_id, pop_label in enumerate(sample_map['pop_id'].unique()):
        pop_dict[pop_label] = pop_id + 1

    return pop_dict
