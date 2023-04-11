import numpy.typing as npt

def to_nucleotides(variants : npt.ArrayLike, ref : npt.ArrayLike, alt : npt.ArrayLike) -> npt.ArrayLike:
    """Encode bi-allelic variants stored as 0 (reference allele) or 1 
    (alternate allele) as nucleotides.
    
    Parameters
    ----------
    variants : npt.ArrayLike
        NumPy array of bi-allelic variants encoded as 0's and 1's.
    ref : npt.ArrayLike
        NumPy array of reference allele nucleotides encoded as strings.
    alt : npt.ArrayLike
        NumPy array of alternate allele nucleotides encoded as strings.
        
    Returns
    -------
    nucs : npt.ArrayLike
        NumPy array of nucleotides encoded as strings.
    """
    is_ref = variants == 0
    is_alt = variants > 0

    nucs = (ref * is_ref) + (alt * is_alt)

    return nucs


def recode_nucleotides(nucs : npt.ArrayLike, miss_encode : int = -9) -> npt.ArrayLike:
    """Recodes nucleotides as integers.
    
    Parameters
    ----------
    nucs : npt.ArrayLike
        NumPy array of nucleotides encoded as strings.
    miss_encode : int
        Integer used to encode missing data. Defaults to -9.
    Returns
    -------
    recoded_nucs : npt.ArrayLike
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
