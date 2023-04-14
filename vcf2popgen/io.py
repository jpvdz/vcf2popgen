import allel as al
import numpy as np
from .formats import bayescan, genepop, structure

def read(vcf_file : str, sample_map_file : str, header : bool = True) -> dict:
    """Read data from a VCF file and sample map file into NumPy arrays. Only data
    from samples listed in the sample map is read.  
    
    Parameters
    ----------
    input_file : str
        Path to a VCF file on the local filesystem containing bi-allelic variant calls.
    sample_map_file : str
        Path to a sample map file on the local filesystem in CSV format. The first 
        column should contain sample IDs, the second population IDs encoded as integers. 

    Returns
    -------
    data : dict
        Dictionary containing arrays holding bi-allelic variant calls from samples and
        populations listed in the sample map.
    """
    sample_map = {}

    try:
        with open(sample_map_file, 'r') as fh:
            for i, line in enumerate(fh):
                if i == 0 and header == True:
                    continue
                else:
                    key, value = line.strip().split(',')
                    sample_map[key] = int(value)
    except:
        raise ValueError(
    """Could not read data from sample map file. Check whether it exists, is 
    properly formatted and/or has the correct file extension."""
            )
    
    try: 
        data = al.read_vcf(vcf_file, alt_number=1, samples=sample_map.keys())
    except:
        raise ValueError(
    """Could not read data from VCF file. Check whether it exists, is properly 
    formatted and/or has the correct file extension."""
            )
    
    data['populations'] = np.array(list(sample_map.values()), dtype=int)
    
    return data


def write(data: dict, format : str, output_file : str) -> None:
    """Write bi-allelic variant calls to an output file in various formats.
    
    Parameters
    ----------
    data : dict
        Dictionary containing arrays holding VCF data and population IDs.
    format : str
        The output file format. Possible values are 'bayescan', 'genepop'
        or 'structure'.
    output_file : str
        Name or path of output file. 
    """
    
    match format:
        case 'bayescan':
            bayescan(data, output_file)
            
        case 'genepop':
            genepop(data, output_file)
            
        case 'structure':
            structure(data, output_file)
            
        case _:
            raise ValueError(
    """Invalid output format. Valid options are: 'bayescan', 'genepop' and 'structure'."""
            )
