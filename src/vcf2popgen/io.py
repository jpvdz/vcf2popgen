import allel as al
import numpy as np
from .classes import Variants, PopGenData

def read(vcf_file : str, sample_map_file : str, header : bool = True) -> PopGenData:
    """Read data from a VCF file and sample map file into NumPy arrays. Only data
    from samples listed in the sample map is read.  
    
    Parameters
    ----------
    input_file : str
        Path to a VCF file on the local filesystem containing bi-allelic SNPs.
    sample_map_file : str
        Path to a sample map file on the local filesystem in CSV format. The first 
        column should contain sample IDs, the second population IDs encoded as 
        integers. 
    
    Returns
    -------
    data : PopGenData
        A PopGenData object that holds bi-allelic SNPs from a set of samples from 
        one or more populations.
    """
    
    sample_map = {}

    try:
        print(f"Reading sample-population mappings from: {sample_map_file}")
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
        print(f"Reading genotypic data from: {vcf_file}")
        vcf = al.read_vcf(vcf_file, alt_number=1, samples=sample_map.keys())
    except:
        raise ValueError(
    """Could not read data from VCF file. Check whether it exists, is properly 
    formatted and/or has the correct file extension."""
            )
    
    data = PopGenData(samples=vcf['samples'], 
                      populations=np.array(list(sample_map.values()), dtype=int), 
                      genotypes=al.GenotypeArray(vcf['calldata/GT']), 
                      ref_variants=Variants(vcf['variants/REF']),
                      alt_variants=Variants(vcf['variants/ALT']))
    
    return data

