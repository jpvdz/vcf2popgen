import allel as al
import numpy as np
from .test_classes import Variants, PopGenData
import os

def read(vcf_file : str, sample_map_file : str, header : bool = True) -> PopGenData:
    sample_map = {}

    try:
        with open(sample_map_file, 'r') as fh:
            print(f"Reading sample-population mappings from: {sample_map_file}")
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
                      variants=Variants(chromosome=vcf['variants/CHROM'],
                                        position=vcf['variants/POS'],
                                        ref_alleles=vcf['variants/REF'],
                                        alt_alleles=vcf['variants/ALT']))
    return data


def test_reading_from_vcf():
    test_vcf = os.path.join(os.path.dirname(__file__), 'test.vcf')
    test_sample_map = os.path.join(os.path.dirname(__file__), 'test.csv')
    data = read(vcf_file=test_vcf,
                sample_map_file=test_sample_map,
                header=True)
    assert isinstance(data, PopGenData)


def test_reading_from_vcf_gz():
    test_vcf_gz = os.path.join(os.path.dirname(__file__), 'test.vcf.gz')
    test_sample_map = os.path.join(os.path.dirname(__file__), 'test.csv')
    data = read(vcf_file=test_vcf_gz,
                sample_map_file=test_sample_map,
                header=True)
    assert isinstance(data, PopGenData)