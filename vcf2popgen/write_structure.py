import allel as al
import numpy as np
from .utils import to_nucleotides, recode_nucleotides

def write_structure(input_file, sample_map_file, output_file, one_row_per_sample = False, header = True):
    """Write bi-allelic variant calls to an output file in GENEPOP format.
    
    Parameters
    ----------
    input_file
        A VCF file containing bi-allelic variant calls.
    sample_map_file
        A sample map file in CSV format. The first column should contain sample
        IDs, the second population IDs encoded as integers. 
    output_file
        The name of the output file.
    one_row_per_sample
        Whether samples should be encoded on a single row, which results in two
        columns per locus (one for each allele). Defaults to False.
    header
        Whether the sample map contains a header. Defaults to True.
    """

    try: 
        vcf = al.read_vcf(input_file)
    except:
        raise ValueError(
            """Could not read data from VCF file. Check whether it exists,
            is properly formatted and/or has the correct file extension."""
            )
    
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
            """Could not read data from sample map file. Check whether
            it exists, is properly formatted and/or has the correct file 
            extension."""
            )
    
    gt = al.GenotypeArray(vcf['calldata/GT'])

    sample_ids = vcf['samples']
    ref_allele = vcf['variants/REF']
    alt_allele = vcf['variants/ALT'].transpose()[0]

    pops = np.array([
        sample_map[sample_id] for sample_id in vcf['samples']
    ])
    
    if one_row_per_sample == False:
        variant_ids = [f"variant_{x+1}" for x in range(len(vcf['variants/ID']))]
         
        with open(output_file, 'w') as fh:
            variant_cols = '\t'.join(str(variant_id) for variant_id in variant_ids)
            fh.write(f'\t\t{variant_cols}\n')

            for i, sample_id in enumerate(sample_ids):
                if sample_id in sample_map.keys():
                    alleles_recoded0 = recode_nucleotides(to_nucleotides(gt[:, i][:, 0], ref_allele, alt_allele))
                    alleles_recoded1 = recode_nucleotides(to_nucleotides(gt[:, i][:, 1], ref_allele, alt_allele))

                    alleles0 = '\t'.join(str(allele) for allele in alleles_recoded0)
                    alleles1 = '\t'.join(str(allele) for allele in alleles_recoded1)

                    fh.write(f'{sample_id}\t{pops[i]}\t{alleles0}\n')
                    fh.write(f'{sample_id}\t{pops[i]}\t{alleles1}\n')

                else:
                    continue
    
    elif one_row_per_sample == True:
        variant_ids0 = [f'variant_{x+1}_1' for x in range(len(vcf['variants/ID']))]
        variant_ids1 = [f'variant_{x+1}_2' for x in range(len(vcf['variants/ID']))]
        variant_ids = np.ravel([variant_ids0, variant_ids1], order = 'F')

        with open(output_file, 'w') as fh:
            variant_cols = '\t'.join(str(variant_id) for variant_id in variant_ids)
            fh.write(f'\t\t{variant_cols}\n')

            for i, sample_id in enumerate(sample_ids):
                if sample_id in sample_map.keys():
                    alleles_recoded0 = recode_nucleotides(to_nucleotides(gt[:, i][:, 0], ref_allele, alt_allele))
                    alleles_recoded1 = recode_nucleotides(to_nucleotides(gt[:, i][:, 1], ref_allele, alt_allele))
                    alleles_recoded = np.ravel([alleles_recoded0, alleles_recoded1], order = 'F')

                    alleles = '\t'.join(str(allele) for allele in alleles_recoded)
                    fh.write(f'{sample_id}\t{pops[i]}\t{alleles}\n')

                else:
                    continue

    else:        
        print('Warning: invalid argument specified.')
