import sys
import allel as al
import numpy as np
import pandas as pd

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
    is_ref = ~np.array(variants, dtype=bool)
    is_alt = np.array(variants, dtype=bool)

    nucs = (ref * is_ref) + (alt * is_alt)

    return nucs


def recode_nucleotides(nucs):
    """Recodes nucleotides as integers. Encodes missing data as -9.
    
    Parameters
    ----------
    nucs
        NumPy array of nucleotides encoded as strings.
    
    Returns
    -------
    recoded_nucs
        NumPy array of nucleotides encoded as integers.
    """
    nuc_codes = {'A': 1, 'T': 2, 'C': 3, 'G': 4, '' : -9}

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


def create_pop_dict(sample_map):
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


def write_structure(input_file, sample_map_file, output_file, one_row_per_sample = False):
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
    """

    vcf = al.read_vcf(input_file)
    sample_map = pd.read_csv(sample_map_file)
    
    gt = al.GenotypeArray(vcf['calldata/GT'])

    sample_ids = vcf['samples']
    ref_allele = vcf['variants/REF']
    alt_allele = vcf['variants/ALT'].transpose()[0]

    sample_selection = np.array(sample_map['sample_id'])

    pop_dict = create_pop_dict(sample_map)

    if one_row_per_sample == False:
        variant_ids = [f"variant_{x+1}" for x in range(len(vcf['variants/ID']))]
         
        with open(output_file, 'w') as fh:
            variant_cols = '\t'.join(str(variant_id) for variant_id in variant_ids)
            fh.write(f'\t\t{variant_cols}\n')

            for i, sample_id in enumerate(sample_ids):
                if sample_id in sample_selection:
                    pop_id = pop_dict[sample_map.loc[sample_map['sample_id'] == sample_id]['pop_id'].values[0]]

                    alleles_recoded0 = recode_nucleotides(to_nucleotides(gt[:, i][:, 0], ref_allele, alt_allele))
                    alleles_recoded1 = recode_nucleotides(to_nucleotides(gt[:, i][:, 1], ref_allele, alt_allele))

                    alleles0 = '\t'.join(str(allele) for allele in alleles_recoded0)
                    alleles1 = '\t'.join(str(allele) for allele in alleles_recoded1)

                    fh.write(f'{sample_id}\t{pop_id}\t{alleles0}\n')
                    fh.write(f'{sample_id}\t{pop_id}\t{alleles1}\n')

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
                if sample_id in sample_selection:
                    pop_id = pop_dict[sample_map.loc[sample_map['sample_id'] == sample_id]['pop_id'].values[0]]

                    alleles_recoded0 = recode_nucleotides(to_nucleotides(gt[:, i][:, 0], ref_allele, alt_allele))
                    alleles_recoded1 = recode_nucleotides(to_nucleotides(gt[:, i][:, 1], ref_allele, alt_allele))
                    alleles_recoded = np.ravel([alleles_recoded0, alleles_recoded1], order = 'F')

                    alleles = '\t'.join(str(allele) for allele in alleles_recoded)
                    fh.write(f'{sample_id}\t{pop_id}\t{alleles}\n')

                else:
                    continue

    else:        
        print('Warning: invalid argument specified.')


def main():
    input_file = sys.argv[1]
    sample_map_file = sys.argv[2]
    output_file = sys.argv[3]
    one_row_per_sample = sys.argv[4]

    write_structure(input_file, sample_map_file, output_file, one_row_per_sample)


if __name__ == "__main__":
    main()