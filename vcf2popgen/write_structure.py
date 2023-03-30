import sys
import allel
import numpy as np
import pandas as pd

def get_nucleotide(variants, ref, alt):

    '''
    Determines the nucleotide of each variant.  
    '''

    # create boolean vectors encoding whether a variant is the reference or the alternate allele
    is_ref = ~np.array(variants, dtype=bool)
    is_alt = np.array(variants, dtype=bool)

    # determine nucleotides
    nucleotides = (ref * is_ref) + (alt * is_alt)

    # return nucleotides
    return(nucleotides)


def recode_nucleotides(nucleotides):
    
    '''
    Recodes nucleotides to integers. Missing data is encoded as -9.
    '''
    
    # create dictionary encoding nucleotides as integers
    nucleotide_codes = {'A': 1, 'T': 2, 'C': 3, 'G': 4, '' : -9}

    # create boolean vectors for nucleotides and missing data
    is_A = nucleotides == 'A'
    is_T = nucleotides == 'T'
    is_C = nucleotides == 'C'
    is_G = nucleotides == 'G'
    is_N = nucleotides == ''

    # recode nucleotides
    nucleotides_recoded = (
        (is_A * nucleotide_codes['A']) +
        (is_T * nucleotide_codes['T']) +
        (is_C * nucleotide_codes['C']) +
        (is_G * nucleotide_codes['G']) +
        (is_N * nucleotide_codes[''])
    )

    return(nucleotides_recoded)


def create_pop_dict(pop_map):

    ''' 
    Creates a population dictionary that maps unique population labels
    to a population ID encoded as integer.
    '''

    # create empty population dictionary
    pop_dict = {}

    # enumerate over population column
    for pop_id, pop_label in enumerate(pop_map['population'].unique()):

        # labels as keys and IDs as values
        pop_dict[pop_label] = pop_id + 1

    # return population dictionary
    return(pop_dict)


def write_structure(vcf_file, pop_map_file, str_file, one_row_per_sample = False):
    
    ''' 
    Writes a STRUCTURE file.
    '''

    # import data
    vcf = allel.read_vcf(vcf_file)
    pop_map = pd.read_csv(pop_map_file)

    # get relevant info
    sample_ids = vcf['samples']
    ref_allele = vcf['variants/REF']
    alt_allele = vcf['variants/ALT'].transpose()[0]

    # convert to genotype array
    gt = allel.GenotypeArray(vcf['calldata/GT'])

    # read sample IDs to be included from population map
    sample_selection = np.array(pop_map['sample_id'])

    # create population dictionary
    pop_dict = create_pop_dict(pop_map)

    # encode samples on two rows (with one column per locus)
    if one_row_per_sample == False:

        # set SNP IDs
        snp_ids = [f'snp_{x+1}' for x in range(len(vcf['variants/ID']))]
         
        # open STRUCTURE output file for writing
        with open(str_file, 'w') as fh:
            
            # join SNP IDs as string
            snp_cols = '\t'.join(str(snp_id) for snp_id in snp_ids)

            # write header line with SNP IDs
            fh.write(f'\t\t{snp_cols}\n')

            # enumerate over sample IDs
            for i, sample_id in enumerate(sample_ids):

                # check if sample ID should be included, else skip sample
                if sample_id in sample_selection:

                    # get population ID from population dictionary
                    pop_id = pop_dict[pop_map.loc[pop_map['sample_id'] == sample_id]['population'].values[0]]

                    # recode alleles
                    alleles_recoded0 = recode_nucleotides(get_nucleotide(gt[:, i][:, 0], ref_allele, alt_allele))
                    alleles_recoded1 = recode_nucleotides(get_nucleotide(gt[:, i][:, 1], ref_allele, alt_allele))

                    # convert to string
                    alleles0 = '\t'.join(str(allele) for allele in alleles_recoded0)
                    alleles1 = '\t'.join(str(allele) for allele in alleles_recoded1)

                    # write to file
                    fh.write(f'{sample_id}\t{pop_id}\t{alleles0}\n')
                    fh.write(f'{sample_id}\t{pop_id}\t{alleles1}\n')

                else:
                    continue
    
    # encode samples on one row (with two columns per locus)
    elif one_row_per_sample == True:

        # set SNP IDs
        snp_ids0 = [f'snp_{x+1}_1' for x in range(len(vcf['variants/ID']))]
        snp_ids1 = [f'snp_{x+1}_2' for x in range(len(vcf['variants/ID']))]
        snp_ids = np.ravel([snp_ids0, snp_ids1], order = 'F')

        # open STRUCTURE output file for writing
        with open(str_file, 'w') as fh:
            
            # join SNP IDs as string
            snp_cols = '\t'.join(str(snp_id) for snp_id in snp_ids)

            # write header line with SNP IDs
            fh.write(f'\t\t{snp_cols}\n')

            # enumerate over sample IDs
            for i, sample_id in enumerate(sample_ids):

                # check if sample ID should be included, else skip sample
                if sample_id in sample_selection:

                    # get population ID from population dictionary
                    pop_id = pop_dict[pop_map.loc[pop_map['sample_id'] == sample_id]['population'].values[0]]

                    # recode alleles
                    alleles_recoded0 = recode_nucleotides(get_nucleotide(gt[:, i][:, 0], ref_allele, alt_allele))
                    alleles_recoded1 = recode_nucleotides(get_nucleotide(gt[:, i][:, 1], ref_allele, alt_allele))
                    alleles_recoded = np.ravel([alleles_recoded0, alleles_recoded1], order = 'F')

                    # convert to string
                    alleles = '\t'.join(str(allele) for allele in alleles_recoded)

                    # write to file
                    fh.write(f'{sample_id}\t{pop_id}\t{alleles}\n')

                else:
                    continue

    else:        
        print('Warning: invalid argument specified.')


def main():

    # read command line args
    vcf_file = sys.argv[1]
    pop_map_file = sys.argv[2]
    str_file = sys.argv[3]
    one_row_per_sample = sys.argv[4]

    # convert vcf to structure
    write_structure(vcf_file, pop_map_file, str_file, one_row_per_sample)


if __name__ == "__main__":
    main()