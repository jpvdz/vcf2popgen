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

    # return recoded nucleotides
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


def write_genepop(vcf_file, pop_map_file, genepop_file):

    '''
    Writes a GENEPOP file. 
    '''

    # import data
    vcf = allel.read_vcf(vcf_file)
    pop_map = pd.read_csv(pop_map_file)

    # convert to genotype array
    gt = allel.GenotypeArray(vcf['calldata/GT'])

    # read sample IDs to be included from population map
    sample_ids = vcf['samples']
    ref_allele = vcf['variants/REF']
    alt_allele = vcf['variants/ALT'].transpose()[0]

    sample_selection = np.array(pop_map['sample_id'])

    # create population dictionary
    pop_dict = create_pop_dict(pop_map)

    # set SNP IDs
    snp_ids = [f'snp_{x+1}' for x in range(len(vcf['variants/ID']))]

    # create empty lists for sample IDs, population IDs and genotypes encoded in genepop format
    samples = []
    populations = []
    genotypes = []

    # loop over sample IDs
    for i, sample_id in enumerate(sample_ids):

        # check if sample ID should be included, else skip sample
        if sample_id in sample_selection:

            # get population ID from population dictionary
            pop_id = pop_dict[pop_map.loc[pop_map['sample_id'] == sample_id]['population'].values[0]]

            # get alleles and convert characters characters
            allele0 = np.array(recode_nucleotides(get_nucleotide(gt[:, i][:, 0], ref_allele, alt_allele))).astype(str)
            allele1 = np.array(recode_nucleotides(get_nucleotide(gt[:, i][:, 1], ref_allele, alt_allele))).astype(str)

            # combine alleles, add padding and convert to string
            genotype = ' '.join(np.char.add(np.char.add(np.zeros(len(allele0), dtype = 'int8').astype(str), allele0), 
                np.char.add(np.zeros(len(allele1), dtype = 'int8').astype(str), allele1)))

            # append samples, populations and genotypes
            samples.append(sample_id)
            populations.append(pop_id)
            genotypes.append(genotype)

        else:
            continue

    # combine sample IDs, population IDs and genotypes
    results = np.vstack((np.asarray(populations), np.asarray(samples), np.asarray(genotypes))).T

    # sort results on population ID
    sorted_results = results[results[:, 0].argsort()]

    # write genepop file
    with open(genepop_file, 'w') as fh:

        # write title line
        fh.write(f'Title line: GENEPOP file created from {vcf_file}\n')

        # write SNP IDs
        for snp_id in snp_ids:
            fh.write(f'{snp_id}\n')
            
        # loop over rows of sorted_results matrix
        for i in range(0, len(sorted_results)):

            pop_id = sorted_results[i][0]
            sample_id = sorted_results[i][1]
            genotype = sorted_results[i][2]

            # on first iteration write POP statement in genepop file
            if i == 0:
                fh.write(f'POP\n')
            
            # for subsequent iterations compare current population ID to previous one
            elif i > 0:
                prev_pop = sorted_results[i-1][0]

                # write POP statement if current and previous population IDs differ                
                if pop_id != prev_pop:
                    fh.write(f'POP\n')
            
            # write sample IDs and genotypes
            fh.write(f'{sample_id}, {genotype}\n')


def main():

    # read command line args
    vcf_file = sys.argv[1]
    pop_map_file = sys.argv[2]
    genepop_file = sys.argv[3]

    # convert vcf to structure
    write_genepop(vcf_file, pop_map_file, genepop_file)


if __name__ == "__main__":
    main()