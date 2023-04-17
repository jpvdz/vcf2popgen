import numpy as np
import allel as al
from .utils import to_nucleotides, recode_nucleotides

def bayescan(data : dict, output_file : str) -> None:
    """Write bi-allelic variant calls to BAYESCAN format.
    
    Parameters
    ----------
    data : dict
        Dictionary containing arrays holding VCF data and population IDs.
    output_file : str
        Name or path of output file. 
    """
    gt = al.GenotypeArray(data['calldata/GT'])
    
    with open(output_file, 'w') as fh:
        fh.write(f"[loci]={len(gt)}\n\n")
        fh.write(f"[populations]={len(set(data['populations']))}\n\n")
        
        for i in set(data['populations']):
            fh.write(f"[pop]={i}\n")
            counts = gt.count_alleles(subpop = np.where((data['populations'] == i) == True)[0])
            for j, count in enumerate(counts):
                fh.write(f"{j+1} {count[0] + count[1]} 2 {count[0]} {count[1]}\n")
            fh.write("\n")


def genepop(data : dict, output_file : str) -> None:
    """Write bi-allelic variant calls to GENEPOP format.
    
    Parameters
    ----------
    data : dict
        Dictionary containing arrays holding VCF data and population IDs.
    output_file : str
        Name or path of output file. 
    """
    gt = al.GenotypeArray(data['calldata/GT'])
    ref_allele = data['variants/REF']
    alt_allele = data['variants/ALT']
    
    variant_ids = [f"variant_{x+1}" for x in range(len(data['variants/ID']))]

    samples = []
    populations = []
    genotypes = []

    for i, sample_id in enumerate(data['samples']):
        allele0 = np.array(recode_nucleotides(to_nucleotides(gt[:, i][:, 0], ref_allele, alt_allele), miss_encode = 0)).astype(str)
        allele1 = np.array(recode_nucleotides(to_nucleotides(gt[:, i][:, 1], ref_allele, alt_allele), miss_encode = 0)).astype(str)

        genotype = ' '.join(np.char.add(np.char.add(np.zeros(len(allele0), dtype = 'int8').astype(str), allele0), 
            np.char.add(np.zeros(len(allele1), dtype = 'int8').astype(str), allele1)))

        samples.append(sample_id)
        populations.append(data['populations'][i])
        genotypes.append(genotype)

    results = np.vstack((np.asarray(populations), np.asarray(samples), np.asarray(genotypes))).T
    sorted_results = results[results[:, 0].argsort()]
    
    with open(output_file, 'w') as fh:
        fh.write(f"Title line: GENEPOP file created with VCF2PopGen\n")

        for variant_id in variant_ids:
            fh.write(f'{variant_id}\n')
            
        for i in range(0, len(sorted_results)):
            pop_id = sorted_results[i][0]
            sample_id = sorted_results[i][1]
            genotype = sorted_results[i][2]

            if i == 0:
                fh.write(f"POP\n")
            
            elif i > 0:
                prev_pop = sorted_results[i-1][0]
                if pop_id != prev_pop:
                    fh.write(f"POP\n")
            
            fh.write(f"{sample_id}, {genotype}\n")


def structure(data : dict, output_file : str, one_row_per_sample : bool = False) -> None:
    """Write bi-allelic variant calls to STRUCTURE format.
    
    Parameters
    ----------
    data : dict
        Dictionary containing arrays holding VCF data and population IDs.
    output_file : str
        Name or path of output file.
    one_row_per_sample : bool
        Whether samples should be encoded on a single row, which results in two
        columns per locus (one for each allele). Defaults to False.
    """        
    gt = al.GenotypeArray(data['calldata/GT'])
    ref_allele = data['variants/REF']
    alt_allele = data['variants/ALT']
    
    if one_row_per_sample == True:
        variant_ids0 = [f'variant_{x+1}_1' for x in range(len(data['variants/ID']))]
        variant_ids1 = [f'variant_{x+1}_2' for x in range(len(data['variants/ID']))]
        variant_ids = np.ravel([variant_ids0, variant_ids1], order = 'F')
        
        with open(output_file, 'w') as fh:
            variant_cols = '\t'.join(str(variant_id) for variant_id in variant_ids)
            fh.write(f'\t\t{variant_cols}\n')

            for i, sample_id in enumerate(data['samples']):
                alleles_recoded0 = recode_nucleotides(to_nucleotides(gt[:, i][:, 0], ref_allele, alt_allele))
                alleles_recoded1 = recode_nucleotides(to_nucleotides(gt[:, i][:, 1], ref_allele, alt_allele))
                alleles_recoded = np.ravel([alleles_recoded0, alleles_recoded1], order = 'F')

                pop = data['populations'][i]
                
                alleles = '\t'.join(str(allele) for allele in alleles_recoded)
                fh.write(f'{sample_id}\t{pop}\t{alleles}\n')
        
    else:
        variant_ids = [f"variant_{x+1}" for x in range(len(data['variants/ID']))]
        
        with open(output_file, 'w') as fh:
            variant_cols = '\t'.join(str(variant_id) for variant_id in variant_ids)
            fh.write(f'\t\t{variant_cols}\n')

            for i, sample_id in enumerate(data['samples']):
                alleles_recoded0 = recode_nucleotides(to_nucleotides(gt[:, i][:, 0], ref_allele, alt_allele))
                alleles_recoded1 = recode_nucleotides(to_nucleotides(gt[:, i][:, 1], ref_allele, alt_allele))

                alleles0 = '\t'.join(str(allele) for allele in alleles_recoded0)
                alleles1 = '\t'.join(str(allele) for allele in alleles_recoded1)

                pop = data['populations'][i]
                
                fh.write(f'{sample_id}\t{pop}\t{alleles0}\n')
                fh.write(f'{sample_id}\t{pop}\t{alleles1}\n')
