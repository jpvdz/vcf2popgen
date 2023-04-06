import sys
import allel as al
import numpy as np
import pandas as pd
from vcf2popgen.utils import to_nucleotides, recode_nucleotides

def write_genepop(input_file, sample_map_file, output_file, header = True):
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

    variant_ids = [f"variant_{x+1}" for x in range(len(vcf['variants/ID']))]

    samples = []
    populations = []
    genotypes = []

    for i, sample_id in enumerate(sample_ids):
        if sample_id in sample_map.keys():
            allele0 = np.array(recode_nucleotides(to_nucleotides(gt[:, i][:, 0], ref_allele, alt_allele), miss_encode = 0)).astype(str)
            allele1 = np.array(recode_nucleotides(to_nucleotides(gt[:, i][:, 1], ref_allele, alt_allele), miss_encode = 0)).astype(str)

            genotype = ' '.join(np.char.add(np.char.add(np.zeros(len(allele0), dtype = 'int8').astype(str), allele0), 
                np.char.add(np.zeros(len(allele1), dtype = 'int8').astype(str), allele1)))

            samples.append(sample_id)
            populations.append(pops[i])
            genotypes.append(genotype)

        else:
            continue

    results = np.vstack((np.asarray(populations), np.asarray(samples), np.asarray(genotypes))).T
    sorted_results = results[results[:, 0].argsort()]

    with open(output_file, 'w') as fh:
        fh.write(f"Title line: GENEPOP file created from {input_file}\n")

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


def main():
    input_file = sys.argv[1]
    sample_map_file = sys.argv[2]
    output_file = sys.argv[3]

    write_genepop(input_file, sample_map_file, output_file)


if __name__ == "__main__":
    main()