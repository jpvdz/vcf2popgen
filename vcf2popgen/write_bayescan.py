import sys
import allel as al
import numpy as np
 
def write_bayescan(input_file, sample_map_file, output_file, header = True):
    """Write bi-allelic variant calls to an output file in BAYESCAN format.
    
    Parameters
    ----------
    input_file
        A VCF file containing bi-allelic variant calls.
    sample_map_file
        A sample map file in CSV format. The first column should contain sample
        IDs, the second population IDs encoded as integers. 
    output_file
        The name of the output file.
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
    
    assert sorted(vcf['samples']) == sorted(sample_map), "Sample IDs in VCF and sample map don't match."
     
    gt = al.GenotypeArray(vcf['calldata/GT'])

    pops = np.array([
        sample_map[sample_id] for sample_id in vcf['samples']
    ])
    
    with open(output_file, 'w') as fh:
        fh.write(f"[loci]={len(gt)}\n\n")
        fh.write(f"[populations]={len(set(sample_map.values()))}\n\n")
        
        for i in set(sample_map.values()):
            fh.write(f"[pop]={i}\n")
            counts = gt.count_alleles(subpop = np.where((pops == i) == True)[0])
            for j, count in enumerate(counts):
                fh.write(f"{j+1} {count[0] + count[1]} 2 {count[0]} {count[1]}\n")
            fh.write("\n")


def main():
    input_file = sys.argv[1]
    sample_map_file = sys.argv[2]
    output_file = sys.argv[3]
    
    write_bayescan(input_file, sample_map_file, output_file)


if __name__ == "__main__":
    main()