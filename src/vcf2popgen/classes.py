import allel as al
import numpy as np
import numpy.typing as npt

class Variants:
    """An object that holds a one-dimensional array of variants. 
    
    Parameters
    ----------
    variants
        A NumPy array containing variants encoded as nucleotides.
    """
    def __init__(self, variants : npt.ArrayLike):
        self.variants = variants
    
    
    def to_matrix(self, n_loci : int, n_samples : int, n_alleles : int) -> npt.ArrayLike:
        """Convert a one-dimensional array containing variants to a three-dimensional matrix, 
        suitable for multiplying with a GenotypeArray. 
        
        Parameters
        ----------
        n_loci : int
            Number of loci (first dimension).
        n_samples : int
            Number of samples (second dimension).
        n_alleles : int
            Ploidy (third dimension).
            
        Returns
        -------
        variants_matrix : npt.ArrayLike
            A three-dimensional matrix containing variants encoded as nucleotides.
        """
        variants_matrix = np.repeat(self.variants, n_samples * n_alleles, axis = 0).reshape(
            n_loci, n_samples, n_alleles)
        return variants_matrix
        

class PopGenData:
    """An object that holds bi-allelic SNPs from a set of samples from one or
    more populations.
    
    Parameters
    ----------
    samples : npt.ArrayLike
        A NumPy array containing sample IDs encoded as strings.
    populations : npt.ArrayLike
        A NumPy array containing population IDs encoded as integers.
    genotypes : al.GenotypeArray
        A GenotypeArray object.
    ref_variants : Variants
        A Variant object holding reference alleles.
    alt_variants : Variants
        A Variant object holding alternate alleles.
        
    """
    def __init__(self, 
                 samples : npt.ArrayLike, 
                 populations : npt.ArrayLike, 
                 genotypes : al.GenotypeArray, 
                 ref_variants : Variants, 
                 alt_variants : Variants):
        self.samples = samples
        self.populations = populations
        self.genotypes = genotypes
        self.ref_variants = ref_variants
        self.alt_variants = alt_variants
        

    def n_loci(self):
        """Return the total number of loci."""
        return self.genotypes.shape[0]
    
    
    def n_samples(self):
        """Return the total number of samples."""
        return self.genotypes.shape[1]
    
    
    def n_alleles(self):
        """Return the number of alleles per locus."""
        return self.genotypes.shape[2]

    
    def to_nucleotides(self):
        """Return the genotype data as a matrix of nucleotides encoded as characters."""
        ref_matrix = self.ref_variants.to_matrix(self.n_loci(), self.n_samples(), self.n_alleles())
        alt_matrix = self.alt_variants.to_matrix(self.n_loci(), self.n_samples(), self.n_alleles())
        return (self.genotypes == 0) * ref_matrix + (self.genotypes == 1) * alt_matrix
    
    
    def recode_nucleotides(self, missing : int = -9) -> npt.ArrayLike:
        """Recode nucleotides to integers.
        
        Parameters
        ----------
        missing : int
            Encoding for missing data. Defaults to -9.
        
        Returns
        -------
        recoded_nucleotides : npt.ArrayLike
            A matrix containing genotypic data with nucleotides encoded as integers.
        """
        encoding = {'A': 1, 'T': 2, 'C': 3, 'G': 4}
        
        if missing not in encoding:
            nucleotides = self.to_nucleotides()
            
            encoding[''] = missing
            
            recoded_nucleotides = (
                ((nucleotides == 'A') * encoding['A']) +
                ((nucleotides == 'T') * encoding['T']) +
                ((nucleotides == 'C') * encoding['C']) +
                ((nucleotides == 'G') * encoding['G']) +
                ((nucleotides == '') * encoding[''])
            )
            
            return recoded_nucleotides
        
        else:
            print("Invalid encoding for missing data.")

    def to_bayescan(self, output_file : str) -> None:
        """Write bi-allelic SNPs to an output file in BAYESCAN format.
        
        Parameters
        ----------
        output_file : str
            Name of output file. 
        """
        with open(output_file, 'w') as fh:
            fh.write(f"[loci]={len(self.genotypes)}\n\n")
            fh.write(f"[populations]={len(set(self.populations))}\n\n")
            
            for i in set(self.populations):
                fh.write(f"[pop]={i}\n")
                counts = self.genotypes.count_alleles(subpop = np.where((self.populations == i) == True)[0])
                for j, count in enumerate(counts):
                    fh.write(f"{j+1} {count[0] + count[1]} {self.n_alleles()} {count[0]} {count[1]}\n")
                fh.write("\n")
                    
    
    def to_genepop(self, output_file : str) -> None:
        """Write bi-allelic SNPs to an output file in GENEPOP format.
    
        Parameters
        ----------
        output_file : str
            Name of output file. 
        """
        variant_ids = [f"locus_{x+1}" for x in range(self.n_loci())]
        recoded_nucs = self.recode_nucleotides(missing=0)

        samples = []
        populations = []
        genotypes = []

        for i, sample_id in enumerate(self.samples):
            allele0 = recoded_nucs[:, i][:, 0].astype(str)
            allele1 = recoded_nucs[:, i][:, 1].astype(str)

            genotype = ' '.join(np.char.add(np.char.add(np.zeros(len(allele0), dtype = 'int8').astype(str), allele0), 
                np.char.add(np.zeros(len(allele1), dtype = 'int8').astype(str), allele1)))
            
            samples.append(sample_id)
            populations.append(self.populations[i])
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
        
        
    def to_structure(self, output_file : str, one_row_per_sample: bool = False) -> None:
        """Write bi-allelic SNPs to an output file in STRUCTURE format.
        
        Parameters
        ----------
        output_file : str
            Name of output file.
        one_row_per_sample : bool
            Whether samples should be encoded on a single row, which results in two
            columns per locus (one for each allele). Defaults to False.
        """
        recoded_nucs = self.recode_nucleotides(missing=-9)
        
        if one_row_per_sample == True:
            variant_ids0 = [f"locus_{x+1}_1" for x in range(self.n_loci())]
            variant_ids1 = [f"locus_{x+1}_2" for x in range(self.n_loci())]
            
            variant_ids = np.ravel([variant_ids0, variant_ids1], order = 'F')
            variant_cols = '\t'.join(str(variant_id) for variant_id in variant_ids)

            with open(output_file, 'w') as fh:
                fh.write(f"\t\t{variant_cols}\n")
                for i, sample_id in enumerate(self.samples):
                    allele0 = recoded_nucs[:, i][:, 0]
                    allele1 = recoded_nucs[:, i][:, 1]
                                    
                    genotype = '\t'.join(str(allele) for allele in np.ravel([allele0, allele1], 
                                                                            order = 'F'))
                    fh.write(f"{sample_id}\t{self.populations[i]}\t{genotype}\n")
            
        else:
            variant_ids = [f"locus_{x+1}" for x in range(self.n_loci())]
            variant_cols = '\t'.join(str(variant_id) for variant_id in variant_ids)
            with open(output_file, 'w') as fh:
                fh.write(f"\t\t{variant_cols}\n")

                for i, sample_id in enumerate(self.samples):
                    allele0 = '\t'.join(str(nucleotide) for nucleotide in recoded_nucs[:, i][:, 0])
                    allele1 = '\t'.join(str(nucleotide) for nucleotide in recoded_nucs[:, i][:, 1])
                    
                    fh.write(f"{sample_id}\t{self.populations[i]}\t{allele0}\n")
                    fh.write(f"{sample_id}\t{self.populations[i]}\t{allele1}\n")
