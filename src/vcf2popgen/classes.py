import allel as al
import numpy as np
import numpy.typing as npt

class Variants:
    """A class that holds the chromosome, position, reference allele and alternate 
    allele of each variant. 
    
    Parameters
    ----------
    chromosome : npt.ArrayLike
        A NumPy array containing the chromosome of each variant.
    position : npt.ArrayLike
        A NumPy array containing the position of each variant.
    ref_alleles : npt.ArrayLike
        A NumPy array containing reference alleles encoded as nucleotides.
    alt_alleles : npt.ArrayLike
        A NumPy array containing alternate alleles encoded as nucleotides.
    """
    def __init__(self, chromosome : npt.ArrayLike, position : npt.ArrayLike, ref_alleles : npt.ArrayLike, alt_alleles : npt.ArrayLike):
        self.chromosome = chromosome
        self.position = position
        self.alleles = {'REF': ref_alleles,
                        'ALT': alt_alleles}
    
    
    def to_ndarray(self, allele_type, n_loci, n_samples, ploidy):
        """Convert a one-dimensional array containing variants to a three-dimensional array, 
        suitable for multiplying with a GenotypeArray. 
        
        Parameters
        ----------
        allele_type : str
            The allele type to convert. Valid options are 'REF', which selects the reference
            alleles, or 'ALT', which selects the alternate alleles. 
        n_loci : int
            Number of loci (first dimension).
        n_samples : int
            Number of samples (second dimension).
        ploidy : int
            Ploidy (third dimension).
            
        Returns
        -------
        variants_ndarray : npt.ArrayLike
            A three-dimensional matrix containing variants encoded as nucleotides.
        """
        variants_ndarray = np.repeat(self.alleles[allele_type], n_samples * ploidy, axis = 0).reshape(
                n_loci,
                n_samples,
                ploidy
            )
        return variants_ndarray
        

class PopGenData:
    """A class that holds bi-allelic SNPs from a set of samples from one or
    more populations.
    
    Parameters
    ----------
    samples : npt.ArrayLike
        A NumPy array containing sample IDs encoded as strings.
    populations : npt.ArrayLike
        A NumPy array containing population IDs encoded as integers.
    genotypes : al.GenotypeArray
        A GenotypeArray object.
    variants : Variants
        A Variants object.
        
    """
    def __init__(self, 
                 samples : npt.ArrayLike, 
                 populations : npt.ArrayLike, 
                 genotypes : al.GenotypeArray, 
                 variants : Variants):
        self.samples = samples
        self.populations = populations
        self.genotypes = genotypes
        self.variants = variants
  

    def n_loci(self):
        """Return the total number of loci."""
        return self.genotypes.shape[0]
    
    
    def n_samples(self):
        """Return the total number of samples."""
        return self.genotypes.shape[1]
    
    
    def ploidy(self):
        """Return the ploidy of a locus."""
        return self.genotypes.shape[2]

    
    def to_nucleotide_array(self):
        """Return the genotype data as an array of nucleotides encoded as characters."""
        ref_ndarray = self.variants.to_ndarray('REF', self.n_loci(), self.n_samples(), self.ploidy())
        alt_ndarray = self.variants.to_ndarray('ALT', self.n_loci(), self.n_samples(), self.ploidy())
        return (self.genotypes == 0) * ref_ndarray + (self.genotypes == 1) * alt_ndarray
    
    
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
        encoding = {'A': 1, 'T': 2, 'C': 3, 'G': 4, '' : missing}

        nuc_array = self.to_nucleotide_array()
        
        return (
            ((nuc_array == 'A') * encoding['A']) +
            ((nuc_array == 'T') * encoding['T']) +
            ((nuc_array == 'C') * encoding['C']) +
            ((nuc_array == 'G') * encoding['G']) +
            ((nuc_array == '') * encoding[''])
        )
        

    def to_bayescan(self, output_file : str) -> None:
        """Write bi-allelic SNPs to an output file in BAYESCAN format.
        
        Parameters
        ----------
        output_file : str
            Name of output file. 
        """
        print(f"Writing genotypic data in BAYESCAN format to: {output_file}")
        
        with open(output_file, 'w') as fh:
            fh.write(f"[loci]={len(self.genotypes)}\n\n")
            fh.write(f"[populations]={len(set(self.populations))}\n\n")
            
            for i in set(self.populations):
                fh.write(f"[pop]={i}\n")
                counts = self.genotypes.count_alleles(subpop = np.where((self.populations == i) == True)[0])
                for j, count in enumerate(counts):
                    fh.write(f"{j+1} {count[0] + count[1]} {self.ploidy()} {count[0]} {count[1]}\n")
                fh.write("\n")
                    
    
    def to_genepop(self, output_file : str) -> None:
        """Write bi-allelic SNPs to an output file in GENEPOP format.
    
        Parameters
        ----------
        output_file : str
            Name of output file. 
        """
        print(f"Writing genotypic data in GENEPOP format to: {output_file}")
        
        variant_ids = [f"snp{x+1}_{self.variants.chromosome[x]}_{self.variants.position[x]}" for x in range(self.n_loci())]
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
            print(f"Writing genotypic data in STRUCTURE format (one row per sample) to: {output_file}")

            variant_ids0 = [f"snp{x+1}a_{self.variants.chromosome[x]}_{self.variants.position[x]}" for x in range(self.n_loci())]
            variant_ids1 = [f"snp{x+1}b_{self.variants.chromosome[x]}_{self.variants.position[x]}" for x in range(self.n_loci())]
            
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
            print(f"Writing genotypic data in STRUCTURE format to: {output_file}")

            variant_ids = [f"snp{x+1}_{self.variants.chromosome[x]}_{self.variants.position[x]}" for x in range(self.n_loci())]
            variant_cols = '\t'.join(str(variant_id) for variant_id in variant_ids)
            with open(output_file, 'w') as fh:
                fh.write(f"\t\t{variant_cols}\n")

                for i, sample_id in enumerate(self.samples):
                    allele0 = '\t'.join(str(nucleotide) for nucleotide in recoded_nucs[:, i][:, 0])
                    allele1 = '\t'.join(str(nucleotide) for nucleotide in recoded_nucs[:, i][:, 1])
                    
                    fh.write(f"{sample_id}\t{self.populations[i]}\t{allele0}\n")
                    fh.write(f"{sample_id}\t{self.populations[i]}\t{allele1}\n")
