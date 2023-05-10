import allel as al
import numpy as np

class Variants:
    def __init__(self, variants):
        self.variants = variants
    
    
    def to_ndarray(self, n_loci, n_samples, ploidy):
        return np.repeat(
            self.variants, 
            n_samples * ploidy,
            axis = 0).reshape(
                n_loci,
                n_samples,
                ploidy
            )
        

class PopGenData:
    def __init__(self, samples, populations, genotypes : al.GenotypeArray, ref_variants : Variants, alt_variants : Variants):
        self.samples = samples
        self.populations = populations
        self.genotypes = genotypes
        self.ref_variants = ref_variants
        self.alt_variants = alt_variants
        

    def n_loci(self):
        return self.genotypes.shape[0]
    
    
    def n_samples(self):
        return self.genotypes.shape[1]
    
    
    def ploidy(self):
        return self.genotypes.shape[2]

    
    def to_nucleotide_array(self):
        ref_ndarray = self.ref_variants.to_ndarray(self.n_loci(), self.n_samples(), self.ploidy())
        alt_ndarray = self.alt_variants.to_ndarray(self.n_loci(), self.n_samples(), self.ploidy())
        return (self.genotypes == 0) * ref_ndarray + (self.genotypes == 1) * alt_ndarray
    
    
    def recode_nucleotides(self, missing = -9):
        
        encoding = {'A': 1, 'T': 2, 'C': 3, 'G': 4, '' : missing}

        nuc_array = self.to_nucleotide_array()
        
        return (
            ((nuc_array == 'A') * encoding['A']) +
            ((nuc_array == 'T') * encoding['T']) +
            ((nuc_array == 'C') * encoding['C']) +
            ((nuc_array == 'G') * encoding['G']) +
            ((nuc_array == '') * encoding[''])
        )


    def test_to_bayescan(self, output_file):
        print(f"Writing genotypic data in BAYESCAN format to: {output_file}")
        print(f"[loci]={len(self.genotypes)}\n")
        print(f"[populations]={len(set(self.populations))}\n")
        
        for i in set(self.populations):
            print(f"[pop]={i}")
            counts = self.genotypes.count_alleles(subpop = np.where((self.populations == i) == True)[0])
            for j, count in enumerate(counts):
                print(f"{j+1} {count[0] + count[1]} {self.ploidy()} {count[0]} {count[1]}")
            print("")


    def test_to_genepop(self, output_file):
        print(f"Writing genotypic data in GENEPOP format to: {output_file}")
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
        
        for variant_id in variant_ids:
            print(f"{variant_id}")
            
        for i in range(0, len(sorted_results)):
            pop_id = sorted_results[i][0]
            sample_id = sorted_results[i][1]
            genotype = sorted_results[i][2]

            if i == 0:
                print(f"POP")
            
            elif i > 0:
                prev_pop = sorted_results[i-1][0]
                if pop_id != prev_pop:
                    print(f"POP")
            
            print(f"{sample_id}, {genotype}")
    
    
    def test_to_structure(self, output_file, one_row_per_sample = False):
        recoded_nucs = self.recode_nucleotides(missing=-9)
        
        if one_row_per_sample == True:
            print(f"Writing genotypic data in STRUCTURE format (one row per sample) to: {output_file}")
            variant_ids0 = [f"locus_{x+1}_1" for x in range(self.n_loci())]
            variant_ids1 = [f"locus_{x+1}_2" for x in range(self.n_loci())]
            
            variant_ids = np.ravel([variant_ids0, variant_ids1], order = 'F')
               
            variant_cols = '\t'.join(str(variant_id) for variant_id in variant_ids)
            print(f"\t\t{variant_cols}")

            for i, sample_id in enumerate(self.samples):
                allele0 = recoded_nucs[:, i][:, 0]
                allele1 = recoded_nucs[:, i][:, 1]
                                
                genotype = '\t'.join(str(allele) for allele in np.ravel([allele0, allele1], 
                                                                        order = 'F'))
                print(f"{sample_id}\t{self.populations[i]}\t{genotype}")
            
        else:
            print(f"Writing genotypic data in STRUCTURE format to: {output_file}")
            variant_ids = [f"locus_{x+1}" for x in range(self.n_loci())]
            
            variant_cols = '\t'.join(str(variant_id) for variant_id in variant_ids)
            print(f"\t\t{variant_cols}")

            for i, sample_id in enumerate(self.samples):
                allele0 = '\t'.join(str(nucleotide) for nucleotide in recoded_nucs[:, i][:, 0])
                allele1 = '\t'.join(str(nucleotide) for nucleotide in recoded_nucs[:, i][:, 1])
                
                print(f"{sample_id}\t{self.populations[i]}\t{allele0}")
                print(f"{sample_id}\t{self.populations[i]}\t{allele1}")


def create_test_data():
    
    test_samples = np.array(['sample1', 'sample2', 'sample3'], dtype=object)
    test_pops = np.array([1, 2, 3], dtype = int)

    test_variants = np.array(
        [[[1,0], [1,1], [0,0]], 
        [[0,0], [1,0], [0,0]],
        [[0,1], [0,0], [1,1]],
        [[0,0], [1,1], [0,0]]], dtype = int
    )
    
    test_genotypes = al.GenotypeArray(test_variants)
    
    test_refs = Variants(np.array(['A', 'C', 'G', 'T'], dtype = object))
    test_alts = Variants(np.array(['T', 'G', 'A', 'C'], dtype = object))
    
    test_data = PopGenData(samples=test_samples,
                      populations=test_pops,
                      genotypes=test_genotypes,
                      ref_variants=test_refs,
                      alt_variants=test_alts)

    assert test_data.n_samples() == 3
    assert test_data.n_loci() == 4
    assert test_data.ploidy() == 2
    
    return test_data


def test_nucleotide_array():
    test_data = create_test_data()
    test_nucs = np.array(
        [[['T', 'A'], 
          ['T', 'T'], 
          ['A', 'A']],
         [['C', 'C'], 
          ['G', 'C'],
          ['C', 'C']],
         [['G', 'A'], 
          ['G', 'G'], 
          ['A', 'A']],
         [['T', 'T'], 
          ['C', 'C'], 
          ['T', 'T']]], dtype = object
    )
    assert (test_data.to_nucleotide_array() == test_nucs).all()


def test_recode_nucleotides():
    test_data = create_test_data()
    test_recoded_nucs = np.array(
        [[[2, 1], 
          [2, 2], 
          [1, 1]],
         [[3, 3], 
          [4, 3],
          [3, 3]],
         [[4, 1], 
          [4, 4], 
          [1, 1]],
         [[2, 2], 
          [3, 3], 
          [2, 2]]], dtype = int
    )
    assert (test_data.recode_nucleotides() == test_recoded_nucs).all()
    

def test_to_bayescan():
    test_data = create_test_data()
    test_data.test_to_bayescan(output_file="testout.bayescan")
    assert True
    

def test_to_genepop():
    test_data = create_test_data()
    test_data.test_to_genepop(output_file="testout.genepop")
    assert True
    

def test_to_structure():
    test_data = create_test_data()
    test_data.test_to_structure(output_file="testout.str")
    test_data.test_to_structure(output_file="testout.str", one_row_per_sample=True)
    assert True