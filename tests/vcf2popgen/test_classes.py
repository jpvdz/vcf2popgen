import allel as al
import numpy as np

class PopGenData:
    def __init__(self, samples, populations, genotypes, ref_allele, alt_allele):
        self.samples = samples
        self.populations = populations
        self.genotypes = genotypes
        
        self.alleles = {
            'REF' : ref_allele,
            'ALT' : alt_allele
        }
        
        self.n_loci = self.genotypes.shape[0]
        self.n_samples = self.genotypes.shape[1]
        self.ploidy = self.genotypes.shape[2]
        
        self.nucleotide_array = self._to_nucleotide_array()
        self.nucleotide_recoded = self._recode_nucleotides()
        

    def summary(self):
        print(f"Number of loci: {self.n_loci}")
        print(f"Number of samples: {self.n_samples}")
        print(f"Ploidy: {self.ploidy}")
    
    
    def _resize_alleles(self, allele_type):
        return np.repeat(self.alleles[allele_type], self.n_samples * self.ploidy, axis = 0).reshape(self.n_loci, self.n_samples, self.ploidy)
    
    
    def _to_nucleotide_array(self):
        ref_array, alt_array = map(self._resize_alleles, ('REF', 'ALT'))
        return (self.genotypes == 0) * ref_array + (self.genotypes == 1) * alt_array
    
    
    def _recode_nucleotides(self, missing = -9):
        
        encoding = {'A': 1, 'T': 2, 'C': 3, 'G': 4, '' : missing}

        return (
            ((self.nucleotide_array == 'A') * encoding['A']) +
            ((self.nucleotide_array == 'T') * encoding['T']) +
            ((self.nucleotide_array == 'C') * encoding['C']) +
            ((self.nucleotide_array == 'G') * encoding['G']) +
            ((self.nucleotide_array == '') * encoding[''])
        )


    def test_to_bayescan(self):
            print(f"[loci]={len(self.genotypes)}\n")
            print(f"[populations]={len(set(self.populations))}\n")
            
            for i in set(self.populations):
                print(f"[pop]={i}")
                counts = self.genotypes.count_alleles(subpop = np.where((self.populations == i) == True)[0])
                for j, count in enumerate(counts):
                    print(f"{j+1} {count[0] + count[1]} {self.ploidy} {count[0]} {count[1]}")
                print("")


    def test_to_genepop(self):
        variant_ids = [f"locus_{x+1}" for x in range(self.n_loci)]

        samples = []
        populations = []
        genotypes = []

        for i, sample_id in enumerate(self.samples):
            allele0 = self.nucleotide_recoded[:, i][:, 0].astype(str)
            allele1 = self.nucleotide_recoded[:, i][:, 1].astype(str)

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
    
    
    def test_to_structure(self, one_row_per_sample = False):
                
        if one_row_per_sample == True:
            variant_ids0 = [f"locus_{x+1}_1" for x in range(self.n_loci)]
            variant_ids1 = [f"locus_{x+1}_2" for x in range(self.n_loci)]
            
            variant_ids = np.ravel([variant_ids0, variant_ids1], order = 'F')
               
            variant_cols = '\t'.join(str(variant_id) for variant_id in variant_ids)
            print(f"\t\t{variant_cols}")

            for i, sample_id in enumerate(self.samples):
                allele0 = self.nucleotide_recoded[:, i][:, 0]
                allele1 = self.nucleotide_recoded[:, i][:, 1]
                                
                genotype = '\t'.join(str(allele) for allele in np.ravel([allele0, allele1], 
                                                                        order = 'F'))
                print(f"{sample_id}\t{self.populations[i]}\t{genotype}")
            
        else:
            variant_ids = [f"locus_{x+1}" for x in range(self.n_loci)]
            
            variant_cols = '\t'.join(str(variant_id) for variant_id in variant_ids)
            print(f"\t\t{variant_cols}")

            for i, sample_id in enumerate(self.samples):
                allele0 = '\t'.join(str(nucleotide) for nucleotide in self.nucleotide_recoded[:, i][:, 0])
                allele1 = '\t'.join(str(nucleotide) for nucleotide in self.nucleotide_recoded[:, i][:, 1])
                
                print(f"{sample_id}\t{self.populations[i]}\t{allele0}")
                print(f"{sample_id}\t{self.populations[i]}\t{allele1}")


    # def to_bayescan(self, output_file):
    #     """Write bi-allelic variant calls to BAYESCAN format.
        
    #     Parameters
    #     ----------
    #     output_file : str
    #         Name or path of output file. 
    #     """
    #     with open(output_file, 'w') as fh:
    #         fh.write(f"[loci]={len(self.genotypes)}\n\n")
    #         fh.write(f"[populations]={len(set(self.populations))}\n\n")
            
    #         for i in set(self.populations):
    #             fh.write(f"[pop]={i}\n")
    #             counts = self.genotypes.count_alleles(subpop = np.where((self.populations == i) == True)[0])
    #             for j, count in enumerate(counts):
    #                 fh.write(f"{j+1} {count[0] + count[1]} {self.ploidy} {count[0]} {count[1]}\n")
    #             fh.write("\n")
                    
    
    # def to_genepop(self, output_file):
    #     """Write bi-allelic variant calls to GENEPOP format.
    
    #     Parameters
    #     ----------
    #     data : dict
    #         Dictionary containing arrays holding VCF data and population IDs.
    #     output_file : str
    #         Name or path of output file. 
    #     """
    #     gt = self.genotypes
    #     ref_allele = self.ref_allele
    #     alt_allele = self.alt_allele
        
    #     variant_ids = [f"variant_{x+1}" for x in range(len(self.ref_allele))] # TODO method that returns number of loci, samples, etc

    #     samples = []
    #     populations = []
    #     genotypes = []

    #     # TODO fix this section so it works with class attributes and methods
    #     for i, sample_id in enumerate(self.samples):
    #         allele0 = np.array(self._recode_nucleotides(gt[:, i][:, 0], missing = 0)).astype(str)
    #         allele1 = np.array(self._recode_nucleotides(gt[:, i][:, 1], missing = 0)).astype(str)

    #         genotype = ' '.join(np.char.add(np.char.add(np.zeros(len(allele0), dtype = 'int8').astype(str), allele0), 
    #             np.char.add(np.zeros(len(allele1), dtype = 'int8').astype(str), allele1)))

    #         samples.append(sample_id)
    #         populations.append(data['populations'][i])
    #         genotypes.append(genotype)

    #     results = np.vstack((np.asarray(populations), np.asarray(samples), np.asarray(genotypes))).T
    #     sorted_results = results[results[:, 0].argsort()]
        
    #     # with open(output_file, 'w') as fh:
    #         # fh.write(f"Title line: GENEPOP file created with VCF2PopGen\n")

    #     for variant_id in variant_ids:
    #         print(f"{variant_id}\n")
    #         # fh.write(f'{variant_id}\n')
            
    #     for i in range(0, len(sorted_results)):
    #         pop_id = sorted_results[i][0]
    #         sample_id = sorted_results[i][1]
    #         genotype = sorted_results[i][2]

    #         if i == 0:
    #             # fh.write(f"POP\n")
    #             print(f"POP\n")
            
    #         elif i > 0:
    #             prev_pop = sorted_results[i-1][0]
    #             if pop_id != prev_pop:
    #                 # fh.write(f"POP\n")
    #                 print(f"POP\n")
            
    #         # fh.write(f"{sample_id}, {genotype}\n")
    #         print(f"{sample_id}, {genotype}\n")
        
        
    # def to_structure(self, output_file, one_row_per_sample):
    #     """Write bi-allelic variant calls to STRUCTURE format.
        
    #     Parameters
    #     ----------
    #     output_file : str
    #         Name or path of output file.
    #     one_row_per_sample : bool
    #         Whether samples should be encoded on a single row, which results in two
    #         columns per locus (one for each allele). Defaults to False.
    #     """        
    #     gt = self.genotypes
    #     ref_allele = self.ref_allele
    #     alt_allele = self.alt_allele
        
    #     if one_row_per_sample == True:
    #         variant_ids0 = [f'variant_{x+1}_1' for x in range(len(self.ref_allele))] # TODO use loci counting method here
    #         variant_ids1 = [f'variant_{x+1}_2' for x in range(len(self.ref_allele))] # TODO use loci counting method here
    #         variant_ids = np.ravel([variant_ids0, variant_ids1], order = 'F')
            
    #         with open(output_file, 'w') as fh:
    #             variant_cols = '\t'.join(str(variant_id) for variant_id in variant_ids)
    #             fh.write(f'\t\t{variant_cols}\n')

    #             for i, sample_id in enumerate(self.samples):
    #                 alleles_recoded0 = recode_nucleotides(to_nucleotides(gt[:, i][:, 0], ref_allele, alt_allele))
    #                 alleles_recoded1 = recode_nucleotides(to_nucleotides(gt[:, i][:, 1], ref_allele, alt_allele))
    #                 alleles_recoded = np.ravel([alleles_recoded0, alleles_recoded1], order = 'F')

    #                 pop = self.populations[i]
                    
    #                 alleles = '\t'.join(str(allele) for allele in alleles_recoded)
    #                 fh.write(f'{sample_id}\t{pop}\t{alleles}\n')
            
    #     else:
    #         variant_ids = [f"variant_{x+1}" for x in range(len(self.ref_allele))]
            
    #         with open(output_file, 'w') as fh:
    #             variant_cols = '\t'.join(str(variant_id) for variant_id in variant_ids)
    #             fh.write(f'\t\t{variant_cols}\n')

    #             for i, sample_id in enumerate(self.samples):
    #                 alleles_recoded0 = recode_nucleotides(to_nucleotides(gt[:, i][:, 0], ref_allele, alt_allele))
    #                 alleles_recoded1 = recode_nucleotides(to_nucleotides(gt[:, i][:, 1], ref_allele, alt_allele))

    #                 alleles0 = '\t'.join(str(allele) for allele in alleles_recoded0)
    #                 alleles1 = '\t'.join(str(allele) for allele in alleles_recoded1)

    #                 pop = self.populations[i]
                    
    #                 fh.write(f'{sample_id}\t{pop}\t{alleles0}\n')
    #                 fh.write(f'{sample_id}\t{pop}\t{alleles1}\n')


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
    
    test_refs = np.array(['A', 'C', 'G', 'T'], dtype = object)
    test_alts = np.array(['T', 'G', 'A', 'C'], dtype = object)
    
    test_data = PopGenData(samples=test_samples,
                      populations=test_pops,
                      genotypes=test_genotypes,
                      ref_allele=test_refs,
                      alt_allele=test_alts)

    assert test_data.n_samples == 3
    assert test_data.n_loci == 4
    assert test_data.ploidy == 2
    
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
    assert (test_data.nucleotide_array == test_nucs).all()


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
    assert (test_data.nucleotide_recoded == test_recoded_nucs).all()
    

def test_to_bayescan():
    test_data = create_test_data()
    test_data.test_to_bayescan()
    assert True
    

def test_to_genepop():
    test_data = create_test_data()
    test_data.test_to_genepop()
    assert True
    

def test_to_structure():
    test_data = create_test_data()
    test_data.test_to_structure()
    test_data.test_to_structure(one_row_per_sample=True)
    assert True