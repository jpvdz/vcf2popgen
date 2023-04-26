# VCF2PopGen
`VCF2PopGen` is a simple Python package for converting bi-allelic single nucleotide polymorphism (SNP) data stored in VCF files to various population genetic data formats. Currently, output data can be generated in STRUCTURE, GENEPOP or BAYESCAN format. `VCF2PopGen` is able to handle missing data, which it encodes according to the default way the desired output format handles missing data, e.g., encoding missing data as -9 in the case of formatting data for STRUCTURE.

## Dependencies
Requires [scikit-allel](https://scikit-allel.readthedocs.io/en/stable/) and [NumPy](https://numpy.org/doc/stable/), which should be installed automatically alongside `VCF2PopGen` if you follow the installation instructions below.

## How to install
Simply download the latest release and install either the source distribution (`.tar.gz`) or built distribution (`.whl`) using `pip`:
```
pip install /path/to/dir/vcf2popgen-0.1.0.tar.gz
```

## How to use
Import the a VCF containing bi-allelic SNPs and a sample map using the `vcf2popgen.read()` function, which returns a `PopGenData` object:

```
import vcf2popgen
data = vcf2popgen.read('examples/example_snps.vcf', 'examples/example_sample_map.csv')
```

Next, simply use one of the available methods associated with `PopGenData` objects to write output in the desired format. For example, to write SNP data to an output file in BAYESCAN format, use the `.to_bayescan()` method:
```
data.to_bayescan('examples/example_snps.bayescan')
```
Or, for example, generate an output file in STRUCTURE format but with samples encoded on a single row (and thus two columns per locus):
```
data.to_structure('examples/example_snps.str', one_row_per_sample=True)
```

## Bugs
Please report bugs in the [issue tracker](https://github.com/jpvdz/vcf2popgen/issues). 

