# VCF2PopGen
`VCF2PopGen` is a tool for converting bi-allelic single nucleotide polymorphism (SNP) data stored in VCF files to various population genetic data formats. Currently, output data can be generated in STRUCTURE, GENEPOP and BAYESCAN formats. VCF2PopGen is able to handle missing data, which it encodes according to the default way the desired output format handles missing data, e.g., encoding missing data as -9 in the case of formatting data for STRUCTURE. 

## Dependencies
Requires the [scikit-allel](https://scikit-allel.readthedocs.io/en/stable/) and [NumPy](https://numpy.org/doc/stable/) packages to be installed. You can easily install all dependencies from the requirements file by running `pip install -r requirements.txt`.

## How to use
The easiest way to use `VCF2PopGen` is to download the GitHub repository and put your own files under the `examples/` directory. Next, you can import `VCF2PopGen` into a Jupyter Notebook or within IPython console (be sure your current working directory is the repository), decide which output format you need and then run the appropriate function. For example, to generate a BAYESCAN-formatted output file using the example data:

```
import vcf2popgen
vcf2popgen.write_bayescan(input_file = 'examples/example_snps.vcf', sample_map_file = 'examples/example_snps.csv', output_file = 'examples/example_snps.bayescan')
```
