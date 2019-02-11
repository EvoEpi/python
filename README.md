# python

A collection of Python scripts for handling various data files.

### Epigen[etics|omics]

`methTable.py`: More often than not you want to categorize loci based on their methylation enrichment. `methTable.py` will do just this. As input `methTable.py` requires a file generated from the [bedtools](https://bedtools.readthedocs.io/en/latest/) intersection of an allc file (a vcf-like file generated from [methylpy](https://github.com/yupenghe/methylpy)) and a bed file of loci chromosomal coordinates and identifier (i.e., four columns: chromosome, start, end, identifier):

```bash
#allc file requires the vcf header for bedtools recognize it
sed -i '1s/^/##fileformat=VCFv4.2\n/' <allc file>
#intersection
bedtools intersect -a <allc file> \
-b <bed file> \
-wa \
-wb \
> <intersection outfile>
```

Other required inputs are a taxonomic group for sequence contexts to parse (`animal`: CG and CH, or `plant`: CG, CHG, and CHH) and the name of an `outfile` to write to. A statistical enrichment test with False Discovery Rate (FDR) correction is performed for each locus following a method very similar to [Takuno and Gaut (2012)](https://www.ncbi.nlm.nih.gov/pubmed/21813466). I have removed filtering by number of sequence contexts (`n_[cg|ch|chg|chh]` column in `outfile`) present in the locus and sequencing coverage (`tor_[cg|ch|chg|chh]` column in `outfile`). I highly recommend you post-filter on both of these variables.
