pandasVCF
=========
VCF parser using the Python pandas library


Update: Nov 10 2014
pdVCFsingle.py can now parse a dataframe with a single individual, either from a multi-sample VCF or a single-sample VCF.  Missing genotype calls maked with '.' are dropped when add_variant_annotations are called.  100,000 variants are parsed in ~10sec.  


