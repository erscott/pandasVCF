pandasVCF
=========
VCF parser using the Python pandas library for interactive analysis

<h3>Update: April 13 2018</h3>
Updates to Python3 print statements and changed default pandas compression to 'infer'

<BR>


<h3>Update: March 9 2016</h3>
VCF header parsing requires the tabix -H command which is broken in tabix version 1.2.X  Please update your tabix version to 1.3. 

<BR>

<h3>Update: August 21 2015</h3>
pandasVCF handles both multi-sample and single-sample VCF files. Please see ipynb/ for usage. pandasVCFmulti and pandasVCFsingle are now depracated.  

<BR>

<h3>Update: February 12 2015</h3>
pandasVCFmulti now handles both multi-sample and single-sample VCF files. Please see http://nbviewer.ipython.org/github/erscott/pandasVCF/blob/master/ipynb/multi_sample_ex.ipynb for usage. pdVCFsingle.py is now depracated and will be removed in the near future.  

Command line support will be added in the near future. 

<h3>Update: Nov 10 2014</h3>
pdVCFsingle.py can now parse a dataframe with a single individual, either from a multi-sample VCF or a single-sample VCF.  Missing genotype calls maked with '.' are dropped when add_variant_annotations are called.  100,000 variants are parsed in ~10sec.  


