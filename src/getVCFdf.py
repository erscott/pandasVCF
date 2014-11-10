import os
import pandas as pd

def get_vcf_df(file_path, compressed=''):
    
    def get_header(vcf_path):
        if vcf_path.endswith('.gz'):
            return os.popen('tabix -H ' + vcf_path).readlines()
        header_lines = os.popen('head -5000 ' + vcf_path).readlines()
        return [l for l in header_lines if l.startswith('#')]
    
    
    
    header = get_header(file_path)
    df = pd.read_table(file_path, sep="\t", compression=compressed, skiprows=(len(header)-1))
    df.set_index(['#CHROM', 'POS', 'REF', 'ALT'], inplace=True, drop=False)

    return df
    
    df['vartype2'] = df.apply(vartype, args=['a2'], axis=1)  #vartype for allele2: snp, mnp, insertion, deletion, sv
    return df
