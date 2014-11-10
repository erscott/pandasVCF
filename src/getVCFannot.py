def get_vcf_annotations(df, sample_name):
    df = get_allele_bases(df, sample_name)
    df = zygosity_fast(df)
    df['vartype1'] = map(vartype_mp, df[['REF','a1']].values)
    df['vartype2'] = map(vartype_mp, df[['REF','a2']].values)
    df['multiallele'] = df.ALT.map(lambda x: 1 if "," in x else 0)  #1 is more than 1 alt allele, 0 else
    return df



def get_genotype_counts(df):
    df['vartype1_2'] = ["_".join(i) for i in map(set, df[['vartype1','vartype2']].values)]   
    df_groups = df.groupby('zygosity')['vartype1_2'].value_counts()
    df_groups.index.set_names(['zygosity','vartype1_2'], inplace=True)
    df_groups.sortlevel(level=['zygosity','vartype1_2'])
    return df_groups

def get_allele_counts(df):
    df1 = df.groupby('zygosity')['vartype1'].value_counts()
    df2 = rtg_nist.groupby('zygosity')['vartype2'].value_counts
    df1_2 = df1 + df2
    return df1_2.unstack().sum()
