def zygosity_fast(df):
    '''
    This function quickly assigns zygosity states to 
    each variant using set logic
    
    '''
    
    df_hom_ref = df[ (df['a1'] == df['REF']) & (df['a2'] == df['REF'])]
    df_hom_ref['zygosity'] = 'hom-ref'
    
    df_hom_miss = df[ (df['a1'] == '.') & (df['a2'] == '.')]
    df_hom_miss['zygosity'] = 'hom-miss'
    
    df_het_miss = df[ (df['a1'] == '.') | (df['a2'] == '.') ]
    df_het_miss['zygosity'] = 'het-miss'
    
    
    df_not_miss = df[~df.index.isin( set(df_hom_miss.index) | set(df_het_miss.index))  ]
    
    
    df_het_alt = df_not_miss[ ((df_not_miss['a1'] != df_not_miss['REF']) & (df_not_miss['a2'] != df_not_miss['REF'])) & (df_not_miss['a1'] != df_not_miss['a2']) ]
    df_het_alt['zygosity'] = 'het-alt'
    
    df_hom_alt = df_not_miss[ (((df_not_miss['a1'] != df_not_miss['REF']) & (df_not_miss['a2'] != df_not_miss['REF']))) & (df_not_miss['a1'] == df_not_miss['a2']) ]
    df_hom_alt['zygosity'] = 'hom-alt'
    
    df_het_ref = df_not_miss[ ((df_not_miss['a1'] == df_not_miss['REF']) & (df_not_miss['a2'] != df_not_miss['REF'])) | ((df_not_miss['a1'] != df_not_miss['REF']) & (df_not_miss['a2'] == df_not_miss['REF']))  ]
    df_het_ref['zygosity'] = 'het-ref'
    

    
    df_zygosity = pd.concat([df_hom_ref, df_hom_miss, df_het_miss, df_het_ref, df_het_alt, df_hom_alt])
    
    assert len(df_zygosity) == len(df)
    return df_zygosity


def zygosity_slow(line):
    allele_set = set([line['a1'], line['a2'], line['REF']])
        
    if len(allele_set) == 1:
        return 'hom-ref'
    if len(allele_set) ==2 and line['a1'] == line['a2']:
        if "." in allele_set:
            return 'hom-miss'
        return 'hom-alt'
    if len(allele_set) == 2 and line['a1'] != line['a2']:
        if "." in allele_set:
            return 'het-miss'
        return 'het-ref'
    if len(allele_set) ==3:
        if "." in allele_set:
            return 'het-miss'
        return 'het-alt'
    assert False

