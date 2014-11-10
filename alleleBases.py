def get_allele_bases(df, sample_col, single_sample_vcf=True):
    '''
    This function is 10X faster than previous iterations
    
    '''
    
    def _get_allele(line, gt_col):
        '''
        Returns allelic base, handles multi-allelic variants
        '''
        alleles = [line['REF']]
        alleles.extend(line['ALT'].split(","))
        a1 = "."
        try:
            a1 = alleles[int(line[gt_col])]  #returns missing if gt_int_call is "."
        except:
            a1 = "."
        return a1
    
    
    def get_GT_multisample_vcf(line, sample_col, gt_index):
        return int( line[sample_col].split(line['phase'])[int(gt_index)])
    
    
    if single_sample_vcf:
        df['phase'] = df[sample_col].str[1]
        df['GT1'] = df[sample_col].str[0]
        df['GT1'] = df['GT1'].astype(int)
        df['GT2'] = df[sample_col].str[2]
        df['GT2'] = df['GT2'].astype(int)
    
    
    if not single_sample_vcf:
        df['phase'] = df.apply(get_phase, args=['GT'], axis=1)  #get phase
        df = df[df.phase != "-"]  #likley occurs at sex chromosome sites
        df['GT1'] = df.apply(get_GT_multisample, args=[sample_col, 0], axis=1)
        df['GT2'] = df.apply(get_GT_multisample, args=[sample_col, 1], axis=1)
        
        
    
    #SLOW PROCESS MULTIPLE ALLELE GENOTYPES
    df_multi = df[(df.GT1>1) | (df.GT2>1)] #select all multi-alleleic variants
    df_multi['a1'] = df_multi.apply(_get_allele, args=['GT1'], axis=1)  #
    df_multi['a2'] = df_multi.apply(_get_allele, args=['GT2'], axis=1)
    
    
    #FAST PROCESS SIMPLE ALLELE GENOTYPES
    df_simple = df[~df.index.isin(df_multi.index)][['REF', 'ALT', 'GT1', 'GT2']]  #dropping multiallele variants, minimize memory usage
    
    df_gt1_ref = df_simple[df_simple.GT1==0][['REF']]  #get a1 ref alleles
    df_gt1_ref.columns = ['a1']
    df_gt2_ref = df_simple[df_simple.GT2==0][['REF']]  #get a2 ref alleles
    df_gt2_ref.columns = ['a2']
    
    
    df_gt1_alt = df_simple[df_simple.GT1==1][['ALT']]  #get a1 alt alleles
    df_gt1_alt.columns = ['a1']
    df_gt2_alt = df_simple[df_simple.GT2==1][['ALT']]  #get a2 alt alleles
    df_gt2_alt.columns = ['a2']
    
    
    gt1_alleles = pd.concat([df_gt1_ref,df_gt1_alt])  #merging GT1 allele bases into a single df
    #del gt1_alleles[0]
    gt2_alleles = pd.concat([df_gt2_ref,df_gt2_alt])  #merging GT2 allele bases into a single df
    #del gt2_alleles[0]
    gt1_2_allele_df = gt1_alleles.join(gt2_alleles, how='outer')  #Joining the GT1 and GT2 simple allele bases 
    
    
    df = df.join(gt1_2_allele_df, how='inner')  #Adding simle allele a1 and a2 columns to original df
    df = df.append(df_multi)  #Adding multi-alleleic bases to original df
    
    #df['vartype1'] = df.apply(vartype, args=['a1'], axis=1)  #vartype for allele1: snp, mnp, insertion, deletion, sv
    #df['vartype2'] = df.apply(vartype, args=['a2'], axis=1)  #vartype for allele2: snp, mnp, insertion, deletion, sv

    return df


