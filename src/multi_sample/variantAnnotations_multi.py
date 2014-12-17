
'''
The following methods generate annotations for each
VCF DNA variant
'''

import os,sys,gzip
import pandas as pd
import numpy as np
import multiprocessing as mp
import gc

def get_multiallelic_bases(df_orig, sample_col, single_sample_vcf=True):
    '''
    This function parses multiallele variants into DNA base representations.
    It currently does not support haploid chromosomes.
    '''
    
    haploid_chromosomes = ['X', 'chrX', 'Y', 'chrY', 'M', 'chrM']

    df = df_orig.copy()
    
    def get_phase(line, sample_id):
        '''
        Returns phase from genotype
        '''
        genotype = str(line[sample_id])
        if "|" in genotype:
            return "|"
        if "/" in genotype:
            return "/"
        else:
            return '-'
    
    
    def _get_allele(line, gt_col):
        '''
        Returns allele base call from multi-allelic variants
        '''
        
        alleles = [line['REF']]
        alleles.extend(list( line['ALT'].split(",")) )
        a1 = "."
        try:
            a1 = alleles[int(line[gt_col])]  #returns missing if gt_int_call is "."
        except:
            a1 = "."
        return a1


    def get_GT_multisample_vcf(line, sample_col, gt_index):
        return line[sample_col].split(':')[0].split(line['phase'])[int(gt_index)]
    
    def get_GT_multisample_vcf_haploid(line, sample_col, gt_index):
        return str(line[sample_col]).split(':')[0]

    
#    if single_sample_vcf:
#        df['phase'] = df[sample_col].str[1]
#        df = df[df['phase']!=':']  #removing haploid variants
#        
#        df['GT1'] = df[sample_col].str[0]
#        df = df[df['GT1']!='.']  #removing variants with missing calls
#        df['GT1'] = df['GT1'].astype(int)
#        
#        df['GT2'] = df[sample_col].str[2]
#        df = df[df['GT2']!='.']  #removing variants with missing calls
#        df['GT2'] = df['GT2'].astype(int)


    if not single_sample_vcf:
        df['phase'] = df.apply(get_phase, args=[sample_col], axis=1)  #get phase
        haploid_df = df[df.phase == "-"]  #likley occurs at sex chromosome sites
        haploid_df = haploid_df[haploid_df.CHROM.isin(haploid_chromosomes)]
        
        if len(haploid_df) > 0:
            haploid_df['GT1'] = df.apply(get_GT_multisample_vcf_haploid, args=[sample_col, 0], axis=1)
            haploid_df = haploid_df[ (haploid_df['GT1']!='.') & (haploid_df['GT1']!=np.NaN)]
            haploid_df['GT1'] = haploid_df['GT1'].astype(int)
            haploid_df['GT2'] = 0
            haploid_df['a1'] = haploid_df.apply(_get_allele, args=['GT1'], axis=1)
            haploid_df['a2'] = haploid_df.apply(_get_allele, args=['GT2'], axis=1)
        
        if len(df) > 0:
        
            df = df[df.phase != "-"]
            
            df['GT1'] = df.apply(get_GT_multisample_vcf, args=[sample_col, 0], axis=1)
            df = df[ (df['GT1']!='.') & (df['GT1']!=np.NaN)]
            df['GT1'] = df['GT1'].astype(int)
            
            
            df['GT2'] = df.apply(get_GT_multisample_vcf, args=[sample_col, 1], axis=1)
            df = df[ (df['GT2']!='.') & (df['GT2']!=np.NaN)]
            df['GT2'] = df['GT2'].astype(int)
            df['a1'] = df.apply(_get_allele, args=['GT1'], axis=1)
            df['a2'] = df.apply(_get_allele, args=['GT2'], axis=1)
    
    #if len(df_multi) > 0:
    #    df = df.append(df_multi)
    
        
    #df['a1'] = df.apply(_get_allele, args=['GT1'], axis=1)
    #df['a2'] = df.apply(_get_allele, args=['GT2'], axis=1)


    if len(df) > 0:
        if len(haploid_df) > 0:
            df = df.append(haploid_df)  #adding haploid variants to dataframe
            return df
        else:
            return df
    else:
        return haploid_df




def get_biallelic_bases(df, sample_col, single_sample_vcf=True):
    '''
    This function is 10X faster than previous iterations

    '''
    haploid_chromosomes = ['X', 'chrX', 'Y', 'chrY', 'M', 'chrM']
    
    
    def get_phase(line, sample_id):
        '''
        Returns phase from genotype
        '''
        genotype = str(line[sample_id])
        if "|" in genotype:
            return "|"
        if "/" in genotype:
            return "/"
        else:
            return '-'
    

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
        return int(line[sample_col].split(line['phase'])[int(gt_index)])
    
    def get_GT_multisample_vcf_haploid(line, sample_col, gt_index):
        return str(line[sample_col]).split(':')[0]


    if single_sample_vcf:
        df['GT_len'] = df[sample_col].str.split(':').str[0].str.len()
        haploid_df = df[df['GT_len'] <= 1]
        haploid_df['phase'] = '-'
        haploid_df = haploid_df[haploid_df['CHROM'].isin(haploid_chromosomes)]
        
        df = df[df['GT_len'] > 1]
        df['phase'] = df[sample_col].str[1]
        df = df[ (df['phase']!='-') ]
        
        del df['GT_len']
        del haploid_df['GT_len']
        
        if len(haploid_df) > 0:
            haploid_df['GT1'] = haploid_df[sample_col].str[0]
            haploid_df = haploid_df[ (haploid_df['GT1']!='.') & (haploid_df['GT1']!=np.NaN)]
            haploid_df['GT2'] = 0
        
        
        if len(df) > 0:
            df['GT1'] = df[sample_col].str[0]
            df = df[(df['GT1']!='.') & (df['GT1']!=np.NaN)]
            df['GT1'] = df['GT1'].astype(int)
            df['GT2'] = df[sample_col].str[2]
            df = df[(df['GT2']!='.') & (df['GT2']!=np.NaN)]
            df['GT2'] = df['GT2'].astype(int)
            
        


    if not single_sample_vcf:
        df['GT_len'] = df[sample_col].str.split(':').str[0].str.len()
        haploid_df = df[df['GT_len'] > 1]
        haploid_df['phase'] = '-'
        haploid_df = haploid_df[haploid_df['CHROM'].isin(haploid_chromosomes)]
        
        df = df[~df.index.isin(haploid_df.index)]
        df['phase'] = df[sample_col].str[1]
        df = df[ (df['phase']!='-') ]
        
        del df['GT_len']
        del haploid_df['GT_len']
        
        if len(df) > 0:
            df['GT1'] = df.apply(get_GT_multisample, args=[sample_col, 0], axis=1)
            df = df[(df['GT1']!='.') & (df['GT1']!=np.NaN)]
            df['GT2'] = df.apply(get_GT_multisample, args=[sample_col, 1], axis=1)
            df = df[(df['GT2']!='.') & (df['GT2']!=np.NaN)]
            
        if len(haploid_df) > 0:
            haploid_df['GT1'] = df.apply(get_GT_multisample_vcf_haploid, args=[sample_col, 0], axis=1)
            haploid_df = haploid_df[ (haploid_df['GT1']!='.') & (haploid_df['GT1']!=np.NaN)]
            haploid_df['GT1'] = haploid_df['GT1'].astype(int)
            haploid_df['GT2'] = 0



    if len(df) > 0:
        if len(haploid_df) > 0:
            df = df.append(haploid_df)
        else: pass
    else:
        df = haploid_df
        


    #FAST PROCESS SIMPLE ALLELE GENOTYPES
    df_simple = df
    
    df_gt1_ref = df_simple[df_simple.GT1.astype(int)==0][['REF']]  #get a1 ref alleles
    df_gt1_ref.columns = ['a1']
    df_gt2_ref = df_simple[df_simple.GT2.astype(int)==0][['REF']]  #get a2 ref alleles
    df_gt2_ref.columns = ['a2']


    df_gt1_alt = df_simple[df_simple.GT1.astype(int)==1][['ALT']]  #get a1 alt alleles
    df_gt1_alt.columns = ['a1']
    df_gt2_alt = df_simple[df_simple.GT2.astype(int)==1][['ALT']]  #get a2 alt alleles
    df_gt2_alt.columns = ['a2']


    gt1_alleles = pd.concat([df_gt1_ref,df_gt1_alt])  #merging GT1 allele bases into a single df
    #del gt1_alleles[0]
    gt2_alleles = pd.concat([df_gt2_ref,df_gt2_alt])  #merging GT2 allele bases into a single df
    #del gt2_alleles[0]
    gt1_2_allele_df = gt1_alleles.join(gt2_alleles)  #Joining the GT1 and GT2 simple allele bases 

    #print len(df)
    
    return df.join(gt1_2_allele_df)
    
    #print len(df)






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
    #return df_zygosity
    assert len(df_zygosity) == len(df)
    return df_zygosity




def vartype_map(ref_alt_bases):
    '''
    This function assigns the following vartypes to the 
    allele specified by allele_base_col: snp, mnp, ins, del, indel or SV

    '''
    ref, alt = ref_alt_bases
    len_diff = len(ref) - len(alt)

    if ref == alt: return 'ref' #Orderd by frequency of the variant to reduce complexity

    if len_diff == 0:
        base_diff = [nt for i,nt in enumerate(alt) if ref[i] != alt[i]]
        if len(base_diff) == 1: return 'snp' 
        else: return 'mnp'

    if len_diff > 0:
        base_diff = [nt for i,nt in enumerate(alt) if ref[i] != alt[i]]
        if len(base_diff) > 0: return 'indel'
        else: return 'del'

    if len_diff < 0:
        #base_diff = [nt for i,nt in enumerate(ref) if ref[i] != alt[i]]
        return 'ins'

    elif is_sv(ref,alt): return 'sv'

    else: return 'indel or SV'




def parse_single_genotype_data(df, sample_id, split_cols=''):
    '''
    This function parses the genotype sample column and left joins to
    the df object.

    split_cols is a dictionary specifying the name of the column and the number
    of values expected for that column, e.g. {'AD':2, 'PL':3}
    '''


    genotypes = df.groupby(by='FORMAT')  #genotypes grouped by FORMAT variant annotations

    #Iterate through genotype groups, dropping missing calls
    master_df = []
    for name,group in genotypes:
        temp_group = group[sample_id].astype(str)  #group of interest
        #del temp_group['FORMAT']  #remove the format column
        temp_group.replace(to_replace='.', value='', inplace=True)  #replace . with none, allows stack to remove null columns, space savings

        temp_group_data = pd.DataFrame.from_records(list(temp_group.str.split(':')))
        temp_group_data.index = temp_group.index
        temp_group_data.columns = name.split(':')
        temp_group_data.replace(to_replace='.', value='', inplace=True)

        master_df.append(temp_group_data)

    #Concatenating all genotype groups
    sample_df = pd.concat(master_df)
    sample_df.index.names = ['CHROM', 'POS', 'REF', 'ALT', 'sample_ids']

    #spliting user-defined columns
    if split_cols != '':
        for col in split_cols:
            for i in range(0, split_cols[col]):
                sample_df[col + '_' + str(i)] = sample_df[col].str.split(',').str[i]
            del sample_df[col]
    return sample_df


def process_variant_annotations(df_vars_split_cols_sample_id):
    '''
    This function stacks a pandas vcf dataframe and adds annotations
    
    '''
    df_vars, split_columns, sample_id = df_vars_split_cols_sample_id
    df_groups = df_vars.groupby('FORMAT')
    
    parsed_df = []
    for i,df_format in df_groups:
        
        df_format = df_format[df_format['ALT'] != '.']
        df_format = df_format[sample_id]
        df_format = df_format.replace(to_replace='.', value=np.NaN)
        df_format = pd.DataFrame( df_format.stack(), columns=['sample_genotypes'] )
        df_format.index.names = ['CHROM', 'POS', 'REF', 'ALT', 'sample_ids']
        df_format.reset_index(inplace=True)
        df_format.set_index(['CHROM', 'POS', 'REF', 'ALT', 'sample_ids'], drop=False, inplace=True)
        df_format['FORMAT'] = i
        df_format.drop_duplicates(inplace=True)
        parsed_df.append( get_vcf_annotations(df_format, 'sample_genotypes', split_columns=split_columns) )
        
    
    df_annot = pd.concat(parsed_df)
    
    return df_annot


def mp_variant_annotations(df_mp, df_split_cols, df_sampleid, n_cores=1):
    pool = mp.Pool(int(n_cores))
    tasks = np.array_split(df_mp.copy(), int(n_cores))
    tasks = [[pd.DataFrame(t), df_split_cols, df_sampleid] for t in tasks]
    results =[]
    del df_mp
    gc.collect()
    r = pool.map_async(process_variant_annotations, tasks, callback=results.append)
    r.wait()
    pool.close()
    pool.join()
    return pd.concat(results[0])


def get_vcf_annotations(df, sample_name, split_columns='', drop_hom_ref=True):
        '''
        This function adds the following annotations for each variant:
        multiallele, phase, a1, a2, GT1, GT2, vartype1, vartype2, zygosity,
        and parsed FORMAT values, see below for additional information.
                        
        Parameters
        --------------
        sample_name: str, required 
                    sample column header id, e.g. NA12878
        
        split_columns: dict, optional
                    key:FORMAT id value:#fields expected
                    e.g. {'AD':2} indicates Allelic Depth should be
                    split into 2 columns.
        
        drop_hom_ref: bool, optional
                    specifies whether to drop all homozygous reference
                    variants from dataframe.
                    FALSE REQUIRES LARGE MEMORY FOOTPRINT
        
        Output
        --------------
        This function adds the following annotations to each variant:
        
        multiallele: {0,1} 0=biallele  1=multiallelic
    
        phase: {'/', '|'} /=unphased, |=phased
        
        a1: DNA base representation of allele1 call, e.g. A
        a2: DNA base representation of allele2 call, e.g. A
        
        GT1: numeric representation of allele1 call, e.g. 0
        GT2: numeric representation of allele2 call, e.g. 1
        
        vartype1: {snp, mnp, ins, del, indel or SV} variant type of first allele
        vartype2: {snp, mnp, ins, del, indel or SV} variant type of second allele
        
        zygosity: {het-ref, hom-ref, alt-ref, het-miss, hom-miss}
        
        FORMAT values: any values associated with the genotype calls are 
                        added as additional columns, split_columns are further
                        split by ',' into individual columns
        
        
        
        '''
     
        
        
        
        
        
        df['multiallele'] = df.ALT.str.count(',')
        multidf = df[df['multiallele'] > 0]
        df = df[~df.index.isin(multidf.index)]
        #print len(multidf), 'multidf rows'
        
        if len(multidf) > 0:
            multidf = get_multiallelic_bases(multidf, sample_name, single_sample_vcf=False)
        
        
        #print 'single alleles', len(df)
        
        df = get_biallelic_bases(df, sample_name)
        
        
        if len(multidf) > 0:
            df = df.append(multidf)
        
        
        df = zygosity_fast(df)
        
        if drop_hom_ref:
            #recording number of homozygous reference calls for each variant
            hom_ref_counts = df.groupby(level=[0,1,2,3])['zygosity'].value_counts()
            hom_ref_counts = hom_ref_counts.unstack(level=4)['hom-ref']
            hom_ref_counts.name = 'hom_ref_counts'
            
            df.reset_index(level=4, inplace=True, drop=True)
            df = df[df['zygosity']!='hom-ref']  #dropping all homozygous reference variants
            df = df.join(hom_ref_counts, how='left')
            df['hom_ref_counts'].fillna(value=0, inplace=True)
        
        df['vartype1'] = map(vartype_map, df[['REF','a1']].values)
        df['vartype2'] = map(vartype_map, df[['REF','a2']].values)
        
        df.set_index(['CHROM', 'POS', 'REF', 'ALT', 'sample_ids'], inplace=True)
        
        #df.sortlevel(level=['CHROM','POS','REF','ALT','sample_ids'],inplace=True)  #sorting biallelic and multiallele variants 
        
        #print 'before parse_single_genotype_data', len(df)
        
        
        df = df.join( parse_single_genotype_data(df, sample_name, split_cols=split_columns), how='left' )
        del df[sample_name]
        if 'FORMAT' in df.columns:
            del df['FORMAT']
        
        
        return df




