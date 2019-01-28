
"""
The following methods generate annotations for each
VCF DNA variant
"""

import pandas as pd
import numpy as np
from functools import partial
import multiprocessing as mp
import gc



def add_allelic_bases(df, sample_col, single_sample_vcf=True):
    """This function returns the base call for each biallelic base
    10X faster than previous iterations
    """
    haploid_chromosomes = ['X', 'chrX', 'Y', 'chrY', 'M', 'chrM']

    def vector_GT_alleles(ref_alt_gt):
        '''Retrieve GT1, GT2, a1, a2'''

        def get_phase(genotype):
            """Returns phase from genotype"""
            if "|" in genotype:
                return "|"
            if "/" in genotype:
                return "/"
            else:
                return '-'

        REF, ALT, GT = ref_alt_gt
        phase = get_phase(GT)
        bases = [REF] + str(ALT).split(',')
        gt1, gt2 = np.NaN, np.NaN
        a1, a2 = '.', '.'
        GT = str(GT).split(phase)

        if len(GT) == 2:  #Diploid
            gt1, gt2 = GT

            if gt1 == ".":
                a1 = "."
            else:
                a1 = bases[int(gt1)]

            if gt2 == ".":
                a2 = "."
            else:
                a2 = bases[int(gt2)]

        if len(GT) == 1:  #Haploid
            gt1 = GT[0]
            a1 = bases[int(gt1)]

        return (gt1, gt2, a1, a2, phase)

    df = df.assign(**{'GT1':-1, 'GT2':-1, 'a1':'', 'a2':'', 'phase':''})
    df.loc[:, ['GT1', 'GT2', 'a1', 'a2', 'phase']] = [i for i in map(vector_GT_alleles, df[['REF','ALT', sample_col]].values)]
    return df


def zygosity_fast(df):
    """
    This function quickly assigns zygosity states to
    each variant using set logic
    """

    def check_empty_df(df, zygosity):

        try:
            df.loc[:, 'zygosity'] = zygosity
            return df
        except ValueError:
            if len(df) == 0:
                return pd.DataFrame()
            else:
                assert False


    df_hom_ref = df[(df['a1'] == df['REF']) & (df['a2'] == df['REF'])].copy()
    if len(df_hom_ref) > 0:
        df_hom_ref.loc[:, 'zygosity'] = 'hom-ref'
    else:
        df_hom_ref = check_empty_df(df_hom_ref, 'hom-ref')

    df_hom_miss = df[(df['a1'] == '.') & (df['a2'] == '.')].copy()
    if len(df_hom_miss) > 0:
        df_hom_miss = check_empty_df(df_hom_miss, 'hom-miss')
        #df_hom_miss.loc[:, 'zygosity'] = 'hom-miss'

    df_het_miss = df[(df['a1'] == '.') | (df['a2'] == '.')].copy()
    if len(df_het_miss) > 0:
        df_het_miss = check_empty_df(df_het_miss, 'het-miss')
        #df_het_miss.loc[:, 'zygosity'] = 'het-miss'

    df_not_miss = df.drop(set(df_hom_miss.index) |
                              set(df_het_miss.index)).copy()

    df_het_alt = df_not_miss[((df_not_miss['a1'] != df_not_miss['REF']) &
                              (df_not_miss['a2'] != df_not_miss['REF'])) &
                             (df_not_miss['a1'] != df_not_miss['a2'])].copy()
    df_het_alt = check_empty_df(df_het_alt, 'het-alt')
    #df_het_alt.loc[:, 'zygosity'] = 'het-alt'

    df_hom_alt = df_not_miss[(((df_not_miss['a1'] != df_not_miss['REF']) &
                               (df_not_miss['a2'] != df_not_miss['REF']))) &
                             (df_not_miss['a1'] == df_not_miss['a2'])].copy()
    df_hom_alt = check_empty_df(df_hom_alt, 'hom-alt')
    #df_hom_alt.loc[:, 'zygosity'] = 'hom-alt'

    df_het_ref = df_not_miss[((df_not_miss['a1'] == df_not_miss['REF']) &
                              (df_not_miss['a2'] != df_not_miss['REF'])) |
                             ((df_not_miss['a1'] != df_not_miss['REF']) &
                              (df_not_miss['a2'] == df_not_miss['REF']))].copy()
    df_het_ref = check_empty_df(df_het_ref, 'het-ref')
    #df_het_ref.loc[:, 'zygosity'] = 'het-ref'

    df_zygosity = pd.concat([df_hom_ref, df_hom_miss,
                             df_het_miss, df_het_ref,
                             df_het_alt, df_hom_alt])

    df_zygosity.loc[:, 'zygosity'] = df_zygosity['zygosity'].astype('category')
    # return df_zygosity
    assert len(df_zygosity) == len(df)
    return df_zygosity


def vartype_map(ref_alt_bases):
    """
    This function assigns the following vartypes to the
    allele specified by allele_base_col: snp, mnp, ins, del, indel or SV
    """
    ref, alt = str(ref_alt_bases[0]), str(ref_alt_bases[-1])
    len_diff = len(ref) - len(alt)

    if ref == alt:
        return 'ref'  # Orderd by frequency of the variant to reduce complexity

    if len_diff == 0:
        base_diff = [nt for i, nt in enumerate(alt) if ref[i] != alt[i]]
        if len(base_diff) == 1:
            return 'snp'
        else:
            return 'mnp'

    if len_diff > 0:
        base_diff = [nt for i, nt in enumerate(alt) if ref[i] != alt[i]]
        if len(base_diff) > 0:
            return 'indel'
        else:
            return 'del'

    if len_diff < 0:
        # base_diff = [nt for i,nt in enumerate(ref) if ref[i] != alt[i]]
        return 'ins'

    # elif is_sv(ref,alt): return 'sv'

    else:
        return 'indel or SV'


def get_hom_ref_counts(df):
    """
    This function calculates the number of homozygous reference variant
    calls in a dataframe assuming the df is indexed on:
    ['CHROM', 'POS', 'REF', 'ALT'] in that order.

    Also assumes the homozygous reference values are ascribed:
        0|0 , 0/0 , 0
    """
    if 'hom-ref' in df.zygosity.value_counts().index:
        hom_ref = df.groupby(['CHROM', 'POS', 'REF', 'ALT'])['zygosity'].value_counts().xs('hom-ref',level=4)
        hom_ref = pd.DataFrame(hom_ref).reset_index()
        hom_ref = hom_ref.rename(columns={'zygosity':'hom_ref_counts'})
        return hom_ref
    else:
        return pd.DataFrame()


def parse_single_genotype_data(df, sample_id, split_cols=''):
    """
    This function parses the genotype sample column and left joins to
    the df object.

    split_cols is a dictionary specifying the name of the column and the number
    of values expected for that column, e.g. {'AD':2, 'PL':3}
    """

    # genotypes grouped by FORMAT variant annotations
    genotypes = df.groupby(by='FORMAT')

    # Iterate through genotype groups, dropping missing calls
    master_df = []
    for name, group in genotypes:
        temp_group = group[sample_id].astype(str)  # group of interest
        # del temp_group['FORMAT']  #remove the format column
        # replace . with none, allows stack to remove null columns, space
        # savings
        temp_group = temp_group.replace(to_replace='.', value='')

        temp_group_data = pd.DataFrame.from_records(
            list(temp_group.str.split(':')))
        temp_group_data.index = temp_group.index
        temp_group_data.columns = name.split(':')
        temp_group_data = temp_group_data.replace(to_replace='.', value='')

        master_df.append(temp_group_data)

    # Concatenating all genotype groups
    sample_df = pd.concat(master_df)
    sample_df.index.names = ['CHROM', 'POS', 'REF', 'ALT', 'sample_ids']

    # spliting user-defined columns
    if split_cols != '':
        for col in split_cols:
            for i in range(0, split_cols[col]):
                sample_df[col + '_' + str(i)] = sample_df[col].str.split(',').str[i]
            del sample_df[col]
    return sample_df



def process_variant_annotations(df_vars_split_cols_sample_id_drop_hom_ref):
    """
    This function stacks a pandas vcf dataframe and adds annotations for
    each genotype

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
    """
    df_vars, split_columns, sample_id, drop_hom_ref = df_vars_split_cols_sample_id_drop_hom_ref

    df_groups = df_vars.groupby('FORMAT')

    parsed_df = []
    # iterate through different FORMAT types
    for format, df_format in df_groups:

        # dropping missing ALT alleles
        df_format = df_format[df_format['ALT'] != '.']
        df_format = df_format[sample_id]  # only consider sample columns
        # replacing missing calls with None
        df_format = df_format.replace(to_replace='.', value=np.NaN)

        # stacks sample calls and drops none calls
        df_format = pd.DataFrame(
            df_format.stack(), columns=['sample_genotypes'])

        if len(df_format) < 1:  # occurs when all calls are empty
            continue

        # SAVE QUALITY INFORMATION SEPARETELY TO AVOID ANNOTATION PROCESSING
        # IDENTICAL GENOTYPE CALLS (DIFFERENT QUALITY DOESNT MATTER)
        if format.count(':') > 0:
            # qual df, setting aside for later joining
            df_qual = pd.DataFrame(list(df_format['sample_genotypes'].str.split(':')),
                                   index=df_format.index)
            # print df_format.head(), format.split(':')
            df_qual.columns = format.split(':')  # setting quality column names
            # setting index names for joining with df_format later
            df_qual.index.names = ['CHROM', 'POS', 'REF', 'ALT', 'sample_ids']
            # setting just the GT calls
            df_format['sample_genotypes'] = df_qual[format.split(':')[0]]
            # removing from df_qual to avoid joining problems with df_format
            # after add_annotations
            del df_qual['GT']

        # DROPPING MISSING CALLS
        df_format = df_format[(df_format['sample_genotypes'] != './.') &
                              (df_format['sample_genotypes'] != '.|.') &
                              (df_format['sample_genotypes'] != '.')]

        # SETTING INDICES
        # setting index names
        df_format.index.names = ['CHROM', 'POS', 'REF', 'ALT', 'sample_ids']
        df_format.reset_index(inplace=True)

        # ONLY NEED TO PASS UNIQUE GENOTYPE CALLS DF TO get_vcf_annotations,
        # then broadcast back to df_format
        df_annotations = df_format.drop_duplicates(subset=['CHROM', 'POS',
                                                           'REF', 'ALT',
                                                           'sample_genotypes'])
        df_annotations['FORMAT'] = format.split(':')[0]  # setting format id
        df_annotations.set_index(['CHROM', 'POS', 'REF',
                                  'ALT', 'sample_genotypes'],
                                 drop=False, inplace=True)
        # getting annotations
        df_annotations = get_vcf_annotations(df_annotations,
                                             'sample_genotypes',
                                             split_columns=split_columns)

        # SETTING INDICES AGAIN
        if len(df_annotations) < 1:
            continue  # continue if no variants within this FORMAT category
        df_format.set_index(['CHROM', 'POS', 'REF',
                             'ALT', 'sample_genotypes'],
                            drop=True, inplace=True)
        df_annotations.index.names = ['CHROM', 'POS', 'REF',
                                      'ALT', 'sample_genotypes']
        df_format = df_format.join(df_annotations)

        # df_format.set_index('sample_ids', drop=True, inplace=True, append=True)
        df_format['FORMAT'] = format
        df_format.reset_index(level=4, inplace=True, drop=False)

        if drop_hom_ref:
            hom_ref_counts = get_hom_ref_counts(df_format)
            hom_ref_counts.name = 'hom_ref_counts'
            # dropping all homozygous reference variants
            df_format = df_format[df_format['zygosity'] != 'hom-ref']
            df_format = df_format.join(hom_ref_counts)
            df_format['hom_ref_counts'].fillna(value=0, inplace=True)

        del df_format['sample_genotypes']
        df_format.set_index('sample_ids', inplace=True, append=True, drop=True)

        # JOINING QUAL INFO BACK TO DF
        if format.count(':') > 0 and len(df_qual) > 0:
            df_format = df_format.join(df_qual, how='left')
            pass

        # SPLITTING GENOTYPE QUALITY COLUMNS
        if split_columns != '':
            for col in split_columns:
                split_col_names = [col + '_' + str(n)
                                   for n in range(0, split_columns[col])]
                df_format = df_format.join(pd.DataFrame(list(df_format[col].str.split(',').str[:len(split_col_names)]),
                                                        index=df_format.index,
                                                        columns=split_col_names))
                del df_format[col]

        parsed_df.append(df_format)

    if len(parsed_df) > 0:
        df_annot = pd.concat(parsed_df)
        # reseting sample_ids from index
        df_annot.reset_index('sample_ids', drop=False, inplace=True)
        return df_annot
    else:
        print('No Annotations generated, please check for excessive missing values')
        return df_vars


def df_split(df, split_level):
    """
    Splits pandas dataframe into roughly
    equal sizes

    Parameters
    ---------------
    df: pandas df, required
        VCF pandas dataframe

    split_level: int, required
        Specifies the number of chunks to split df into

    """
    row_count = len(df)
    split_size = int(row_count / split_level)
    split_df = []
    for n, i in enumerate(range(0, row_count, split_size)):
        if n + 1 == split_level:
            split_df.append(df.ix[df.index[i:]])
            break
        else:
            split_df.append(df.ix[df.index[i: i + split_size]])
    return split_df


def mp_variant_annotations(df_mp, df_split_cols='', df_sampleid='all',
                           drop_hom_ref=True, n_cores=1):
    """
    Multiprocessing variant annotations

    see variantAnnotations.process_variant_annotations for description of annotations


    This function coordinates the annotation of variants using the
    multiprocessing library.

    Parameters
    ---------------
    df_mp: pandas df, required
        VCF DataFrame

    df_split_cols: dict, optional
        key:FORMAT id value:#fields expected
        e.g. {'AD':2} indicates Allelic Depth should be
        split into 2 columns.

    df_sampleid: list, required
        list of sample_ids, can be 'all'

    drop_hom_ref: bool, optional
        specifies whether to drop all homozygous reference
        variants from dataframe.
        FALSE REQUIRES LARGE MEMORY FOOTPRINT

    n_cores: int, optional
        Number of multiprocessing jobs to start.
        Be careful as memory is copied to each process, RAM intensive
    """
    from functools import partial
    import multiprocessing as mp
    import gc

    print('starting multiprocessing')
    pool = mp.Pool(int(n_cores))
    # tasks = np.array_split(df_mp.copy(), int(n_cores))  #breaks with older
    # pandas/numpy
    dfs = df_split(df_mp.copy(), int(n_cores))

    mp_process =  partial(process_variant_annotations, sample_id=df_sampleid,
                 split_columns=df_split_cols, drop_hom_ref=drop_hom_ref)

    results = []
    del df_mp
    gc.collect()
    r = pool.map_async(mp_process, \
                       dfs, callback=results.append)
    r.wait()
    pool.close()
    pool.join()
    pool.terminate()

    print('multiprocessing complete')
    res_df = pd.concat([df for df in results[0] if len(df) > 0])

    cat_cols = ['vartype1', 'vartype2', 'a1', 'a2', \
                 'GT1', 'GT2', 'GT','sample_ids', 'zygosity']
    res_df.loc[:, cat_cols] = res_df[cat_cols].astype('category')
    return res_df


def get_vcf_annotations(df, sample_name, split_columns='', drop_hom_ref=True):
    """
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



    """

    df.loc[:, 'multiallele'] = df.ALT.str.count(',')
    multidf = df[df['multiallele'] > 0]

    while len(df) > 0:


        df = add_allelic_bases(df, sample_name)

        df = zygosity_fast(df)


        df.loc[:, 'vartype1'] = [vtype for vtype in map(vartype_map, df[['REF', 'a1']].values)]
        df.loc[:, 'vartype2'] = [vtype for vtype in map(vartype_map, df[['REF', 'a2']].values)]

        cat_cols = ['vartype1', 'vartype2', 'a1', 'a2', 'GT1', 'GT2', 'sample_genotypes', 'phase']
        for c in cat_cols:
            df.loc[:, c] = df[c].astype('category')


        df.loc[:, 'GT'] = df['sample_genotypes'].astype('category')
        del df[sample_name]
        if 'FORMAT' in df.columns:
            del df['FORMAT']


        return df

    return pd.DataFrame()


def process_variant_annotations(df_vars, sample_id='all', split_columns='', drop_hom_ref=False):
    """
    This function stacks a pandas vcf dataframe and adds annotations for
    each genotype

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
    """
    #df_vars, split_columns, sample_id, drop_hom_ref = df_vars_split_cols_sample_id_drop_hom_ref


    def _format_preprocess(df, sample_id):
        """
        Handles stacking the wide dataframe into a long dataframe of sample genotypes
        """
         # dropping missing ALT alleles
        df = df[df['ALT'] != '.']
        df = df[['CHROM', 'POS', 'REF', 'ALT'] + sample_id]  # only consider sample columns
        # replacing missing calls with None
        df = df.replace(to_replace='.', value=np.NaN)

        # stacks sample calls and drops none calls

        df = df.set_index(['CHROM', 'POS', 'REF', 'ALT'])

        df = pd.DataFrame(df.stack(),
                          columns=['sample_genotypes'])
        df.index.names = ['CHROM', 'POS', 'REF', 'ALT', 'sample_ids']

        return df.reset_index()


    def _sampleid_preprocess(df):
        """
        Identifies sample columns
        """
        s_ids = set(df.columns) - set(['CHROM', 'POS', 'REF', 'ALT', \
                                           'ID', 'QUAL', 'FILTER', 'INFO','FORMAT'])
        return list(s_ids)


    def _qual_preprocess(df, form):
        """
        Creates dataframe with non-GT data for each genotype call,
        often quality data
        """
         # qual df, setting aside for later joining
        df_qual = pd.DataFrame(list(df['sample_genotypes'].str.split(':')),
                               df.index)
        df_qual.columns = form.split(':')  # setting quality column names
        df_qual = df[['CHROM', 'POS', 'REF', 'ALT', 'sample_ids']].join(df_qual)
        # print df_format.head(), format.split(':')

        # setting index names for joining with df_format later
        # setting just the GT calls
        df.loc[:, 'sample_genotypes'] = df_qual[form.split(':')[0]]
        # removing from df_qual to avoid joining problems with df_format
        # after add_annotations
        #del df_qual['GT']
        return df_qual


    def _missing_preprocess(df):
        """
        Filters dataframe for missing values
        """
        df_nonmissing = df[(df['sample_genotypes'] != './.') &
                           (df['sample_genotypes'] != '.|.') &
                           (df['sample_genotypes'] != '.')]

        return df_nonmissing


    def _coordinate_variant_annotation(df_format, sample_id):
        """
        Coordinates variant annotations, hom-ref counting and dropping,
        and stacking into a tidy df

       Parameters
        --------------
        df_format: pandas DataFrame, required
                    pd.DataFrame containing CHROM, POS, REF, ALT, sample genotype cols

        sample_id: list, required
                    list of sample genotype columns


        Output
        --------------
        This function produces the df_annot dataframe containing the following columns:
        CHROM, POS, REF, ALT, GT, GT1, GT2, a1, a2, multiallele, phase, zygosity,
        vartype1, vartype2, FORMAT, hom_ref_counts
        """

        df_format = _format_preprocess(df_format, sample_id)

        if len(df_format) < 1:  # occurs when all calls are empty
            return pd.DataFrame()

        # SAVE QUALITY INFORMATION SEPARETELY TO AVOID ANNOTATION PROCESSING
        # IDENTICAL GENOTYPE CALLS (DIFFERENT QUALITY DOESNT MATTER)
        if form.count(':') > 0:
            df_qual = _qual_preprocess(df_format, form)

        # DROPPING MISSING CALLS
        df_format = _missing_preprocess(df_format)


        # SETTING INDICES
        # setting index names
        # df_format.index.names = ['CHROM', 'POS', 'REF', 'ALT', 'sample_ids']
        # df_format = df_format.reset_index()

        # ONLY NEED TO PASS UNIQUE GENOTYPE CALLS DF TO get_vcf_annotations,
        # then broadcast back to df_format
        annot_cols = ['CHROM', 'POS','REF', 'ALT','sample_genotypes']
        df_annotations = df_format[annot_cols].drop_duplicates(subset=annot_cols)

        df_annotations.loc[:, 'FORMAT'] = form.split(':')[0]  # setting format id

        # get variant annotations
        df_annotations = get_vcf_annotations(df_annotations,
                                             'sample_genotypes',
                                             split_columns=split_columns)

        # BROADCASTING VARIANT ANNOTATIONS BACK TO SAMPLE GENOTYPE DF
        if len(df_annotations) < 1:
            return pd.DataFrame()  # continue if no variants within this FORMAT category

        df_format.rename(columns={'sample_genotypes':'GT'}, inplace=True)
        df_format.loc[:, 'GT'] = df_format['GT'].astype('category')
        df_format = df_format.merge(df_annotations, how='left',
                                  left_on = ['CHROM', 'POS', 'REF', 'ALT','GT'],
                                  right_on = ['CHROM', 'POS', 'REF', 'ALT','GT'])

        del df_annotations
        gc.collect()
        # df_format.set_index('sample_ids', drop=True, inplace=True, append=True)
        df_format.loc[:, 'FORMAT'] = form


        # DROPPING HOMOZYGOUS REFERENCE VARIANTS IF SPECIFIED BY USER
        hom_ref_counts = get_hom_ref_counts(df_format)
        if len(hom_ref_counts) > 0:
            df_format = df_format.merge(hom_ref_counts, how='left', \
                                            left_on=['CHROM', 'POS','REF','ALT'], \
                                            right_on=['CHROM', 'POS','REF','ALT'])
            df_format['hom_ref_counts'].fillna(value=0, inplace=True)
            df_format.loc[:, 'hom_ref_counts'] = df_format['hom_ref_counts'].astype(np.uint8)
        else:
            df_format.loc[:, 'hom_ref_counts'] = -1

        if drop_hom_ref:
            # dropping all homozygous reference variants
            df_format = df_format[df_format['zygosity'] != 'hom-ref']

        # JOINING QUAL INFO BACK TO DF, IF NON-GT FIELDS IN SAMPLE COLUMNS
        if form.count(':') > 0 and len(df_qual) > 0:
            df_format = df_format.merge(df_qual, how='left',
                                      left_on = ['CHROM', 'POS', 'REF', 'ALT', 'GT', 'sample_ids'],
                                      right_on = ['CHROM', 'POS', 'REF', 'ALT', 'GT', 'sample_ids'])
            del df_qual
            gc.collect()
            pass

        # SPLITTING GENOTYPE QUALITY COLUMNS, IF SPECIFIED BY USER
        if split_columns != '':
            for col in split_columns:
                if split_columns[col] > 1:  #only parse split_columns with more than expected 1 column
                    split_col_names = [col + '_' + str(n) for n in range(0, split_columns[col])]
                    try:
                        split_col_df = pd.DataFrame(list(df_format[col].str.split(',') \
                                                                       .str[:len(split_col_names)]),
                                                            index=df_format.index,
                                                            columns=split_col_names)
                        df_format = df_format.join(split_col_df)
                        del df_format[col]
                    except AssertionError:
                        print('{} has incorrect column number, '.format(col) + \
                               'please check split_cols value. ' + \
                               'Leaving {} unparsed.'.format(col))
                        print()
                else:
                    continue
        return df_format


    if '#CHROM' in df_vars.columns:
        df_vars = df_vars.rename(columns={'#CHROM':'CHROM'}).set_index(['CHROM', 'POS', 'REF', 'ALT'],drop=False)

    if sample_id == 'all':
        sample_id = _sampleid_preprocess(df_vars)

    df_groups = df_vars.groupby('FORMAT')

    parsed_df = []
    # iterate through different FORMAT types
    for form, df_format in df_groups:

        df_format = _coordinate_variant_annotation(df_format, sample_id)

        parsed_df.append(df_format)

    if len(parsed_df) > 0:
        df_annot = pd.concat(parsed_df)
        
        for c in ['sample_ids', 'FORMAT']:
            df_annot.loc[:, c] = df_annot[c].astype('category')
        df_annot.loc[:, ['multiallele']] = df_annot[['multiallele']].astype(np.uint8)
        if 'hom_ref_counts' not in df_annot.columns:
            df_annot.loc[:, 'hom_ref_counts'] = -1
        df_annot.loc[:, 'hom_ref_counts'] = df_annot['hom_ref_counts'].astype(np.uint8)

        return df_annot
    else:
        print('No Annotations generated, please check for excessive missing values')
        return pd.DataFrame()
