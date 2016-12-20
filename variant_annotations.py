
"""
The following methods generate annotations for each
VCF DNA variant
"""

import pandas as pd
import numpy as np
import multiprocessing as mp
import gc


def get_multiallelic_bases(df_orig, sample_col, single_sample_vcf=True):
    """
    This function parses multiallele variants into DNA base representations.
    It currently does not support haploid chromosomes.
    """

    haploid_chromosomes = ['X', 'chrX', 'Y', 'chrY', 'M', 'chrM']

    df = df_orig.copy()

    def get_phase(line, sample_id):
        """Returns phase from genotype"""
        genotype = str(line[sample_id])
        if "|" in genotype:
            return "|"
        if "/" in genotype:
            return "/"
        else:
            return '-'

    def _get_allele(line, gt_col):
        """Returns allele base call from multi-allelic variants"""

        alleles = [line['REF']]
        alleles.extend(list(line['ALT'].split(",")))
        a1 = "."
        try:
            # returns missing if gt_int_call is "."
            a1 = alleles[int(line[gt_col])]
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
        df['phase'] = df.apply(
            get_phase, args=[sample_col], axis=1)  # get phase
        # likley occurs at sex chromosome sites
        haploid_df = df[df.phase == "-"]
        haploid_df = haploid_df[haploid_df['CHROM'].isin(haploid_chromosomes)]

        if len(haploid_df) > 0:
            haploid_df['GT1'] = df.apply(get_GT_multisample_vcf_haploid,
                                         args=[sample_col, 0], axis=1)
            haploid_df = haploid_df[(haploid_df['GT1'] != '.') &
                                    (haploid_df['GT1'] != np.NaN)]
            haploid_df['GT1'] = haploid_df['GT1'].astype(int)
            haploid_df['GT2'] = 0
            haploid_df['a1'] = haploid_df.apply(_get_allele,
                                                args=['GT1'], axis=1)
            haploid_df['a2'] = haploid_df.apply(_get_allele,
                                                args=['GT2'], axis=1)
        df = df[df.phase != "-"]
        if len(df) > 0:
            df['GT1'] = df.apply(get_GT_multisample_vcf,
                                 args=[sample_col, 0], axis=1)
            df = df[(df['GT1'] != '.') & (df['GT1'] != np.NaN)]
            df['GT1'] = df['GT1'].astype(int)

            df['GT2'] = df.apply(get_GT_multisample_vcf,
                                 args=[sample_col, 1], axis=1)
            df = df[(df['GT2'] != '.') & (df['GT2'] != np.NaN)]
            df['GT2'] = df['GT2'].astype(int)
            df['a1'] = df.apply(_get_allele, args=['GT1'], axis=1)
            df['a2'] = df.apply(_get_allele, args=['GT2'], axis=1)

    # if len(df_multi) > 0:
    #    df = df.append(df_multi)

    # df['a1'] = df.apply(_get_allele, args=['GT1'], axis=1)
    # df['a2'] = df.apply(_get_allele, args=['GT2'], axis=1)

    if len(df) > 0:
        if len(haploid_df) > 0:
            df = df.append(haploid_df)  # adding haploid variants to dataframe
            return df
        else:
            return df
    if len(haploid_df) > 0:
        return haploid_df
    else:
        return pd.DataFrame()


def get_biallelic_bases(df, sample_col, single_sample_vcf=True):
    """This function returns the base call for each biallelic base
    10X faster than previous iterations
    """
    haploid_chromosomes = ['X', 'chrX', 'Y', 'chrY', 'M', 'chrM']

    def get_phase(line, sample_id):
        """Returns phase from genotype"""
        genotype = str(line[sample_id])
        if "|" in genotype:
            return "|"
        if "/" in genotype:
            return "/"
        else:
            return '-'

    def _get_allele(line, gt_col):
        """Returns allelic base, handles multi-allelic variants"""
        alleles = [line['REF']]
        alleles.extend(line['ALT'].split(","))
        a1 = "."
        try:
            # returns missing if gt_int_call is "."
            a1 = alleles[int(line[gt_col])]
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
        df = df[(df['phase'] != '-')]

        del df['GT_len']
        del haploid_df['GT_len']

        if len(haploid_df) > 0:
            haploid_df['GT1'] = haploid_df[sample_col].str[0]
            haploid_df = haploid_df[(haploid_df['GT1'] != '.') &
                                    (haploid_df['GT1'] != np.NaN)]
            haploid_df['GT2'] = 0

        if len(df) > 0:
            df['GT1'] = df[sample_col].str[0]
            df = df[(df['GT1'] != '.') & (df['GT1'] != np.NaN)]
            df['GT1'] = df['GT1'].astype(int)
            df['GT2'] = df[sample_col].str[2]
            df = df[(df['GT2'] != '.') & (df['GT2'] != np.NaN)]
            df['GT2'] = df['GT2'].astype(int)

    # 16th December 2014 not sure this is needed now that
    # get_multiallelic_bases is separate function
    if not single_sample_vcf:
        df['GT_len'] = df[sample_col].str.split(':').str[0].str.len()
        haploid_df = df[df['GT_len'] <= 1]
        haploid_df['phase'] = '-'
        haploid_df = haploid_df[haploid_df['CHROM'].isin(haploid_chromosomes)]

        df = df[~df.index.isin(haploid_df.index)]
        df['phase'] = df[sample_col].str[1]
        df = df[(df['phase'] != '-')]

        del df['GT_len']
        del haploid_df['GT_len']

        if len(df) > 0:
            df['GT1'] = df.apply(get_GT_multisample,
                                 args=[sample_col, 0], axis=1)
            df = df[(df['GT1'] != '.') & (df['GT1'] != np.NaN)]
            df['GT2'] = df.apply(get_GT_multisample,
                                 args=[sample_col, 1], axis=1)
            df = df[(df['GT2'] != '.') & (df['GT2'] != np.NaN)]

        if len(haploid_df) > 0:
            haploid_df['GT1'] = df.apply(get_GT_multisample_vcf_haploid,
                                         args=[sample_col, 0], axis=1)
            haploid_df = haploid_df[(haploid_df['GT1'] != '.') & (haploid_df['GT1'] != np.NaN)]
            haploid_df['GT1'] = haploid_df['GT1'].astype(int)
            haploid_df['GT2'] = 0

    if len(df) > 0:
        if len(haploid_df) > 0:
            df = df.append(haploid_df)
        else:
            pass
    else:
        df = haploid_df

    # FAST PROCESS SIMPLE ALLELE GENOTYPES
    df_simple = df
    if len(df_simple) > 0:
        # get a1 ref alleles
        df_gt1_ref = df_simple[df_simple.GT1.astype(int) == 0][['REF']]
        df_gt1_ref.columns = ['a1']
        # get a2 ref alleles
        df_gt2_ref = df_simple[df_simple.GT2.astype(int) == 0][['REF']]
        df_gt2_ref.columns = ['a2']

        # get a1 alt alleles
        df_gt1_alt = df_simple[df_simple.GT1.astype(int) == 1][['ALT']]
        df_gt1_alt.columns = ['a1']
        # get a2 alt alleles
        df_gt2_alt = df_simple[df_simple.GT2.astype(int) == 1][['ALT']]
        df_gt2_alt.columns = ['a2']

        # merging GT1 allele bases into a single df
        gt1_alleles = pd.concat([df_gt1_ref, df_gt1_alt])
        # del gt1_alleles[0]
        # merging GT2 allele bases into a single df
        gt2_alleles = pd.concat([df_gt2_ref, df_gt2_alt])
        # del gt2_alleles[0]
        # Joining the GT1 and GT2 simple allele bases
        gt1_2_allele_df = gt1_alleles.join(gt2_alleles)

        # print len(df)

        return df.join(gt1_2_allele_df)

    else:
        return pd.DataFrame()
    # print len(df)


def zygosity_fast(df):
    """
    This function quickly assigns zygosity states to
    each variant using set logic
    """

    df_hom_ref = df[(df['a1'] == df['REF']) & (df['a2'] == df['REF'])]
    df_hom_ref['zygosity'] = 'hom-ref'

    df_hom_miss = df[(df['a1'] == '.') & (df['a2'] == '.')]
    df_hom_miss['zygosity'] = 'hom-miss'

    df_het_miss = df[(df['a1'] == '.') | (df['a2'] == '.')]
    df_het_miss['zygosity'] = 'het-miss'

    df_not_miss = df[~df.index.isin(set(df_hom_miss.index) |
                                    set(df_het_miss.index))]

    df_het_alt = df_not_miss[((df_not_miss['a1'] != df_not_miss['REF']) &
                              (df_not_miss['a2'] != df_not_miss['REF'])) &
                             (df_not_miss['a1'] != df_not_miss['a2'])]
    df_het_alt['zygosity'] = 'het-alt'

    df_hom_alt = df_not_miss[(((df_not_miss['a1'] != df_not_miss['REF']) &
                               (df_not_miss['a2'] != df_not_miss['REF']))) &
                             (df_not_miss['a1'] == df_not_miss['a2'])]
    df_hom_alt['zygosity'] = 'hom-alt'

    df_het_ref = df_not_miss[((df_not_miss['a1'] == df_not_miss['REF']) &
                              (df_not_miss['a2'] != df_not_miss['REF'])) |
                             ((df_not_miss['a1'] != df_not_miss['REF']) &
                              (df_not_miss['a2'] == df_not_miss['REF']))]
    df_het_ref['zygosity'] = 'het-ref'

    df_zygosity = pd.concat([df_hom_ref, df_hom_miss,
                             df_het_miss, df_het_ref,
                             df_het_alt, df_hom_alt])
    # return df_zygosity
    assert len(df_zygosity) == len(df)
    return df_zygosity


def vartype_map(ref_alt_bases):
    """
    This function assigns the following vartypes to the
    allele specified by allele_base_col: snp, mnp, ins, del, indel or SV
    """
    ref, alt = ref_alt_bases
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
    hom_ref = df[df['GT'].isin(['0|0', '0/0', '0'])]
    hom_ref['hom_ref'] = 1
    return hom_ref.groupby(level=[0, 1, 2, 3])['hom_ref'].aggregate(np.sum)


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
        temp_group.replace(to_replace='.', value='', inplace=True)

        temp_group_data = pd.DataFrame.from_records(
            list(temp_group.str.split(':')))
        temp_group_data.index = temp_group.index
        temp_group_data.columns = name.split(':')
        temp_group_data.replace(to_replace='.', value='', inplace=True)

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
        print 'No Annotations generated, please check for excessive missing values'
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
    split_size = row_count / split_level
    split_df = []
    for n, i in enumerate(range(0, row_count, split_size)):
        if n + 1 == split_level:
            split_df.append(df.ix[df.index[i:]])
            break
        else:
            split_df.append(df.ix[df.index[i: i + split_size]])
    return split_df


def mp_variant_annotations(df_mp, df_split_cols, df_sampleid,
                           drop_hom_ref, n_cores=1):
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

    pool = mp.Pool(int(n_cores))
    # tasks = np.array_split(df_mp.copy(), int(n_cores))  #breaks with older
    # pandas/numpy
    tasks = df_split(df_mp.copy(), int(n_cores))
    tasks = [
        [pd.DataFrame(t), df_split_cols, df_sampleid, drop_hom_ref] for t in tasks]
    results = []
    del df_mp
    gc.collect()
    r = pool.map_async(
        process_variant_annotations, tasks, callback=results.append)
    r.wait()
    pool.close()
    pool.join()
    pool.terminate()

    return pd.concat(results[0])


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

    df['multiallele'] = df.ALT.str.count(',')
    multidf = df[df['multiallele'] > 0]

    while len(df) + len(multidf) > 0:

        df = df[~df.index.isin(multidf.index)]
        # print len(multidf), 'multidf rows'

        if len(multidf) > 0:
            multidf = get_multiallelic_bases(
                multidf, sample_name, single_sample_vcf=False)

        # print 'single alleles', len(df)

        df = get_biallelic_bases(df, sample_name)

        if len(multidf) > 0:
            df = df.append(multidf)

        df = zygosity_fast(df)

        df['vartype1'] = map(vartype_map, df[['REF', 'a1']].values)
        df['vartype2'] = map(vartype_map, df[['REF', 'a2']].values)

        df.set_index(['CHROM', 'POS', 'REF', 'ALT', 'sample_ids'],
                     inplace=True)

        # df.sortlevel(level=['CHROM','POS','REF','ALT','sample_ids'],inplace=True)
        # #sorting biallelic and multiallele variants

        # print 'before parse_single_genotype_data', len(df)

        # df = df.join( parse_single_genotype_data(df, sample_name, split_cols=split_columns), how='left' )
        df['GT'] = df['sample_genotypes']
        del df[sample_name]
        if 'FORMAT' in df.columns:
            del df['FORMAT']

        df.reset_index(level=4, inplace=True, drop=True)
        df.set_index('GT', inplace=True, drop=False, append=True)
        # print df
        return df

    return pd.DataFrame()
