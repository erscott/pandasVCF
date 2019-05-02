import pandas as pd

from variant_annotations import process_variant_annotations, mp_variant_annotations
from vcf_metadata import VCFMetadata


class VCF(object):
    """Loads in a vcf file, aware of gzipped files.


    Parameters
    --------------------------------------
    filename: str, required
        path to vcf file

    sample_id: str or list, default='all'
        specifies the sample column ids to read and parse

        'all' means all sample columns

         can use a str (e.g. 'NA12878')
         or
         can use a list (e.g. ['NA12878', 'NA12877']


    cols: list, default ['#CHROM', 'POS', 'REF', 'ALT', 'FORMAT']
        specifies the VCF column names, EXCEPT SAMPLE COLS, to read and parse

        Must include ['#CHROM', 'POS', 'REF', 'ALT', 'FORMAT']

        Additional columns such as QUAL, FILTER, INFO will be accepted
            e.g. ['#CHROM', 'POS', 'REF', 'ALT', 'FORMAT', 'INFO', 'QUAL']

    chunksize: int, default=5000
        specifies the number of VCF lines to read and parse in 1 chunk

        Note using a large chunksize with large n_cores requires LOTS OF RAM

        requires ~40 seconds to parse 1000 rows with 2500 samples


    Methods
    -----------------------------------------
    get_vcf_df_chunk
        returns VCF pandasDF with chunksize


    add_variant_annotations
        Annotates each variant
        See docstring for details



    Returns VCF Obj with following attributes
    -----------------------------------------
    header_df: pandas df
        VCF header as a pandas df

    samples: list
        sample column IDs

    all_columns: list
        all sample column IDs in VCF

    vcf_chunks: pandas.io.parsers.TextFileReader chunk
        VCF chunk
        Access to chunk provided by get_vcf_df_chunk()

    df: pandas DF
        Index: CHROM, POS, REF, ALT
        Columns: CHROM, POS, REF, ALT, SAMPLE(S) +/- {QUAL, FILTER, INFO if specified}


    """

    def __init__(self, filename, sample_id='all',
                 cols=['#CHROM', 'POS', 'REF', 'ALT', 'FORMAT'],
                 chunksize=5000):

        # Header
        header_parsed = VCFMetadata(filename)
        # header parsed into key/values dataframe
        self.header_df = self.get_header_df(header_parsed.header)
        # Sample IDs
        self.samples = list(self.header_df.loc['SampleIDs'])[0]
        self.sample_id = self.get_sample_ids(sample_id)

        self.set_cols(cols)

        self.set_dtypes()
        
        # Open pandas chunk object (TextReader)
        self.chunksize = chunksize
        self.vcf_chunks = pd.read_csv(filename, sep="\t",
                                        compression=header_parsed.compression,
                                        skiprows=(len(self.header_df) - 2),
                                        usecols=self.usecols,
                                        chunksize=chunksize,
                                        dtype=self.vcf_dtypes)


    def get_header_df(self, header_txt):
        """Parses header into pandas DataFrame"""
        try:
            key_value_header = [i.replace('##', '').replace(
                '\n', '').split('=', 1) for i in header_txt if '##' in i]
            key_value_header.append(
                ['SampleIDs', header_txt[-1].rstrip('\n').split('\t')[9:]])
            key_value_header.append(
                ['ColumnHeader', header_txt[-1].rstrip('\n').split('\t')])
            header_df = pd.DataFrame.from_records(key_value_header)
            header_df.set_index(0, inplace=True)
            header_df.index.name = 'header_keys'
            header_df.columns = ['header_values']
            return header_df
        except IndexError:
            print("VCF header parsing failed, "
                  "this may be due to the use of "
                  "tabix version 1.2.x, please upgrade to tabix 1.3 or greater")
            return

    def get_sample_ids(self, sample_id):
        """
        Identifies and stores sample_id(s)
        """
        if sample_id == 'all':
            return self.samples[:]
        else:
            if type(sample_id) == str:
                return [sample_id]
            else:
                return sample_id

    def set_cols(self, cols):
        # Columns
        self.all_columns = list(self.header_df.ix['ColumnHeader'])[0]
        self.FORMAT = self.all_columns[8]

        assert len(set(cols) & set(['#CHROM', 'POS', 'REF', 'ALT', 'FORMAT'])) > 4, "cols requires the following columns: ['#CHROM', 'POS', 'REF', 'ALT', 'FORMAT']"
        self.cols = cols
        if len(cols) > 0:  # columns specified
            self.usecols = [c for c in self.all_columns if c in cols]

            if len(self.sample_id) > 0:
                self.usecols.extend(self.sample_id)
                # print self.usecols
            else:
                assert False, 'no sample IDs'
        else:  # columns not specified
            self.usecols = [s for s in self.cols if s not in self.samples]
            self.usecols.extend(self.sample_id)

    def set_dtypes(self):
        self.vcf_dtypes = {'CHROM':'category',
                           'POS':'int32',
                           'REF':'category',
                           'ALT':'category',
                           'FORMAT':'category',
                           'QUAL':'int8',
                           'FILTER':'category'}

    def get_vcf_df_chunk(self):
        """
        This function iterates through the VCF files using the user-defined
        chunksize (default = 5000 lines).
        """
        try:
            self.df = self.vcf_chunks.get_chunk()
            self.stopIteration = False
        except StopIteration:
            self.stopIteration = True
            print("End of File Reached")
            # self.df = None
            return 1
        self.df.drop_duplicates(inplace=True)  # dropping duplicate rows
        self.df.columns = [c.replace('#', '') for c in self.usecols]
        self.df['CHROM'] = self.df['CHROM'].astype(str).str.replace('chr', '').astype('category')
        self.df.set_index(
            ['CHROM', 'POS', 'REF', 'ALT'], inplace=True, drop=False)

        self.df_bytes = self.df.values.nbytes + \
            self.df.index.nbytes + self.df.columns.nbytes

        return 0

    def add_variant_annotations(self, split_columns='', verbose=False,
                                inplace=False, drop_hom_ref=True,
                                n_cores=1):
        """
        This function adds the following annotations for each variant:
        multiallele, phase, a1, a2, GT1, GT2, vartype1, vartype2, zygosity,
        and parsed FORMAT values, see below for additional information.

        Parameters
        --------------

        split_columns: dict, optional
            key:FORMAT id value:#fields expected
            e.g. {'AD':2} indicates Allelic Depth should be
            split into 2 columns.

        drop_hom_ref: bool, default=True
            This will drop homozygous reference genotype calls from
            the long dataframe.  As most calls in a multisample vcf
            are homozygous reference, this will reduce memory requirements
            dramatically.

        verbose: bool, default=False
            This will describe how many missing variants were dropped

        inplace: bool, default=False
            This will replace the sample_id column with parsed columns,
            and drop the FORMAT field.  If True, this will create an
            additional dataframe, df_annot, to the VCF object composed of
            the parsed columns (memory intensive)

        n_cores: int, default=1
            specifies the number of cpus to use during variantAnnotation

            Note using a large chunksize with large n_cores requires LOTS OF RAM

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

        if self.stopIteration:
            print('End of File Reached')
            return 1

        self.drop_hom_ref = drop_hom_ref

        df_vcf_cols = self.df[list(set(self.df.columns)
                                   - {'CHROM', 'POS', 'REF', 'ALT', 'FORMAT'}
                                   - set(self.sample_id))]

        self.df = self.df.reset_index(drop=True)

        if n_cores==1:
            if inplace:
                self.df = process_variant_annotations(self.df,
                                                       split_columns=split_columns,
                                                       sample_id=self.sample_id,
                                                       drop_hom_ref=drop_hom_ref)
                # joining QUAL, FILTER, and/or INFO columns
            else:
                self.df_annot = process_variant_annotations(self.df,
                                                             split_columns=split_columns,
                                                             sample_id=self.sample_id,
                                                             drop_hom_ref=drop_hom_ref)
        else:
            if inplace:
                self.df = mp_variant_annotations(self.df, 
                                                 n_cores=n_cores,
                                                 df_split_cols=split_columns,
                                                 df_sampleid=self.sample_id,
                                                 drop_hom_ref=drop_hom_ref)
            else:
                self.df_annot = mp_variant_annotations(self.df, 
                                                 n_cores=n_cores,
                                                 df_split_cols=split_columns,
                                                 df_sampleid=self.sample_id,
                                                 drop_hom_ref=drop_hom_ref)
        if inplace:
            self.df = self.df.set_index(['CHROM', 'POS', 'REF', 'ALT'])
        else:
            self.df_annot = self.df_annot.set_index(['CHROM', 'POS', 'REF', 'ALT'])
        return 0








