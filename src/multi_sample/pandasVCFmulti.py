import os,sys,gzip
import pandas as pd
from variantAnnotations import *
from Vcf_metadata import *



class Vcf(object):
    '''
        Loads in a vcf file, aware of gzipped files.
        testing
        '''
    
    def __init__(self, filename, sample_id='', cols='', chunksize=5000):
        
        #Header
        header_parsed = Vcf_metadata(filename)
        self.header_df = self.get_header_df(header_parsed.header)  #header parsed into key/values dataframe
        
        
        
        #Sample IDs
        self.samples = list(self.header_df.ix['SampleIDs'])[0]
        self.sample_id = sample_id
        
        
        #Columns
        self.all_columns = list(self.header_df.ix['ColumnHeader'])[0]
        self.FORMAT = self.all_columns[8]
        
        
        assert len(set(cols) & set(['#CHROM', 'POS', 'REF', 'ALT', 'FORMAT'])), "cols requires the following columns: ['#CHROM', 'POS', 'REF', 'ALT', 'FORMAT']"
        self.cols = cols
        if len(cols) > 0:
            self.usecols = [c for c in self.all_columns if c in cols]
            if self.sample_id == 'all':
                self.usecols.extend(self.samples)
                self.sample_id = self.samples
            elif: len(sample_id) > 1:
                self.usecols.extend(self.sample_id)
            else:
                self.usecols.extend([self.sample_id])
        else:
            self.usecols = [s for s in self.cols if s not in self.samples]
            if sample_id == 'all':
                self.usecols.extend(self.samples)
                self.sample_id = self.samples
            else:
                self.usecols.extend(self.sample_id)
        
        #Open pandas chunk object (TextReader)
        self.chunksize = chunksize
        print self.usecols
        assert False
        self.vcf_chunks = pd.read_table(filename, sep="\t", compression=header_parsed.compression, skiprows=(len(self.header_df)-2), usecols=self.usecols, chunksize=chunksize)
    
    
    
    
    def get_header_df(self, header_txt):
        '''
            Parses header into pandas DataFrame
            '''
        key_value_header = [i.replace('##','').replace('\n','').split('=',1) for i in header_txt if '##' in i]
        key_value_header.append(['SampleIDs',header_txt[-1].rstrip('\n').split('\t')[9:]])
        key_value_header.append(['ColumnHeader', header_txt[-1].rstrip('\n').split('\t')])
        header_df =  pd.DataFrame.from_records(key_value_header)
        header_df.set_index(0,inplace=True)
        header_df.index.name = 'header_keys'
        header_df.columns = ['header_values']
        return header_df
    
    
    
    
    def get_vcf_df_chunk(self):
        '''
            This function iterates through the VCF files using the user-defined
            chunksize (default = 5000 lines).
            '''
        
        self.df = self.vcf_chunks.get_chunk()
        self.df.columns = [c.replace('#', '') for c in self.usecols]
        self.df.set_index(['CHROM', 'POS', 'REF', 'ALT'], inplace=True, drop=False)
        
        self.df_bytes = self.df.values.nbytes + self.df.index.nbytes + self.df.columns.nbytes
        
        return 0
    
    
    def add_variant_annotations(self, split_columns='', verbose=False, inplace=False):
        
        '''
            This function adds the following annotations for each variant:
            multiallele, phase, a1, a2, GT1, GT2, vartype1, vartype2, zygosity,
            and parsed FORMAT values, see below for additional information.
            
            Parameters
            --------------
            
            split_columns: dict, optional
            key:FORMAT id value:#fields expected
            e.g. {'AD':2} indicates Allelic Depth should be
            split into 2 columns.
            
            verbose: bool, default=False
            This will describe how many missing variants were dropped
            
            inplace: bool, default=False
            This will replace the sample_id column with parsed columns,
            and drop the FORMAT field.  If True, this will create an
            additional dataframe, df_annot, to the Vcf object composed of
            the parsed columns (memory intensive)
            
            
            
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
        
        if self.sample_id not in self.df.columns:
            print 'Sample genotype column not found, add_variant_annotations can only be called if the FORMAT and sample genotype columns are available?'
            return 1
        
        if verbose:
            print 'input variant rows:', len(self.df)
        
        self.df.drop_duplicates(inplace=True)
        
        var_counts = len(self.df)
        
        self.df = self.df[self.df[self.sample_id]!='.']  #dropping missing genotype calls
        
        if verbose:
            
            print 'dropping',var_counts - len(self.df), 'variants with genotype call == "." '
            print 'current variant rows:', len(self.df)
        
        
        if inplace:
            if 'POS' not in self.df.columns:
                self.df.reset_index(inplace=True)
                self.df.set_index(['CHROM', 'POS', 'REF', 'ALT'],inplace=True, drop=False)
            self.df = get_vcf_annotations(self.df, self.sample_id, split_columns)
        else:
            self.df_annot = get_vcf_annotations(self.df, self.sample_id, split_columns)
            self.df.set_index(['CHROM', 'POS', 'REF', 'ALT'],inplace=True)
        
        return 0


