import os


class VCFMetadata(object):
    """This class parses a VCF header into a pandas dataframe object.

    It recognizes gzip and uncompressed file formats.
    This function assumes the header does not extent past 5000 lines.
    """

    def __init__(self, filename):
        if filename.endswith('.gz'):
            self.compression = 'gzip'
            if filename + '.tbi' in os.listdir(os.path.split(filename)[0]):
                header_lines = os.popen('tabix -H ' + filename).readlines()
                self.header = [l.replace('#CHROM', 'CHROM')
                               for l in header_lines if l.startswith('#')]
            os.system('tabix -p vcf ' + filename)
            header_lines = os.popen('tabix -H ' + filename).readlines()
            self.header = [l for l in header_lines if l.startswith('#')]

        else:
            self.compression = 'infer'
            header_lines = os.popen('head -5000 ' + filename).readlines()
            self.header = [l for l in header_lines if l.startswith('#')]
