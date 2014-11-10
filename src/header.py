import os
def get_header(vcf_path):
    if vcf_path.endswith('.gz'):
        return os.popen('tabix -H ' + vcf_path).readlines()
    header_lines = os.popen('head -5000 ' + vcf_path).readlines()
    return [l for l in header_lines if l.startswith('#')]

