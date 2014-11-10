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





