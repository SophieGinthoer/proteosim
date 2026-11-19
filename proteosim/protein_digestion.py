import re

enzyme_cleavage_patterns = {
    'LysC': r'(?<=K)',
    'LysN': r'(?=K)',
    'ArgC': r'(?<=R)',
    'Trypsin': r'(?<=[KR])(?!P)',
}

enzyme_cleavage_patterns



def digest_protein_sequence(protein_seq, cleave_pattern, min_pep_len=5, max_pep_len=30):
    """
    Digest a protein in silico & filter peptides by length.

    Parameters
    ----------
    protein_seq : str
        Amino-acid sequence of the protein.
    cleave_pattern : str
        Regular expression defining protease cleavage sites.
    min_pep_len : int, optional
        Minimum peptide length (default 5).
    max_pep_len : int, optional
        Maximum peptide length (default 30).

    Returns
    -------
    list of str
        Peptides within the length window.
    """
    if min_pep_len > max_pep_len:
        raise ValueError("minimum must be smaller than maximum")

    peptides = re.split(cleave_pattern, protein_seq)
    filtered_peptides = [pep for pep in peptides if min_pep_len <= len(pep) <= max_pep_len]
    return filtered_peptides



def digest_protein_collection(protein_map, cleave_pattern, min_pep_len=5, max_pep_len=30):
    """
    Digest every protein in a collection in silico.

    Parameters
    ----------
    protein_map : dict of str to str
    cleave_pattern : str
    min & max peptide length default

    Returns
    -------
    dict of str to list of str
    """
    digested = {} # empty dictionary to collect results

    for protein_id, sequence in protein_map.items(): # loop over every prot, get ID & seq as pair
        digested[protein_id] = digest_protein_sequence(sequence, cleave_pattern, min_pep_len, max_pep_len) # digest this seq & store peptide list under same ID
        # turn into two lines, makes this more comprehensible
    
    return digested # return new dictionary, protein_id -> list of peptides



def compute_sequence_coverage(protein_seq, peptides):
    """
    Compute percentage of protein seq covered by peptides.
    Each residue counted once, even if several peptides overlap it.

    Parameters
    ----------
    protein_seq : str
        AA seq of protein
    peptides : list of str
        peptides to map onto prot seq

    Returns
    -------
    float
    """
    
    if len(protein_seq) == 0:
        return 0.0 # empty seq has no residue to cover & stops division by zero at end

    covered = [False] * len(protein_seq) # one entry per residue, not yet covered by peptide

    for pep in peptides: # map every peptide onto prot seq
        start = protein_seq.find(pep) # find first position where peptide occurs (negative means not found)

        while start != -1: # repeat as long as another occurence found
            for i in range(start, start + len(pep)): # set all positions of this peptide to True
                covered[i] = True

            start = protein_seq.find(pep, start+1) # look again after pos of last hit -> catch repeated peptides

    return sum(covered) / len(protein_seq) * 100.0 # True counts as 1 -> number of coverd residues; percentage
