def read_fasta(filepath):
    """
    Read a protein FASTA file into a dictionary of UniProt IDs and sequences.

    Parameters
    ----------
    filepath : str
        Path to the FASTA file.

    Returns
    -------
    dict of str to str
        Maps each UniProt protein ID to its amino-acid sequence.
    """
    protein_map = {}
    current_id = None
    current_sequence = []

    with open(filepath, 'r') as fasta_handle:
        for line in fasta_handle:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith('>'):
                if current_id is not None:
                    protein_map[current_id] = ''.join(current_sequence)
                    current_sequence = []
                current_id = stripped.split('|')[1]
            else:
                current_sequence.append(stripped)

    if current_id is not None:
        protein_map[current_id] = ''.join(current_sequence)

    return protein_map

