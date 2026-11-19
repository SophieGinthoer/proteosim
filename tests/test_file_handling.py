from proteosim import read_fasta

def test_read_fasta(tmp_path):
    # Create a small dummy FASTA with a wrapped sequence and blank lines
    dummy_fasta = """>sp|W|Dummy protein Tryptophan
Trypto
phan

>sp|A|Dummy protein Alanin
Alanin


"""
    fasta_path = tmp_path / "dummy_proteins.fasta"
    fasta_path.write_text(dummy_fasta)

    protein_map = read_fasta(fasta_path)

    # Keys are the UniProt IDs, values the joined sequences
    assert protein_map["W"] == "Tryptophan"
    assert protein_map["A"] == "Alanin"
    assert len(protein_map) == 2