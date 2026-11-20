from proteosim import (
    enzyme_cleavage_patterns,
    digest_protein_sequence,
    digest_protein_collection,
    compute_sequence_coverage,
)


def test_digest_protein_sequence():
    """Test digest_protein_sequence with different proteases and length filters."""
    trypsin = enzyme_cleavage_patterns['Trypsin']
    lysc = enzyme_cleavage_patterns['LysC']
    argc = enzyme_cleavage_patterns['ArgC']

    dummy_seq_1 = "GGGKLLALLLRMMM"
    dummy_seq_2 = "GGGGKPLLLLLLR"
    dummy_seq_3 = "GKLLLLLLK"
    dummy_seq_4 = "RKLLLLGGG"

    # trypsin cuts after K and R
    assert digest_protein_sequence(dummy_seq_1, trypsin, min_pep_len=3) == ["GGGK", "LLALLLR", "MMM"]
    # default filter (5-30): 'GGGK' & 'MMM' too short
    assert digest_protein_sequence(dummy_seq_1, trypsin) == ["LLALLLR"]
    # proline -> no cut between K & P
    assert digest_protein_sequence(dummy_seq_2, trypsin) == ["GGGGKPLLLLLLR"]
    # LysC cuts after K only
    assert digest_protein_sequence(dummy_seq_3, lysc, min_pep_len=2) == ["GK", "LLLLLLK"]
    # ArgC cuts after R only
    assert digest_protein_sequence(dummy_seq_4, argc, min_pep_len=1) == ["R", "KLLLLGGG"]


def test_digest_protein_collection():
    trypsin = enzyme_cleavage_patterns['Trypsin']

    dummy_proteins = {
        "protA" : "KLWWWWWAAAAAR",
        "protB" : "RPLLLLLAAA",
        "protC" : "AAAAA",
        "protD" : "KRAAAAWLM"
    }
    expected_output = {
        "protA" : ["LWWWWWAAAAAR"],
        "protB" : ["RPLLLLLAAA"],
        "protC" : ["AAAAA"],
        "protD" : ["AAAAWLM"],
    }

    assert digest_protein_collection (dummy_proteins, trypsin) == expected_output


def test_compute_sequence_coverage():
    dummy_prot_seq = "WWRKLLLLAAGGG"
    dummy_peps = ["AAG", "WRKLL", "LLLLLLLL", "LAAAG"]

    coverage = compute_sequence_coverage(dummy_prot_seq, dummy_peps)
    assert coverage == 8/13 * 100.0