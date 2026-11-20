from proteosim import (
    calculate_mol_mass,
    calculate_mol_mass_collection,
    calculate_mz_collection,
    fragment_peptide,
)


def test_calculate_mol_mass():
    aa_mass_dict = {"A" : 12, "D" : 3, "R" : 5}
    A, D, R = aa_mass_dict["A"], aa_mass_dict["D"], aa_mass_dict["R"]

    test_peptide = "AADR"
    expected_mass = {"AADR" : 2*A+D+R}
    calculated_mass = calculate_mol_mass(test_peptide, aa_mass_dict)
    assert calculated_mass == expected_mass

    assert calculate_mol_mass("DAR", aa_mass_dict)["DAR"] == calculate_mol_mass("RAD", aa_mass_dict)["RAD"]

    error_peptide = "CCC"
    try: # run code that should fail
        calculate_mol_mass(error_peptide, aa_mass_dict) # this should raise KeyError
        assert False, "expected KeyError for unexpected aa" # only reached if no error detected -> fail
    except KeyError: # as expected: KeyError
        pass # no action, test passes



def test_calculate_mol_mass_collection():
    aa_mass_dict = {"A" : 3, "R" : 55, "L" : 2}
    A, L, R = aa_mass_dict["A"], aa_mass_dict["L"], aa_mass_dict["R"]
    peptides = ["ARA", "ALL"]
    expected = {"ARA" : 2*A + R, 
                "ALL" : A + 2*L}
    actual = calculate_mol_mass_collection(peptides, aa_mass_dict)

    assert actual == expected



def test_calculate_mz_collection():
    peptide_mass_map = {"ARL":3, "WLR":55}
    actual = calculate_mz_collection(peptide_mass_map, charge=2)
    expected = {
        "ARL" : (3 + 2*1.007) / 2,
        "WLR" : (55 + 2*1.007) / 2,
    }

    assert actual == expected



def test_fragment_peptide():
    peptide = "ABCDE"
    expected = ["A", "AB", "ABC", "ABCD", "E", "DE", "CDE", "BCDE"]
    actual = fragment_peptide(peptide)

    assert set(actual) == set(expected)
    assert peptide not in actual  # which string should not be in the expected?