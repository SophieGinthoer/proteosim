import numpy as np
import matplotlib.pyplot as plt



amino_acid_mass_dalton = {
    'A': 71.08, 'R': 156.19, 'N': 114.10, 'D': 115.09,
    'C': 103.15, 'E': 129.12, 'Q': 128.13, 'G': 57.05,
    'H': 137.14, 'I': 113.16, 'L': 113.16, 'K': 128.17,
    'M': 131.19, 'F': 147.18, 'P': 97.12, 'S': 87.08,
    'T': 101.11, 'W': 186.21, 'Y': 163.18, 'V': 99.13,
}



def calculate_mol_mass(peptide_seq, amino_acid_mass_dict):
    """
    Calculate molecular mass of peptide by summing its aa masses.

    Parameters
    ----------
    peptide_seq : str
        Amino-acid sequence of the peptide.
    amino_acid_mass_dict : dict of str to float
        Maps each amino acid to its mass in Dalton.

    Returns
    -------
    dict of str to float
    
    """
    mass = 0

    for aa in peptide_seq:
        mass += amino_acid_mass_dict[aa]

    return {peptide_seq : mass}



def calculate_mol_mass_collection(peptides, amino_acid_mass_dict):
    """
    Calculate molecular masses of list of peptides.
    
    Parameters
    ----------
    peptides : list of str
    amino_acid_mass_dict : dict of str to float
    Returns
    -------
    dict of str to float
    """
    mass_map = {}
    for pep in peptides:
        mass_map.update(calculate_mol_mass(pep, amino_acid_mass_dict))
    return mass_map



def calculate_mz_collection(peptide_mass_map, charge, proton_mass=1.007):
    """
    Convert peptide masses to m/z values for given charge.

    Parameters
    ----------
    peptide_mass_map : dict of str to float
    charge : int
    proton_mass : float, optional

    Returns
    -------
    dict of str to float
    """
    mz_map = {}
    for pep, mass in peptide_mass_map.items():
        mz_map[pep] = (mass + charge * proton_mass)  / charge
    return mz_map



def plot_spectrum(mz_values, random_count_range=(0, 30000), seed=42):
    """
    Plot simulated MS with random intensities.

    Parameters
    ----------
    mz_values : list of float
        m/z values.
    random_count_range : tuple of int, optional
        range for random intensities (default (0, 30000).
    seed : int, optional
        random seed for reproducible results (default 42).

    Returns
    -------
    None 
        bar Chart.
    """
    mz_values = list(mz_values)
    np.random.seed(seed)
    intensities = np.random.randint(random_count_range[0], random_count_range[1], size=len(mz_values))

    plt.bar(mz_values, intensities, width=1)
    plt.xlabel("m/z")
    plt.ylabel("intensity [counts]")
    plt.show()



def fragment_peptide(peptide):
    """
    Generate b- & y-ion fragments of peptide.

    Parameters
    ----------
    peptide : str
        peptide's aa seq.
    Returns
    -------
    list of str
        b-ions (N-term), y-ions (C-term).
    """
    b_ions = [peptide[:i] for i in range(1, len(peptide))]
    y_ions = [peptide[i:] for i in range(len(peptide)-1, 0, -1)]
    return b_ions + y_ions