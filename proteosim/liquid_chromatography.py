import matplotlib.pyplot as plt
from pyteomics import achrom



def predict_lc_retention_times(peptides):
    """
    Predict relative RT for list of peptides.
    Uses achrom.calculate_RT with data from Guo et al. at pH 7.

    Parameters
    ----------
    peptides : list of str
        Peptide sequences.

    Returns
    -------
    dict of str to float
        Maps each peptide to its corresponding relative RT (rounded to 2 decimals).
    """

    rt_map = {} # empty dict

    for pep in peptides: # predict RT for every pep
        rt = achrom.calculate_RT(pep, achrom.RCs_guo_ph7_0) # sum of aa retention coefficients
        rt_map[pep] = round(rt, 2) # store RT under pep as key

    return rt_map # dictionary peptide : RT



def plot_retention_time(retention_times, resolution=30):
    """
    Histogram where x is retention time [min] and y is approx. I by counting # of peptides within a bin.

    Parameters
    ----------
    retention_times : list of float
        Predicted rt of peps
    resolution : int, optional
        # of histogram bins (default=30)
    Returns
    -------
    None
        Shows histogram
    """
    plt.hist(list(retention_times), bins=resolution)
    plt.xlabel("retention time [min]")
    plt.ylabel("intensity [counted peptides]")
    plt.show()



def select_retention_time_window(peptide_rt_map, lower_ret_time, upper_ret_time):
    """
    Select peptides that elute within given RT.

    Parameters
    ----------
    peptide_rt_map : dict or str to float
    lower_ret_time : float
    upper_ret_time : float
    Returns
    -------
    list of str
    """
    selected_peptides = []
    for pep, rt in peptide_rt_map.items():
        if lower_ret_time <= rt <= upper_ret_time:
            selected_peptides.append(pep)
    return selected_peptides