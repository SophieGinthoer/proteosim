from pyteomics import achrom

from proteosim import (
    predict_lc_retention_times,
    select_retention_time_window,
)


def test_predict_lc_retention_times():

    rc = achrom.RCs_guo_ph7_0["aa"]
    A, W, K, L, M, G = rc["A"], rc["W"], rc["K"], rc["L"], rc["M"], rc["G"]

    peptides = ["AWWKLMMMM", "WWWWLLLLAAA", "AAAA", "WWWWWWALGGG", "KKKLKL"]
    
    expected = {
        "AWWKLMMMM" : round(1*A + 2*W + 1*K + 1*L + 4*M, 2),
        "WWWWLLLLAAA" : round(4*W + 4*L + 3*A, 2),
        "AAAA" : round(4*A, 2),
        "WWWWWWALGGG" : round(6*W + 1*A + 1*L + 3*G, 2),
        "KKKLKL" : round(4*K + 2*L, 2),
    }

    actual = predict_lc_retention_times(peptides)
    assert actual == expected


def test_select_retention_time_window():
    peptide_rt_map = {"AAA":1, "LAA":2, "WWW":3, "AWW":4, "A":5}
    selected = select_retention_time_window(peptide_rt_map, lower_ret_time=3, upper_ret_time=4)

    assert selected == ["WWW", "AWW"]