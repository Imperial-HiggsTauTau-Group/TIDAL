import json


def reweighting_string(era, cat, ind):

    dm = {
        "sel_mupi": 0,
        "sel_murho": 1,
        "sel_mua11pr": 2,
        "sel_mua1": 10,
    }

    with open("Draw/residual_correcting/corrections.json", "r") as f:
        corrections = json.load(f)
    corr = corrections[era][cat]
    
    reweighting_str = (
        "( 1.0*(pt_X < 20)"
        + f" + {corr['20-25']:.4f}*(pt_X >= 20 && pt_X < 25)"
        + f" + {corr['25-30']:.4f}*(pt_X >= 25 && pt_X < 30)"
        + f" + {corr['30-35']:.4f}*(pt_X >= 30 && pt_X < 35)"
        + f" + {corr['35-40']:.4f}*(pt_X >= 35 && pt_X < 40)"
        + f" + {corr['40-45']:.4f}*(pt_X >= 40 && pt_X < 45)"
        + f" + {corr['45-50']:.4f}*(pt_X >= 45 && pt_X < 50)"
        + " + 1.0*(pt_X >= 50) )"
    )

    full_mask = f"*( (genPartFlav_X != 5) + (genPartFlav_X == 5)*(decayModePNet_X != {dm[cat]}) + (genPartFlav_X == 5)*(decayModePNet_X == {dm[cat]})*{reweighting_str} )"
    
    return full_mask.replace("X", ind)


def test():
    print(reweighting_string("I", "sel_mupi", "1"))


if __name__ == "__main__":
    test()
