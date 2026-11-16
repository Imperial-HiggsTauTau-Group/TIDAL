import json


def reweighting_string(era, ind):

    dms = {
        "sel_mupi": 0,
        "sel_murho": 1,
        "sel_mua11pr": 2,
        "sel_mua1": 10,
    }

    with open("Draw/residual_correcting/corrections.json", "r") as f:
        corrections = json.load(f)

    string = "*( 1.0 + (genPartFlav_X == 5)*( "

    for cat, dm in dms.items():
        string += f"(decayModePNet_X == {dm})*( "

        # Putting the corrections in cumulative form
        corr = corrections[era][cat]
        c20 = f"{(corr['20-25'] - 1.0):.4f}"
        c25 = f"{(corr['25-30'] - corr['20-25']):.4f}"
        c30 = f"{(corr['30-35'] - corr['25-30']):.4f}"
        c35 = f"{(corr['35-40'] - corr['30-35']):.4f}"
        c40 = f"{(corr['40-45'] - corr['35-40']):.4f}"
        c45 = f"{(corr['45-50'] - corr['40-45']):.4f}"
        c50 = f"{(1.0 - corr['45-50']):.4f}"

        string += (
            f"({c20})*(pt_X >= 20)"
            f" + ({c25})*(pt_X >= 25)"
            f" + ({c30})*(pt_X >= 30)"
            f" + ({c35})*(pt_X >= 35)"
            f" + ({c40})*(pt_X >= 40)"
            f" + ({c45})*(pt_X >= 45)"
            f" + ({c50})*(pt_X >= 50)"
        )

        string += " ) + "

    string = string[:-3] + " ) )"
    
    return string.replace("X", ind)


def test():
    print(reweighting_string("I", "1"))


if __name__ == "__main__":
    test()
