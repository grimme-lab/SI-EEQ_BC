"""
Python script to plot the partial charge of an atom in a
two-atom-molecule depending on the distance between the atoms.
"""

import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns  # type: ignore
import pandas as pd

# Constants
CM_TO_INCH = 1 / 2.54
VERBOSITY = 1
PSE = {
    0: "X",
    1: "H",
    2: "He",
    3: "Li",
    4: "Be",
    5: "B",
    6: "C",
    7: "N",
    8: "O",
    9: "F",
    10: "Ne",
    11: "Na",
    12: "Mg",
    13: "Al",
    14: "Si",
    15: "P",
    16: "S",
    17: "Cl",
    18: "Ar",
    19: "K",
    20: "Ca",
    21: "Sc",
    22: "Ti",
    23: "V",
    24: "Cr",
    25: "Mn",
    26: "Fe",
    27: "Co",
    28: "Ni",
    29: "Cu",
    30: "Zn",
    31: "Ga",
    32: "Ge",
    33: "As",
    34: "Se",
    35: "Br",
    36: "Kr",
    37: "Rb",
    38: "Sr",
    39: "Y",
    40: "Zr",
    41: "Nb",
    42: "Mo",
    43: "Tc",
    44: "Ru",
    45: "Rh",
    46: "Pd",
    47: "Ag",
    48: "Cd",
    49: "In",
    50: "Sn",
    51: "Sb",
    52: "Te",
    53: "I",
    54: "Xe",
    55: "Cs",
    56: "Ba",
    57: "La",
    58: "Ce",
    59: "Pr",
    60: "Nd",
    61: "Pm",
    62: "Sm",
    63: "Eu",
    64: "Gd",
    65: "Tb",
    66: "Dy",
    67: "Ho",
    68: "Er",
    69: "Tm",
    70: "Yb",
    71: "Lu",
    72: "Hf",
    73: "Ta",
    74: "W",
    75: "Re",
    76: "Os",
    77: "Ir",
    78: "Pt",
    79: "Au",
    80: "Hg",
    81: "Tl",
    82: "Pb",
    83: "Bi",
    84: "Po",
    85: "At",
    86: "Rn",
    87: "Fr",
    88: "Ra",
    89: "Ac",
    90: "Th",
    91: "Pa",
    92: "U",
    93: "Np",
    94: "Pu",
    95: "Am",
    96: "Cm",
    97: "Bk",
    98: "Cf",
    99: "Es",
    100: "Fm",
    101: "Md",
    102: "No",
    103: "Lr",
    104: "Rf",
    105: "Db",
    106: "Sg",
    107: "Bh",
    108: "Hs",
    109: "Mt",
    110: "Ds",
    111: "Rg",
    112: "Cn",
    113: "Nh",
    114: "Fl",
    115: "Mc",
    116: "Lv",
    117: "Ts",
    118: "Og",
}
PSE_NUMBERS: dict[str, int] = {k.lower(): v for v, k in PSE.items()}
PSE_SYMBOLS: dict[int, str] = {v: k.lower() for v, k in PSE.items()}
DESIRED_ATOM = 8
SVGFILE = Path("plot_methaneo2.svg")


def main(
    args: argparse.Namespace,
    relevant_atoms: dict[str, list[int]] | None = None,
) -> None:
    df = pd.read_csv(args.file, comment="#")

    # get all unique Method values
    methods = list(df["Method"].unique())

    # pivot the dataframe
    df = df.pivot(
        index=["CID", "Atom number", "Atom type"], columns="Method", values="Charge"
    ).reset_index()
    df.to_csv("charges_pivoted.csv", index=False)

    if args.gas:
        # exclude all methods that end with "_DIELECTRIC" or "_CPCM"
        methods = [
            method
            for method in methods
            if not method.endswith("_DIELECTRIC") and not method.endswith("_CPCM")
        ]
        print("Methods:", methods)
    if args.solvation:
        # include only methods that end with "_DIELECTRIC" or "_CPCM"
        methods = [
            method
            for method in methods
            if method.endswith("_DIELECTRIC") or method.endswith("_CPCM")
        ]
        print("Methods:", methods)

    # keep only the data point with "Atom type" == 8
    df = df[df["Atom type"] == DESIRED_ATOM]
    # check if per (i) CID, only one atom of the desired "Atom type" (here DESIRED_ATOM) is present
    if df["CID"].groupby(df["CID"]).count().max() > 1:
        cumulative_df = df.copy(deep=True)
        for method in methods:
            # calculate sum of values (atomic charges) for each CID and discard all other entries
            cumulative_df[method] = cumulative_df.groupby("CID")[method].transform(
                "sum"
            )
            # only keep the first row per CID
        cumulative_df = cumulative_df.drop_duplicates(subset="CID", keep="first")
        df = cumulative_df

    if VERBOSITY > 0:
        print(df)
        if VERBOSITY > 1:
            df.to_csv("tmp.csv", index=False)

    # currently, the CID column is a string (1.5A, 2.0A, ...), convert it to a float, after removing the "A"
    df["CID"] = df["CID"].str.replace("A", "").astype(float)
    # plot the charges
    plot_charges(df, methods)


def plot_charges(data: pd.DataFrame, methods: list[str]) -> None:
    """
    Actual plotting using seaborn.
    """
    # fig, ax = plt.subplots(figsize=(8.5 * CM_TO_INCH, 6 * CM_TO_INCH))
    plt.figure(figsize=(8.5 * CM_TO_INCH, 6 * CM_TO_INCH))
    sns.set(style="darkgrid")

    # define Roboto Condensed as the default font
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = "Roboto Condensed"
    plt.rcParams["font.weight"] = "regular"

    # plot a vertical dashed line at 3.477 Å
    plt.axvline(x=3.477, color="#909085", linestyle=":", linewidth=0.5)

    for method in methods:
        if method in ["EEQ", "EEQ_DIELECTRIC"]:
            color = "#fcba00"
            marker = "p"
            label = method
            if "_DIELECTRIC" in method:
                label = r"EEQ$_\mathrm{dielec}$"
        elif method in ["EEQ_BC", "EEQ_BC_DIELECTRIC"]:
            color = "#D773F0"
            marker = "*"
            label = r"EEQ$_\mathrm{BC}$"
            if "_DIELECTRIC" in method:
                label = r"EEQ$_\mathrm{BC,dielec}$"
        elif method == "CEH-v2":
            color = "#07529a"
            marker = "o"
            label = "CEH"
        elif method in ["GFN1-xTB", "GFN1-xTB_CPCM"]:
            color = "#7AC284"
            marker = "s"
            label = method
            if "_CPCM" in method:
                label = r"GFN1-xTB$_\mathrm{CPCM}$"
        elif method in ["GFN2-xTB", "GFN2-xTB_CPCM"]:
            color = "#C73C5F"
            marker = "^"
            label = method
            if "_CPCM" in method:
                label = r"GFN2-xTB$_\mathrm{CPCM}$"
        elif method in ["wB97M-V", "wB97M-V_CPCM"]:
            color = "#909085"
            marker = "<"
            label = "ωB97M-V"
            if "_CPCM" in method:
                label = r"ωB97M-V$_\mathrm{CPCM}$"
        else:
            raise ValueError(f"Unknown method: {method}")

        p = sns.lineplot(
            data=data,
            x="CID",
            y=method,
            marker=marker,
            label=label,
            color=color,
            markeredgecolor="black",
            markeredgewidth=0.25,
            markersize=2,
            linewidth=1,
            linestyle="-",
        )

    # set y axis limits from 0.05 to -1.05
    # plt.ylim(-1.05, 0.05)
    # include legend and put it above the plot
    plt.legend(fontsize=8, loc="best", ncols=2)
    # set y axis label
    plt.ylabel(r"cumulated atomic charge ($q$) on O2 / $e^{-}$", fontsize=8)
    # set x axis label
    plt.xlabel(r"C$_{\mathrm{CH4}}$–O2$^*$ distance / Å", fontsize=8)
    # set xticks font size
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=8)
    plt.savefig(SVGFILE, dpi=600, bbox_inches="tight", format="svg")


def get_args():
    # Argument parser
    parser = argparse.ArgumentParser(
        description="Plot the partial charge of an atom in a two-atom-molecule depending on the distance between the atoms."
    )
    parser.add_argument(
        "--file",
        "-f",
        type=str,
        required=True,
        help="Input file with the partial charges.",
    )
    parser.add_argument(
        "--gas",
        "-g",
        action="store_true",
        default=False,
        help="Plot the partial charges of the gas phase.",
    )
    parser.add_argument(
        "--solvation",
        "-s",
        action="store_true",
        default=False,
        help="Plot the partial charges of the solvated phase.",
    )

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = get_args()
    main(args)
