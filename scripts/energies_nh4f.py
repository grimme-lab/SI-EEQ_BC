"""
Plot the total energy of the NH4F molecule depending on the F–H distance.
"""

import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns  # type: ignore
import pandas as pd


CM_TO_INCH = 1 / 2.54
VERBOSITY = 1
HARTREE_TO_KCAL = 627.50947428


def main(
    args: argparse.Namespace,
) -> None:
    df = pd.read_csv(args.file, comment="#")

    # get all unique Method values
    methods = list(df["Method"].unique())

    # pivot the dataframe
    df = df.pivot(index=["CID"], columns="Method", values="Energy").reset_index()
    df.to_csv("energies_pivoted.csv", index=False)

    # subtract the energy of the last (longest) F–H distance with wB97M-V_dSCF from all energies
    kcal = False
    if args.set_to_zero:
        # reset the values to the dissociation limit
        for method in methods:
            df[method] -= df[args.set_to_zero].iloc[-1]
            # convert all values to kcal/mol
            df[method] *= HARTREE_TO_KCAL
        kcal = True

    if VERBOSITY > 0:
        print(df)
        if VERBOSITY > 1:
            df.to_csv("tmp.csv", index=False)
    # plot the charges
    plot_charges(df, methods, kcal)


def plot_charges(data: pd.DataFrame, methods: list[str], kcal: bool = False) -> None:
    """
    Actual plotting using seaborn.
    """
    # fig, ax = plt.subplots(figsize=(8.5 * CM_TO_INCH, 6 * CM_TO_INCH))
    plt.figure(figsize=(8.5 * CM_TO_INCH, 5 * CM_TO_INCH))
    sns.set(style="darkgrid")

    # define Roboto Condensed as the default font
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = "Roboto Condensed"
    plt.rcParams["font.weight"] = "regular"

    for method in methods:
        if method in ["GFN1-xTB", "GFN1-xTB_CPCM"]:
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
            label = r"ωB97M-V$_\mathrm{RKS}$"
            if "_CPCM" in method:
                label = r"ωB97M-V$_\mathrm{CPCM}$"
        elif method in ["wB97M-V_dSCF", "wB97M-V_dSCF_CPCM"]:
            color = "#06768d"
            marker = ">"
            label = r"ωB97M-V$_\mathrm{UKS\ singlet}$"
            if "_CPCM" in method:
                label = r"ωB97M-V$_\mathrm{dSCF,CPCM}$"
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
            markersize=3,
            linewidth=1,
            linestyle="-",
        )

    # include legend and put it above the plot
    plt.legend(fontsize=8, loc="best")
    # set y axis label
    if kcal:
        plt.ylabel(r"relative energy / kcal mol$^{-1}$", fontsize=8)
    else:
        plt.ylabel(r"total energy / $E_{\mathrm{H}}$", fontsize=8)
    # set x axis label
    plt.xlabel(r"F–H$_\mathrm{NH4}$ distance / Å", fontsize=8)
    # set xticks font size
    plt.xticks(np.arange(min(data["CID"]), max(data["CID"]) + 1, 3.0), fontsize=8)
    plt.yticks(fontsize=8)
    # Generate output file name from the method names
    filename = "plot_nh4f_energies.svg"
    # Saving the plot
    plt.savefig(filename, dpi=600, bbox_inches="tight")


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
        "--set-to-zero",
        "-stz",
        type=str,
        default=None,
        required=False,
        help="Reset the values to the dissociation limit. Option: <method>.",
    )

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = get_args()
    main(args)
