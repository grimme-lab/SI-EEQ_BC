"""
Python script that calls the MindlessGen API with inputs over the whole actinides series.
"""

import warnings

from itertools import combinations

from mindlessgen.generator import generator  # type: ignore
from mindlessgen.prog import ConfigManager  # type: ignore
from mindlessgen.molecules import PSE_SYMBOLS  # type: ignore


def main():
    """
    Script for execution of MindlessGen via Python API.
    """
    config = ConfigManager()

    # General settings
    config.general.max_cycles = 5000
    config.general.parallel = 28
    config.general.verbosity = 0
    config.general.num_molecules = 1
    config.general.postprocess = True
    config.general.write_xyz = False

    # Settings for the random molecule generation
    config.generate.min_num_atoms = 8
    config.generate.max_num_atoms = 20
    config.generate.forbidden_elements = "21-30,39-48,57-80" # avoid transition metals and lanthanides

    # xtb-related settings
    config.xtb.level = 1 # GFN1-xTB

    # ORCA-related settings
    config.orca.functional = "PBE"
    config.orca.basis = "def2-mTZVPP"
    config.orca.scf_cycles = 150

    # iterate over all two-element combinations of actinides (89-103)
    actinide_combinations = list(combinations(range(89, 104), 2))
    print(f"Number of actinide combinations: {len(actinide_combinations)}")
    for actinide_pair in actinide_combinations:
        actinide_symbols = [PSE_SYMBOLS[actinide] for actinide in actinide_pair]
        element_composition_string = f"{actinide_symbols[0]}:1-1,{actinide_symbols[1]}:1-1"
        config.generate.element_composition = (element_composition_string)
        print(f"Element composition string: {element_composition_string}")
        try:
            molecules, exitcode = generator(config)
        except RuntimeError as e:
            print(f"\nGeneration failed: {e}")
            raise RuntimeError("Generation failed.") from e
        if exitcode != 0:
            warnings.warn("Generation completed with errors for parts of the generation.")
        for molecule in molecules:
            molecule.write_xyz_to_file()
            print(
                "\n###############\nProperties of molecule "
                + f"'{molecule.name}' with formula {molecule.sum_formula()}:"
            )
            print(molecule)

if __name__ == "__main__":
    main()
