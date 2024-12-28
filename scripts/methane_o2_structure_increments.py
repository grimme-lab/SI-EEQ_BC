from pathlib import Path
from shutil import copy


def increment_and_write_xyz(
    input_file: Path, output_dir: Path, start: float, stop: float, step: float
) -> None:
    # Read the input file
    with input_file.open("r") as f:
        lines = f.readlines()

    # Extract the header and coordinates
    header = lines[:2]
    atoms = lines[2:]

    # set up a list with all directories
    all_directories: list[Path] = []

    # Process each increment
    current_y = start
    while current_y <= stop:
        modified_atoms = []
        for atom in atoms:
            parts = atom.split()
            # Modify only the last two oxygen atoms' y-coordinates
            if (
                parts[0] == "O" and float(parts[2]) == 2.0
            ):  # Check if y-coordinate is 2.0
                parts[2] = f"{current_y:.6f}"
            # Reformat the line with consistent spacing and precision
            formatted_atom = f"{parts[0]:<2} {float(parts[1]):>14.8f} {float(parts[2]):>14.8f} {float(parts[3]):>14.8f}"
            modified_atoms.append(formatted_atom)

        # Prepare the output directory and file
        increment_dir = output_dir / f"{current_y:.1f}A"
        increment_dir.mkdir(parents=True, exist_ok=True)
        all_directories.append(increment_dir)
        increment_file = increment_dir / f"{current_y:.1f}A.xyz"

        # Write the new xyz file
        with increment_file.open("w") as f:
            f.writelines(header)
            f.writelines("\n".join(modified_atoms) + "\n")

        # copy the xyz file to the name "struc.xyz" in the same directory
        copy(increment_file, increment_dir / "struc.xyz")

        # Increment y-coordinate
        current_y += step

    # write all the directories created to a file called .list.dirs
    with open(output_dir / ".list.dirs", "w") as f:
        for directory in all_directories:
            f.write(f"{directory.name}\n")


# Define paths and parameters
input_xyz = Path("ch4o2.xyz")
output_base_dir = Path("output_xyz")
output_base_dir.mkdir(exist_ok=True)
start_increment = 1.5
stop_increment = 8.0
step_increment = 0.1

# Run the function
increment_and_write_xyz(
    input_xyz, output_base_dir, start_increment, stop_increment, step_increment
)
