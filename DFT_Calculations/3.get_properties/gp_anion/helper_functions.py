import io
import re
import os
import sys
import glob
import subprocess
import pandas as pd
from PIL import Image
import matplotlib.pyplot as plt
from matplotlib.colors import ColorConverter
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from IPython.display import display
import ipywidgets as widgets


def setup_and_clean_folder(folder):
    """Create the folder if it doesn't exist and remove all its contents."""
    os.makedirs(folder, exist_ok=True)
    for file in os.listdir(folder):
        os.remove(os.path.join(folder, file))


def load_structures_from_cdxml(root_dir: str) -> dict:
    """
    Load structures from .cdxml files in the given directory.

    Args:
        root_dir (str): The directory containing .cdxml files.

    Returns:
        dict: A dictionary mapping file names (without extension) to RDKit Mol objects.
    """
    structures = {}
    for file in glob.glob("*.cdxml", root_dir=root_dir):
        # get the filename without the extension, we will use it as the prefix
        key = file.split(".")[0]
        substructure = Chem.MolsFromCDXMLFile(os.path.join(root_dir, file))
        # output to Smarts format
        temp_smarts = Chem.MolToSmarts(substructure[0])
        # reparse the Smarts, this is in case the structure from the cdxml file is not the same as from the SMILES string
        substructure = Chem.MolFromSmarts(temp_smarts)
        # add the substructure to the dictionary
        structures[key] = substructure

    return structures


def load_files_grouped_by_prefix(
    file_type: str, root_dir: str, available_prefixes: list
) -> dict:
    """
    Load files from the given directory and group them by prefix.

    Args:
        file_type (str): The type of files to load (e.g., 'sdf', 'com').
        root_dir (str): The directory containing the files.
        available_prefixes (list): List of prefixes to match.

    Returns:
        dict: A dictionary mapping prefixes to lists of file paths.
    """
    files_dict = {}
    files = glob.glob(f"*.{file_type}", root_dir=root_dir)
    if not files:
        print(f"No {file_type.upper()} files found in the input folder: {root_dir}")
        raise ValueError(
            f"No {file_type.upper()} files found in the input folder: {root_dir}"
        )

    for file in files:
        print(f"file: {file}", end=", ")
        prefix = re.search(r"^(\D+)\d+\D*_", file).group(1)
        print(f'prefix: "{prefix}"', end=", ")
        for common_prefix in available_prefixes:
            # we want exact match
            if prefix == common_prefix:
                print(f'Matched template: "{common_prefix}"')
                # create a list of files for each prefix
                if common_prefix not in files_dict:
                    files_dict[common_prefix] = []
                files_dict[common_prefix].append(file)

    return files_dict


def print_stat_on_files_by_prefix(files_dict: dict, file_type: str) -> None:
    """
    Print statistics on the number of files for each prefix.

    Args:
        files_dict (dict): A dictionary mapping prefixes to lists of file paths.
        file_type (str): The type of files (e.g., 'SDF', 'COM').
        prefix (str): The prefix to match.
    """
    print(f"Stat on the {file_type.upper()} files by prefix: ")
    total_count = 0
    for key, value in files_dict.items():
        total_count += len(value)
        print(f"{key}: {len(value)}")
    print(f"Total count: {total_count}")


def sort_files_by_prefix(files_dict: dict) -> None:
    """
    Sort files in the dictionary by the number after the prefix.

    Args:
        files_dict (dict): A dictionary mapping prefixes to lists of file paths.
    """
    for key, value in files_dict.items():
        for file in value:
            match = re.search(r"\D+(\d+)_{1}", file)
            files_dict[key] = sorted(
                value, key=lambda x: int(re.search(r"^\D+(\d+)_{1}", x).group(1))
            )


# we will optimized the coordinates of the Gaussian input file using result from the DFT calculation
def update_gaussian_com(log_file_path, com_file_path, output_folder):
    """Updates the Gaussian input .com file with the last coordinates sections from the log file.
    Parameters:
        log_file_path (str): Path to the Gaussian log file.
        com_file_path (str): Path to the Gaussian .com file.
        output_folder (str): Path to the folder where the modified .com file will be saved.
    Returns:
        None
    """
    # Read the log file
    with open(log_file_path, "r") as log_file:
        print(f"Reading log file: {log_file_path}")
        log_file_content = log_file.read()

    # keyword to search for the atomic coordinates
    coordinates_pattern = re.compile(
        r"Standard\sorientation:\s*[-]+\s+[^-]+\s+[-]+\s+"
        r"((?:\s*\d+\s+\d+\s+\d+\s+-?\d+\.?\d*\s+-?\d+\.?\d*\s+-?\d+\.?\d*\s*?)+)"
        r"\s*[-]+",
        re.MULTILINE,
    )

    # find all occurrences of the input orientation
    coordinates_matches = coordinates_pattern.findall(log_file_content)
    if coordinates_matches:
        print(f"Found {len(coordinates_matches)} standard coordinates")
    else:
        print("No standard orientation found. Try input orientation.")
        coordinates_pattern = re.compile(
            r"Input\sorientation:\s*[-]+\s+[^-]+\s+[-]+\s+"
            r"((?:\s*\d+\s+\d+\s+\d+\s+-?\d+\.?\d*\s+-?\d+\.?\d*\s+-?\d+\.?\d*\s*?)+)"
            r"\s*[-]+",
            re.MULTILINE,
        )
        coordinates_matches = coordinates_pattern.findall(log_file_content)
        if coordinates_matches:
            print(f"Found {len(coordinates_matches)} input coordinates")
        else:
            raise ValueError(
                f"No standard or input orientation found in the log file: {log_file_path.split(os.sep)[-1]}"
            )

    # Extract atomic coordinates from the last occurrence
    last_orientation = coordinates_matches[-1]
    atom_coordinates = {}  # build a dictionary of atom index and its coordinates
    for line in last_orientation.strip().split("\n"):
        parts = re.split(r"\s+", line.strip())
        atom_index = parts[0]
        x, y, z = parts[3:6]
        atom_coordinates[atom_index] = (x, y, z)
    print(f"Extracted the last atom coordinates: {atom_coordinates}")

    # Read the .com file
    with open(com_file_path, "r") as com_file:
        com_file_content = com_file.read()

    # Regex to match the atomic coordinates section in the .com file
    com_coordinates_pattern = re.compile(
        r"0\s+1\s+((?:\S+\s+-?\d+\.?\d*\s+-?\d+\.?\d*\s+-?\d+\.?\d*\s)+)", re.MULTILINE
    )
    coordinates_match = com_coordinates_pattern.search(com_file_content)
    if not coordinates_match:
        raise ValueError(
            f"No coordinates section found in the .com file: {com_file_path.split(os.sep)[-1]}"
        )

    # Update the coordinates in the .com file
    original_coordinates = coordinates_match.group(1)
    updated_coordinates = []
    for i, line in enumerate(original_coordinates.strip().split("\n"), start=1):
        # print(f"Processing line {i}: {line}")
        parts = re.split(r"\s+", line.strip())
        x, y, z = atom_coordinates.get(
            str(i), parts[1:4]
        )  # Default to original values if index not found
        updated_coordinates.append(f"{parts[0]}    {x}    {y}    {z}")

    updated_coordinates_str = "\n".join(updated_coordinates) + "\n"

    # Replace the coordinates in the .com file
    updated_com_file_content = com_coordinates_pattern.sub(
        f"0 1\n{updated_coordinates_str}", com_file_content
    )

    # print(f'Original coordinates: "{original_coordinates}"')
    # print(f'Updated coordinates: "{updated_coordinates_str}"')

    # Write the updated .com file
    output_file_path = os.path.join(output_folder, os.path.basename(com_file_path))
    os.makedirs(output_folder, exist_ok=True)
    with open(output_file_path, "w") as output_file:
        output_file.write(updated_com_file_content)

    print(
        f"Updated coordinates in the Gaussian input file {com_file_path.split(os.sep)[-1]}, saved to path: {output_file_path}"
    )


def find_benzylic_hydrogens_index(
    sdf_files_dict: dict,
    input_sdf_folder: str,
    common_structures: dict,
    output_close_shell_sdf_folder: str,
    output_open_shell_sdf_folder: str,
    output_anion_sdf_folder: str,
    display_images: bool = True,
) -> tuple[dict, dict]:
    """
    Find the benzylic hydrogens in the given SDF files and return a dictionary of the found hydrogen atom index.
    Parameters:
        sdf_files_dict (dict): A dictionary of SDF files.
        input_sdf_folder (str): The folder containing the input SDF files.
        common_structures (dict): A dictionary of common structures.
        output_close_shell_sdf_folder (str): The folder to save the output SDF files.
        output_open_shell_sdf_folder (str): The folder to save the open shell SDF files.
        output_anion_sdf_folder (str): The folder to save the anion SDF files.
        display_images (bool): Whether to display the images or not to the jupyter notebook. Default is True.
    Returns:
        tuple: A tuple of two dictionaries, the first one is the found hydrogen atom index, the second one is the warnings.
    """
    # go through the sdf_files_dict, load the sdf file using Chem.SDMolSupplier
    mols_dict = {}
    warnings = {}
    for template, files_list in sdf_files_dict.items():
        print("" + "=" * 100)
        print(f"Processing prefix: {template}")
        for sdf_file in files_list:
            print("" + "-" * 80)
            print(
                f"Processing sdf file: {sdf_file}, base filename: {sdf_file.split('.')[0]}"
            )
            supplier = Chem.SDMolSupplier(
                os.path.join(input_sdf_folder, sdf_file), removeHs=False
            )

            suffix = 1
            for mol in supplier:
                mol_without_H = Chem.RemoveHs(mol)

                key = f"{sdf_file.split('.')[0]}-{suffix}"
                if suffix == 1:
                    print(
                        f"Conformer: {key} (This should match the based filename of the com files)"
                    )
                else:
                    print(f"Continue finding conformer: {key}")
                mols_dict[key] = mol

                # write the closed shell molecule to sdf file
                mol_closed_shell = Chem.RWMol(mol)
                mol_closed_shell.SetProp("_Name", f"{key}.sdf")
                sd_writer_closed_shell = Chem.SDWriter(
                    os.path.join(output_close_shell_sdf_folder, f"{key}.sdf")
                )
                sd_writer_closed_shell.SetKekulize(True)
                sd_writer_closed_shell.write(mol_closed_shell.GetMol())
                sd_writer_closed_shell.close()

                # there is a possibility that there are multiple matches
                submatches = mol.GetSubstructMatches(common_structures[template])
                # combine the submatches into a single list
                matcheslist = list([item for sublist in submatches for item in sublist])

                print(
                    f"All common structure match: {submatches}, Combined match atom: {matcheslist}"
                )
                # if we find multiple matches, print a warning
                if len(submatches) > 1:
                    print(
                        f"Warning: multiple matches found in conformer {key} under file {sdf_file} using common structure {template}"
                    )
                    warnings[key] = (
                        f"Warning: multiple matches found in conformer {key} under file {sdf_file} using common structure {template}"
                    )
                if len(submatches) == 0:
                    print(
                        f"Error: no match found in conformer {key} under file {sdf_file} using common structure {template}"
                    )
                    warnings[key] = (
                        f"Error: no match found in conformer {key} under file {sdf_file} using common structure {template}"
                    )

                if display_images:
                    # Compute 2D coordinates and redraw the molecule
                    drawer = Chem.Draw.rdMolDraw2D.MolDraw2DCairo(500, 500)
                    drawer.drawOptions().comicMode = True
                    drawer.drawOptions().addAtomIndices = True
                    drawer.drawOptions().continuousHighlight = True
                    drawer.DrawMolecule(
                        mol_without_H,
                        highlightAtoms=matcheslist,
                        highlightAtomColors={
                            atom: ColorConverter().to_rgba("Salmon", alpha=0.75)
                            for atom in matcheslist
                        },
                    )
                    drawer.FinishDrawing()
                    img1 = drawer.GetDrawingText()

                for atom in matcheslist:
                    if (
                        mol.GetAtomWithIdx(atom).GetSymbol() == "C"
                    ):  # check if the atom is carbon
                        for neighbor in mol.GetAtomWithIdx(atom).GetNeighbors():
                            # limit to C neighbors, ignore those that are already in the matcheslist
                            if (
                                neighbor.GetSymbol() == "C"
                                and neighbor.GetIdx() not in matcheslist
                            ):
                                # all bonds on this neighbor should be single bond, and contains at least one hydrogen
                                if (
                                    neighbor.GetBonds()
                                    is not None  # check if there is a bond
                                    and neighbor.GetTotalNumHs(includeNeighbors=True)
                                    >= 1  # check if there is at least one hydrogen
                                    and all(
                                        [
                                            bond.GetBondType()
                                            == Chem.rdchem.BondType.SINGLE
                                            for bond in neighbor.GetBonds()
                                        ]
                                    )  # check if all bonds are single bond
                                ):
                                    # get all hydrogens atom index of this neighbor
                                    hydrogens = [
                                        neighbor_neighbors.GetIdx()
                                        for neighbor_neighbors in neighbor.GetNeighbors()
                                        if neighbor_neighbors.GetSymbol() == "H"
                                    ]
                                    hydrogens.sort()

                                    # remove the lowest index hydrogen atom from the mol to create the open shell species and anion
                                    H_to_remove = hydrogens[0]
                                    # remove the H and the bond from the molecule
                                    mol_open_shell = Chem.RWMol(mol)
                                    mol_open_shell.RemoveBond(
                                        neighbor.GetIdx(), H_to_remove
                                    )
                                    mol_open_shell.RemoveAtom(H_to_remove)
                                    # set the radical electron on the neighbor atom
                                    mol_open_shell.GetAtomWithIdx(
                                        neighbor.GetIdx()
                                    ).SetNumRadicalElectrons(1)
                                    mol_open_shell.UpdatePropertyCache()
                                    # write to sdf file
                                    mol_open_shell.SetProp("_Name", f"{key}_openshell")
                                    sd_writer_open_shell = Chem.SDWriter(
                                        os.path.join(
                                            output_open_shell_sdf_folder,
                                            f"{key}_openshell.sdf",
                                        )
                                    )
                                    sd_writer_open_shell.SetKekulize(True)
                                    sd_writer_open_shell.write(mol_open_shell.GetMol())
                                    sd_writer_open_shell.close()

                                    # remove the H and the bond from the molecule
                                    mol_anion = Chem.RWMol(mol)
                                    mol_anion.RemoveBond(neighbor.GetIdx(), H_to_remove)
                                    mol_anion.RemoveAtom(H_to_remove)
                                    # hack to get the valence to 3
                                    mol_anion.GetAtomWithIdx(
                                        neighbor.GetIdx()
                                    ).SetFormalCharge(-1)
                                    mol_anion.UpdatePropertyCache(strict=False)
                                    # write to sdf file
                                    mol_anion.SetProp("_Name", f"{key}_anion")
                                    sd_writer_anion = Chem.SDWriter(
                                        os.path.join(
                                            output_anion_sdf_folder, f"{key}_anion.sdf"
                                        )
                                    )
                                    sd_writer_anion.SetKekulize(True)
                                    sd_writer_anion.write(mol_anion.GetMol())
                                    sd_writer_anion.close()

                                    # store the found hydrogen atom index at value, filename as key
                                    #! increase all value in the hydrogen index by 1 since GaussianView is 1-based index
                                    hydrogens_index_GaussianView = [
                                        hydrogen + 1 for hydrogen in hydrogens
                                    ]
                                    mols_dict[key] = hydrogens_index_GaussianView

                                    print(
                                        f"Valid neighbor atom index: {neighbor.GetIdx()}, neighbor's H neighbors index(0-based): {hydrogens}, index at GaussianView(1-based): {hydrogens_index_GaussianView}"
                                    )
                                    if display_images:
                                        # highlight the neighbor in gold, highlight the hydrogens in yellow, high light the substructure in red
                                        highlightAtoms = (
                                            matcheslist
                                            + [neighbor.GetIdx()]
                                            + hydrogens
                                        )
                                        highlightAtomColors = {}
                                        for atom in highlightAtoms:
                                            if atom in matcheslist:
                                                highlightAtomColors[atom] = (
                                                    ColorConverter().to_rgba(
                                                        "Salmon", alpha=0.75
                                                    )
                                                )
                                            elif atom == neighbor.GetIdx():
                                                highlightAtomColors[atom] = (
                                                    ColorConverter().to_rgba(
                                                        "gold", alpha=0.75
                                                    )
                                                )
                                            elif atom in hydrogens:
                                                highlightAtomColors[atom] = (
                                                    ColorConverter().to_rgba(
                                                        "DeepSkyBlue", alpha=0.75
                                                    )
                                                )

                                        # set explicit H count to the len of hydrogens for the neighbor atom
                                        drawer = Chem.Draw.rdMolDraw2D.MolDraw2DCairo(
                                            600, 600
                                        )
                                        Chem.rdDepictor.GenerateDepictionMatching2DStructure(
                                            mol, mol_without_H
                                        )
                                        drawer.drawOptions().comicMode = True
                                        drawer.drawOptions().addAtomIndices = True
                                        drawer.drawOptions().continuousHighlight = True
                                        drawer.DrawMolecule(
                                            mol,
                                            highlightAtoms=highlightAtoms,
                                            highlightAtomColors=highlightAtomColors,
                                        )
                                        drawer.FinishDrawing()
                                        img2 = drawer.GetDrawingText()

                                        # Compute 2D coordinates and redraw the molecule
                                        Chem.rdDepictor.Compute2DCoords(mol)
                                        Chem.rdDepictor.NormalizeDepiction(mol)
                                        drawer = Chem.Draw.rdMolDraw2D.MolDraw2DCairo(
                                            600, 600
                                        )
                                        drawer.drawOptions().comicMode = True
                                        drawer.drawOptions().addAtomIndices = True
                                        drawer.drawOptions().continuousHighlight = True
                                        drawer.DrawMolecule(
                                            mol,
                                            highlightAtoms=highlightAtoms,
                                            highlightAtomColors=highlightAtomColors,
                                        )
                                        drawer.FinishDrawing()
                                        img3 = drawer.GetDrawingText()

                if display_images:
                    # display both images in a grid
                    fig, ax = plt.subplots(1, 3, figsize=(12, 4))
                    ax[0].imshow(Image.open(io.BytesIO(img1)))
                    ax[0].axis("off")
                    ax[0].set_title("Substructure (in red)")
                    ax[1].imshow(Image.open(io.BytesIO(img2)))
                    ax[1].axis("off")
                    ax[1].set_title("Benzylic hydrogen (in blue)")
                    ax[2].imshow(Image.open(io.BytesIO(img3)))
                    ax[2].axis("off")
                    ax[2].set_title("Drawing in 2D")
                    plt.show()

                # !IF YOU HAVE MULTIPLE CONFORMER, THERE WILL BE MULTIPLE MOLECULES IN ONE SDF FILE, thus the suffix will be increased
                # increase the suffix to differentiate different conformers
                suffix += 1

    return mols_dict, warnings


def parse_goodvibes_output(output: str, temp: float = 298.15) -> dict:
    """
    Parse the output of GoodVibes to extract thermodynamic properties.
    Args:
        output (str): The output string from GoodVibes.
        temp (float): The temperature in Kelvin.
    Returns:
        dict: A dictionary containing the extracted properties.
    """
    # This function extracts the desired values from the GoodVibes output
    lines = output.splitlines()
    data = {}
    column_mapping = {
        "E_SPC": "E_spc (Hartree)",
        "E": "E (Hartree)",
        "ZPE": "ZPE(Hartree)",
        "H_SPC": "H_spc(Hartree)",
        "T.S": "T*S",
        "T.qh-S": "T*qh_S",
        "G(T)_SPC": "G(T)_spc(Hartree)",
        "qh-G(T)_SPC": "qh_G(T)_spc(Hartree)",
    }

    # Find the index positions of the two lines of asterisks
    start_index = None
    end_index = None
    header_line = None
    for i, line in enumerate(lines):
        if re.match(
            r"^\s*\*{12,}\s*$", line
        ):  # Matches lines with 12 or more asterisks
            if start_index is None:
                start_index = i
                header_line = lines[
                    i - 1
                ]  # The header line is the one before the first line of asterisks
            else:
                end_index = i
                break  # We only need the first two lines of asterisks

    # Parse the header line to determine the order of properties
    headers = []
    if header_line:
        headers = re.split(r"\s+", header_line.strip())[
            1:
        ]  # get rid of the first column, which is the structure column

    # Extract relevant lines between the two asterisk lines
    if start_index is not None and end_index is not None and end_index > start_index:
        for line in lines[start_index + 1 : end_index]:
            if re.match(
                r"^\s*o", line
            ):  # Matches lines starting with 'o' (with any amount of whitespace before)
                parts = re.split(r"\s+", line.strip())  # Split the line by whitespace
                structure_name = parts[1]  # The structure name is in parts[1]

                # Initialize a dictionary for this structure
                structure_data = {}
                # Populate the structure's data dictionary using headers and corresponding values
                for i, header in enumerate(headers):
                    if header in column_mapping:
                        structure_data[column_mapping[header]] = float(
                            parts[i + 2]
                        )  # Offset by 2 for correct column indexing
                    else:
                        structure_data[header] = float(parts[i + 2])

                structure_data["T"] = temp  # Add temperature to each structure's data
                data[structure_name] = structure_data  # Add to the main data dictionary

    return data


def get_goodvibes_e_batch(input_log_folder: str, temp: float = 298.15) -> dict:
    """
    Run GoodVibes on a batch of log files and parse the output.
    Args:
        input_log_folder (str): The folder containing the log files.
        temp (float): The temperature in Kelvin.
    Returns:
        dict: A dictionary containing the parsed data.
    """
    # Construct command-line arguments for GoodVibes
    cmd_args = [
        sys.executable,
        "-m",
        "goodvibes",
        os.path.join(input_log_folder, "*.log"),
        "--spc",
        "link",
        "-t",
        str(temp),
    ]
    # Run the GoodVibes command and capture the output
    result = subprocess.run(
        cmd_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    # Parse the output
    parsed_data = parse_goodvibes_output(result.stdout, temp)
    return parsed_data


def group_results_by_substrates(result: dict) -> dict:
    # Initialize a new dictionary to
    result_grouped_by_substrates = {}
    # Regex pattern to capture the substrate name before "_", for example pyrd1_conf-1 -> pyrd1
    substrate_pattern = re.compile(r"^(.+?)_")
    for conformer_name, properties in result.items():
        # Extract substrate name using regex
        match = substrate_pattern.match(conformer_name)
        if match:
            substrate_name = match.group(1)

            # Append conformer data to the list for this substrate
            if substrate_name not in result_grouped_by_substrates:
                result_grouped_by_substrates[substrate_name] = []
            result_grouped_by_substrates[substrate_name].append(
                {conformer_name: properties}
            )
    return result_grouped_by_substrates


# Function to select the conformer with the lowest specified energy property for each substrate
def select_lowest_energy_conformer(
    substrates: dict, column: str = "G(T)_spc(Hartree)"
) -> dict:
    print(f"filtering based on {column}")
    result_filtered = {}
    for substrate, conformers in substrates.items():
        print("-" * 80)
        print(f"for substrate {substrate}, we have {len(conformers)} conformers")
        for conformer in conformers:
            print(
                f"{list(conformer.keys())[0]} have {column}: {list(conformer.values())[0].get(column, None)}"
            )
        # Find the conformer with the minimum specified energy value
        min_conformer = min(
            conformers, key=lambda x: list(x.values())[0].get(column, float("inf"))
        )
        print(f"lowest energy conformer: {list(min_conformer.keys())[0]}")
        # Add the lowest energy conformer to the result dictionary
        result_filtered[substrate] = min_conformer
    return result_filtered


def draw_grid_image(
    temp_folder: str,
    img_list: list,
    title_list: list,
    num_Cols: int,
    save_image_prefix: str,
    item_suffix: str,
):
    """Draws a grid image of the given images and saves it to the temp folder.
    Parameters:
        temp_folder (str): The folder to save the image.
        img_list (list): A list of images to be drawn.
        title_list (list): A list of titles for each image.
        num_Cols (int): The number of columns in the grid.
        save_image_prefix (str): The prefix for the image file name.
        item_suffix (str): The suffix for the image file name.
    Returns:
        None
    """
    num_images = len(img_list)
    num_rows = (
        num_images + num_Cols - 1
    ) // num_Cols  # Calculate number of rows needed
    # calculate the combined image size
    plt.figure(figsize=(num_Cols, num_rows), layout="constrained", dpi=600)
    # get the number of pixes in the figure
    dpi = plt.gcf().get_dpi()
    fig_width, fig_height = plt.gcf().get_size_inches()
    fig_width *= dpi
    fig_height *= dpi
    # print(f"fig_width: {fig_width}, fig_height: {fig_height}")

    # if the width exceeds 2^16, throw an error tell the user to shrink the number of columns
    if fig_width > 2**16:
        raise ValueError(
            "The combined width exceeds the maximum allowed size. Please reduce the number of columns."
        )
    # if the height exceeds 2^16, we will perform cutoff generate multiple images
    num_combined_images = int(fig_height // (2**16) + 1)
    # print(f"num_combined_images: {num_combined_images}")

    if num_combined_images > 1:
        # calculate the maximum images for each combined image
        max_images_per_combined_image = (
            int((2**16) // (fig_height // num_rows)) * num_rows
        )
        # print(f"max_images_per_combined_image: {max_images_per_combined_image}")
        for i in range(num_combined_images):
            start_index = i * max_images_per_combined_image
            end_index = min((i + 1) * max_images_per_combined_image, num_images)
            plt.figure(
                figsize=(num_Cols, int((end_index - start_index) / num_Cols)),
                layout="constrained",
                dpi=600,
            )
            for j in range(start_index, end_index):
                plt.subplot(num_rows, num_Cols, j + 1)
                plt.imshow(img_list[j])
                plt.axis("on")
                plt.grid(True)  # Enable grid lines
                plt.title(title_list[j], fontsize=6)
                plt.xticks([])  # Remove x-ticks
                plt.yticks([])  # Remove y-ticks
            plt.savefig(
                temp_folder + os.sep + save_image_prefix + item_suffix + f"_{i + 1}.png"
            )
            plt.close()
    else:
        plt.figure(figsize=(num_Cols, num_rows), layout="constrained", dpi=600)
        for i in range(num_images):
            plt.subplot(num_rows, num_Cols, i + 1)
            plt.imshow(img_list[i])
            plt.axis("on")  # Turn on the axis to show grid lines
            plt.grid(True)  # Enable grid lines
            plt.title(title_list[i], fontsize=6)
            plt.xticks([])  # Remove x-ticks
            plt.yticks([])  # Remove y-ticks
        plt.savefig(
            temp_folder + os.sep + save_image_prefix + item_suffix + ".png"
        )  # save the grid image to temp folder
        plt.close()


def search_for_substructure_close_shell(
    temp_folder: str,
    all_compounds: Chem.rdmolfiles.SDMolSupplier,
    common_structure: Chem.rdchem.Mol,
    prefix: str,
    num_Cols: int = 5,
) -> pd.DataFrame:
    """
    Search for the substructure in the given compounds and return a list of atom indices.
    Parameters:
        temp_folder (str): The folder to save the image.
        all_compounds (Chem.rdmolfiles.SDMolSupplier): A list of compounds load via Chem.SDMolSupplier.
        common_structure (Chem.rdchem.Mol): The substructure to search for.
        prefix (str): The prefix for the image file name.
        num_Cols (int): The number of columns in the grid.
    Returns:
        pd.DataFrame: A dataframe of the atom indices.
    """
    # uses RDKit to search for the substructure in each compound you will analyze
    atoms = []
    img_list = []
    for molecule in all_compounds:
        if molecule is not None:
            submatches = molecule.GetSubstructMatches(
                common_structure
            )  # find substructure
            matchlist = []

            global_found = False
            for submatch in submatches:
                found = False
                temp_matchlist = [
                    item for item in submatch
                ]  # list of zero-indexed atom numbers

                # !this is specific to this project which is chlorination on the alkyl side chain on the substrucutre, where we care about the alkyl side chain
                # search atoms that are connected to the found substructure, which are also not part of the current substructure
                # they must be a carbon atom and it have to have at least one hydrogen neighbors and all bond must be single bond
                connected_atoms = []
                for atom in temp_matchlist:
                    if (
                        molecule.GetAtomWithIdx(atom).GetSymbol() == "C"
                    ):  # check if the atom is carbon
                        for neighbor in molecule.GetAtomWithIdx(
                            atom
                        ).GetNeighbors():  # iterate through all neighbors
                            if (
                                neighbor.GetSymbol() == "C"
                                and neighbor.GetIdx() not in temp_matchlist
                            ):  # filter neighbor that is not carbon and is not in the temp_matchlist
                                if (
                                    neighbor.GetBonds() is not None
                                ):  # check if there is any bond going to the neighbor
                                    if (
                                        neighbor.GetTotalNumHs(includeNeighbors=True)
                                        >= 1
                                    ):  # check if there is at least one hydrogen, !specific for closed shell molecules
                                        if all(
                                            [
                                                bond.GetBondType()
                                                == Chem.rdchem.BondType.SINGLE
                                                for bond in neighbor.GetBonds()
                                            ]
                                        ):  # check if all bonds are single bond
                                            found = True
                                            if not global_found:
                                                global_found = True
                                            else:
                                                print("Warning: double match found")
                                                display(molecule)
                                            # the atom now is the ipso carbon, the neighbor is the benzylic carbon
                                            connected_atoms.append(atom)
                                            connected_atoms.append(neighbor.GetIdx())
                                            # we don't add they right now to the temp_matchlist, in case the top loop is not finished

                                            temp_matchlist.remove(
                                                atom
                                            )  # drop the ipso carbon from temp_matchlist for duplicate
                                            # get all hydrogens atom index of this neighbor
                                            hydrogens = [
                                                atom.GetIdx()
                                                for atom in neighbor.GetNeighbors()
                                                if atom.GetSymbol() == "H"
                                            ]
                                            hydrogens.sort()
                                            # add only the first hydrogen atoms index to the connected_atoms list
                                            connected_atoms.append(hydrogens[0])

                if found:
                    # append connected_atoms to the matchlist, which will fix the position of the ipso carbon, benzylic carbon, and the hydrogen atom at the end of the matchlist
                    matchlist.extend(temp_matchlist)
                    matchlist.extend(connected_atoms)

                    match_atom_symbol = [
                        molecule.GetAtomWithIdx(x).GetSymbol() for x in matchlist
                    ]  # find atom symbols
                    match_idx = [
                        x + 1 for x in matchlist
                    ]  # changes from 0-indexed to 1-indexed (for Gaussian view)
                    match_combined = [
                        str(match_atom_symbol[i]) + str(match_idx[i])
                        for i in range(len(match_atom_symbol))
                    ]  # combine atom symbol and number

                    atoms.append(
                        match_combined
                    )  # append 1-indexed list to atoms (a list of lists)
            if not global_found:
                print("Warning: no match found")
                display(molecule)
            # add a label to the atom that is being matched
            for atom in matchlist:
                molecule.GetAtomWithIdx(atom).SetProp(
                    "atomLabel", match_combined[matchlist.index(atom)]
                )  # label the atom being matched

            # now create a grid image of all the molecules, label the atom being matched
            # the atom that had the property name GaussianMap added to it will be labeled
            # we will draw each molecule with the atom number labeled and the substructure highlighted and then combine them into a grid image with captions of the file name
            Chem.rdDepictor.Compute2DCoords(molecule)

            # Prepare the highlight color dictionary
            # last item in the matchlist is the hydrogen atom, the second last is the benzylic carbon
            # the third last is the ipso carbon, the rest are the ring atoms
            highlight_colors = {}
            for atom in matchlist[:-3]:
                highlight_colors[atom] = ColorConverter().to_rgba("Pink", alpha=0.8)
            for atom in matchlist[-3:-2]:
                highlight_colors[atom] = ColorConverter().to_rgba("Gold", alpha=0.8)
            for atom in matchlist[-2:-1]:
                highlight_colors[atom] = ColorConverter().to_rgba("Green", alpha=0.8)
            for atom in matchlist[-1:]:
                highlight_colors[atom] = ColorConverter().to_rgba(
                    "DeepSkyBlue", alpha=0.8
                )

            # Create a drawer
            drawer = Chem.Draw.rdMolDraw2D.MolDraw2DCairo(
                600, 600
            )  # Use Cairo backend for drawing
            drawer.drawOptions().continuousHighlight = True
            # Prepare highlight dictionary
            drawer.DrawMolecule(
                molecule, highlightAtoms=matchlist, highlightAtomColors=highlight_colors
            )
            drawer.FinishDrawing()

            png_data = drawer.GetDrawingText()
            img = Image.open(io.BytesIO(png_data))
            img_list.append(img)

    # this loop extracts log names from log_ids and splits them to the desired format
    filenames = open(temp_folder + os.sep + "log_ids_" + prefix + ".txt", "r")
    # it is a text file that contains the file name for every molecule you will analyze
    list_of_filenames = [
        (line.strip()).split() for line in filenames
    ]  # list of the file names (each of which includes all conformers)
    list_of_files = []
    for filename in list_of_filenames:
        file = filename[0].split(".")
        list_of_files.append(file[0])
    filenames.close()

    draw_grid_image(
        temp_folder=temp_folder,
        img_list=img_list,
        title_list=list_of_files,
        num_Cols=num_Cols,
        save_image_prefix="common_structure_",
        item_suffix=prefix,
    )

    # put the atom numbers for the substructure for each log file into a dataframe
    prelim_df = pd.DataFrame(atoms)
    prelim_df.insert(0, column="log_name", value=list_of_files)

    return prelim_df


def search_for_substructure_anion(
    temp_folder: str,
    all_compounds: Chem.rdmolfiles.SDMolSupplier,
    common_structure: Chem.rdchem.Mol,
    prefix: str,
    num_Cols: int = 5,
) -> pd.DataFrame:
    """
    Search for the substructure in the given compounds and return a list of atom indices.
    Parameters:
        temp_folder (str): The folder to save the image.
        all_compounds (Chem.rdmolfiles.SDMolSupplier): A list of compounds load via Chem.SDMolSupplier.
        common_structure (Chem.rdchem.Mol): The substructure to search for.
        prefix (str): The prefix for the image file name.
        num_Cols (int): The number of columns in the grid.
    Returns:
        pd.DataFrame: A dataframe of the atom indices.
    """
    # uses RDKit to search for the substructure in each compound you will analyze
    pt = Chem.GetPeriodicTable()
    atoms = []
    img_list = []
    for molecule in all_compounds:
        if molecule is not None:
            submatches = molecule.GetSubstructMatches(
                common_structure
            )  # find substructure
            matchlist = []

            global_found = False
            for submatch in submatches:
                found = False
                temp_matchlist = [
                    item for item in submatch
                ]  # list of zero-indexed atom numbers

                # !this is specific to this project which is chlorination on the alkyl side chain on the substrucutre, where we care about the alkyl side chain
                # search atoms that are connected to the found substructure, which are also not part of the current substructure
                # they must be a carbon atom and it have to have at least one hydrogen neighbors and all bond must be single bond
                connected_atoms = []
                for atom in temp_matchlist:
                    if (
                        molecule.GetAtomWithIdx(atom).GetSymbol() == "C"
                    ):  # check if the atom is carbon
                        for neighbor in molecule.GetAtomWithIdx(
                            atom
                        ).GetNeighbors():  # iterate through all neighbors
                            if (
                                neighbor.GetSymbol() == "C"
                                and neighbor.GetIdx() not in temp_matchlist
                            ):  # filter neighbor that is not carbon and is not in the temp_matchlist
                                if (
                                    neighbor.GetBonds() is not None
                                ):  # check if there is any bond going to the neighbor
                                    if (
                                        pt.GetDefaultValence(neighbor.GetSymbol())
                                        != neighbor.GetTotalValence()
                                        or neighbor.GetFormalCharge() == -1
                                    ):  # check for either valence mismatch or a formal charge of -1 on the neighbor atom
                                        if all(
                                            [
                                                bond.GetBondType()
                                                == Chem.rdchem.BondType.SINGLE
                                                for bond in neighbor.GetBonds()
                                            ]
                                        ):  # check if all bonds are single bond
                                            found = True
                                            if not global_found:
                                                global_found = True
                                            else:
                                                print("Warning: double match found")
                                                display(molecule)
                                            # the atom now is the ipso carbon, the neighbor is the benzylic carbon
                                            connected_atoms.append(atom)
                                            connected_atoms.append(neighbor.GetIdx())
                                            # we don't add they right now to the temp_matchlist, in case the top loop is not finished
                                            # drop the ipso carbon from temp_matchlist for duplicate
                                            temp_matchlist.remove(atom)

                if found:
                    # append connected_atoms to the matchlist, which will fix the position of the ipso carbon, benzylic carbon, and the hydrogen atom at the end of the matchlist
                    matchlist.extend(temp_matchlist)
                    matchlist.extend(connected_atoms)

                    match_atom_symbol = [
                        molecule.GetAtomWithIdx(x).GetSymbol() for x in matchlist
                    ]  # find atom symbols
                    match_idx = [
                        x + 1 for x in matchlist
                    ]  # changes from 0-indexed to 1-indexed (for Gaussian view)
                    match_combined = [
                        str(match_atom_symbol[i]) + str(match_idx[i])
                        for i in range(len(match_atom_symbol))
                    ]  # combine atom symbol and number

                    # append 1-indexed list to atoms (a list of lists)
                    atoms.append(match_combined)
            if not global_found:
                print("Warning: no match found")
                display(molecule)
            # add a label to the atom that is being matched
            for atom in matchlist:
                molecule.GetAtomWithIdx(atom).SetProp(
                    "atomLabel", match_combined[matchlist.index(atom)]
                )  # label the atom being matched

            # now create a grid image of all the molecules, label the atom being matched
            # the atom that had the property name GaussianMap added to it will be labeled
            # we will draw each molecule with the atom number labeled and the substructure highlighted and then combine them into a grid image with captions of the file name
            Chem.rdDepictor.Compute2DCoords(molecule)

            # Prepare the highlight color dictionary
            # last item in the matchlist is the hydrogen atom, the second last is the benzylic carbon
            # the third last is the ipso carbon, the rest are the ring atoms
            highlight_colors = {}
            for atom in matchlist[:-2]:
                highlight_colors[atom] = ColorConverter().to_rgba("Pink", alpha=0.8)
            for atom in matchlist[-2:-1]:
                highlight_colors[atom] = ColorConverter().to_rgba("Gold", alpha=0.8)
            for atom in matchlist[-1:]:
                highlight_colors[atom] = ColorConverter().to_rgba("Green", alpha=0.8)

            # Create a drawer
            drawer = Chem.Draw.rdMolDraw2D.MolDraw2DCairo(
                600, 600
            )  # Use Cairo backend for drawing
            drawer.drawOptions().continuousHighlight = True
            # Prepare highlight dictionary
            drawer.DrawMolecule(
                molecule, highlightAtoms=matchlist, highlightAtomColors=highlight_colors
            )
            drawer.FinishDrawing()

            png_data = drawer.GetDrawingText()
            img = Image.open(io.BytesIO(png_data))
            img_list.append(img)

    # this loop extracts log names from log_ids and splits them to the desired format
    filenames = open(temp_folder + os.sep + "log_ids_" + prefix + ".txt", "r")
    # it is a text file that contains the file name for every molecule you will analyze
    list_of_filenames = [
        (line.strip()).split() for line in filenames
    ]  # list of the file names (each of which includes all conformers)
    list_of_files = []
    for filename in list_of_filenames:
        file = filename[0].split(".")
        list_of_files.append(file[0])
    filenames.close()

    draw_grid_image(
        temp_folder=temp_folder,
        img_list=img_list,
        title_list=list_of_files,
        num_Cols=num_Cols,
        save_image_prefix="common_structure_",
        item_suffix=prefix,
    )

    # put the atom numbers for the substructure for each log file into a dataframe
    prelim_df = pd.DataFrame(atoms)
    prelim_df.insert(0, column="log_name", value=list_of_files)

    return prelim_df
