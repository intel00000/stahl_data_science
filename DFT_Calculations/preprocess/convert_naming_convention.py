import os
import json
import pandas as pd
import argparse

# This script converts the com file from Tom naming convention (filename by common structure, for example pyrd1, pyrd2...)
# to Leah naming convention Het001, Het002, Het003...
# The mapping is defined in the excel file named smiles_with_mapping.xlsx

# for each file, split the filename by '_' and check if the prefix has a match in the mapping dictionary
# if there is, replace the filename prefix with the value in the mapping dictionary.
# Also read the original file content, replace ANY old prefix with the new prefix (to change the output job filenames)
# write to the new file


def main(original_folder, new_folder, file_extension):
    # Read in the mapping from smiles_with_mapping.xlsx
    df = pd.read_excel("smiles_with_mapping.xlsx", header=0)
    df = df[["id", "mapping"]]

    # Convert the mapping column to a dictionary with id being key, mapping being value
    mapping_dict_from_tom_to_leah = dict(zip(df["id"], df["mapping"]))

    # Dump the dictionary to a json file
    with open("tom_to_leah_mapping.json", "w") as f:
        json.dump(mapping_dict_from_tom_to_leah, f, indent=2)

    # Grep all the files in the current directory, limited to the desired file extension
    files = [
        f
        for f in os.listdir(original_folder)
        if os.path.isfile(os.path.join(original_folder, f))
        and f.endswith(file_extension)
    ]
    print(f"Found {len(files)} files with {file_extension} extension.")

    # Make a new directory to store the renamed files
    os.makedirs(new_folder, exist_ok=True)

    # Process and rename the files based on the mapping
    for file in files:
        if "_" in file:
            prefix = file.split("_")[0]
            # the rest stays the same
            rest = "_".join(file.split("_")[1:])

            print(f"---------------------------------------------")
            print(f"prefix: {prefix}, rest: {rest}")
            if prefix in mapping_dict_from_tom_to_leah:
                new_filename = mapping_dict_from_tom_to_leah[prefix] + "_" + rest
                with open(os.path.join(original_folder, file), "r") as f:
                    with open(os.path.join(new_folder, new_filename), "w") as nf:
                        for line in f:
                            nf.write(
                                line.replace(
                                    prefix, mapping_dict_from_tom_to_leah[prefix]
                                )
                            )
                print(f"{file} -> {new_filename}")
            else:
                print(f"{file} not found in mapping dictionary.")
        else:
            print(
                f"{file} is not in the expected naming convention (should contain an underscore)."
            )


if __name__ == "__main__":
    # Argument parser for handling command-line inputs
    parser = argparse.ArgumentParser(
        description="Convert and rename files based on a predefined mapping."
    )
    parser.add_argument(
        "original_folder", type=str, help="The folder containing the original files."
    )
    parser.add_argument(
        "new_folder", type=str, help="The folder to save renamed files."
    )
    parser.add_argument(
        "file_extension", type=str, help="The file extension to filter (e.g., '.com')."
    )

    args = parser.parse_args()

    main(args.original_folder, args.new_folder, args.file_extension)
