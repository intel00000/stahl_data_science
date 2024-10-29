import os
import filecmp
from difflib import unified_diff


def compare_files_in_subfolders(subfolder1, subfolder2):
    # Get the current working directory
    current_dir = os.getcwd()

    # Get the absolute paths of the two subfolders
    folder1_path = os.path.join(current_dir, subfolder1)
    folder2_path = os.path.join(current_dir, subfolder2)

    # Get the list of files in both subfolders
    files_in_folder1 = set(os.listdir(folder1_path))
    files_in_folder2 = set(os.listdir(folder2_path))

    # Find common files between both subfolders
    common_files = files_in_folder1.intersection(files_in_folder2)
    unique_to_folder1 = files_in_folder1 - files_in_folder2
    unique_to_folder2 = files_in_folder2 - files_in_folder1

    # Show files unique to each folder
    print("Files unique to", subfolder1, ":")
    for file_name in unique_to_folder1:
        print(file_name)

    print("Files unique to", subfolder2, ":")
    for file_name in unique_to_folder2:
        print(file_name)

    # Compare the content of each file with the same name
    for file_name in common_files:
        file1_path = os.path.join(folder1_path, file_name)
        file2_path = os.path.join(folder2_path, file_name)

        # Compare the contents of the files
        if filecmp.cmp(file1_path, file2_path, shallow=False):
            print(f"{file_name}: The files are identical.")
        else:
            print(f"{file_name}: The files are different.")
            show_file_differences(file1_path, file2_path)


def show_file_differences(file1_path, file2_path):
    # Open both files and compare their contents line by line
    with open(file1_path, "r") as file1, open(file2_path, "r") as file2:
        file1_lines = file1.readlines()
        file2_lines = file2.readlines()

        diff = unified_diff(
            file1_lines, file2_lines, fromfile=file1_path, tofile=file2_path
        )
        print("".join(diff))


if __name__ == "__main__":
    subfolder1 = "closed_shell_tom_old"
    subfolder2 = "closed_shell_tom_new"

    compare_files_in_subfolders(subfolder1, subfolder2)
