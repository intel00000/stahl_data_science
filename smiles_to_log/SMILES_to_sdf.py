# How to use on the CHPC:
# install any missing package like RDKit and openpyxl with:
#      run "pip install rdkit-pypi"
#      run "pip install openpyxl"
# load Python 3.7:
#      run "ml python/3.7"
# load OpenBabel:
#      run "ml openbabel"
# run this script by navigating into the directory with the Excel file "smiles.xlsx" that includes the SMILES codes and substrate IDs you want to make into .sdf files
#      run "python SMILES_to_sdf.py" or "python ~/bin/SMILES_to_sdf.py"

import pandas as pd
import shutil
import uuid
import getpass
import socket
import os
import sys
import numpy as np
import scipy.spatial as scsp
import copy
import time
import subprocess
import shlex
import yaml


import scipy.spatial as scsp
import scipy.linalg as scli
from rdkit import Chem
from rdkit.Chem import AllChem

from rdkit.Chem import MolFromSmiles as smi2mol
from rdkit.Chem import MolToSmiles as mol2smi

# os.system("module load openbabel")


def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def input_ext(argv, ext):
    # for a given extension, loads argv specified files, else all in folder
    files = []
    if len(argv) > 1:
        for i in range(1, len(argv)):
            try:
                if argv[i].split(".")[1] == ext:  # split extension
                    files.append(argv[i])
            except:  # if argument not of file format
                pass

    if (
        files == []
    ):  # if list remains empty, take all files in directory with given extension
        files = [file for file in glob("*.%s" % ext)]
        print("Examining " + str(len(files)) + " file(s). \n")
    assert files != [], "No files with %s extension in directory" % ext
    return files


def readXYZ(filename):
    infile = open(filename, "r")
    coords = []
    elements = []
    lines = infile.readlines()
    if len(lines) < 3:
        exit("ERROR: no coordinates found in %s/%s" % (os.getcwd(), filename))
    for line in lines[2:]:
        elements.append(line.split()[0].capitalize())
        coords.append(
            [float(line.split()[1]), float(line.split()[2]), float(line.split()[3])]
        )
    infile.close()
    coords = np.array(coords)
    return coords, elements


def get_coords_from_smiles(smiles, conversion_method, id):
    if conversion_method == "any":
        to_try = ["rdkit", "obabel"]
    elif conversion_method == "rdkit":
        to_try = ["rdkit"]
    elif conversion_method == "molconvert":
        to_try = ["molconvert"]
    elif conversion_method == "obabel":
        to_try = ["obabel"]

    error = ""
    for m in to_try:
        # print("   ---   try to convert %s to 3D using %s (of %s)"%(smiles, m, str(to_try)))
        if m == "molconvert":

            if which("molconvert") != None:

                coords, elements = get_coords_from_smiles_marvin(smiles, id)
                if coords is None or elements is None:
                    error += " molconvert_failed "
                    pass
                else:
                    if abs(np.max(coords.T[2]) - np.min(coords.T[2])) > 0.01:
                        # print("   ---   conversion done with molconvert")
                        return (coords, elements)
                    else:
                        error += " molconvert_mol_flat "
                        pass
                        # print("WARNING: molconvert produced a flat molecule. proceed with other methods (obabel or rdkit)")
            else:
                error += "molconvert_not_available "

        if m == "obabel":
            if which("obabel") != None:
                # print("use obabel")
                coords, elements = get_coords_from_smiles_obabel(smiles, id)
                if coords is None or elements is None:
                    error += " obabel_failed "
                    pass
                else:
                    if abs(np.max(coords.T[2]) - np.min(coords.T[2])) > 0.01:
                        # print("   ---   conversion done with obabel")
                        return (coords, elements)
                    else:
                        error += " obabel_failed "
                        pass

            else:
                error += " obabel_not_available "

        if m == "rdkit":
            # print("use rdkit")
            coords, elements = get_coords_from_smiles_rdkit(smiles, id)
            if coords is None or elements is None:
                error += " rdkit_failed "
                pass
            else:
                if abs(np.max(coords.T[2]) - np.min(coords.T[2])) > 0.01:
                    # print("   ---   conversion done with rdkit")
                    return (coords, elements)
                else:
                    error += " rdkit_failed "
                    pass
    print("ERROR: No 3D conversion worked: %s, %s" % (id, error))
    conversion_error.write("No 3D conversion worked: %s, %s \n" % (id, error))
    return (None, None)


def get_coords_from_smiles_obabel(smiles, id):
    name = uuid.uuid4()

    if not os.path.exists("input_structures"):
        try:
            os.makedirs("input_structures")
        except:
            pass

    filename = "input_structures/%s.xyz" % (name)
    os.system('obabel -:"%s" --gen3D -oxyz > %s' % (smiles, filename))
    if not os.path.exists(filename):
        return (None, None)
        # print("ERROR: could not convert %s to 3D using obabel. Exit!"%(smiles))
        # exit()

    coords, elements = readXYZ(filename)
    if len(coords) == 0:
        return (None, None)
        # print("ERROR: could not convert %s to 3D using obabel. Exit!"%(smiles))
        # exit()
    os.system("rm %s" % (filename))
    xyzfilename = id + ".xyz"
    exportXYZ(coords, elements, xyzfilename)
    os.system("obabel " + id + ".xyz -O " + id + ".sdf")

    # sdf_stem = sdf.split(".")[0]
    mult = 0
    with open(id + ".sdf", "r") as inp:
        lines = inp.readlines()
        for line in lines:
            if "RAD" in line:
                mult += 1
    if mult >= 1:
        os.system("rm " + id + ".sdf")
        os.system("rm " + id + ".xyz")
        return (None, None)

    return (coords, elements)


def get_coords_from_smiles_rdkit(smiles, id):
    try:
        m = Chem.MolFromSmiles(smiles)
    except:
        return (None, None)
        print(id)
        # print("could not convert %s to rdkit molecule. Exit!"%(smiles))
        # exit()
    try:
        m = Chem.AddHs(m)
    except:
        return (None, None)
        # print("ERROR: could not add hydrogen to rdkit molecule of %s. Exit!"%(smiles))
        # exit()
    try:
        AllChem.EmbedMolecule(m)
    except:
        return (None, None)
        # print("ERROR: could not calculate 3D coordinates from rdkit molecule %s. Exit!"%(smiles))
        # exit()
    try:
        block = Chem.MolToMolBlock(m)
        blocklines = block.split("\n")
        coords = []
        elements = []
        for line in blocklines[4:]:
            if len(line.split()) == 4:
                break
            elements.append(line.split()[3])
            coords.append(
                [float(line.split()[0]), float(line.split()[1]), float(line.split()[2])]
            )
        coords = np.array(coords)
        mean = np.mean(coords, axis=0)
        distances = scsp.distance.cdist([mean], coords)[0]
        if np.max(distances) < 0.1:
            return (None, None)
            # print("ERROR: something is wrong with rdkit molecule %s. Exit!"%(smiles))
            # print("%i\n"%(len(coords)))
            # for atomidx, atom in enumerate(coords):
            #    print("%s %f %f %f"%(elements[atomidx], atom[0], atom[1], atom[2]))
            # exit()
    except:
        return (None, None)
        # print("ERROR: could not read xyz coordinates from rdkit molecule %s. Exit!"%(smiles))
        # exit()

    xyzfilename = id + ".xyz"
    exportXYZ(coords, elements, xyzfilename)
    os.system("obabel " + id + ".xyz -O " + id + ".sdf")

    # sdf_stem = sdf.split(".")[0]
    mult = 0
    with open(id + ".sdf", "r") as inp:
        lines = inp.readlines()
        for line in lines:
            if "RAD" in line:
                mult += 1
    if mult >= 1:
        os.system("rm " + id + ".sdf")
        os.system("rm " + id + ".xyz")
        return (None, None)

    return (coords, elements)


def get_coords_from_smiles_marvin(smiles, id):

    name = uuid.uuid4()

    if not os.path.exists("tempfiles"):
        try:
            os.makedirs("tempfiles")
        except:
            pass
    if not os.path.exists("input_structures"):
        try:
            os.makedirs("input_structures")
        except:
            pass

    outfile = open("tempfiles/%s.smi" % (name), "w")
    outfile.write("%s\n" % (smiles))
    outfile.close()

    path_here = os.getcwd()
    os.system(
        "molconvert -2 mrv:+H %s/tempfiles/%s.smi > tempfiles/%s.mrv"
        % (path_here, name, name)
    )
    filename = "tempfiles/%s.mrv" % (name)
    if not os.path.exists(filename):
        os.system("rm tempfiles/%s.smi" % (name))
        return (None, None)
        # print("ERROR: could not convert %s to 2D (mrv) using marvin. Exit!"%(smiles))
        # exit()

    os.system(
        "molconvert -3 xyz %s/tempfiles/%s.mrv > input_structures/%s.xyz"
        % (path_here, name, name)
    )
    filename = "input_structures/%s.xyz" % (name)
    if not os.path.exists(filename):
        os.system("rm tempfiles/%s.smi tempfiles/%s.mrv" % (name, name))
        return (None, None)
        # print("ERROR: could not convert %s to 3D (xyz) using marvin. Exit!"%(smiles))
        # exit()

    coords, elements = readXYZ(filename)
    if len(coords) == 0:
        os.system(
            "rm tempfiles/%s.smi tempfiles/%s.mrv input_structures/%s.xyz"
            % (name, name, name)
        )
        # print("ERROR: could not convert %s to 3D (coords in empty) using marvin. Exit!"%(smiles))
        # return(None, None)
        # exit()
    os.system(
        "rm tempfiles/%s.smi tempfiles/%s.mrv input_structures/%s.xyz"
        % (name, name, name)
    )
    return (coords, elements)


def exportXYZ(coords, elements, filename, mask=[]):
    outfile = open(filename, "w")

    if len(mask) == 0:
        outfile.write("%i\n\n" % (len(elements)))
        for atomidx, atom in enumerate(coords):
            outfile.write(
                "%s %f %f %f\n"
                % (elements[atomidx].capitalize(), atom[0], atom[1], atom[2])
            )
    else:
        outfile.write("%i\n\n" % (len(mask)))
        for atomidx in mask:
            atom = coords[atomidx]
            outfile.write(
                "%s %f %f %f\n"
                % (elements[atomidx].capitalize(), atom[0], atom[1], atom[2])
            )
    outfile.close()


filename = "smiles.xlsx"
conversion_method = "any"

df = pd.read_excel(filename, engine="openpyxl")

with open("no_conversion.txt", "a") as conversion_error:
    for i, row in df.iterrows():
        smi = row["SMILES"]
        id = row["id"]
        print("Trying to convert " + str(id) + ":")
        get_coords_from_smiles(smi, conversion_method, id)

# if under linux, remove all .xyz files and empty directories
if sys.platform == "linux" or sys.platform == "linux2":
    os.system("rm *.xyz")
elif sys.platform == "win32":
    os.system("del *.xyz")
os.system("find . -type d -empty -delete")
