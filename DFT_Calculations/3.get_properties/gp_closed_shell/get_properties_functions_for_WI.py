import pandas as pd
import numpy as np
import re
import math
from morfeus import Sterimol
from morfeus import BuriedVolume
from morfeus import Pyramidalization
from morfeus import SASA
from morfeus import ConeAngle
from morfeus import Sterimol

import dbstep.Dbstep as db
from matplotlib import rcParams

import sys
import subprocess

# updated by Tom Tan at 2024.8.30
# 1. change all pandas dataframe append to using pd.concat
# see https://pandas.pydata.org/pandas-docs/version/1.4/whatsnew/v1.4.0.html#deprecated-dataframe-append-and-series-append
# Pandas ditch dataframe1.append(dataframe2) method bacause every you call it, it create a brand new dataframe take up the same amount of memory.
# Using pd.concat will append two dataframe "in-place" which take up less memory.
# 2. move the get_one_lp_energy from notebook to here.
# 3. update get_goodvibes_e to work with the newest version of goodvibes
# Goodvibes change significantly, now they don't support calling it from a function anymore.
# Instead what they do is to call goovibes in the command line with the log file and it outputs all the property.
# 4. adapt the code to work when the cell in the atom map is atom_label + atom_number (e.g. C1, N1) instead of just number.

homo_pattern = re.compile("Alpha  occ. eigenvalues")
npa_pattern = re.compile("Summary of Natural Population Analysis:")
nbo_os_pattern = re.compile("beta spin orbitals")
nmrstart_pattern = " SCF GIAO Magnetic shielding tensor (ppm):\n"
nmrend_pattern = re.compile("End of Minotr F.D.")
nmrend_pattern_os = re.compile("g value of the free electron")
zero_pattern = re.compile("zero-point Energies")
cputime_pattern = re.compile("Job cpu time:")
walltime_pattern = re.compile("Elapsed time:")
volume_pattern = re.compile("Molar volume =")
polarizability_pattern = re.compile("Dipole polarizability, Alpha")
dipole_pattern = "Dipole moment (field-independent basis, Debye)"
frqs_pattern = re.compile("Red. masses")
frqsend_pattern = re.compile("Thermochemistry")
chelpg1_pattern = re.compile("(CHELPG)")
chelpg2_pattern = re.compile("Charges from ESP fit")
hirshfeld_pattern = re.compile("Hirshfeld charges, spin densities, dipoles, and CM5 charges")
nborbs_pattern = "NATURAL BOND ORBITALS (Summary):" # "Natural Bond Orbitals (Summary)"
nborbs2_pattern = re.compile("NATURAL BOND ORBITALS (Summary):")
nrt_pattern = re.compile("NATURAL RESONANCE THEORY ANALYSIS:")
nav_pattern = re.compile("Natural Atomic Valencies, Electron Counts, and Charges:") # Natural Atomic Valencies keyword
nborder_pattern = re.compile("Natural Bond Order:") # Natural Bond Order keyword

def get_geom(streams): #extracts the geometry from the compressed stream
    geom = []
    for item in streams[-1][16:]:
        if item == "":
            break
        geom.append([item.split(",")[0],float(item.split(",")[-3]),float(item.split(",")[-2]),float(item.split(",")[-1])])
    return(geom)

def get_outstreams(log): #gets the compressed stream information at the end of a Gaussian job
    streams = []
    starts,ends = [],[]
    error = ""
    an_error = True
    try:
        with open(log+".log") as f:
            loglines = f.readlines()
    except:
        with open(log+".LOG") as f:
            loglines = f.readlines()

    # line in loglines[::-12]:
        #if "Normal termination" in line:
            #an_error = False
        #if an_error:
            #error = "****Failed or incomplete jobs for " + log + ".log"

    for i in range(len(loglines)):
        if "1\\1\\" in loglines[i]:
            starts.append(i)
        if "@" in loglines[i]:
            ends.append(i)
    #    if "Normal termination" in loglines[i]:
    #        error = ""


    #if len(starts) != len(ends) or len(starts) == 0: #probably redundant
        #error = "****Failed or incomplete jobs for " + log + ".log"
        #return(streams,error)
    for i in range(len(starts)):
        tmp = ""
        for j in range(starts[i],ends[i]+1,1):
            tmp = tmp + loglines[j][1:-1]
        streams.append(tmp.split("\\"))
    return(streams,error)

def get_filecont(log): #gets the entire job output
    error = "" #default unless "normal termination" is in file
    an_error = True
    with open(log+".log") as f:
        loglines = f.readlines()
    #for line in loglines[::-1]:
        #if "Normal termination" in line:
            #an_error = False
        #if an_error:
            #error = "****Failed or incomplete jobs for " + log + ".log"
    return(loglines, error)

def get_sterimol_morfeus(dataframe, sterimol_list): #uses morfeus to calculate sterimol L, B1, B5 for two input atoms for every entry in df
    sterimol_dataframe = pd.DataFrame(columns=[])

    for index, row in dataframe.iterrows():
        try:
            # Parsing the Sterimol axis defined in the list from input line
            sterimolnums_list = []
            for sterimol in sterimol_list:
                atomnum_list = [] # The atom numbers used to collect sterimol values (i.e. [18 16 17 15]) are collected from the df using the input list (i.e. [["O2", "C1"], ["O3", "H5"]])
                for atom in sterimol:
                    atomnum = row[str(atom)]
                    atomnum = re.findall(r'\d+', atomnum)[0] # Extract only the number from the atom label
                    atomnum_list.append(str(atomnum))
                sterimolnums_list.append(atomnum_list) # Append atomnum_list for each sterimol axis defined in the input to make a list of the form [['18', '16'], ['16', '15']]

            # This makes column headers based on Sterimol axis defined in the input line
            sterimoltitle_list = []
            for sterimol in sterimol_list:
                sterimoltitle = str(sterimol[0]) + "_" + str(sterimol[1])
                sterimoltitle_list.append(sterimoltitle)

            log_file = row['log_name']
            streams, error = get_outstreams(log_file) # Need to add file path if you're running from a different directory than file
            if error != "":
                print(error)
                row_i = {}
                for a in range(len(sterimolnums_list)):
                    entry = {'Sterimol_L_' + str(sterimoltitle_list[a]) + '(Å)_morfeus': "no data",
                             'Sterimol_B1_' + str(sterimoltitle_list[a]) + '(Å)_morfeus': "no data",
                             'Sterimol_B5_' + str(sterimoltitle_list[a]) + '(Å)_morfeus': "no data"}
                    row_i.update(entry)
                sterimol_dataframe = pd.concat([sterimol_dataframe, pd.DataFrame([row_i])], ignore_index=True)
                continue

            geom = get_geom(streams)

            # Checks for if the wrong number of atoms are input, input is not of the correct form, or calls atom numbers that do not exist in the molecule
            error = ""
            for sterimol in sterimolnums_list:
                if len(sterimol) % 2 != 0:
                    error = f"Number of atom inputs given for Sterimol is not divisible by two. {len(sterimol)} atoms were given."
                for atom in sterimol:
                    if not atom.isdigit():
                        error += f" {atom}: Only numbers accepted as input for Sterimol."
                    if int(atom) > len(geom):
                        error += f" {atom} is out of range. Maximum valid atom number: {len(geom)}."
                if error:
                    print(error)

            elements = np.array([geom[i][0] for i in range(len(geom))])
            coordinates = np.array([np.array(geom[i][1:]) for i in range(len(geom))])

            # This collects Sterimol values for each pair of inputs
            sterimolout = []
            for sterimol in sterimolnums_list:
                sterimol_values = Sterimol(elements, coordinates, int(sterimol[0]), int(sterimol[1])) # Calls morfeus
                sterimolout.append(sterimol_values)

            # This adds the data from sterimolout into the new property df
            row_i = {}
            for a in range(len(sterimolnums_list)):
                entry = {'Sterimol_L_' + str(sterimoltitle_list[a]) + '(Å)_morfeus': sterimolout[a].L_value,
                         'Sterimol_B1_' + str(sterimoltitle_list[a]) + '(Å)_morfeus': sterimolout[a].B_1_value,
                         'Sterimol_B5_' + str(sterimoltitle_list[a]) + '(Å)_morfeus': sterimolout[a].B_5_value}
                row_i.update(entry)
            sterimol_dataframe = pd.concat([sterimol_dataframe, pd.DataFrame([row_i])], ignore_index=True)
        except:
            print(f'****Unable to acquire Morfeus Sterimol parameters for: {row["log_name"]}.log')
            row_i = {}
            try:
                for a in range(len(sterimolnums_list)):
                    entry = {'Sterimol_L_' + str(sterimoltitle_list[a]) + '(Å)_morfeus': "no data",
                             'Sterimol_B1_' + str(sterimoltitle_list[a]) + '(Å)_morfeus': "no data",
                             'Sterimol_B5_' + str(sterimoltitle_list[a]) + '(Å)_morfeus': "no data"}
                    row_i.update(entry)
                sterimol_dataframe = pd.concat([sterimol_dataframe, pd.DataFrame([row_i])], ignore_index=True)
            except:
                print("****Ope, there's a problem with your atom inputs.")
    print(f"Morfeus Sterimol function has completed for {sterimol_list}")
    return pd.concat([dataframe, sterimol_dataframe], axis=1)

def get_sterimol_dbstep(dataframe, sterimol_list): #uses DBSTEP to calculate sterimol L, B1, B5 for two input atoms for every entry in df
    sterimol_dataframe = pd.DataFrame(columns=[])

    for index, row in dataframe.iterrows():
        try:
            log_file = row['log_name']

            # Parsing the Sterimol axis defined in the list from input line
            sterimolnums_list = []
            for sterimol in sterimol_list:
                atomnum_list = [] # The atom numbers used to collect sterimol values (i.e. [18 16 17 15]) are collected from the df using the input list (i.e. [["O2", "C1"], ["O3", "H5"]])
                for atom in sterimol:
                    atomnum = row[str(atom)]
                    atomnum = re.findall(r'\d+', atomnum)[0] # Extract only the number from the atom label
                    atomnum_list.append(str(atomnum))
                sterimolnums_list.append(atomnum_list) # Append atomnum_list for each sterimol axis defined in the input to make a list of the form [['18', '16'], ['16', '15']]

            # Checks for if the wrong number of atoms are input or input is not of the correct form
            error = ""
            for sterimol in sterimolnums_list:
                if len(sterimol) % 2 != 0:
                    error = f"****Number of atom inputs given for Sterimol is not divisible by two. {len(sterimol)} atoms were given."
                for atom in sterimol:
                    if not atom.isdigit():
                        error += f"**** {atom}: Only numbers accepted as input for Sterimol."
                if error:
                    print(error)

            # This collects Sterimol values for each pair of inputs
            sterimol_out = []
            fp = log_file + ".log"
            for sterimol in sterimolnums_list:
                sterimol_values = db.dbstep(fp, atom1=int(sterimol[0]), atom2=int(sterimol[1]), commandline=True, verbose=False, sterimol=True, measure='grid')
                sterimol_out.append(sterimol_values)

            # This makes column headers based on Sterimol axis defined in the input line
            sterimoltitle_list = []
            for sterimol in sterimol_list:
                sterimoltitle = str(sterimol[0]) + "_" + str(sterimol[1])
                sterimoltitle_list.append(sterimoltitle)

            # This adds the data from sterimol_out into the new property df
            row_i = {}
            for a in range(len(sterimolnums_list)):
                entry = {'Sterimol_B1_' + str(sterimoltitle_list[a]) + "(Å)_dbstep": sterimol_out[a].Bmin,
                         'Sterimol_B5_' + str(sterimoltitle_list[a]) + "(Å)_dbstep": sterimol_out[a].Bmax,
                         'Sterimol_L_' + str(sterimoltitle_list[a]) + "(Å)_dbstep": sterimol_out[a].L}
                row_i.update(entry)
            sterimol_dataframe = pd.concat([sterimol_dataframe, pd.DataFrame([row_i])], ignore_index=True)
        except:
            print(f'****Unable to acquire DBSTEP Sterimol parameters for: {row["log_name"]}.log')
            row_i = {}
            try:
                for a in range(len(sterimolnums_list)):
                    entry = {'Sterimol_L_' + str(sterimoltitle_list[a]) + '(Å)_dbstep': "no data",
                             'Sterimol_B1_' + str(sterimoltitle_list[a]) + '(Å)_dbstep': "no data",
                             'Sterimol_B5_' + str(sterimoltitle_list[a]) + '(Å)_dbstep': "no data"}
                    row_i.update(entry)
                sterimol_dataframe = pd.concat([sterimol_dataframe, pd.DataFrame([row_i])], ignore_index=True)
            except:
                print("****Ope, there's a problem with your atom inputs.")
    print(f"DBSTEP Sterimol function has completed for {sterimol_list}")
    return pd.concat([dataframe, sterimol_dataframe], axis=1)

def get_sterimol2vec(dataframe, sterimol_list, end_r, step_size): #uses DBSTEP to calculate sterimol Bmin and Bmax for two input atoms at intervals from 0 to end_r at step_size
    sterimol_dataframe = pd.DataFrame(columns=[])
    num_steps = int((end_r) / step_size + 1)
    radii_list = [0 + step_size * i for i in range(num_steps)]

    for index, row in dataframe.iterrows():
        try:
            log_file = row['log_name']

            # Parsing the Sterimol axis defined in the list from input line
            sterimolnums_list = []
            for sterimol in sterimol_list:
                atomnum_list = [] # The atom numbers used to collect sterimol values (i.e. [18 16 17 15]) are collected from the df using the input list (i.e. [["O2", "C1"], ["O3", "H5"]])
                for atom in sterimol:
                    atomnum = row[str(atom)]
                    atomnum = re.findall(r'\d+', atomnum)[0] # Extract only the number from the atom label
                    atomnum_list.append(str(atomnum))
                sterimolnums_list.append(atomnum_list) # Append atomnum_list for each sterimol axis defined in the input to make a list of the form [['18', '16'], ['16', '15']]

            # Checks for if the wrong number of atoms are input or input is not of the correct form
            error = ""
            for sterimol in sterimolnums_list:
                if len(sterimol) % 2 != 0:
                    error = f"Number of atom inputs given for Sterimol is not divisible by two. {len(sterimol)} atoms were given."
                for atom in sterimol:
                    if not atom.isdigit():
                        error += f" {atom}: Only numbers accepted as input for Sterimol."
                if error:
                    print(error)

            # This collects Sterimol values for each pair of inputs
            sterimol2vec_out = []
            fp = log_file + ".log"
            for sterimol in sterimolnums_list:
                sterimol2vec_values = db.dbstep(fp, atom1=int(sterimol[0]), atom2=int(sterimol[1]), scan='0.0:{}:{}'.format(end_r, step_size), commandline=True, verbose=False, sterimol=True, measure='grid')
                sterimol2vec_out.append(sterimol2vec_values)

            # This makes column headers based on Sterimol axis defined in the input line
            sterimoltitle_list = []
            for sterimol in sterimol_list:
                sterimoltitle = str(sterimol[0]) + "_" + str(sterimol[1])
                sterimoltitle_list.append(sterimoltitle)

            scans = radii_list
            # This adds the data from sterimol2vec_out into the new property df
            row_i = {}
            for a in range(len(sterimolnums_list)):
                for i in range(len(scans)):
                    entry = {'Sterimol_Bmin_' + str(sterimoltitle_list[a]) + "_" + str(scans[i]) + "Å(Å)": sterimol2vec_out[a].Bmin[i],
                             'Sterimol_Bmax_' + str(sterimoltitle_list[a]) + "_" + str(scans[i]) + "Å(Å)": sterimol2vec_out[a].Bmax[i]}
                    row_i.update(entry)
            sterimol_dataframe = pd.concat([sterimol_dataframe, pd.DataFrame([row_i])], ignore_index=True)
        except:
            print(f'****Unable to acquire DBSTEP Sterimol2Vec parameters for: {row["log_name"]}.log')
            row_i = {}
            try:
                for a in range(len(sterimolnums_list)):
                    for i in range(len(scans)):
                        entry = {'Sterimol_Bmin_' + str(sterimoltitle_list[a]) + "_" + str(scans[i]) + "Å(Å)": "no data",
                                 'Sterimol_Bmax_' + str(sterimoltitle_list[a]) + "_" + str(scans[i]) + "Å(Å)": "no data"}
                        row_i.update(entry)
                sterimol_dataframe = pd.concat([sterimol_dataframe, pd.DataFrame([row_i])], ignore_index=True)
            except:
                print("****Ope, there's a problem with your atom inputs.")
    print(f"DBSTEP Sterimol2Vec function has completed for {sterimol_list}")
    return pd.concat([dataframe, sterimol_dataframe], axis=1)

def get_vbur_one_radius(dataframe, a1, radius): #uses morfeus to calculate vbur at a single radius for an atom (a1) in df
    atom = str(a1)
    vbur_dataframe = pd.DataFrame(columns=[])

    for index, row in dataframe.iterrows():
        try:
            log_file = row['log_name']
            atom1 = row[str(a1)] #gets numerical value (e.g. 16) for a1 (e.g. C1)
            atom1 = re.findall(r'\d+', atom1)[0] # extract only the number from the atom label
            streams, error = get_outstreams(log_file) #need to add file path if you're running from a different directory than file
            if error != "":
                print(error)
                row_i = {'%Vbur_'+str(atom)+"_"+str(radius)+"Å": "no data"}
                vbur_dataframe = pd.concat([vbur_dataframe, pd.DataFrame([row_i])], ignore_index=True)
                continue

            log_coordinates = get_geom(streams)
            elements = np.array([log_coordinates[i][0] for i in range(len(log_coordinates))])
            coordinates = np.array([np.array(log_coordinates[i][1:]) for i in range(len(log_coordinates))])
            vbur = BuriedVolume(elements, coordinates, int(atom1), include_hs=True, radius=radius) #calls morfeus
            row_i = {'%Vbur_'+str(atom)+"_"+str(radius)+"Å": vbur.percent_buried_volume * 100}
            vbur_dataframe = pd.concat([vbur_dataframe, pd.DataFrame([row_i])], ignore_index=True)
        except:
            print(f'****Unable to acquire Vbur parameters for: {row["log_name"]}.log')
            row_i = {'%Vbur_'+str(atom)+"_"+str(radius)+"Å": "no data"}
            vbur_dataframe = pd.concat([vbur_dataframe, pd.DataFrame([row_i])], ignore_index=True)
    return vbur_dataframe

def get_vbur_scan(dataframe, a_list, start_r, end_r, step_size): #uses morfeus via get_vbur_one_radius to scan vbur across a range of radii
    num_steps = int((end_r-start_r)/step_size + 1)
    radii = [start_r + step_size*i for i in range(num_steps)]
    frames = []
    for radius in radii:
        for a in a_list:
            frames.append(get_vbur_one_radius(dataframe, a, radius))
    vbur_scan_dataframe = pd.concat(frames, axis=1)
    print(f"Vbur scan function has completed for {a_list} from {start_r} to {end_r}")
    return pd.concat([dataframe, vbur_scan_dataframe], axis=1)

def get_pyramidalization(dataframe, a_list): #uses morfeus to calculate pyramidalization for all atoms in a_list (form ["C1", "C4", "O2"]) in a dataframe that contains file name and atom number
    pyr_dataframe = pd.DataFrame(columns=[])

    for index, row in dataframe.iterrows():
        try:
            atom_list = []
            for label in a_list:
                atom = row[str(label)] # the atom number (i.e. 16) to add to the list is the df entry of this row for the labeled atom (i.e. "C1")
                atom = re.findall(r'\d+', atom)[0] # extract only the number from the atom label
                atom_list.append(str(atom)) # append that to atom_list to make a list of the form [16, 17, 29]

            log_file = row['log_name']
            streams, error = get_outstreams(log_file) # need to add file path if you're running from a different directory than file
            if error != "":
                print(error)
                row_i = {}
                for a in range(len(atom_list)):
                    entry = {'pyramidalization_Gavrish_' + str(a_list[a]) + '(°)': "no data",
                             'pyramidalization_Agranat-Radhakrishnan_' + str(a_list[a]): "no data"} # details on these values can be found here: https://kjelljorner.github.io/morfeus/pyramidalization.html
                    row_i.update(entry)
                pyr_dataframe = pd.concat([pyr_dataframe, pd.DataFrame([row_i])], ignore_index=True)
                continue

            log_coordinates = get_geom(streams)
            elements = np.array([log_coordinates[i][0] for i in range(len(log_coordinates))])
            coordinates = np.array([np.array(log_coordinates[i][1:]) for i in range(len(log_coordinates))])

            pyrout = []
            for atom in atom_list:
                pyr = Pyramidalization(coordinates, int(atom)) # calls morfeus
                pyrout.append(pyr)

            row_i = {}
            for a in range(len(atom_list)):
                entry = {'pyramidalization_Gavrish_' + str(a_list[a]) + '(°)': pyrout[a].P_angle,
                         'pyramidalization_Agranat-Radhakrishnan_' + str(a_list[a]): pyrout[a].P} # details on these values can be found here: https://kjelljorner.github.io/morfeus/pyramidalization.html
                row_i.update(entry)
            pyr_dataframe = pd.concat([pyr_dataframe, pd.DataFrame([row_i])], ignore_index=True)
        except:
            print(f'****Unable to acquire pyramidalization parameters for: {row["log_name"]}.log')
            row_i = {}
            for a in range(len(atom_list)):
                entry = {'pyramidalization_Gavrish_' + str(a_list[a]) + '(°)': "no data",
                         'pyramidalization_Agranat-Radhakrishnan_' + str(a_list[a]): "no data"} # details on these values can be found here: https://kjelljorner.github.io/morfeus/pyramidalization.html
                row_i.update(entry)
            pyr_dataframe = pd.concat([pyr_dataframe, pd.DataFrame([row_i])], ignore_index=True)
    print(f"Pyramidalization function has completed for {a_list}")
    return pd.concat([dataframe, pyr_dataframe], axis=1)

def get_specdata(atoms, prop): #input a list of atom numbers of interest and a list of pairs of all atom numbers and property of interest for use with NMR, NBO, possibly others with similar output structures
    propout = []
    for atom in atoms:
        if atom.isdigit():
            a = int(atom)-1
            if a <= len(prop):
                propout.append(float(prop[a][1]))
            else: continue
        else: continue
    return(propout)

def get_nbo(dataframe, a_list): #a function to get the nbo for all atoms (a_list, form ["C1", "C4", "O2"]) in a dataframe that contains file name and atom number
    nbo_dataframe = pd.DataFrame(columns=[]) #define an empty df to place results in

    for index, row in dataframe.iterrows(): #iterate over the dataframe
        try: #try to get the data
            atomnum_list = []
            for atom in a_list:
                atomnum = row[str(atom)] # the atom number (i.e. 16) to add to the list is the df entry of this row for the labeled atom (i.e. "C1")
                atomnum = re.findall(r'\d+', atomnum)[0] # extract only the number from the atom label
                atomnum_list.append(str(atomnum)) # append that to atomnum_list to make a list of the form [16, 17, 29]

            log_file = row['log_name'] #read file name from df
            filecont, error = get_filecont(log_file) #read the contents of the log file
            if error != "":
                print(error)
                row_i = {}
                for a in range(len(a_list)):
                    entry = {'NBO_charge_'+str(a_list[a]): "no data"}
                    row_i.update(entry)
                nbo_dataframe = pd.concat([nbo_dataframe, pd.DataFrame([row_i])], ignore_index=True)
                continue

            nbo, nbostart, nboout, skip = [], 0, "", 0
            #this section finds the line (nbostart) where the nbo data is located
            for i in range(len(filecont)-1, 0, -1): #search the file contents for the phrase "beta spin orbitals" to check for open shell molecules
                if re.search(nbo_os_pattern, filecont[i]) and skip == 0:
                    skip = 2 # retrieve only combined orbitals NPA in open shell molecules
                if npa_pattern.search(filecont[i]): #search the file content for the phrase which indicates the start of the NBO section
                    if skip != 0:
                        skip = skip-1
                        continue
                    nbostart = i + 6 #skips the set number of lines between the search key and the start of the table
                    break
            if nbostart == 0:
                error = f"****no Natural Population Analysis found in: {row['log_name']}.log"
                print(error)
                row_i = {}
                for a in range(len(a_list)):
                    entry = {'NBO_charge_'+str(a_list[a]): "no data"}
                    row_i.update(entry)
                nbo_dataframe = pd.concat([nbo_dataframe, pd.DataFrame([row_i])], ignore_index=True)
                continue

            #this section splits the table where nbo data is located into just the atom number and charge to generate a list of lists (nbo)
            ls = []
            for line in filecont[nbostart:]:
                if "==" in line: break
                ls = [str.split(line)[1], str.split(line)[2]]
                nbo.append(ls)

            #this uses the nbo list to return only the charges for only the atoms of interest as a list (nboout)
            nboout = get_specdata(atomnum_list, nbo)

            #this adds the data from the nboout into the new property df
            row_i = {}
            for a in range(len(a_list)):
                entry = {'NBO_charge_'+str(a_list[a]): nboout[a]}
                row_i.update(entry)
            nbo_dataframe = pd.concat([nbo_dataframe, pd.DataFrame([row_i])], ignore_index=True)
        except:
            print(f'****Unable to acquire NBO charges for: {row["log_name"]}.log')
            row_i = {}
            for a in range(len(a_list)):
                entry = {'NBO_charge_'+str(a_list[a]): "no data"}
                row_i.update(entry)
            nbo_dataframe = pd.concat([nbo_dataframe, pd.DataFrame([row_i])], ignore_index=True)
    print(f"NBO function has completed for {a_list}")
    return pd.concat([dataframe, nbo_dataframe], axis=1)

def get_nmr(dataframe, a_list): # a function to get the nmr for all atoms (a_list, form ["C1", "C4", "O2"]) in a dataframe that contains file name and atom number
    nmr_dataframe = pd.DataFrame(columns=[]) #define an empty df to place results in

    for index, row in dataframe.iterrows(): #iterate over the dataframe

        try: #try to get the data
            atom_list = []
            for new_a in a_list:
                new_atom = row[str(new_a)] #the atom number (i.e. 16) to add to the list is the df entry of this row for the labeled atom (i.e.) "C1")
                new_atom = re.findall(r'\d+', new_atom)[0] # extract only the number from the atom label
                atom_list.append(str(new_atom)) #append that to atom_list to make a list of the form [16, 17, 29]
            log_file = row['log_name'] #read file name from df
            filecont, error = get_filecont(log_file) #read the contents of the log file
            if error != "":
                print(error)
                row_i = {}
                for a in range(0, len(a_list)):
                    entry = {'NMR_shift_'+str(a_list[a]): "no data"}
                    row_i.update(entry)
                nmr_dataframe = pd.concat([nmr_dataframe, pd.DataFrame([row_i])], ignore_index=True)
                continue

            #determining the locations/values for start and end of NMR section
            start, end, i = 0, 0, 0
            if nmrstart_pattern in filecont:
                start = filecont.index(nmrstart_pattern)+1
                for i in range(start, len(filecont), 1):
                    if nmrend_pattern.search(filecont[i]) or nmrend_pattern_os.search(filecont[i]):
                        end = i
                        break
            if start == 0:
                error = f"****no NMR data found in file: {row['log_name']}.log"
                print(error)
                row_i = {}
                for a in range(0, len(a_list)):
                    entry = {'NMR_shift_'+str(a_list[a]): "no data"}
                    row_i.update(entry)
                nmr_dataframe = pd.concat([nmr_dataframe, pd.DataFrame([row_i])], ignore_index=True)
                continue

            atoms = int((end - start) / 5) #total number of atoms in molecule (there are 5 lines generated per atom)
            nmr = []
            for atom in range(atoms):
                element = str.split(filecont[start+5*atom])[1]
                shift_s = str.split(filecont[start+5*atom])[4]
                nmr.append([element, shift_s])
            nmrout = get_specdata(atom_list, nmr) #revisit
            #print(nmrout)

            #this adds the data from the nmrout into the new property df
            row_i = {}
            for a in range(0, len(a_list)):
                entry = {'NMR_shift_'+str(a_list[a]): nmrout[a]}
                row_i.update(entry)
            nmr_dataframe = pd.concat([nmr_dataframe, pd.DataFrame([row_i])], ignore_index=True)
        except:
            print(f'****Unable to acquire NMR shifts for: {row["log_name"]}.log')
            row_i = {}
            for a in range(0, len(a_list)):
                entry = {'NMR_shift_'+str(a_list[a]): "no data"}
                row_i.update(entry)
            nmr_dataframe = pd.concat([nmr_dataframe, pd.DataFrame([row_i])], ignore_index=True)
    print("NMR function has completed for", a_list)
    return pd.concat([dataframe, nmr_dataframe], axis = 1)

def get_angles(dataframe, angle_list): # a function to get the angles for all atoms (angle_list, form [[O3, C1, O2], [C4, C1, O3]]) in a dataframe that contains file name and atom number
    angle_dataframe = pd.DataFrame(columns=[]) #define an empty df to place results in

    for index, row in dataframe.iterrows(): #iterate over the dataframe
        try:
            #parsing the angle list from input line
            anglenums_list = []
            for angle in angle_list:
                atomnum_list = [] #the atom numbers for an angle (i.e. 17 16 18) are collected from the df using the input list (i.e.["O3", "C1", "O2"])
                for atom in angle:
                    atomnum = row[str(atom)]
                    atomnum = re.findall(r'\d+', atomnum)[0] # extract only the number from the atom label
                    atomnum_list.append(str(atomnum))
                anglenums_list.append(atomnum_list) #append atomnum_list for each angle to make a list of the form [['17', '16', '18'], ['15', '16', '17']]

            angletitle_list = []
            for angle in angle_list:
                angletitle = str(angle[0]) + "_" + str(angle[1]) + "_" + str(angle[2])
                angletitle_list.append(angletitle)

            log_file = row['log_name'] #read file name from df
            streams, error = get_outstreams(log_file)
            if error != "":
                print(error)
                row_i = {}
                for a in range(len(anglenums_list)):
                    entry = {'angle_'+str(angletitle_list[a]) + '(°)': "no data"}
                    row_i.update(entry)
                angle_dataframe = pd.concat([angle_dataframe, pd.DataFrame([row_i])], ignore_index=True)
                continue

            geom = get_geom(streams)

            # checks for if the wrong number of atoms are input, input is not of the correct form, or calls atom numbers that do not exist in the molecule.
            error = ""
            for angle in anglenums_list:
                if len(angle) % 3 != 0:
                    error = f"****Number of atom inputs given for angle is not divisible by three. {len(angle)} atoms were given."
                for atom in angle:
                    if not atom.isdigit():
                        error += f"**** {atom}: Only numbers accepted as input for angles"
                    if int(atom) > len(geom):
                        error += f"**** {atom} is out of range. Maximum valid atom number: {len(geom)}"
                if error:
                    print(error)

            anglesout = []
            for angle in anglenums_list:
                a = geom[int(angle[0])-1][:4] # atom coords
                b = geom[int(angle[1])-1][:4]
                c = geom[int(angle[2])-1][:4]
                ba = np.array(a[1:]) - np.array(b[1:])
                bc = np.array(c[1:]) - np.array(b[1:])
                cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
                anglevalue = np.arccos(cosine_angle)

                anglesout.append(float(round(np.degrees(anglevalue), 3)))

            #this adds the data from the anglesout into the new property df
            row_i = {}
            for a in range(len(anglenums_list)):
                entry = {'angle_'+str(angletitle_list[a]) + '(°)': anglesout[a]}
                row_i.update(entry)
            angle_dataframe = pd.concat([angle_dataframe, pd.DataFrame([row_i])], ignore_index=True)
        except Exception as e:
            print(f"E: {e}")
            print(f'****Unable to acquire angles for: {row["log_name"]}.log')
            row_i = {}
            try:
                for a in range(len(anglenums_list)):
                    entry = {'angle_'+str(angletitle_list[a]) + '(°)': "no data"}
                    row_i.update(entry)
                angle_dataframe = pd.concat([angle_dataframe, pd.DataFrame([row_i])], ignore_index=True)
            except Exception as e_inner:
                print(f"****Ope, there's a problem with your atom inputs: {e_inner}")
    print(f"Angles function has completed for {angle_list}")
    return pd.concat([dataframe, angle_dataframe], axis=1)

def get_dihedral(dataframe, dihedral_list): # a function to get the dihedrals for all atoms (dihedral_list, form [[O2, C1, O3, H5], [C4, C1, O3, H5]]) in a dataframe that contains file name and atom number
    dihedral_dataframe = pd.DataFrame(columns=[]) #define an empty df to place results in

    for index, row in dataframe.iterrows(): #iterate over the dataframe
        try:
            #parsing the dihedral list from input line
            dihedralnums_list = []
            for dihedral in dihedral_list:
                atomnum_list = [] #the atom numbers for a dihedral (i.e. 18 16 17 50) are collected from the df using the input list (i.e.["O2", "C1", "O3", "H5"])
                for atom in dihedral:
                    atomnum = row[str(atom)]
                    atomnum = re.findall(r'\d+', atomnum)[0] # extract only the number from the atom label
                    atomnum_list.append(str(atomnum))
                dihedralnums_list.append(atomnum_list) #append atomnum_list for each dihedral to make a list of the form [['18', '16', '17', '50'], ['18', '16', '17', '50']]
            dihedraltitle_list = []
            for dihedral in dihedral_list:
                dihedraltitle = str(dihedral[0]) + "_" + str(dihedral[1]) + "_" + str(dihedral[2]) + "_" + str(dihedral[3])
                dihedraltitle_list.append(dihedraltitle)

            log_file = row['log_name'] #read file name from df
            streams, error = get_outstreams(log_file)
            if error != "":
                print(error)
                row_i = {}
                for a in range(len(dihedralnums_list)):
                    entry = {'dihedral_'+str(dihedraltitle_list[a]) + '(°)': "no data"}
                    row_i.update(entry)
                dihedral_dataframe = pd.concat([dihedral_dataframe, pd.DataFrame([row_i])], ignore_index=True)
                continue
            geom = get_geom(streams)

            #checks for if the wrong number of atoms are input, input is not of the correct form, or calls atom numbers that do not exist in the molecule.
            error = ""
            for dihedral in dihedralnums_list:
                if len(dihedral) % 4 != 0:
                    error = f"****Number of atom inputs given for dihedral angle is not divisible by four. {len(dihedral)} atoms were given."
                for atom in dihedral:
                    if not atom.isdigit():
                        error += f"**** {atom}: Only numbers accepted as input for dihedral angles."
                    if int(atom) > len(geom):
                        error += f"**** {atom} is out of range. Maximum valid atom number: {len(geom)}."
                if error:
                    print(error)

            dihedralsout = []
            for dihedral in dihedralnums_list:
                a = geom[int(dihedral[0])-1][:4] # atom coords
                b = geom[int(dihedral[1])-1][:4]
                c = geom[int(dihedral[2])-1][:4]
                d = geom[int(dihedral[3])-1][:4]

                ab = np.array([a[1]-b[1], a[2]-b[2], a[3]-b[3]]) # vectors
                bc = np.array([b[1]-c[1], b[2]-c[2], b[3]-c[3]])
                cd = np.array([c[1]-d[1], c[2]-d[2], c[3]-d[3]])

                n1 = np.cross(ab, bc) # normal vectors
                n2 = np.cross(bc, cd)

                dihedral = round(np.degrees(np.arccos(np.dot(n1, n2) / (np.linalg.norm(n1) * np.linalg.norm(n2)))), 3)
                dihedralsout.append(float(dihedral))

            #this adds the data from the dihedralsout into the new property df
            row_i = {}
            for a in range(len(dihedralnums_list)):
                entry = {'dihedral_'+str(dihedraltitle_list[a]) + '(°)': dihedralsout[a]}
                row_i.update(entry)
            dihedral_dataframe = pd.concat([dihedral_dataframe, pd.DataFrame([row_i])], ignore_index=True)
        except Exception as e:
            print(f'****Unable to acquire dihedral angles for: {row["log_name"]}.log with error: {e}')
            row_i = {}
            try:
                for a in range(len(dihedralnums_list)):
                    entry = {'dihedral_'+str(dihedraltitle_list[a]) + '(°)': "no data"}
                    row_i.update(entry)
                dihedral_dataframe = pd.concat([dihedral_dataframe, pd.DataFrame([row_i])], ignore_index=True)
            except Exception as e_inner:
                print(f"****Ope, there's a problem with your atom inputs: {e_inner}")
    print(f"Dihedral function has completed for {dihedral_list}")
    return pd.concat([dataframe, dihedral_dataframe], axis=1)

def get_distance(dataframe, dist_list): # a function to get the distances for all atoms (dist_list, form [[C1, O2], [C4, C1]]) in a dataframe that contains file name and atom number
    dist_dataframe = pd.DataFrame(columns=[]) #define an empty df to place results in

    for index, row in dataframe.iterrows(): #iterate over the dataframe
        try:
            #parsing the distances list from input line
            distnums_list = []
            for dist in dist_list:
                atomnum_list = [] #the atom numbers for a distance (i.e. 18 16 16 15) are collected from the df using the input list (i.e.["O2", "C1", "O3", "H5"])
                for atom in dist:
                    atomnum = row[str(atom)]
                    atomnum = re.findall(r'\d+', atomnum)[0] # extract only the number from the atom label
                    atomnum_list.append(str(atomnum))
                distnums_list.append(atomnum_list) #append atomnum_list for each distance to make a list of the form [['18', '16'], ['16', '15']]

            disttitle_list = []
            for dist in dist_list:
                disttitle = str(dist[0]) + "_" + str(dist[1])
                disttitle_list.append(disttitle)

            log_file = row['log_name'] #read file name from df
            streams, error = get_outstreams(log_file)
            if error != "":
                print(error)
                row_i = {}
                for a in range(0, len(distnums_list)):
                    entry = {'distance_' + str(disttitle_list[a]) + '(Å)': "no data"}
                    row_i.update(entry)
                dist_dataframe = pd.concat([dist_dataframe, pd.DataFrame([row_i])], ignore_index=True)
                continue
            geom = get_geom(streams)

            #checks for if the wrong number of atoms are input, input is not of the correct form, or calls atom numbers that do not exist in the molecule.
            error = ""
            for dist in distnums_list:
                if len(dist) % 2 != 0:
                    error = f"****Number of atom inputs given for distance is not divisible by two. {len(dist)} atoms were given."
                for atom in dist:
                    if not atom.isdigit():
                        error += f"**** {atom}: Only numbers accepted as input for distances"
                    if int(atom) > len(geom):
                        error += f"**** {atom} is out of range. Maximum valid atom number: {len(geom)} "
                if error != "":
                    print(error)

            distout = []
            for dist in distnums_list:
                a = geom[int(dist[0])-1][:4] # Atomcoords
                b = geom[int(dist[1])-1][:4]
                ba = np.array(a[1:]) - np.array(b[1:])
                dist = round(np.linalg.norm(ba), 5)
                distout.append(float(dist))

            #this adds the data from the distout into the new property df
            row_i = {}
            for a in range(0, len(distnums_list)):
                entry = {'distance_' + str(disttitle_list[a]) + '(Å)': distout[a]}
                row_i.update(entry)
            dist_dataframe = pd.concat([dist_dataframe, pd.DataFrame([row_i])], ignore_index=True)
        except:
            print(f'****Unable to acquire distance for: {row["log_name"]}.log')
            row_i = {}
            try:
                for a in range(0, len(distnums_list)):
                    entry = {'distance_' + str(disttitle_list[a]) + '(Å)': "no data"}
                    row_i.update(entry)
                dist_dataframe = pd.concat([dist_dataframe, pd.DataFrame([row_i])], ignore_index=True)
            except:
                print("****Ope, there's a problem with your atom inputs.")
    print("Distance function has completed for", dist_list)
    return pd.concat([dataframe, dist_dataframe], axis = 1)

def get_enthalpies(dataframe): # gets thermochemical data from freq jobs
    enthalpy_dataframe = pd.DataFrame(columns=[]) #define an empty df to place results in

    for index, row in dataframe.iterrows(): #iterate over the dataframe
        try: #try to get the data
            log_file = row['log_name'] #read file name from df
            filecont = get_filecont(log_file) #read the contents of the log file

            evals = []
            error = "no thermochemical data found;;"
            e_hf, ezpe, h, g = 0, 0, 0, 0
            for i in range(len(filecont) - 1): #uses the zero_pattern that denotes this section to gather relevant energy terms
                if zero_pattern.search(filecont[i]):
                    e_hf = round(-eval(str.split(filecont[i-4])[-2]) + ezpe, 6)
                    evals.append(e_hf)
                    ezpe = eval(str.split(filecont[i])[-1])
                    evals.append(ezpe)
                    h = eval(str.split(filecont[i+2])[-1])
                    evals.append(h)
                    g = eval(str.split(filecont[i+3])[-1])
                    evals.append(g)
                    error = ""

            # This adds the data from the energy_values list (evals) into the new property df
            row_i = {'ZP_correction(Hartree)': evals[0], 'E_ZPE(Hartree)': evals[1], 'H(Hartree)': evals[2], 'G(Hartree)': evals[3]}
            #print(row_i)

            enthalpy_dataframe = pd.concat([enthalpy_dataframe, pd.DataFrame([row_i])], ignore_index=True)
        except:
            print(f'Unable to acquire enthalpies for: {row["log_name"]}.log')
    print("Enthalpies function has completed")
    return pd.concat([dataframe, enthalpy_dataframe], axis=1)

def get_time(dataframe): # gets wall time and CPU for all jobs
    time_dataframe = pd.DataFrame(columns=[]) #define an empty df to place results in

    for index, row in dataframe.iterrows(): #iterate over the dataframe
        try: #try to get the data
            log_file = row['log_name'] #read file name from df
            filecont, error = get_filecont(log_file) #read the contents of the log file
            if error != "":
                print(error)
                row_i = {'CPU_time_total(hours)': "no data", 'Wall_time_total(hours)': "no data"}
                time_dataframe = pd.concat([time_dataframe, pd.DataFrame([row_i])], ignore_index=True)
                continue

            cputime, walltime = 0, 0
            timeout = []
            for line in filecont:
                if cputime_pattern.search(line):
                    lsplt = str.split(line)
                    cputime = float(lsplt[-2])/3600 + float(lsplt[-4])/60 + float(lsplt[-6]) + float(lsplt[-8])*24
                    timeout.append(round(cputime, 5))
                if walltime_pattern.search(line):
                    lsplt = str.split(line)
                    walltime = float(lsplt[-2])/3600 + float(lsplt[-4])/60 + float(lsplt[-6]) + float(lsplt[-8])*24
                    timeout.append(walltime)
            CPU_time = sum(timeout[i] for i in range(len(timeout)) if i % 2 == 0)
            Wall_time = sum(timeout[i] for i in range(len(timeout)) if i % 2 != 0)

            #this adds the data from the CPU_time and Wall_time into the property df
            row_i = {'CPU_time_total(hours)': CPU_time, 'Wall_time_total(hours)': Wall_time}
            time_dataframe = pd.concat([time_dataframe, pd.DataFrame([row_i])], ignore_index=True)
        except:
            print(f'****Unable to acquire CPU time and wall time for: {row["log_name"]}.log')
            row_i = {'CPU_time_total(hours)': "no data", 'Wall_time_total(hours)': "no data"}
            time_dataframe = pd.concat([time_dataframe, pd.DataFrame([row_i])], ignore_index=True)
    print("Time function has completed")
    return pd.concat([dataframe, time_dataframe], axis=1)

def get_frontierorbs(dataframe): # HOMO, LUMO energies and derived values of last job in file
    frontierorbs_dataframe = pd.DataFrame(columns=[]) #define an empty df to place results in

    for index, row in dataframe.iterrows(): #iterate over the dataframe
        try: #try to get the data
            log_file = row['log_name'] #read file name from df
            filecont, error = get_filecont(log_file) #read the contents of the log file
            if error != "":
                print(error)
                row_i = {'HOMO': "no data", 'LUMO': "no data", "μ": "no data", "η": "no data", "ω": "no data"}
                frontierorbs_dataframe = pd.concat([frontierorbs_dataframe, pd.DataFrame([row_i])], ignore_index=True)
                continue

            frontierout = []
            index = 0
            for line in filecont[::-1]:
                if homo_pattern.search(line):
                    index += 1 #index ensures only the first entry is included
                    if index == 1:
                        homo = float(str.split(line)[-1])
                        lumo = float(str.split(filecont[filecont.index(line)+1])[4])
                        mu = (homo + lumo) / 2 # chemical potential or negative of molecular electronegativity
                        eta = lumo - homo # hardness/softness
                        omega = round(mu ** 2 / (2 * eta), 5) # electrophilicity index
                        frontierout.append(homo)
                        frontierout.append(lumo)
                        frontierout.append(mu)
                        frontierout.append(eta)
                        frontierout.append(omega)

            #this adds the data from the frontierout into the new property df
            row_i = {'HOMO': frontierout[0], 'LUMO': frontierout[1], "μ": frontierout[2], "η": frontierout[3], "ω": frontierout[4]}
            frontierorbs_dataframe = pd.concat([frontierorbs_dataframe, pd.DataFrame([row_i])], ignore_index=True)
        except:
            print(f'****Unable to acquire frontier orbitals for: {row["log_name"]}.log')
            row_i = {'HOMO': "no data", 'LUMO': "no data", "μ": "no data", "η": "no data", "ω": "no data"}
            frontierorbs_dataframe = pd.concat([frontierorbs_dataframe, pd.DataFrame([row_i])], ignore_index=True)
    print("Frontier orbitals function has completed")
    return pd.concat([dataframe, frontierorbs_dataframe], axis=1)

def get_volume(dataframe): #gets the molar volume of the molecule
    volume_dataframe = pd.DataFrame(columns=[]) #define an empty df to place results in

    for index, row in dataframe.iterrows(): #iterate over the dataframe
        try: #try to get the data
            log_file = row['log_name'] #read file name from df
            filecont, error = get_filecont(log_file) #read the contents of the log file
            if error != "":
                print(error)
                row_i = {'volume(Bohr_radius³/mol)': "no data"}
                volume_dataframe = pd.concat([volume_dataframe, pd.DataFrame([row_i])], ignore_index=True)
                continue

            volume = []
            for line in filecont:
                if volume_pattern.search(line):
                    volume.append(line.split()[3])
            #this adds the data into the new property df
            row_i = {'volume(Bohr_radius³/mol)': float(volume[0])}
            volume_dataframe = pd.concat([volume_dataframe, pd.DataFrame([row_i])], ignore_index=True)

        except:
            print(f'****Unable to acquire volume for: {row["log_name"]}.log')
            row_i = {'volume(Bohr_radius³/mol)': "no data"}
            volume_dataframe = pd.concat([volume_dataframe, pd.DataFrame([row_i])], ignore_index=True)
    print("Volume function has completed")
    return pd.concat([dataframe, volume_dataframe], axis=1)

def get_polarizability(dataframe): # polarizability isotropic and anisotropic
    polarizability_dataframe = pd.DataFrame(columns=[]) #define an empty df to place results in

    for index, row in dataframe.iterrows(): #iterate over the dataframe
        try: #try to get the data
            log_file = row['log_name'] #read file name from df
            filecont, error = get_filecont(log_file) #read the contents of the log file
            if error != "":
                print(error)
                row_i = {'polar_iso(Debye)': "no data", 'polar_aniso(Debye)': "no data"}
                polarizability_dataframe = pd.concat([polarizability_dataframe, pd.DataFrame([row_i])], ignore_index=True)
                continue

            polarout = []
            for i in range(len(filecont)-1, 1, -1):
                if polarizability_pattern.search(filecont[i]):
                    alpha_iso = float(filecont[i+4].split()[1].replace("D","E"))
                    alpha_aniso = float(filecont[i+5].split()[1].replace("D","E"))
                    polarout.append(alpha_iso)
                    polarout.append(alpha_aniso)

            #this adds the data from the polarout into the new property df
            row_i = {'polar_iso(Debye)': polarout[0], 'polar_aniso(Debye)': polarout[1]}
            polarizability_dataframe = pd.concat([polarizability_dataframe, pd.DataFrame([row_i])], ignore_index=True)

        except:
            print('****Unable to acquire polarizability for:' + row['log_name'] + ".log")
            row_i = {'polar_iso(Debye)': "no data", 'polar_aniso(Debye)': "no data"}
            polarizability_dataframe = pd.concat([polarizability_dataframe, pd.DataFrame([row_i])], ignore_index=True)
    print("Polarizability function has completed")
    return(pd.concat([dataframe, polarizability_dataframe], axis = 1))

def get_planeangle(dataframe, planeangle_list): # a function to get the plane angles for all atoms (planeangle_list, form [[O2, C1, O3, H5], [C4, C1, O3, H5]]) in a dataframe that contains file name and atom number
    planeangle_dataframe = pd.DataFrame(columns=[]) #define an empty df to place results in

    for index, row in dataframe.iterrows(): #iterate over the dataframe
        try:
            #parsing the plane angle list from input line
            planeanglenums_list = []
            for planeangle in planeangle_list:
                atomnum_list = [] # the atom numbers for a plane angle (i.e. 18 16 17 50) are collected from the df using the input list (i.e.["O2", "C1", "O3", "H5"])
                for atom in planeangle:
                    atomnum = row[str(atom)]
                    atomnum = re.findall(r'\d+', atomnum)[0] # extract only the number from the atom label
                    atomnum_list.append(str(atomnum))
                planeanglenums_list.append(atomnum_list) # append atomnum_list for each plane angle to make a list of the form [['18', '16', '17', '50'], ['18', '16', '17', '50']]
            
            planeangletitle_list = []
            for planeangle in planeangle_list:
                planeangletitle = str(planeangle[0]) + "_" + str(planeangle[1]) + "_" + str(planeangle[2]) + "_&_" + str(planeangle[3]) + "_" + str(planeangle[4]) + "_" + str(planeangle[5])
                planeangletitle_list.append(planeangletitle)
            
            log_file = row['log_name'] #read file name from df
            streams, error = get_outstreams(log_file)
            if error != "":
                print(error)
                row_i = {}
                for a in range(len(planeanglenums_list)):
                    entry = {'planeangle_' + str(planeangletitle_list[a]) + '(°)': "no data"}
                    row_i.update(entry)
                planeangle_dataframe = pd.concat([planeangle_dataframe, pd.DataFrame([row_i])], ignore_index=True)
                continue

            geom = get_geom(streams)

            #checks for if the wrong number of atoms are input, input is not of the correct form, or calls atom numbers that do not exist in the molecule.
            error = ""
            for planeangle in planeanglenums_list:
                if len(planeangle) % 6 != 0:
                    error = f"****Number of atom inputs given for plane angle is not divisible by six. {len(planeangle)} atoms were given."
                for atom in planeangle:
                    if not atom.isdigit():
                        error += f"**** {atom}: Only numbers accepted as input for plane angles"
                    if int(atom) > len(geom):
                        error += f"**** {atom} is out of range. Maximum valid atom number: {len(geom) + 1}."
                if error:
                    print(error)

            planeanglesout = []
            for planeangle in planeanglenums_list:
                a = geom[int(planeangle[0]) - 1][:4]
                b = geom[int(planeangle[1]) - 1][:4]
                c = geom[int(planeangle[2]) - 1][:4]
                d = geom[int(planeangle[3]) - 1][:4]
                e = geom[int(planeangle[4]) - 1][:4]
                f = geom[int(planeangle[5]) - 1][:4]

                ab = np.array([a[1] - b[1], a[2] - b[2], a[3] - b[3]]) # Vectors
                bc = np.array([b[1] - c[1], b[2] - c[2], b[3] - c[3]])
                de = np.array([d[1] - e[1], d[2] - e[2], d[3] - e[3]])
                ef = np.array([e[1] - f[1], e[2] - f[2], e[3] - f[3]])

                n1 = np.cross(ab, bc) # Normal vectors
                n2 = np.cross(de, ef)

                planeangle_value = round(np.degrees(np.arccos(np.dot(n1, n2) / (np.linalg.norm(n1) * np.linalg.norm(n2)))), 3)
                planeangle_value = min(abs(planeangle_value), abs(180 - planeangle_value))
                planeanglesout.append(planeangle_value)

            #this adds the data from the planeanglesout into the new property df
            row_i = {}
            for a in range(len(planeanglenums_list)):
                entry = {'planeangle_' + str(planeangletitle_list[a]) + '(°)': planeanglesout[a]}
                row_i.update(entry)
            planeangle_dataframe = pd.concat([planeangle_dataframe, pd.DataFrame([row_i])], ignore_index=True)
        except Exception as e:
            print(f'****Unable to acquire plane angle for: {row["log_name"]}.log with error: {e}')
            row_i = {}
            try:
                for a in range(len(planeanglenums_list)):
                    entry = {'planeangle_' + str(planeangletitle_list[a]) + '(°)': "no data"}
                    row_i.update(entry)
                planeangle_dataframe = pd.concat([planeangle_dataframe, pd.DataFrame([row_i])], ignore_index=True)
            except:
                print("****Ope, there's a problem with your atom inputs.")
    print(f"Plane angle function has completed for {planeangle_list}")
    return pd.concat([dataframe, planeangle_dataframe], axis=1)

def get_dipole(dataframe):
    dipole_dataframe = pd.DataFrame(columns=[]) #define an empty df to place results in

    for index, row in dataframe.iterrows(): #iterate over the dataframe
        try: #try to get the data
            log_file = row['log_name'] #read file name from df
            filecont, error = get_filecont(log_file) #read the contents of the log file
            if error != "":
                print(error)
                row_i = {'dipole(Debye)': "no data"}
                dipole_dataframe = pd.concat([dipole_dataframe, pd.DataFrame([row_i])], ignore_index=True)
                continue

            dipole = []
            for i in range(len(filecont) - 1, 0, -1): #search filecont in backwards direction
                if dipole_pattern in filecont[i]:
                    dipole.append(float(str.split(filecont[i + 1])[-1]))
            #this adds the data from the first dipole entry (corresponding to the last job in the file) into the new property df
            row_i = {'dipole(Debye)': dipole[0] if dipole else "no data"}
            dipole_dataframe = pd.concat([dipole_dataframe, pd.DataFrame([row_i])], ignore_index=True)
        except:
            print(f"****Unable to acquire dipole for: {row['log_name']}.log")
            row_i = {'dipole(Debye)': "no data"}
            dipole_dataframe = pd.concat([dipole_dataframe, pd.DataFrame([row_i])], ignore_index=True)
    print("Dipole function has completed")
    return pd.concat([dataframe, dipole_dataframe], axis=1)

def get_SASA(dataframe): #uses morfeus to calculate solvent accessible surface area in a dataframe that contains file name
    #if you want to calculate SASA with different probe radii, morfeus has this functionality, but it has not been implemented here
    sasa_dataframe = pd.DataFrame(columns=[])

    for index, row in dataframe.iterrows():
        try:
            log_file = row['log_name']
            streams, error = get_outstreams(log_file) #need to add file path if you're running from a different directory than the file
            if error != "":
                print(error)
                row_i = {'SASA_surface_area(Å²)': "no data",
                         'SASA_volume(Å³)': "no data",
                         'SASA_sphericity': "no data"}
                sasa_dataframe = pd.concat([sasa_dataframe, pd.DataFrame([row_i])], ignore_index=True)
                continue

            log_coordinates = get_geom(streams)
            elements = np.array([log_coordinates[i][0] for i in range(len(log_coordinates))])
            coordinates = np.array([np.array(log_coordinates[i][1:]) for i in range(len(log_coordinates))])

            sasa = SASA(elements, coordinates) #calls morfeus

            sphericity = np.cbrt((36 * math.pi * sasa.volume ** 2)) / sasa.area

            row_i = {'SASA_surface_area(Å²)': sasa.area,
                     'SASA_volume(Å³)': sasa.volume, #volume inside the solvent accessible surface area
                     'SASA_sphericity': sphericity}
            sasa_dataframe = pd.concat([sasa_dataframe, pd.DataFrame([row_i])], ignore_index=True)
        except:
            print(f'****Unable to acquire SASA parameters for: {row["log_name"]}.log')
            row_i = {'SASA_surface_area(Å²)': "no data",
                     'SASA_volume(Å³)': "no data",
                     'SASA_sphericity': "no data"}
            sasa_dataframe = pd.concat([sasa_dataframe, pd.DataFrame([row_i])], ignore_index=True)
    print("SASA function has completed")
    return pd.concat([dataframe, sasa_dataframe], axis=1)

def parse_goodvibes_output(output, temp=298.15):
    # This function extracts the desired values from the goodvibes output
    lines = output.splitlines()
    data = {}
    column_mapping = {
        "E_SPC": "E_spc (Hartree)",
        "ZPE": "ZPE(Hartree)",
        "H_SPC": "H_spc(Hartree)",
        "T.S": "T*S",
        "T.qh-S": "T*qh_S",
        "G(T)_SPC": "G(T)_spc(Hartree)",
        "qh-G(T)_SPC": "qh_G(T)_spc(Hartree)"
    }
    
    # Find the index positions of the two lines of asterisks
    start_index = None
    end_index = None
    header_line = None
    for i, line in enumerate(lines):
        if re.match(r"^\s*\*{12,}\s*$", line):  # Matches lines with 12 or more asterisks
            if start_index is None:
                start_index = i
                header_line = lines[i - 1]  # The header line is the one before the first line of asterisks
            else:
                end_index = i
                break  # We only need the first two lines of asterisks

    # Parse the header line to determine the order of properties
    if header_line:
        headers = re.split(r'\s+', header_line.strip())[1:]  # get rid of the first column, which is the structure column
    
    # Extract the relevant line between the two asterisk lines
    if start_index is not None and end_index is not None and end_index > start_index:
        for line in lines[start_index+1:end_index]:
            if re.match(r"^\s*o", line):  # Matches lines starting with 'o' (with any amount of whitespace before)
                parts = re.split(r'\s+', line.strip())  # Split the line by whitespace
                for i, header in enumerate(headers):
                    if header in column_mapping:
                        data[column_mapping[header]] = float(parts[i + 2])  # Offset by 2 because parts[0] is 'o', parts[1] is the structure name
                data['T'] = temp  # Assuming the temperature is consistent
                break

    if not data:
        # throw an error if the data is not found
        raise ValueError("No data found in the GoodVibes output.")
    
    return data

def get_goodvibes_e(dataframe, temp):
    e_dataframe = pd.DataFrame(columns=[])
    
    for index, row in dataframe.iterrows():
        try:
            log_file = row['log_name'] + ".log"
            
            # Construct command-line arguments for GoodVibes
            cmd_args = [
                sys.executable, "-m",
                "goodvibes", 
                log_file,
                "--spc", "link",
                "-t", str(temp)
            ]
            
            # Run the GoodVibes command and capture the output
            result = subprocess.run(cmd_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            # Parse the output
            parsed_data = parse_goodvibes_output(result.stdout, temp)
            # Append the parsed data to the DataFrame
            e_dataframe = pd.concat([e_dataframe, pd.DataFrame([parsed_data])], ignore_index=True)
        
        except Exception as e:
            print(f'****Unable to acquire goodvibes energies for: {row["log_name"]}.log with error: {e}')
            row_i = {
                'E_spc (Hartree)': "no data",
                'ZPE(Hartree)': "no data",
                'H_spc(Hartree)': "no data",
                'T*S': "no data",
                'T*qh_S': "no data",
                'G(T)_spc(Hartree)': "no data",
                'qh_G(T)_spc(Hartree)': "no data",
                'T': "no data"
            }
            e_dataframe = pd.concat([e_dataframe, pd.DataFrame([row_i])], ignore_index=True)

    print("Goodvibes function has completed")
    return pd.concat([dataframe, e_dataframe], axis=1)

def get_IR(dataframe, a1, a2, freqmin, freqmax, intmin, intmax, threshold): # a function to get IR values for a pair of atoms at a certain freq and intensity
    IR_dataframe = pd.DataFrame(columns=[]) #define an empty df to place results in
    pair_label = str(a1) + "_" + str(a2)

    for index, row in dataframe.iterrows(): #iterate over the dataframe
        try:
            log_file = row['log_name'] #read file name from df
            filecont, error = get_filecont(log_file)
            if error != "":
                print(error)
                row_i = {'IR_freq_' + str(pair_label): "no data"}
                IR_dataframe = pd.concat([IR_dataframe, pd.DataFrame([row_i])], ignore_index=True)
                continue

            #this changes a1 and a2 (of the form "C1" and "O3") to atomnum_pair (of the form [17, 18])
            atom1 = row[str(a1)]
            atom1 = re.findall(r'\d+', atom1)[0]
            atom2 = row[str(a2)]
            atom2 = re.findall(r'\d+', atom2)[0]

            #this section finds where all IR frequencies are located in the log file
            frq_len = 0
            frq_end = 0
            for i in range(len(filecont)):
                if frqs_pattern.search(filecont[i]) and frq_len == 1: #subsequent times it finds the pattern, it recognizes the frq_len
                    frq_len = i - 3 - frq_start
                if frqs_pattern.search(filecont[i]) and frq_len == 0: #first time it finds the pattern it will set frq_start
                    frq_start = i - 3
                    frq_len = 1
                if frqsend_pattern.search(filecont[i]): #finds the end pattern
                    frq_end = i - 3

            nfrq = filecont[frq_end - frq_len + 1].split()[-1]
            blocks = int((frq_end + 1 - frq_start) / frq_len)
            irdata = []   # list of objects. IR contains: IR.freq, IR.int, IR.deltas = []

            for i in range(0, blocks):
                for j in range(len(filecont[i * frq_len + frq_start].split())):
                    irdata.append(IR(filecont, i * frq_len + frq_start, j, frq_len))

            irout = []
            for i in range(len(irdata)):
                if (irdata[i].freq < freqmax and irdata[i].freq > freqmin and 
                    irdata[i].int > intmin and irdata[i].int < intmax and 
                    irdata[i].deltas[int(atom1)] >= threshold and 
                    irdata[i].deltas[int(atom2)] >= threshold):
                        irout = [irdata[i].freq, irdata[i].int]

            #this adds the frequency data from the irout into the new property df
            row_i = {'IR_freq_' + str(pair_label): irout[0] if irout else "no data"}
            IR_dataframe = pd.concat([IR_dataframe, pd.DataFrame([row_i])], ignore_index=True)
        except:
            print(f'****Unable to acquire IR frequencies for: {row["log_name"]}.log')
            row_i = {'IR_freq_' + str(pair_label): "no data"}
            IR_dataframe = pd.concat([IR_dataframe, pd.DataFrame([row_i])], ignore_index=True)
    print(f"IR function has completed for {a1} and {a2}")
    return pd.concat([dataframe, IR_dataframe], axis=1)

def get_buried_sterimol(dataframe, sterimol_list, r_buried): #uses morfeus to calculate sterimol L, B1, B5 for two input atoms for every entry in df
    sterimol_dataframe = pd.DataFrame(columns=[])
    r_buried -= 0.5 #the function adds

    for index, row in dataframe.iterrows():
        try:
            # Parsing the Sterimol axis defined in the list from input line
            sterimolnums_list = []
            for sterimol in sterimol_list:
                atomnum_list = [] # The atom numbers used to collect sterimol values (i.e. [18 16 17 15]) are collected from the df using the input list (i.e. [["O2", "C1"], ["O3", "H5"]])
                for atom in sterimol:
                    atomnum = row[str(atom)]
                    atomnum = re.findall(r'\d+', atomnum)[0] # Extracts the atom number from the string
                    atomnum_list.append(str(atomnum))
                sterimolnums_list.append(atomnum_list) # Append atomnum_list for each sterimol axis defined in the input to make a list of the form [['18', '16'], ['16', '15']]

            # This makes column headers based on Sterimol axis defined in the input line
            sterimoltitle_list = []
            for sterimol in sterimol_list:
                sterimoltitle = str(sterimol[0]) + "_" + str(sterimol[1])
                sterimoltitle_list.append(sterimoltitle)

            log_file = row['log_name']
            streams, error = get_outstreams(log_file) # Need to add file path if you're running from a different directory than file
            if error != "":
                print(error)
                row_i = {}
                for a in range(len(sterimolnums_list)):
                    entry = {'Buried_Sterimol_L_' + str(sterimoltitle_list[a]) + '_' + str(r_buried) + '(Å)': "no data",
                             'Buried_Sterimol_B1_' + str(sterimoltitle_list[a]) + '_' + str(r_buried) + '(Å)': "no data",
                             'Buried_Sterimol_B5_' + str(sterimoltitle_list[a]) + '_' + str(r_buried) + '(Å)': "no data"}
                    row_i.update(entry)
                sterimol_dataframe = pd.concat([sterimol_dataframe, pd.DataFrame([row_i])], ignore_index=True)
                continue

            geom = get_geom(streams)

            # Checks for if the wrong number of atoms are input, input is not of the correct form, or calls atom numbers that do not exist in the molecule
            error = ""
            for sterimol in sterimolnums_list:
                if len(sterimol) % 2 != 0:
                    error = f"Number of atom inputs given for Sterimol is not divisible by two. {len(sterimol)} atoms were given."
                for atom in sterimol:
                    if not atom.isdigit():
                        error += f" {atom}: Only numbers accepted as input for Sterimol."
                    if int(atom) > len(geom):
                        error += f" {atom} is out of range. Maximum valid atom number: {len(geom)}."
                if error:
                    print(error)

            elements = np.array([geom[i][0] for i in range(len(geom))])
            coordinates = np.array([np.array(geom[i][1:]) for i in range(len(geom))])

            # This collects Sterimol values for each pair of inputs
            sterimolout = []
            for sterimol in sterimolnums_list:
                sterimol_values = Sterimol(elements, coordinates, int(sterimol[0]), int(sterimol[1])) # Calls morfeus
                sterimol_values.bury(method="delete", sphere_radius=float(r_buried))
                sterimolout.append(sterimol_values)

            # This adds the data from sterimolout into the new property df
            row_i = {}
            for a in range(len(sterimolnums_list)):
                entry = {'Buried_Sterimol_L_' + str(sterimoltitle_list[a]) + '_' + str(r_buried) + '(Å)': sterimolout[a].L_value,
                         'Buried_Sterimol_B1_' + str(sterimoltitle_list[a]) + '_' + str(r_buried) + '(Å)': sterimolout[a].B_1_value,
                         'Buried_Sterimol_B5_' + str(sterimoltitle_list[a]) + '_' + str(r_buried) + '(Å)': sterimolout[a].B_5_value}
                row_i.update(entry)
            sterimol_dataframe = pd.concat([sterimol_dataframe, pd.DataFrame([row_i])], ignore_index=True)
        except:
            print(f'****Unable to acquire Morfeus Buried Sterimol parameters for: {row["log_name"]}.log')
            row_i = {}
            try:
                for a in range(len(sterimolnums_list)):
                    entry = {'Buried_Sterimol_L_' + str(sterimoltitle_list[a]) + '_' + str(r_buried) + '(Å)': "no data",
                             'Buried_Sterimol_B1_' + str(sterimoltitle_list[a]) + '_' + str(r_buried) + '(Å)': "no data",
                             'Buried_Sterimol_B5_' + str(sterimoltitle_list[a]) + '_' + str(r_buried) + '(Å)': "no data"}
                    row_i.update(entry)
                sterimol_dataframe = pd.concat([sterimol_dataframe, pd.DataFrame([row_i])], ignore_index=True)
            except:
                print("****Ope, there's a problem with your atom inputs.")
    print(f"Morfeus Buried Sterimol function has completed for {sterimol_list}")
    return pd.concat([dataframe, sterimol_dataframe], axis=1)

def get_chelpg(dataframe, a_list): #a function to get the ChelpG ESP charges for all atoms (a_list, form ["C1", "C4", "O2"]) in a dataframe that contains file name and atom number
    chelpg_dataframe = pd.DataFrame(columns=[]) #define an empty df to place results in

    for index, row in dataframe.iterrows(): #iterate over the dataframe
        try:
            atomnum_list = []
            for atom in a_list:
                atomnum = row[str(atom)] #the atom number (i.e. 16) to add to the list is the df entry of this row for the labeled atom (i.e. "C1")
                atomnum = re.findall(r'\d+', atomnum)[0] #extract only the number from the atom label
                atomnum_list.append(str(atomnum)) #append that to atomnum_list to make a list of the form [16, 17, 29]
            log_file = row['log_name'] #read file name from df
            filecont, error = get_filecont(log_file) #read the contents of the log file
            if error != "":
                print(error)
                row_i = {}
                for a in range(len(a_list)):
                    entry = {'ChelpG_charge_' + str(a_list[a]): "no data"}
                    row_i.update(entry)
                chelpg_dataframe = pd.concat([chelpg_dataframe, pd.DataFrame([row_i])], ignore_index=True)
                continue

            chelpgstart, chelpg, error, chelpgout = 0, False, "", []

            #this section finds the line (chelpgstart) where the ChelpG data is located
            for i in range(len(filecont) - 1, 0, -1):
                if chelpg2_pattern.search(filecont[i]):
                    chelpgstart = i
                if chelpg1_pattern.search(filecont[i]):
                    chelpg = True
                    break
            if chelpgstart != 0 and not chelpg:
                error = f"****Other ESP scheme than ChelpG used in: {log_file}.log"
            if chelpgstart == 0:
                error = f"****no ChelpG ESP charge analysis found in: {log_file}.log"
            if error:
                print(error)
                row_i = {}
                for a in range(len(a_list)):
                    entry = {'ChelpG_charge_' + str(a_list[a]): "no data"}
                    row_i.update(entry)
                chelpg_dataframe = pd.concat([chelpg_dataframe, pd.DataFrame([row_i])], ignore_index=True)
                continue

            for atom in atomnum_list:
                if atom.isnumeric():
                    chelpgout.append(filecont[chelpgstart + int(atom) + 2].split()[-1])

            #this adds the data from the chelpgout into the new property df
            row_i = {}
            for a in range(len(a_list)):
                entry = {'ChelpG_charge_' + str(a_list[a]): chelpgout[a]}
                row_i.update(entry)
            chelpg_dataframe = pd.concat([chelpg_dataframe, pd.DataFrame([row_i])], ignore_index=True)
        except:
            print(f'****Unable to acquire ChelpG charges for: {row["log_name"]}.log')
            row_i = {}
            for a in range(len(a_list)):
                entry = {'ChelpG_charge_' + str(a_list[a]): "no data"}
                row_i.update(entry)
            chelpg_dataframe = pd.concat([chelpg_dataframe, pd.DataFrame([row_i])], ignore_index=True)
    print(f"ChelpG function has completed for {a_list}")
    return pd.concat([dataframe, chelpg_dataframe], axis=1)

def get_hirshfeld(dataframe,a_list): #a function to get the Hirshfeld charge, CM5 charge, and atomic dipole for all atoms (a_list, form ["C1", "C4", "O2"]) in a dataframe that contains file name and atom number
    hirsh_dataframe = pd.DataFrame(columns=[]) #define an empty df to place results in

    for index, row in dataframe.iterrows(): #iterate over the dataframe
        try:#try to get the data
            atomnum_list = []
            for atom in a_list:
                atomnum = row[str(atom)] #the atom number (i.e. 16) to add to the list is the df entry of this row for the labeled atom (i.e. "C1")
                atomnum = re.findall(r'\d+', atomnum)[0] #extract only the number from the atom label
                atomnum_list.append(str(atomnum)) #append that to atomnum_list to make a list of the form [16, 17, 29]

            log_file = row['log_name'] #read file name from df
            filecont, error = get_filecont(log_file) #read the contents of the log file
            if error != "":
                print(error)
                row_i = {}
                for a in range(len(a_list)):
                    entry = {'Hirsh_charge_'+str(a_list[a]): "no data",
                             'Hirsh_CM5_charge_'+str(a_list[a]): "no data",
                             'Hirsh_atom_dipole_'+str(a_list[a]): "no data"}
                    row_i.update(entry)
                hirsh_dataframe = pd.concat([hirsh_dataframe, pd.DataFrame([row_i])], ignore_index=True)
                continue

            hirshstart, error, hirshout = 0, False, []

            #this section finds the line (chelpgstart) where the ChelpG data is located
            for i in range(len(filecont) - 1, 0, -1):
                if hirshfeld_pattern.search(filecont[i]):
                    hirshstart = i
                    break
            if hirshstart == 0:
                error = f"****no Hirshfeld Population Analysis found in: {log_file}.log"
                print(error)
                row_i = {}
                for a in range(len(a_list)):
                    entry = {'Hirsh_charge_'+str(a_list[a]): "no data",
                             'Hirsh_CM5_charge_'+str(a_list[a]): "no data",
                             'Hirsh_atom_dipole_'+str(a_list[a]): "no data"}
                    row_i.update(entry)
                hirsh_dataframe = pd.concat([hirsh_dataframe, pd.DataFrame([row_i])], ignore_index=True)
                continue

            for atom in atomnum_list:
                if atom.isnumeric():
                    cont = filecont[hirshstart + int(atom) + 1].split()
                    qh = cont[2]  # using 0-indexing, this gets the value for Hirshfeld charge from the 2nd column
                    qcm5 = cont[7]  # using 0-indexing, this gets the value for CM5 charge from the 7th column
                    d = np.linalg.norm(np.array((cont[4:8])))
                    hirshout.append([str(qh), str(qcm5), str(d)])

            #this adds the data from the hirshout into the new property df
            row_i = {}
            for a in range(len(a_list)):
                entry = {'Hirsh_charge_'+str(a_list[a]): hirshout[a][0],
                         'Hirsh_CM5_charge_'+str(a_list[a]): hirshout[a][1],
                         'Hirsh_atom_dipole_'+str(a_list[a]): hirshout[a][2]}
                row_i.update(entry)
            hirsh_dataframe = pd.concat([hirsh_dataframe, pd.DataFrame([row_i])], ignore_index=True)
        except:
            print(f'****Unable to acquire Hirshfeld properties for: {row["log_name"]}.log')
            row_i = {}
            for a in range(len(a_list)):
                entry = {'Hirsh_charge_'+str(a_list[a]): "no data",
                         'Hirsh_CM5_charge_'+str(a_list[a]): "no data",
                         'Hirsh_atom_dipole_'+str(a_list[a]): "no data"}
                row_i.update(entry)
            hirsh_dataframe = pd.concat([hirsh_dataframe, pd.DataFrame([row_i])], ignore_index=True)
    print(f"Hirshfeld function has completed for {a_list}")
    return pd.concat([dataframe, hirsh_dataframe], axis=1)

def get_cone_angle(dataframe, a_list): #DOES NOT MATCH VALUES FROM LITERATURE, WORK IN PROGRESS
    cone_angle_dataframe = pd.DataFrame(columns=[])

    for index, row in dataframe.iterrows():
        if True:
        #try:
            atom_list = []
            for label in a_list:
                atom = row[str(label)] #the atom number (i.e. 16) to add to the list is the df entry of this row for the labeled atom (i.e. "C1")
                atom = re.findall(r'\d+', atom)[0] #extract only the number from the atom label
                atom_list.append(str(atom)) #append that to atom_list to make a list of the form [16, 17, 29]

            log_file = row['log_name']
            streams, errors = get_outstreams(log_file) #need to add file path if you're running from a different directory than file
            log_coordinates = get_geom(streams)
            elements = np.array([log_coordinates[i][0] for i in range(len(log_coordinates))])
            coordinates = np.array([np.array(log_coordinates[i][1:]) for i in range(len(log_coordinates))])

            cone_angle_out = []
            for atom in atom_list:
                cone_angle = ConeAngle(elements, coordinates, int(atom)) #calls morfeus
                cone_angle_out.append(cone_angle)
            cone_angle.print_report()

            row_i = {}
            for a in range(0, len(atom_list)):
                entry = {'cone_angle' + str(a_list[a]) + '(°)': cone_angle_out[a].cone_angle} #details on these values can be found here: https://kjelljorner.github.io/morfeus/pyramidalization.html
                row_i.update(entry)
            cone_angle_dataframe = pd.concat([cone_angle_dataframe, pd.DataFrame([row_i])], ignore_index=True)
        #except:
        #    print('Unable to acquire cone_angle parameters for:', row['log_name'], ".log")
    print("cone_angle function has completed for", a_list)
    return(pd.concat([dataframe, cone_angle_dataframe], axis = 1))

def get_one_lp_energy(dataframe, a_list): #a function to get the NB orbitals for all atoms (a_list, form ["C1", "C4", "O2"]) in a dataframe that contains file name and atom number
    nborbs_dataframe = pd.DataFrame(columns=[]) #define an empty df to place results in

    for index, row in dataframe.iterrows(): #iterate over the dataframe 
        try: #try to get the data
            atomnum_list = [] 
            for atom in a_list: 
                atomnum = row[str(atom)] #the atom number (i.e. 16) to add to the list is the df entry of this row for the labeled atom (i.e. "C1")
                atomnum = re.findall(r'\d+', atomnum)[0] #extract only the number from the atom label
                atomnum_list.append(str(atomnum)) #append that to atomnum_list to make a list of the form [16, 17, 29]

            log_file = row['log_name'] #read file name from df
            filecont, error = get_filecont(log_file) #read the contents of the log file
            if error != "":
                print(error)
                row_i = {}
                for a in range(0, len(a_list)):
                    entry = {'NBO_charge_'+str(a_list[a]): "no data"}
                    row_i.update(entry)
                nborbs_dataframe = pd.concat([nborbs_dataframe, pd.DataFrame([row_i])], ignore_index=True)
                continue

            nborbsstart = 0
            # this section finds the line (nborbsstart) where the nbo data is located
            for i in range(len(filecont)-1,0,-1):
                if nborbs_pattern in filecont[i]:#search the file content for the phrase which indicates the start of the NB orbitals section 
                    nborbsstart = i   
            if nborbsstart == 0: 
                error = "****no Natural Bond Orbitals found in: " + str(row['log_name']) + ".log"
                print(error)
                row_i = {}
                for a in range(0, len(a_list)):
                    entry = {'NBO_charge_'+str(a_list[a]): "no data"}
                    row_i.update(entry)
                nborbs_dataframe = pd.concat([nborbs_dataframe, pd.DataFrame([row_i])], ignore_index=True)
                continue

            for atom in a_list: 
                k = 0
                atom_num = row[str(atom)]
                atom_num = re.findall(r'\d+', atom_num)[0]
                for j in range(nborbsstart,len(filecont)):
                    if str(atom_num) in " ".join(re.findall("([A-Z][a-z]? *[0-9]+)",filecont[j])).split() and ("LP" in filecont[j]):
                        orbital_section = re.search("[0-9]+\.[A-Z\*(0-9 ]+\)",filecont[j]).group(0) #type of MO
                        orbital = orbital_section.split(". ")
                        orb = orbital[1]
                        des = orb.split(" ")
                        orb_type = des[0]
                        occ_energy = [x for x in re.findall(r"[-+]?\d*\.\d+",filecont[j])]
                        occ = occ_energy[0]
                        energy = occ_energy[1]
                        k += 1
                        # print(k)
                if k == 0: 
                    error = "****no LPs for atom " + str(atom)+ " in: " + str(row['log_name']) + ".log"
                    print(error)
                    row_i = {}
                    for atom in a_list:
                        entry = {'NBO_LP_occupancy_' + str(atom): "no data", 'NBO_LP_energy_' + str(atom): "no data"}
                        row_i.update(entry)
                    nborbs_dataframe = pd.concat([nborbs_dataframe, pd.DataFrame([row_i])], ignore_index=True)
                    pass
                if k == 2: 
                    error = "****more than one LP for atom " + str(atom)+ " in: " + str(row['log_name']) + ".log"
                    print(error)
                    row_i = {}
                    for atom in a_list:
                        entry = {'NBO_LP_occupancy_' + str(atom): "no data", 'NBO_LP_energy_' + str(atom): "no data"}
                        row_i.update(entry)
                    nborbs_dataframe = pd.concat([nborbs_dataframe, pd.DataFrame([row_i])], ignore_index=True)
                    continue

            # this adds the data from the nboout into the new property df
            row_i = {}
            for atom in a_list:
                entry = {'NBO_LP_occupancy_' + str(atom): occ, 'NBO_LP_energy_' + str(atom): energy}
                row_i.update(entry)
            print(row_i)
            nborbs_dataframe = pd.concat([nborbs_dataframe, pd.DataFrame([row_i])], ignore_index=True)
        except Exception as e:
            print(f'****Unable to acquire NBO orbitals for: {row["log_name"]}.log with error: {e}')
            row_i = {}
            for a in range(0, len(a_list)):
                entry = {'NBO_charge_'+str(a_list[a]): "no data"}
                row_i.update(entry)
            nborbs_dataframe = pd.concat([nborbs_dataframe, pd.DataFrame([row_i])], ignore_index=True)
    print("NBOrbs function has completed for", a_list)
    return(pd.concat([dataframe, nborbs_dataframe], axis = 1))


def get_natural_atomic_valencies(dataframe, a_list):
    """
    Extracts Natural Atomic Valencies for specified atoms from a Gaussian log file.

    Parameters:
        dataframe (pd.DataFrame): DataFrame containing 'log_name' and atom mappings.
        a_list (list): List of atoms to extract valencies for (e.g., ["C1", "C2"]).

    Returns:
        pd.DataFrame: Updated DataFrame with extracted valencies.
    """
    valency_dataframe = pd.DataFrame(columns=[])

    for index, row in dataframe.iterrows():
        try:
            # mappings between a_list and actual labels in the row
            row_to_a_list = {row[atom]: atom for atom in a_list}
            # print(f"Mapped Atoms: {row_to_a_list}")
            # print look something like Mapped Atoms: {'C1': 'C1', 'C2': 'C2'}
            # the key is from the each row cell and the value is the input a_list(reassigned atom label)
            # we will use atom number from the ROW CELL to find the corresponding atom in the log file

            log_file = row["log_name"]  # Load log file content
            filecont, error = get_filecont(log_file)
            if error:
                print(error)
                row_i = {f"Natural_Valency_{atom}": "no data" for atom in a_list}
                valency_dataframe = pd.concat(
                    [valency_dataframe, pd.DataFrame([row_i])], ignore_index=True
                )
                continue

            # Locate the section with Natural Atomic Valencies using regex
            start_line = None
            for i, line in enumerate(filecont):
                if nav_pattern.search(line):
                    start_line = i + 4  # Data starts four lines below the header
                    break
            if start_line is None:
                print(
                    f"****No Natural Atomic Valencies section found in: {row['log_name']}.log"
                )
                row_i = {f"Natural_Valency_{atom}": "no data" for atom in a_list}
                valency_dataframe = pd.concat(
                    [valency_dataframe, pd.DataFrame([row_i])], ignore_index=True
                )
                continue

            # Split and combine headers from both lines
            headers_first = re.split(
                r"\s{1,}", filecont[start_line - 1].strip()
            )  # Split by one or more spaces
            headers_second = filecont[start_line - 2].split()  # Split by single spaces
            # Append headers_second backward to headers_first traversing backward, starting from the last element
            for index, header in enumerate(headers_second[::-1]):
                headers_first[-(index + 1)] = header + headers_first[-(index + 1)]

            # Parse the valency data
            valency_data = {}
            for line in filecont[start_line + 1 :]:
                if (
                    line.strip() == "" or "$" in line
                ):  # Skip empty lines and lines with "$"
                    break  # Stop parsing when the section ends
                parts = re.split(r"\s+", line.strip())
                atom_label = f"{parts[1]}{parts[0].replace('.', '')}"
                if atom_label in row_to_a_list:
                    valency_data[row_to_a_list[atom_label]] = {}
                    for i, valency in enumerate(parts[2:]):
                        valency_data[row_to_a_list[atom_label]][
                            headers_first[i + 1]
                        ] = valency

            exclude_properties = ["Valency", "ElectronCount"]  # a list of properties we don't want to add to the dataframe
            # Add extracted values to the DataFrame
            row_i = {}
            for atom in a_list:
                atom_data = valency_data.get(atom, None)
                if not atom_data:  # Handle missing data
                    row_i[f"Natural_Valency_{atom}"] = "no data"
                    continue
                for property_name, value in atom_data.items():
                    if property_name not in exclude_properties:
                        row_i[f"{atom}_{property_name}"] = value
            valency_dataframe = pd.concat(
                [valency_dataframe, pd.DataFrame([row_i])], ignore_index=True
            )

        except Exception as e:
            print(
                f"****Unable to acquire Natural Atomic Valencies for: {row['log_name']}.log with error: {e}"
            )
            row_i = {f"Natural_Valency_{atom}": "no data" for atom in a_list}
            valency_dataframe = pd.concat(
                [valency_dataframe, pd.DataFrame([row_i])], ignore_index=True
            )

    print(f"Natural Atomic Valencies function has completed for {a_list}")
    return pd.concat([dataframe, valency_dataframe], axis=1)


def get_natural_bond_order(dataframe, bond_list):
    """
    Extracts Natural Bond Order (total/covalent/ionic) for specified bonds from a Gaussian log file.

    Parameters:
        dataframe (pd.DataFrame): DataFrame containing 'log_name' and atom mappings.
        bond_list (list of list): List of bonds to extract (e.g., [["C1", "C2"]]).

    Returns:
        pd.DataFrame: Updated DataFrame with extracted bond orders.
    """
    bond_order_dataframe = pd.DataFrame(columns=[])

    for _, row in dataframe.iterrows():
        try:
            # Map atom labels from the row to the bond list
            bond_mappings = [
                [row[bond[0]], row[bond[1]]]
                for bond in bond_list
                if bond[0] in row and bond[1] in row
            ]
            # if the list is empty, raise an error
            if not bond_mappings:
                print(f"****Invalid bond lists, cannot find atoms in the corresponding row.")
                raise ValueError("No valid atom mappings found in the row.")

            filecont, error = get_filecont(row["log_name"])  # Load log file content
            if error:
                print(error)
                row_data = {f"{bond[0]}_{bond[1]}_Bond_Order": "no data" for bond in bond_list}
                bond_order_dataframe = pd.concat([bond_order_dataframe, pd.DataFrame([row_data])], ignore_index=True)
                continue

            # Locate the section with Natural Bond Order using regex
            start_line = None
            for i, line in enumerate(filecont):
                if nborder_pattern.search(line):
                    start_line = i + 2  # Data starts three lines below the header
                    break
            if start_line is None:
                print(f"****No Natural Bond Order section found in: {row['log_name']}.log")
                row_data = {f"{bond[0]}_{bond[1]}_Bond_Order": "no data" for bond in bond_list}
                bond_order_dataframe = pd.concat([bond_order_dataframe, pd.DataFrame([row_data])], ignore_index=True)
                continue

            # Parse the section until two consecutive newlines
            section_data = []
            blank_line_count = 0
            for line in filecont[start_line:]:
                if line.strip() == "":
                    blank_line_count += 1
                    if blank_line_count == 2:
                        break
                else:
                    blank_line_count = 0
                section_data.append(line)

            # Split section data by atom title lines
            sections = []
            current_section = []
            for line in section_data:
                if line.strip().startswith("Atom"):
                    if current_section:
                        sections.append(current_section)
                    current_section = [line]
                else:
                    current_section.append(line)
            if current_section:
                sections.append(current_section)

            # Parse each section
            bond_orders = {}
            for section in sections:
                title_line = section[0].strip()
                # print(f"****Title Line: {title_line}")
                atom_indices = re.split(r"\s{2,}", title_line)[
                    1:
                ]  # Extract atom indices
                # print(f"****Atom Indices: {atom_indices}")

                for i in range(2, len(section), 4):  # Process 4 lines at a time
                    # print(f"****Processing lines {i} to {i + 2}")
                    if i + 2 >= len(section):
                        break
                    total_line = section[i].strip()
                    covalent_line = section[i + 1].strip()
                    ionic_line = section[i + 2].strip()

                    atom_parts = re.split(r"\s+", total_line)
                    atom_label = (
                        atom_parts[0].strip().replace(".", "")
                    )  # we will only use the number part of the atom label
                    # print(f"****Atom Label: {atom_label}")

                    total_values = atom_parts[3:]  # Skip atom label parts
                    covalent_values = re.split(r"\s+", covalent_line.strip())[1:]
                    ionic_values = re.split(r"\s+", ionic_line.strip())[1:]

                    if atom_label in bond_orders:
                        # print(f"****Updating existing atom label: {atom_label}")
                        bond_orders[atom_label].update(
                            {
                                atom_indices[j]: {
                                    "t": total_values[j],
                                    "c": covalent_values[j],
                                    "i": ionic_values[j],
                                }
                                for j in range(len(atom_indices))
                            }
                        )
                    else:
                        bond_orders[atom_label] = {
                            atom_indices[j]: {
                                "t": total_values[j],
                                "c": covalent_values[j],
                                "i": ionic_values[j],
                            }
                            for j in range(len(atom_indices))
                        }

            # Extract bond orders for the specified bonds
            row_data = {}
            for index, bond in enumerate(bond_mappings):
                atom1, atom2 = bond
                # find the actual atom numbers from the row, eg, C1 -> 1
                atom1_index = str(re.findall(r"\d+", atom1)[0])
                atom2_index = str(re.findall(r"\d+", atom2)[0])
                # print(f"****Processing bond {atom1} - {atom2}, indices: {atom1_index}, {atom2_index}")
                bond_order_data = bond_orders.get(atom1_index, {}).get(atom2_index, None)
                # print(f"****Bond Order Data: {bond_order_data}")
                if bond_order_data:
                    row_data[f"{bond_list[index][0]}_{bond_list[index][1]}_bond_order_total"] = bond_order_data["t"]
                    row_data[f"{bond_list[index][0]}_{bond_list[index][1]}_bond_order_covalent"] = bond_order_data["c"]
                    row_data[f"{bond_list[index][0]}_{bond_list[index][1]}_bond_order_ionic"] = bond_order_data["i"]
                else:
                    row_data[f"{bond_list[index][0]}_{bond_list[index][1]}_bond_order_total"] = "no data"
                    row_data[f"{bond_list[index][0]}_{bond_list[index][1]}_bond_order_covalent"] = "no data"
                    row_data[f"{bond_list[index][0]}_{bond_list[index][1]}_bond_order_ionic"] = "no data"
            bond_order_dataframe = pd.concat(
                [bond_order_dataframe, pd.DataFrame([row_data])], ignore_index=True
            )

        except Exception as e:
            print(f"****Error processing log file {row['log_name']}: {e}")
            row_data = {
                f"Bond_Order_{bond[0]}_{bond[1]}": "no data" for bond in bond_list
            }
            bond_order_dataframe = pd.concat(
                [bond_order_dataframe, pd.DataFrame([row_data])], ignore_index=True
            )

    print(f"Natural Bond Order function has completed for {bond_list}")
    return pd.concat([dataframe, bond_order_dataframe], axis=1)
