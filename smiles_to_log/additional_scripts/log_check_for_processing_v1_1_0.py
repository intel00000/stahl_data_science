# Version 1.1.0:
# 	- .com files now move to the imag or term_error folders when the respective .log file does
# 	- identification of termination errors has been updated to more precisely identify the number of normal terminations in the file (!= rather than > the correct number of jobs run)

# Call with log file as argument. Default all in directory

import sys
from glob import glob
import shutil
import os

# creates a list of files with given extension, specified in argv. If no files specified, returns a list of all files
# with extension in the directory
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

    if files == []:  # if list remains empty, take all files in directory with given extension
        files = [file for file in glob("*.%s" % ext)]
        print("Examining " + str(len(files)) + " file(s). \n")
    assert files != [], "No files with %s extension in directory" % ext
    return files


def normaltermination(log): #checks if "Normal termination" appears in the file.
    log_stem = log.split(".")[0]
    normalt = 0
    with open(log, "r") as inp:
            lines = inp.readlines()
            for line in lines:
                if "Normal termination" in line:
                       normalt += 1
    if normalt != 4: #CHANGE to total number of "Normal termination" that should appear in the .log file 
    #note that you should have one "Normal termination" for each job you run, keep in mind and "Opt + Freq" job is actually two jobs
            print(log + " did not terminate normally. \n")
            source_path = str(os.getcwd()) + "/" + str(log)
            destination_path = str(os.getcwd()) + "/../term_error/" + str(log)
            shutil.move(source_path, destination_path)
            source_path_com = str(os.getcwd()) + "/../" + str(log_stem) + ".com"
            destination_path_com = str(os.getcwd()) + "/../term_error/" + str(log_stem) + ".com"
            shutil.move(source_path_com, destination_path_com)
            
def lowestfrequencies(log): #prints lowest energy frequency in each log file (for files with freq calculations)
    log_stem = log.split(".")[0]
    with open(log, "r") as inp:
            lines = inp.readlines()
            for i in range(0, len(lines)):
                line = lines[i]
                if " and normal coordinates:" in line:
                    lowfreq = lines[i + 3]
                    if float(lowfreq[18:27]) <= 0:
                        print(log + ": " + lowfreq[18:27])
                        source_path = str(os.getcwd()) + "/" + str(log)
                        destination_path = str(os.getcwd()) + "/../imag/" + str(log)
                        shutil.move(source_path, destination_path)
                        source_path_com = str(os.getcwd()) + "/../" + str(log_stem) + ".com"
                        destination_path_com = str(os.getcwd()) + "/../imag/" + str(log_stem) + ".com"
                        shutil.move(source_path_com, destination_path_com)
                        
#dirname = os.path.dirname(__file__)
log_files = input_ext(sys.argv, "log")  # Load specified log files. Default all log files in directory.

for log_file in log_files:
    normaltermination(log_file)

log_files = input_ext(sys.argv, "log")  # reload log files after some have moved out
print("Lowest energy frequency in files with imaginary frequency:")
for log_file in log_files:
    lowestfrequencies(log_file)