# automated_workflow_for_SMILES_to_logs
Version 1.1.0 Updates to post_processing.bash and log_check_for_processing.py
- updated to use Python/3.10 (from Python/3.7, which is no longer available on the CHPC)
- .com files now move to the imag or term_error folders when the respective .log file does
- identification of termination errors has been updated to more precisely identify the number of normal terminations in the file (!= rather than > the correct number of jobs run)

Brittany Haas and Melissa Hardy's workflow scripts for automated generation of sdfs from SMILES, MacroModel conformational searching and clustering, and submitting/checking Gaussian jobs

Automated SMILES to logs Workflow

General Notes: 

Make sure all scripts have the correct permissions (chmod 744)

All scripts can go in bin folder.

The SMILES_to_log_workflow.jpg also shows the order in which to use the scripts. 

1.	SMILES_to_sdf.py 


    a.	Contributors: Guilian, Ellie (kraken)
  
    b.	Uses Rdkit or OpenBabel to convert SMILES strings to .sdf files (using .xyz files as an intermediate step)
  
    c.	File Format Example: “smiles.xlsx”
  
        We strongly recommend your ID contains a prefix + a number or a letter (this will make post-processing easier)
        
    d.	***Cannot be used with radicals. The script currently checks for radicals and deletes any files that contain radicals. Could be modified, if needed. ***
  
    e.	Before running the script you must load OpenBabel (ml openbabel)
  
    f.	python ~/bin/SMILES_to_sdf.py
  
    g.	install any missing packages like RDKit and openpyxl with:
  
        #      run "pip install rdkit-pypi"
        
        #      run "pip install openpyxl"
        
        # load Python 3.7:
        
        #      run "ml python/3.7"
        
        # load OpenBabel:
        
        #      run "ml openbabel"
        
        # run this script by navigating into the directory with the Excel file "smiles.xlsx" that includes the SMILES codes and substrate IDs you want to make into .sdf files
        
        #       run "python SMILES_to_sdf.py" or "python ~/bin/SMILES_to_sdf.py"
        
2.	prep_MM_npsh.bash


    a.	Contributors: Beck
    
    b.	Run in folder with .sdf files that you want to perform conformational searches for
  
    c.	Files Needed: "mm_reference.com" (needs to be in the same folder as the .sdf files and have proper permissions)
  
        Using OnDemand, run a manual MacroModel job and save the .com file after removing the first two lines with input/output file names as “mm_reference.com"
        
    d.	If you need to add a metal or do a substructure search, this is possible but not implemented.
  
    e.	Be cognizant of your Schrodinger token usage when running this. If you have a lot of jobs to run, we suggest running overnight.
  
    f.	Can be run from the OnDemand terminal or from a notchpeak terminal. Uses notchpeak shared short nodes.
  
    g.	# How to use on the CHPC:
  
        # run this script by navigating into the directory with the .sdf files and “mm_reference.com” file 
        
        #       run "sbatch prep_MM_npsh.bash"
        
3.	MM_clustering_and_converting.bash

    a.	Contributors: Beck
  
    b.	Run from the same directly immediately following “prep_MM_npsh.bash”
  
    c.	Clusters according to the minimum Kelley Penalty value (can change the “n” in line 42 if you want a specific number of clusters), if there are ≤20 conformers for a structure (can be changed in line 35).
  
    d.	Can be run from the OnDemand terminal or from a notchpeak terminal. Uses notchpeak shared short nodes.
  
    e.	# How to use on the CHPC:
  
        # run this script in the same directory as you ran “prep_MM_npsh.bash”
        
        #       run "MM_clustering_and_converting.bash"
        
4.	check_manual_clustering.bash


    a.	Run from the “output” folder of “MM_clustering_and_converting.bash” containing .sdf files with conformational ensembles
  
    b.	Converts .sdf files to .com files, and if any were clustered incorrectly, it prints the name of the structure to the terminal for you to manually cluster and deletes the incorrect files
  
    c.	# How to use on the CHPC:
  
        # run this script from the “output” folder of “MM_clustering_and_converting.bash” containing .sdf files with conformational ensembles
        
        #       run "check_manual_clustering.bash"
        
    d.	Manual clustering: use sdfgin.bash script to make .com files
  
5.	Use group scripts to add job information to the generated com files, and submit them to the nodes on the CHPC 


6.	post_processing.bash


    a.	Checks .log files for normal termination, no imaginary frequencies, and that every .com file has a corresponding .log file. Makes intuitive directories to sort various errors.
  
    b.	Supplemental Scripts Needed: "log_check_for_post_processing.py " (needs to be in your bin folder)
  
    c.	Line 36 of “log_check_for_post_processing.py” needs to be modified depending on the number of linked jobs in your .log files. This is the number of optimizations +  + single points. 
  
    d.	# How to use on the CHPC:
  
        # Calls “log_check_for_post_processing.py” from your bin (ensure the correct permissions). 
        
        # Modify line 36 in “log_check_for_post_processing.py” to be the number of "Normal termination" in the file. Note that you should have one "Normal termination" for each job you run, keep in mind and "Opt + Freq" job is actually two jobs. Default is 4.
        
        # run this script from the directory containing Gaussian output files
        
        #       run "bash ~/bin/post_processing.bash"

