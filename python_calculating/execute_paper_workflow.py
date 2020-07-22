__author__ = 'Tom Waschischeck'

import os
from pathlib import Path
import subprocess
from lxml import etree
from lxml import objectify
import numpy as np

def minProb4FDR(root, FDR):        
    for e in root.analysis_summary.getchildren():
        if 'peptideprophet_summary' in e.tag:
            prophet_summary = e.getchildren()

            for error_data in prophet_summary:
                if "roc_error_data" not in error_data.tag:
                    continue
                
                if error_data.attrib['charge'] == 'all':
                    previous_prob = 1
                    for data_point in error_data.getchildren():
                        correct = int(data_point.attrib['num_corr'])
                        incorrect = int(data_point.attrib['num_incorr'])
                        current_FDR =  float(incorrect / float(correct+incorrect))

                        # If current FDR is now worse than wanted, the previous FDR was as close
                        # as we can get. -> return the previous probability
                        # If this happens at the first error data point '1' is returned. That's
                        # obviously a bs FDR and should be handled.
                        if (current_FDR > FDR):
                            return previous_prob
                        
                        previous_prob = float(data_point.attrib['min_prob'])

                    # We get here if the wanted FDR is very high and PeptideProphet doesn't
                    # even consider a probability for this many incorrect hits.
                    # The worst probability PeptideProphet exports is returned.
                    return previous_prob
                
            print("No data for 'all' charges found!\n")
            return 2
        
    print("No 'peptideprophet_summary' found in xml! Did PeptideProphet run?\n")
    return 3

def onlyNovo(proteins):
    for protein in proteins:
        # one db protein is enough
        if 'concatenated_peptides' not in protein:
            return False
    return True

def calculateDBSuit(root, MinProb, consider_alternatives):
    num_db = 0
    num_novo = 0
    for e in root.msms_run_summary.getchildren():
        if 'spectrum_query' not in e.tag:
            continue

        # look for first search hit
        top_hit = -1
        for f in e.search_result.getchildren():
            if 'search_hit' in f.tag:
                top_hit = f
                break

        if top_hit == -1:
            continue

        # check if hit probability is at least MinProb (calculated from FDR)
        prob = float(top_hit.analysis_result.peptideprophet_result.attrib['probability'])        
        if prob < MinProb:
            continue

        possible_proteins = []
        # check alternative protein hits
        if(consider_alternatives):
            for p in top_hit.getchildren():
                if 'alternative_protein' not in p.tag:
                    continue
                if 'DECOY_' not in p.attrib['protein']:
                    possible_proteins.append(p.attrib['protein'])

        if 'DECOY_' not in top_hit.attrib['protein']:
            possible_proteins.append(top_hit.attrib['protein'])

        # only decoy hits found
        if len(possible_proteins) == 0:
            continue
        
        if(onlyNovo(possible_proteins)):
           num_novo += 1
           continue
        num_db += 1
        
    return num_db / float(num_novo + num_db)

# -----------------------------------------------------------------------------
# input

comet_path = r"C:\Development\analysis_executables\comet\comet.2019015.win64"
comet_params = r"C:\Development\python_scripts\comet.params.new"

PeptideProphet_path = r"C:\TPP\bin\PeptideProphetParser.exe"

# file which should be searched
mzML = r"C:\Development\TOPPAS_Outputs\TOPPAS_out\007-FileConverter-out\QEP2_2016_0512_RJ_13_human.mzML"

# file containing the novo protein
novo_fasta = r"C:\Development\TOPPAS_Outputs\TOPPAS_out\011-IDFileConverter-out\QEP2_2016_0512_RJ_13_human.FASTA"

# file containing a list of all databases which should be tested
db_file_list = "databases_human_search.txt"

# file with a crap database
crap_db = r"C:\Development\data\db\crap.FASTA"

# wanted FDR
FDR = 0.01

# Set to 'True' if you want to keep the temporary files produced by this script
keep_files = False

# output file
out = "paper_workflow_out.txt"

# -----------------------------------------------------------------------------
# calculations

with open(db_file_list, 'r') as DBs:
    lines = DBs.readlines()
    db_names = [x.strip() for x in lines]

output = open(out, 'w')
output.write("mzML\tdatabase\tsuitability_no_alternative_proteins\tsuitability_with_alternative_proteins\n")
for database in db_names:
    print(f"Starting with {database}:")
    
    # combine db with novo protein and build decoy database
    
    decoy_db = f"{Path.cwd()}{os.path.sep}{Path(database).stem}_novo_decoy.FASTA"
    print("Calculating decoy database ...")
    os.system(f"DecoyDatabase -in {database} {novo_fasta} {crap_db} -out {decoy_db}")
    print("Done!\n")
  
    # run comet search
    
    pepXML = f"{Path.cwd()}{os.path.sep}{Path(mzML).stem}"
    print("Running Comet ...")
    process = subprocess.Popen([f"{comet_path}", f"-P{comet_params}", f"-D{decoy_db}", f"-N{pepXML}", f"{mzML}"],
                     stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE,
                     shell=True)
    stdout, stderr = process.communicate() # could handle some errors with this ...
    pepXML = pepXML + ".pep.xml"
    print("Done!\n")
    
    pepXMLs = open("pep_xml_list.txt", 'w')
    pepXMLs.write(pepXML)
    pepXMLs.close()
    
    # modify output

    print("Re-ranking Comet output ...")
    with open("modify_comet_pepxml.py", "rb") as source_file:
        code = compile(source_file.read(), "modify_comet_pepxml.py", "exec")
    exec(code)
    print("Done!\n")
    
    re_ranked_pepXML = f"{Path(Path(pepXML).stem).stem}_unNovored.pep.xml"
    
    # establishe FDR with PeptideProphet
    print("Running PeptideProphet ...")
    # couldn't get subprocess to work with PeptideProphet, so this is done the 'old', 'unsecure' way
    os.system(f"{PeptideProphet_path} {re_ranked_pepXML} ZERO")
    print("Done!\n")
    re_ranked_pepXML = re_ranked_pepXML.lower() # PeptideProphet changes the filename to be all lower case for some reason..
    
    with open(re_ranked_pepXML) as InFile:
        root = objectify.parse(InFile).getroot()

    # find MinProb that translates to wanted FDR
    MinProb = minProb4FDR(root, FDR)
    if MinProb == 1:
        print("FDR to low! All hits would be filtered. Choose a higher FDR! Aborted!")
        exit()
    if MinProb > 1:
        exit()    
    
    # calculate db suitability

    print("Calculating database suitability ... ")
    db_suit_no_alternatives = calculateDBSuit(root, MinProb, False)
    db_suit_with_alternatives = calculateDBSuit(root, MinProb, True)
    print("Done!\n")

    # write output

    output.write(f"{mzML}\t{database}\t{db_suit_no_alternatives}\t{db_suit_with_alternatives}\n")

    # delete files to save storage
    if not keep_files:
        print("Removing files ...")
        os.remove(decoy_db)
        os.remove(pepXML)
        os.remove(re_ranked_pepXML)
        print("Done!\n")

output.close()
print(f"Output written to {os.path.abspath(out)}.")

