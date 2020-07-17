import os
from pathlib import Path
import fileinput
import subprocess

comet_path = "C:\\Development\\analysis_executables\\comet\\comet.2019015.win64"
comet_params = "C:\\Development\\python_scripts\\comet.params.new"

PeptideProphet_path = ""

# file which should be searched
in_mzML = "C:\\Development\\TOPPAS_Outputs\\TOPPAS_out\\007-FileConverter-out\\QEP2_2016_0512_RJ_13_human.mzML"

# file containing the novo protein
novo_fasta = "C:\\Development\\TOPPAS_Outputs\\TOPPAS_out\\011-IDFileConverter-out\\QEP2_2016_0512_RJ_13_human.FASTA"

# file containing all databases which should be tested
db_file_list = "one_db.txt"

with open(db_file_list, 'r') as DBs:
    lines = DBs.readlines()
    db_names = [x.strip() for x in lines]
    
for database in db_names:
    # combine db with novo protein and build decoy database
    decoy_db_name = f"{Path.cwd()}{os.path.sep}{Path(database).stem}_novo_decoy.FASTA"
    os.system(f"DecoyDatabase -in {database} {novo_fasta} -out {decoy_db_name}")
  
    # run comet search
    pep_xml = f"{Path.cwd()}{os.path.sep}{Path(in_mzML).stem}.pep.xml"
    # print(f"{comet_path} -P{comet_params} -D{decoy_db_name} -N{pep_xml} {in_mzML}")
    process = subprocess.Popen([f"{comet_path}", f"-P{comet_params}", f"-D{decoy_db_name}", f"-N{pep_xml}", f"{in_mzML}"],
                     stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE,
                     shell=True)
    stdout, stderr = process.communicate() # could handle some errors with this ...
    
    pepXMLs.open("pep_xml_list.txt", 'w')
    pepXMLs.write(pep_xml)
    pepXMLs.close()
    
    # modify output

    execfile("modify_comet_pepxml.py")
    re_ranked_pepXML = f"{Path(pepXML).stem}_unNovored.pep.xml"
    
    # establishe FDR with PeptideProphet

    # calculate db suitability
