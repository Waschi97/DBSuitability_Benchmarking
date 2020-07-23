import os
from pathlib import Path
import csv

# -----------------------------------------------------------------------------
# input

comet_exe = r"C:\Development\analysis_executables\comet\comet.2019015.win64.exe"

# file which should be searched
mzML = r"C:\Development\TOPPAS_Outputs\TOPPAS_out\007-FileConverter-out\QEP2_2016_0512_RJ_13_human.mzML"

# file containing de novo sequences derived from given mzML
novo_seqs = r"C:\Development\TOPPAS_Outputs\TOPPAS_out\006-NovorAdapter-out\QEP2_2016_0512_RJ_13_human.idXML"

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
out = "openMS_workflow_out.txt"

# -----------------------------------------------------------------------------
# calculations

with open(db_file_list, 'r') as DBs:
    lines = DBs.readlines()
    db_names = [x.strip() for x in lines]

output = open(out, 'w')
output.write("mzML\tdatabase\tsuitability\t#db_hits\n")
for database in db_names:
    print(f"Starting with {database}:")
    
    # combine db with novo protein and build decoy database
    
    decoy_db = f"{Path.cwd()}{os.path.sep}{Path(database).stem}_novo_decoy.FASTA"
    print("Calculating decoy database ...")
    os.system(f"DecoyDatabase -in {database} {novo_fasta} {crap_db} -out {decoy_db}")
    print("Done!\n")
  
    # run CometAdapter

    idXML = f"{Path.cwd()}{os.path.sep}{Path(mzML).stem}.idXML"
    print("Running CometAdapter ...")
    os.system(f"CometAdapter -in {mzML} -out {idXML} -database {decoy_db} -comet_executable {comet_exe} -precursor_mass_tolerance 20.0 -isotope_error 0/1/2/3 -allowed_missed_cleavages 2 -num_hits 10 -max_precursor_charge 6 -spectrum_batch_size 15000")
    print("Done!\n")

    # run PeptideIndexer

    indexed_idXML = f"{Path(idXML).stem}_indexed.idXML"
    print("Running PeptideIndexer ...")
    os.system(f"PeptideIndexer -in {idXML} -out {indexed_idXML} -fasta {decoy_db}")
    print("Done!\n")

    # run FDR

    fdr_idXML = f"{Path(idXML).stem}_fdr.idXML"
    print("Running FalseDiscoveryRate ...")
    os.system(f"FalseDiscoveryRate -in {indexed_idXML} -out {fdr_idXML} -protein false -algorithm:use_all_hits -algorithm:add_decoy_peptides -algorithm:add_decoy_proteins")
    print("Done!\n")

    # run DatabaseSuitability
    
    DBSuit_out = "db_suit_out.tsv"
    print("Running DatabaseSuitability ...")
    os.system(f"DatabaseSuitability -in_id {fdr_idXML} -in_novo {novo_seqs} -in_spec {mzML} -out {DBSuit_out} -novor_fract 0.99 -FDR 0.01")
    print("Done!\n")

    # extract wanted data from output
    
    data_file = open(DBSuit_out, 'r')
    data_table = csv.reader(data_file, delimiter="\t")

    db_hits = -1
    suitability = -1
    for row in data_table:
        if row[0] == "#top_db_hits":
            db_hits = int(row[1])
        if row[0] == "db_suitability":
            suitability = float(row[1])

    if db_hits < 0 or suitability < 0:
        print("Couldn't find data in DatabaseSuitability output .. Aborting!")
        exit()

    data_file.close()
    
    output.write(f"{mzML}\t{database}\t{suitability}\t{db_hits}\n")

    # delete files to save storage
    if not keep_files:
        print("Removing files ...")
        os.remove(decoy_db)
        os.remove(idXML)
        os.remove(indexed_idXML)
        os.remove(fdr_idXML)
        os.remove(DBSuit_out)
        print("Done!\n")

output.close()
print(f"Output written to {os.path.abspath(out)}.")

    
