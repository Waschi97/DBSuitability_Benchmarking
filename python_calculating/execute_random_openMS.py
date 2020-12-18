import os
from pathlib import Path
import csv
import random_fasta_portion
import numpy as np
import time

# -----------------------------------------------------------------------------
# input

comet_exe = r"comet.exe"

comet_threads = 5

# file which should be searched
mzML = r"C:\Development\TOPPAS_Outputs\TOPPAS_out\007-FileConverter-out\QEP2_2016_0512_RJ_13_human.mzML"

# file containing de novo sequences derived from given mzML
novo_seqs = r"C:\Development\TOPPAS_Outputs\TOPPAS_out\006-NovorAdapter-out\QEP2_2016_0512_RJ_13_human.idXML"

# file containing the novo protein
novo_fasta = r"C:\Development\TOPPAS_Outputs\TOPPAS_out\011-IDFileConverter-out\QEP2_2016_0512_RJ_13_human.FASTA"

# the database of which the sampled series should be created
db = r"C:\Development\data\db\for_human_analysis\uniprot_human.FASTA"

# file with a crap database
crap_db = r"C:\Development\data\db\crap.FASTA"

# output file
out = "corrected_suitability_human.txt"

# wanted FDR(s)
FDRs = [0.01]

# set the number of samples you wish to be done
# percentages are calculated to be equidistant
num_test = 10

# Set to 'True' if you want to keep the temporary files produced by this script
keep_files = False

# -----------------------------------------------------------------------------
# calculations

ratios = np.linspace(0,1,num_test+1)
ratios = ratios[1:] # delete the 0
ratios = ratios[::-1] # reverse the list

start = time.perf_counter()
ratio_to_suit = {}
for ratio in ratios:
    num_iterations = round(1/ratio) * 2

    print(f"Starting with ratio: {ratio}\n")

    i = 0
    while i < num_iterations:
        i += 1
        print("-----------------------------------------------------------------")
        print(f"Iteration number {i} of {num_iterations} beginning ...")    
        # calculate portioned fasta file
        print("Calculating portioned database ...")
        data = random_fasta_portion.get_random_data(db, ratio)
        database = f"random_database_{ratio}.FASTA"
        with open(database, 'w') as File:
            File.write(data)
        print("Done!\n")
          
        # combine db with novo protein and build decoy database
        
        combined_db = f"{Path.cwd()}{os.path.sep}{Path(database).stem}_novo_combined.FASTA"
        print("Calculating decoy database ...")
        os.system(f"DecoyDatabase -in {database} {crap_db} {novo_fasta} -out {combined_db}")
        print("Done!\n")
  
        # run CometAdapter

        idXML = f"{Path.cwd()}{os.path.sep}{Path(mzML).stem}.idXML"
        print("Running CometAdapter ...")
        os.system(f"CometAdapter -in {mzML} -out {idXML} -database {combined_db} -comet_executable {comet_exe} -precursor_mass_tolerance 20.0 -isotope_error 0/1/2/3 -missed_cleavages 2 -num_hits 10 -max_precursor_charge 6 -spectrum_batch_size 15000 -threads {comet_threads}")
        print("Done!\n")

        # run PeptideIndexer

        indexed_idXML = f"{Path(idXML).stem}_indexed.idXML"
        print("Running PeptideIndexer ...")
        os.system(f"PeptideIndexer -in {idXML} -out {indexed_idXML} -fasta {combined_db}")
        print("Done!\n")

        # run DatabaseSuitability

        for FDR in FDRs:
            DBSuit_out = f"db_suit_out_{FDR}.tsv"
            print(f"Running DatabaseSuitability with FDR {FDR}...")
            os.system(f"DatabaseSuitability -in_id {indexed_idXML} -in_novo {novo_seqs} -in_spec {mzML} -database {database} -novo_database {novo_fasta} -out {DBSuit_out} -algorithm:reranking_cutoff_percentile 0.01 -algorithm:FDR {FDR}")
            print("Done!\n")

            # extract wanted data from output

            # if DatabaseSuitability throws any error, we don't want to lose all the other data
            try:
                data_file = open(DBSuit_out, 'r')
            except FileNotFoundError:
                print(f"Suitability couldn't be calculated for ratio {ratio} with FDR {FDR}.")
                continue
            
            data_table = csv.reader(data_file, delimiter="\t")

            db_hits = " "
            novo_hits = " "
            corr_novo_hits = " "
            suitability = " "
            corr_suitability = " "
            for row in data_table:
                if row[0] == "#top_db_hits":
                    db_hits = int(row[1])
                if row[0] == "#top_novo_hits":
                    novo_hits = int(row[1])
                if row[0] == "#corrected_novo_hits":
                    corr_novo_hits = float(row[1])
                if row[0] == "db_suitability":
                    suitability = float(row[1])
                if row[0] == "corrected_suitability":
                    corr_suitability = float(row[1])

            if db_hits == " " or novo_hits == " " or corr_novo_hits == " " or suitability == " " or corr_suitability == " ":
                print("Couldn't find data in DatabaseSuitability output .. Aborting!")
                exit()

            data_file.close()

            key = (ratio, FDR)
            
            if key in ratio_to_suit:
                ratio_to_suit[key].append((suitability, corr_suitability, db_hits, novo_hits, corr_novo_hits))
            else:
                ratio_to_suit[key] = [(suitability, corr_suitability, db_hits, novo_hits, corr_novo_hits)]

            if not keep_files:
                os.remove(DBSuit_out)

        # delete files to save storage
        if not keep_files:
            print("Removing files ...")
            os.remove(database)
            os.remove(combined_db)
            os.remove(idXML)
            os.remove(indexed_idXML)
            print("Done!\n")

output = open(out, 'w')
output.write("ratio,FDR\t[suitability, corrected_suitability, #db_hits, #novo_hits, #corrected_novo_hits]*N\n")
for key in ratio_to_suit:
    result_list = ratio_to_suit[key]
    line = f"{key}"
    avg_suit = 0
    avg_db_hits = 0
    avg_novo_hits = 0
    for elem in result_list:
        suitability = elem[0]
        corr_suitability = elem[1]
        db_hits = elem[2]
        novo_hits = elem[3]
        corr_novo_hits = elem[4]
        line += f"\t[{suitability},{corr_suitability},{db_hits},{novo_hits},{corr_novo_hits}]"
    output.write(f"{line}\n")
output.close()
print(f"Output written to {os.path.abspath(out)}.")
end = time.perf_counter()
run_time = end - start
h = 0
m = 0
s = run_time
if s > 60:
    m = int(run_time/60)
    s = run_time % 60
if m > 60:
    h = int(m/60)
    m = m % 60
time_str = f"{h} h {m} min {s} sec"
print(f"Execution done in {time_str}.")

    
