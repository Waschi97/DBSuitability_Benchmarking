import os
from pathlib import Path
import csv
import random_fasta_portion
import numpy as np
import time

# -----------------------------------------------------------------------------
# input

comet_exe = r"C:\Development\analysis_executables\comet\comet.2019015.win64.exe"

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
out = "openMS_random_portion_out.txt"

# wanted FDR(s)
FDRs = [0.01, 0.02]

# set the number of samples you wish to be done
# percentages are calculated to be equidistant
num_test = 2

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
    num_iterations = round(1/ratio)
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
    
        decoy_db = f"{Path.cwd()}{os.path.sep}{Path(database).stem}_novo_decoy.FASTA"
        print("Calculating decoy database ...")
        os.system(f"DecoyDatabase -in {database} {novo_fasta} {crap_db} -out {decoy_db}")
        print("Done!\n")
  
        # run CometAdapter

        idXML = f"{Path.cwd()}{os.path.sep}{Path(mzML).stem}.idXML"
        print("Running CometAdapter ...")
        os.system(f"CometAdapter -in {mzML} -out {idXML} -database {decoy_db} -comet_executable {comet_exe} -precursor_mass_tolerance 20.0 -isotope_error 0/1/2/3 -allowed_missed_cleavages 2 -num_hits 10 -max_precursor_charge 6 -spectrum_batch_size 15000 -threads {comet_threads}")
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

        for FDR in FDRs:
            DBSuit_out = f"db_suit_out_{FDR}.tsv"
            print(f"Running DatabaseSuitability with FDR {FDR}...")
            os.system(f"DatabaseSuitability -in_id {fdr_idXML} -in_novo {novo_seqs} -in_spec {mzML} -out {DBSuit_out} -novor_fract 0.99 -FDR {FDR}")
            print("Done!\n")

            # extract wanted data from output

            # if DatabaseSuitability throws any error, we don't want to lose all the other data
            try:
                data_file = open(DBSuit_out, 'r')
            except FileNotFoundError:
                print(f"Suitability couldn't be calculated for ratio {ratio} with FDR {FDR}.")
                continue
            
            data_table = csv.reader(data_file, delimiter="\t")

            db_hits = -1
            novo_hits = -1
            suitability = -1
            for row in data_table:
                if row[0] == "#top_db_hits":
                    db_hits = int(row[1])
                if row[0] == "#top_novo_hits":
                    novo_hits = int(row[1])
                if row[0] == "db_suitability":
                    suitability = float(row[1])

            if db_hits < 0 or suitability < 0:
                print("Couldn't find data in DatabaseSuitability output .. Aborting!")
                exit()

            data_file.close()

            key = (ratio, FDR)
            
            if key in ratio_to_suit:
                ratio_to_suit[key].append((suitability, db_hits, novo_hits))
            else:
                ratio_to_suit[key] = [(suitability, db_hits, novo_hits)]

            if not keep_files:
                os.remove(DBSuit_out)

        # delete files to save storage
        if not keep_files:
            print("Removing files ...")
            os.remove(database)
            os.remove(decoy_db)
            os.remove(idXML)
            os.remove(indexed_idXML)
            os.remove(fdr_idXML)
            print("Done!\n")

output = open(out, 'w')
output.write("ratio,FDR\taverage_suitability\taverage_db_hits\taverage_novo_hits\t[suitability, #db_hits, #novo_hits]*N\n")
for key in ratio_to_suit:
    result_list = ratio_to_suit[key]
    line = ""
    avg_suit = 0
    avg_db_hits = 0
    avg_novo_hits = 0
    for elem in result_list:
        suitability = elem[0]
        db_hits = elem[1]
        novo_hits = elem[2]
        avg_suit += suitability
        avg_db_hits += db_hits
        avg_novo_hits += novo_hits
        line += f"\t[{suitability},{db_hits},{novo_hits}]"
    avg_suit = avg_suit/len(result_list)
    avg_db_hits = avg_db_hits/len(result_list)
    avg_novo_hits = avg_novo_hits/len(result_list)
    line = f"{key}\t{avg_suit}\t{avg_db_hits}\t{avg_novo_hits}" + line
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

    
