import math
from pathlib import Path

FASTA1 = "C:\\Development\\TOPPAS_OutPuts\\TOPPAS_out\\011-IDFileConverter-out\\FU_2016_0606_RJ_28_pooled.fasta"
FASTA2 = "C:\\Development\\TOPPAS_OutPuts\\TOPPAS_out\\011-IDFileConverter-out\\QEP1_2017_0810_RJ_26_ratfishMCX.fasta"
ratio = 0.5

def number_of_fasta_entries_needed(lines, ratio):
    count = 0
    for line in lines:
        if line[0] == '>':
            count += 1
    return math.ceil(count * ratio)

def extract_N_entries(lines, N):
    data = ""
    count = 0
    for line in lines:
        if line[0] == '<':
            count += 1
        if count > N:
            return data
        data += line
    return data

file1 = open(FASTA1,'r')
file2 = open(FASTA2,'r')

lines1 = file1.readlines()
lines2 = file2.readlines()

file1.close()
file2.close()

N1 = number_of_fasta_entries_needed(lines1, ratio)
N2 = number_of_fasta_entries_needed(lines2, 1-ratio)

data = extract_N_entries(lines1, N1) + extract_N_entries(lines2, N2)

output = open("C:\\Development\\data\\db\\" + Path(FASTA1).stem + "_" + Path(FASTA2).stem + "_merged.fasta", 'w')
output.write(data)
output.close()
