import math
from pathlib import Path
import random
import sys

def number_of_fasta_entries_needed(lines, ratio):
    count = 0
    for line in lines:
        if line[0] == '>':
            count += 1
    return count, math.ceil(count * ratio)

def extract_N_random_entries(lines, N, total):
    data = ""
    index = 0
    entry_numbers = random.sample(list(range(1,total+1)),N)
    for line in lines:
        if line[0] == '>':
            index += 1
        if index not in entry_numbers:
            continue
        data += line
    return data

def get_random_data(FASTA, ratio):
    file = open(FASTA,'r')

    lines = file1.readlines()

    file1.close()

    total, N = number_of_fasta_entries_needed(lines1, ratio)

    return extract_N_random_entries(lines, N, total)

if __name__ == '__main__':
    FASTA = sys.argv[0]

    ratio = sys.argv[1]

    data = get_random_data(FASTA, ratio)

    output = open("C:\\Development\\data\\db\\" + Path(FASTA).stem + "_" + str(ratio) + "_random.fasta", 'w')
    output.write(data)
    output.close()
