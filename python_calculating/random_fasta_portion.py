import math
from pathlib import Path
import random
import sys
import os

def build_map_of_entries(lines):
    count = 0
    mapping = {}
    identifier = ""
    sequence = ""
    for line in lines:
        if line[0] == '>':
            mapping[identifier] = sequence
            count += len(sequence)
            identifier = ""
            sequence = ""
            identifier = line[:-1]
            continue
        sequence += line[:-1]
    mapping.pop("")
    print(f"Found {count} total amino acids in FASTA file.")
    return mapping

def number_of_amino_acids_needed(mapping, ratio):
    count = 0
    for identifier in mapping:
        count += len(mapping[identifier])
    return math.ceil(count * ratio)

def extract_random_entries_with_total_N_AS(mapping, N):
    data = ""
    cur_AS = 0
    while cur_AS < N:
        identifier = random.choice(list(mapping.keys())) # random entry
        sequence = mapping[identifier]
        data += identifier + "\n" + sequence + "\n"
        cur_AS += len(sequence)
        mapping.pop(identifier) # avoid double selection
    print(f"Entries with {cur_AS} total amino acids exported.")
    return data

def get_random_data(FASTA, ratio):
    file = open(FASTA,'r')

    lines = file.readlines()

    file.close()

    mapping = build_map_of_entries(lines)
    
    N = number_of_amino_acids_needed(mapping, ratio)

    print(f"Trying to export random entries with total number of amino acids {N}.")
    
    return extract_random_entries_with_total_N_AS(mapping, N)

if __name__ == '__main__':
    FASTA = sys.argv[1]

    ratio = float(sys.argv[2])

    data = get_random_data(FASTA, ratio)

    output = open(f"{Path(FASTA).parent}{os.path.sep}{Path(FASTA).stem}_{ratio}_random.fasta", 'w')
    output.write(data)
    output.close()
