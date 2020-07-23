import csv
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

paper_data = r"C:\Development\DBSuitability_Benchmarking\python_calculating\paper_workflow_out.txt"

InFile = open(paper_data, 'r')
data_table = csv.reader(InFile, delimiter="\t")

db_names = []
db_suit_1 = []
db_suit_2 = []

header = True
for row in data_table:
    if header:
        header = False
        continue
    db_names.append(Path(row[1]).stem)
    db_suit_1.append(float(row[2]))
    db_suit_2.append(float(row[3]))

db_suit, num_db_hits, db_names = (list(t) for t in zip(*sorted(zip(db_suit_1, db_suit_2, db_names), reverse=True)))


pos = list(range(len(db_suit_1)))
width = 0.25

fig, ax = plt.subplots()

plt.bar(pos, db_suit_1, width, alpha=0.5, color='b')
#plt.bar([p + width for p in pos] , db_suit_2, width, alpha=0.5, color='r')

ax.set_ylabel('Suitability')

ax.set_xticks([p + 0.5 * width for p in pos])
ax.set_xticklabels([name for name in db_names], rotation='vertical')

plt.ylim([0, 1])

plt.legend(['alternative\npeptides\nnot considered', 'alternative\npeptides\nconsidered'], loc='upper right')

plt.show()

InFile.close()
