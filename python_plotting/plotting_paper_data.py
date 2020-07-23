import os
import csv
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def readCSVdata(File):
    InFile = open(File, 'r')
    data_table = csv.reader(InFile, delimiter="\t")

    db_names = []
    db_suit = []
    db_hits = []

    header = True
    for row in data_table:
        if header:
            header = False
            continue
        db_names.append(Path(row[1]).stem)
        db_suit.append(float(row[2]))
        db_hits.append(int(row[3]))

    InFile.close()
    
    db_suit, db_hits, db_names = (list(t) for t in zip(*sorted(zip(db_suit, db_hits, db_names), reverse=True)))
    
    return db_names, db_suit, db_hits

#-------------------------------------------------------------------------------------
# input

paper_data = f"..{os.path.sep}python_calculating{os.path.sep}paper_workflow_out.txt"
openMS_data = f"..{os.path.sep}python_calculating{os.path.sep}openMS_workflow_out.txt"

paper_names, paper_suit, paper_hits = readCSVdata(paper_data)
openMS_names, openMS_suit, openMS_hits = readCSVdata(openMS_data)

pos = list(range(len(paper_suit)))
width = 0.25

fig, ax = plt.subplots()

plt.bar(pos, paper_suit, width, alpha=0.5, color='b')
plt.bar([p + width for p in pos] , openMS_suit, width, alpha=0.5, color='r')

ax.set_ylabel('Suitability')

ax.set_xticks([p + 0.5 * width for p in pos])
ax.set_xticklabels([name for name in paper_names], rotation='vertical')

plt.ylim([0, 1])

plt.legend(['Paper', 'OpenMS'], loc='upper right')

plt.show()

InFile.close()
