import os
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from pathlib import Path

if os.path.sep == "\\":
    WIN = True
else:
    WIN = False

def readCSVdata(File):
    InFile = open(File, 'r')
    data_table = csv.reader(InFile, delimiter="\t")

    suit_map = {}
    
    header = True
    for row in data_table:
        if header:
            header = False
            continue
        if WIN:
            suit_map[Path(row[1]).stem] = (float(row[2]), int(row[3]))
        else:
            suit_map[Path(row[1].replace("\\","/")).stem] = (float(row[2]), int(row[3]))

    InFile.close()
    
    return suit_map

#-------------------------------------------------------------------------------------
# input

paper_data = f"..{os.path.sep}python_calculating{os.path.sep}paper_workflow_out.txt"
openMS_data = f"..{os.path.sep}python_calculating{os.path.sep}openMS_workflow_out.txt"

paper_map = readCSVdata(paper_data)
openMS_map = readCSVdata(openMS_data)

suit_paper = []
suit_openMS = []
hits_paper = []
hits_openMS = []
names = []
for key in paper_map:
    if key in openMS_map:
        names.append(key[8:])
        suit_paper.append(paper_map[key][0])
        suit_openMS.append(openMS_map[key][0])
        hits_paper.append(paper_map[key][1])
        hits_openMS.append(openMS_map[key][1])
    else:
        print(f"Skipped {key} from paper output because it wasn't found in the OpenMS output.")


suit_paper, suit_openMS, hits_paper, hits_openMS, names = (list(t) for t in zip(*sorted(zip(suit_paper, suit_openMS, hits_paper, hits_openMS, names), reverse=True)))

pos = list(range(len(suit_paper)))
width = 0.25

params = {'legend.fontsize': 'medium',
          'figure.figsize': (10, 6),
         'axes.labelsize': 'medium',
         'axes.titlesize':'large',
         'xtick.labelsize':'small',
         'ytick.labelsize':'medium'}
pylab.rcParams.update(params)

fig, ax = plt.subplots()

plt.bar(pos, suit_paper, width, alpha=0.5, color='b')
plt.bar([p + width for p in pos] , suit_openMS, width, alpha=0.5, color='r')

ax.set_yticks(list(np.arange(0,1.05,0.05)))
ax.set_ylabel('Suitability')
ax.yaxis.grid(linestyle='dotted')

ax.set_xticks([p + 0.5 * width for p in pos])
ax.set_xticklabels([name for name in names], rotation=35)

plt.ylim([0, 1])

plt.legend(['Paper', 'OpenMS'], loc='upper right')

plt.title("Suitability of Various FASTA Files for Analysis\nof LC-MS/MS Data from Human Tryptic Peptides")

plt.tight_layout()

plt.savefig("paper_vs_openms.png")

