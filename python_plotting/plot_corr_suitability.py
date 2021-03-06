import csv
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab

InFile = open(f"..{os.path.sep}python_calculating{os.path.sep}corrected_suitability_human.txt", 'r')
 
data_table = csv.reader(InFile, delimiter="\t")
FDR_to_data = {}
header = True
for row in data_table:
    if header:
        header = False
        continue
    # handle key
    key = row[0]
    key = key[:-1]
    key = key[1:]
    ratio, FDR = key.split(',')
    ratio = float(ratio)
    FDR = float(FDR)

    #handle values
    data = row[1:]
    suits = []
    corr_suits = []
    db_hits = []
    novo_hits = []
    corr_novo_hits = []
    for data_point in data:
        data_point = data_point[1:-1]
        suit, corr_suit, db_hit, novo_hit, corr_novo_hit = data_point.split(',')
        suits.append(float(suit))
        corr_suits.append(float(corr_suit))
        db_hits.append(int(db_hit))
        novo_hits.append(int(novo_hit))
        corr_novo_hits.append(float(corr_novo_hit))

    if FDR not in FDR_to_data:
        FDR_to_data[FDR] = [(ratio, (suits, corr_suits, db_hits, novo_hits, corr_novo_hits))]
    else:
        FDR_to_data[FDR].append((ratio, (suits, corr_suits, db_hits, novo_hits, corr_novo_hits)))

num_FDRs = len(FDR_to_data.keys())
fig, axs = plt.subplots(2,1)

# get data for FDR 0.01
data = FDR_to_data[0.01]

# Suitability vs Ratio, FDR 0.01
X = []
Y = []
Y_corr = []
for d in data:
    y = d[1][0]
    y_corr = d[1][1]
    x = [d[0]]*len(y)
    axs[0].scatter(x, y, color='g', s=4)
    axs[0].scatter(x, y_corr, color='g', s=4)
    X.append(d[0])
    Y.append(sum(y)/len(y))
    Y_corr.append(sum(y_corr)/len(y_corr))

axs[0].plot(X,Y, ':r', label="Suitability")
axs[0].plot(X,Y_corr, ':b', label="Corrected Suitability")
axs[0].plot(X,X, '-', color="lightgrey")
axs[0].set_ylabel("Suitability")
axs[0].yaxis.grid(linestyle='dotted')
axs[0].yaxis.set_major_locator(plt.MultipleLocator(0.05))
axs[0].set_xlabel("ratio")
axs[0].set_xticks(X)
axs[0].set_title("Suitability against Database Ratio, FDR 0.01")
axs[0].legend()

# novo_hits vs Ratio, FDR 0.01
X = []
Y = []
Y_corr = []
tick_labels = []
for d in data:
    y = d[1][3]
    y_corr = d[1][4]
    x = [d[0]]*len(y)
    axs[1].scatter(x, y, color='g', s=4)
    axs[1].scatter(x, y_corr, color='g', s=4)
    X.append(d[0])
    Y.append(sum(y)/len(y))
    Y_corr.append(sum(y_corr)/len(y_corr))

axs[1].plot(X,Y, ':r', label="deNovo hits before correction")
axs[1].plot(X,Y_corr, ':b', label="deNovo hits after correction")
axs[1].set_ylabel("number of\ndeNovo hits")
axs[1].yaxis.set_major_locator(plt.MultipleLocator(15000))
axs[1].yaxis.grid(linestyle='dotted')
axs[1].set_xlabel("ratio")
axs[1].set_xticks(X)
axs[1].set_title("DeNovo Hits against Database Ratio, FDR 0.01")
axs[1].legend()

fig.set_size_inches(10,8)
fig.tight_layout()
fig.savefig("corr_suitability.png")
