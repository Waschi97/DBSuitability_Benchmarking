import csv
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab

InFile = open(f"..{os.path.sep}python_calculating{os.path.sep}FDR_sampling_result.txt", 'r')
 
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
    data = row[4:]
    suits = []
    db_hits = []
    novo_hits = []
    for data_point in data:
        data_point = data_point[1:-1]
        suit, db_hit, novo_hit = data_point.split(',')
        suits.append(float(suit))
        db_hits.append(int(db_hit))
        novo_hits.append(int(novo_hit))

    if FDR not in FDR_to_data:
        FDR_to_data[FDR] = [(ratio, (suits, db_hits, novo_hits))]
    else:
        FDR_to_data[FDR].append((ratio, (suits, db_hits, novo_hits)))

num_FDRs = len(FDR_to_data.keys())
fig, axs = plt.subplots(2,2)

# Suitability vs FDR
leg = []
for FDR in FDR_to_data:
    x = []
    y = []
    for d in FDR_to_data[FDR]:
        x.append([FDR]*len(d[1][0]))
        y.append(d[1][0])
    for j in range(len(y)):
        axs[0,0].scatter(x[j], y[j], s=4)
        ratio = round(FDR_to_data[FDR][j][0],1)
        if ratio not in leg:
            leg.append(ratio)

axs[0,0].set_ylabel("Suitability")
axs[0,0].yaxis.set_major_locator(plt.MultipleLocator(0.05))
axs[0,0].yaxis.grid(linestyle='dotted')
axs[0,0].set_xlabel("FDR")
axs[0,0].set_xticks([f for f in FDR_to_data])
axs[0,0].legend(leg, title="ratio", loc="center left", bbox_to_anchor=(1,0.5))
axs[0,0].set_title("Suitability against used FDR")

# Suitability vs Ratio, FDR 0.01
data = FDR_to_data[0.01]
X = []
Y = []
for d in data:
    y = d[1][0]
    x = [d[0]]*len(y)
    axs[0,1].scatter(x, y, color='b', s=4)
    X.append(d[0])
    Y.append(sum(y)/len(y))

axs[0,1].plot(X,Y, ':r')
axs[0,1].set_ylabel("Suitability")
axs[0,1].yaxis.grid(linestyle='dotted')
axs[0,1].yaxis.set_major_locator(plt.MultipleLocator(0.05))
axs[0,1].set_xlabel("ratio")
axs[0,1].set_xticks(X)
axs[0,1].set_title("Suitability against Database Ratio, FDR 0.01")

# db_hits vs Ratio, FDR 0.01
X = []
Y = []
tick_labels = []
for d in data:
    y = d[1][1]
    x = [d[0]]*len(y)
    axs[1,0].scatter(x, y, color='b', s=4)
    X += x
    Y += y
    tick_labels.append(x[0])

coef = np.polyfit(X,Y,1)
function = np.poly1d(coef)    
line, = axs[1,0].plot(X,function(X), ':r')
line.set_label(f"slope:\n{round(coef[0]*0.1,1)} hits per 10 % increase")
axs[1,0].set_ylabel("number of\ndatabase hits")
axs[1,0].yaxis.set_major_locator(plt.MultipleLocator(2500))
axs[1,0].yaxis.grid(linestyle='dotted')
axs[1,0].set_xlabel("ratio")
axs[1,0].set_xticks(tick_labels)
axs[1,0].set_title("Database Hits against Database Ratio, FDR 0.01")
axs[1,0].legend()

# novo_hits vs Ratio, FDR 0.01
X = []
Y = []
tick_labels = []
for d in data:
    y = d[1][2]
    x = [d[0]]*len(y)
    axs[1,1].scatter(x, y, color='b', s=4)
    X += x
    Y += y
    tick_labels.append(x[0])

coef = np.polyfit(X,Y,1)
function = np.poly1d(coef)     
line, = axs[1,1].plot(X,function(X), ':r')
line.set_label(f"slope:\n{round(coef[0]*0.1,1)} hits per 10 % increase")
axs[1,1].set_ylabel("number of\ndeNovo hits")
axs[1,1].yaxis.set_major_locator(plt.MultipleLocator(500))
axs[1,1].yaxis.grid(linestyle='dotted')
axs[1,1].set_xlabel("ratio")
axs[1,1].set_xticks(tick_labels)
axs[1,1].set_title("DeNovo Hits against Database Ratio, FDR 0.01")
axs[1,1].legend()

fig.set_size_inches(15,8)
fig.tight_layout()
fig.savefig("ratio_sampling_with_FDR.png")
