import sys
import numpy as np

precision = []
sensitive = []
result_type = []
file = "score.txt"
with open(file, 'r') as f:
    for line in f:
        linestr = line.strip()
        linestrlist = linestr.split('\t')
        if float(linestrlist[0]) > 0.2:
            print(float(linestrlist[0]), float(linestrlist[1]))
            precision.append(float(linestrlist[0]))
            sensitive.append(float(linestrlist[1]))

precision = np.array(precision)
sensitive = np.array(sensitive)
with open('mean.txt', 'w') as f:
    f.write(str(np.mean(precision)) + '\t' + str(np.mean(sensitive)))
print(len(precision), np.mean(precision), np.mean(sensitive))
