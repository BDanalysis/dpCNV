import sys
import numpy as np

result_start = []
result_end = []
result_type = []
file = sys.argv[1]
with open(file, 'r') as f:
    for line in f:
        linestr = line.strip()
        linestrlist = linestr.split('\t')
        result_start.append(int(linestrlist[1]))
        result_end.append(int(linestrlist[2]))
        result_type.append(linestrlist[3])


truth_start = []
truth_end = []
truth_type = []
with open("../../GroundTruthCNV", 'r') as f:
    line = f.readline()
    for line in f:
        linestr = line.strip('\n')
        linestrlist = linestr.split('\t')
        truth_start.append(int(linestrlist[0]))
        truth_end.append(int(linestrlist[1]))
        if linestrlist[3] == 'gain':
            truth_type.append("gain")
        else:
            truth_type.append("loss")


count = 0
for i in range(len(result_type)):
    for j in range(len(truth_type)):
        if truth_start[j] <= result_start[i] <= truth_end[j] and truth_type[j] == result_type[i]:
            if result_end[i] <= truth_end[j]:
                print(result_start[i], result_end[i], truth_start[j], truth_end[j], 1)
                count += (result_end[i] - result_start[i] + 1)
                print(count)
            elif result_end[i] >= truth_end[j]:
                print(result_start[i], result_end[i], truth_start[j], truth_end[j], 2)
                count += (truth_end[j] - result_start[i] + 1)
                print(count)
            break
        elif truth_start[j] >= result_start[i] and truth_type[j] == result_type[i]:
            if truth_start[j] <= result_end[i] <= truth_end[j]:
                print(result_start[i], result_end[i], truth_start[j], truth_end[j],  3)
                count += (result_end[i] - truth_start[j] + 1)
                print(count)
            elif result_end[i] >= truth_end[j]:
                print(result_start[i], result_end[i], truth_start[j], truth_end[j], 4)
                count += (truth_end[j] - truth_start[j] + 1)
                print(count)
            break

result_count = 0
for i in range(len(result_start)):
    result_count += (result_end[i] - result_start[i] + 1)

truth_count = 0
for i in range(len(truth_start)):
    truth_count += (truth_end[i] - truth_start[i] + 1)

print(count, result_count, truth_count)
output = open("score.txt", "a")
output.write(str(count/result_count) + '\t' + str(count/truth_count) + '\n')
print(str(count/result_count) + '\t' + str(count/truth_count) + '\n')
