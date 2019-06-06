import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


def read_sec(path):
    sequences = []
    with open(path) as f:
        line_list = f.readlines()
        for i in range(len(line_list)):
            if line_list[i][0] == ">":
                sequences.append((line_list[i][1:].strip(), line_list[i+1].strip()))

    return sequences

def wd_kernel(seq1, seq2, d):
    result = 0
    beta_k = lambda d,k: (2*(d-k + 1))/(d*(d+1))
    for k in range(1,d+1):
        coeff = beta_k(d,k)
        for l in range(len(seq1) - k + 1):
            if seq1[l:l+k] == seq2[l:l+k]:
                result += coeff
    return result


sequs = read_sec("sequencesMSAfasta.sec")
results = np.empty((len(sequs), len(sequs)))
labels = [seq[0] for seq in sequs]
for i in range(len(sequs)):
    for j in range(len(sequs)):
        results[i][j] = wd_kernel(sequs[i][1], sequs[j][1], 3)

sns.heatmap(results, xticklabels=labels, yticklabels=labels)
plt.show()
plt.savefig("heatmap.png")
print(wd_kernel(sequs[0][1], sequs[5][1], 3))