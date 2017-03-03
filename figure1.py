hun = load("./ABcompartmentsFromHiC/Final/dic_hun_MCF250.pkl")
expressions = load("./Data/gene_expression/rna_seq/dic_FPKM.pkl")

def generate_heatmap():
    import numpy as np
    import math
    import scipy.stats
    keys = ['HUVEC', 'NHEK', 'K562', 'GM12878', 'HMEC', 'MCF7','T47D','IMR90']
    data = np.zeros((len(keys), len(keys)))
    for i in range(data.shape[0]):
        for j in range(i):
            print;print;print keys[i],keys[j]
            ratios = expressionCompartmentAnalysis(keys[i],keys[j],ret=1)
            real_ratio = [[x[0]/x[1] for x in r] for r in ratios]
            data[i,j] = scipy.stats.ks_2samp(real_ratio[1],real_ratio[2])[1]
            data[i,j] = -1*math.log(data[i,j],10)
    return data


import numpy as NP
from matplotlib import pyplot as PLT
from matplotlib import cm as CM

A = NP.random.randint(10, 100, 100).reshape(10, 10)
mask =  NP.tri(A.shape[0], k=-1)
A = NP.ma.array(A, mask=mask) # mask out the lower triangle
fig = PLT.figure()
ax1 = fig.add_subplot(111)
cmap = CM.get_cmap('jet', 10) # jet doesn't have white color
cmap.set_bad('w') # default value is 'k'
ax1.imshow(A, interpolation="nearest", cmap=cmap)
ax1.grid(True)
PLT.show()


maxValueInMap = data.max()
x = np.arange(maxValueInMap+1)
ys = [x+i for i in x]
line_segments = LineCollection([zip(x, y) for y in ys],
                linewidths=(0.5, 3, 1.5, 2),
                linestyles='solid', cmap = matplotlib.cm.Blues)

line_segments.set_array(x)

fig, ax = plt.subplots()
heatmap = ax.pcolor(data, cmap=plt.cm.Blues)
a = ax.set_xticks(np.arange(data.shape[0])+0.5, minor=False)
b = ax.set_yticks(np.arange(data.shape[0])+0.5, minor=False)
c = ax.set_xticklabels(keys, minor=False)
d = ax.set_yticklabels(keys, minor=False)
axcb = fig.colorbar(line_segments)