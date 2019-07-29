import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
sns.set(style="white", color_codes=True)


f = open('data/f_xv.dat', 'r')
X = []
Y = []
for line in f:
    if len(line.split()) == 1:
        continue
    values = line.split()
    X.append(values[1])
    Y.append(values[2])

# Use JointGrid directly to draw a custom plot
grid = sns.JointGrid(X, Y, space=0, height=6, ratio=50)
grid.plot_joint(plt.scatter, color="g")
grid.plot_marginals(sns.rugplot, height=1, color="g")
