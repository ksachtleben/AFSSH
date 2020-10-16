
# coding: utf-8

# In[1]:

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import pandas as pd
from sys import argv

file = argv[1]

x,y,z = np.loadtxt(file, unpack=True)
df = pd.DataFrame({'x': x, 'y': y, 'z': z})


# In[8]:

fig = plt.figure()
ax = Axes3D(fig)
surf = ax.plot_trisurf(df.x, df.y, df.z, cmap=cm.jet, linewidth=0.1)
#fig.colorbar(surf, shrink=0.5, aspect=5)
plt.savefig('teste.svg')
plt.show()

