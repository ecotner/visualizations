"""
Author: Eric Cotner
Date: 2019-05-14

Tests the time series heatmap implementation.
"""

import numpy as np
import matplotlib.pyplot as plt
from timeseries_heatmap import ts_heatmap

x = np.linspace(0, 4*np.pi, 10)
m1 = 1000
m2 = 300
Y1 = np.ones((m1, len(x))) * np.sin(x) + 0.6*np.random.randn(m1, len(x))
Y2 = np.ones((m2, len(x))) * 1.5*np.sin(x+np.pi) + 0.6*np.random.randn(m2, len(x))
Y = np.concatenate([Y1,Y2], axis=0)

fig, (ax1, ax2) = plt.subplots(ncols=2)
plt.sca(ax1)
ts_heatmap(Y, colorscale='linear', cmap='plasma')
plt.sca(ax2)
plt.plot(Y.T, color='black', alpha=1/256)
plt.show()
