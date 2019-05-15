"""
Author: Eric Cotner
Date: 2019-05-014

Matplotlib does not have a very useful way of visualizing large numbers of
time series. It is possible to plot the series with partial transparency
and a monochrome color scheme so that the intensity of the plot gives an
idea of the density of the time series, but this tends to overemphasize
regions of the plot with high slopes.

Perhaps it would be a better idea to plot a heatmap or histogram of a
collection of time series. This would not suffer from the aforementioned
problem, and could be more colorful.
"""

from matplotlib.pyplot import hist2d
from matplotlib.colors import LogNorm, Normalize
import numpy as np

def ts_heatmap(Y, t=None, pts_per_bin=None, n_interp=None,
				colorscale=None, **kwargs):
	"""
	Produces a heatmap of a collection of time series.

	Arguments:
		Y : array
			Array of shape (m,n) containing time series values, where m is
			number of time series, and n is number of points per time
			series. Expects the length of each time series to be the same.
		t : array
			Array containing the values of the time steps. Should be length
			m, corresponding to the length of the 1st axis of Y.

	Returns:
		None
	"""
	aspect_ratio = (1+np.sqrt(5))/2
	default_n_bins_y = 300
	default_n_bins_x = default_n_bins_y * aspect_ratio
	default_pts_per_bin = 10

	Y = np.array(Y)
	assert len(Y.shape) == 2, "Y should have 2 dimensions"
	if t is None:
		t = np.arange(Y.shape[-1])
	t = np.array(t)
	if n_interp is None:
		n_interp = round(default_n_bins_x/(len(t)-1) - 1)
	if pts_per_bin is None:
		pts_per_bin = default_pts_per_bin
	if (colorscale == 'linear') or (colorscale is None):
		norm = Normalize()
	elif colorscale == 'log':
		norm = LogNorm()

	# Combine domain and range
	Y0 = np.ones(Y.shape) * t[np.newaxis,:]
	Y = np.stack([Y0, Y], axis=2)

	# Interpolate between points
	Y1 = Y[:,:-1,:]; Y2 = Y[:,1:,:]
	Y = list()
	for x in np.linspace(0, 1, n_interp+2)[:-1]:
		Y.append(x*Y1 + (1-x)*Y2)
	Y = np.concatenate(Y, axis=0)

	# Convert into collection of points
	Y = Y.reshape(-1, 2)

	# Make 2d histogram
	bins_x = (n_interp+1)*(len(t)-1)
	bins_y = min(default_n_bins_y, max(len(Y)//pts_per_bin, 10))
	bins = [bins_x, bins_y]
	return hist2d(Y[:,0], Y[:,1], bins=bins, norm=norm, **kwargs)






