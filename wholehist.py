# Astronomical Image Processing
# David Christopher Ragusa + Aidan Boxall

from astropy.io import fits
import matplotlib.pyplot as plt

with fits.open("mosaic.fits") as infile:  # with, for automatic closing
    data = infile[0].data.flatten()


# we have to do all the calculations manually as np.hist crashes

rangemax = 65000
rangemin = 0
numbins = 1000

width = (rangemax-rangemin) / numbins
histdata = [0 for i in xrange(numbins)]

for value in data:
    if value <= rangemax:
        histdata[int(value/width)] += 1    # select appropriate bin and add 1 to it

bins = [width*i for i in range(numbins)]

plt.bar(bins, histdata)  # have to do a bar instead of a histogram as array is too large

ax = plt.gca()
ax.get_yaxis().set_tick_params(direction='out')  # ticks outwards
ax.get_xaxis().set_tick_params(direction='out')

plt.xlabel('Counts')
plt.ylabel('Number of Pixels')

plt.show()
