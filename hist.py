# Astronomical Image Processing
# David Christopher Ragusa + Aidan Boxall

from astropy.io import fits  # http://www.astropy.org/
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy as np
import cPickle

with fits.open("mosaic.fits") as infile:  # with, for automatic closing
    data = infile[0].data.flatten()

counts = []

rangemax = 3700
rangemin = 3200  # area of interest

print 'Getting values in range...'
counts = [value for value in data if rangemin <= value <= rangemax]
print 'Done.\n'

def getmean(indata):  # need custom functions as array is too large for numpy
    mean = sum(indata) / float(len(indata))
    return mean

def getstdev(indata, mean):
    stdevlist = []
    for value in indata:
        stdevlist.append(value**2 - mean**2)
    stdev = np.sqrt((1 / float(len(indata) - 1)) * sum(stdevlist))
    return stdev

print 'Getting mean...'
mean = getmean(counts)
print 'Done.\n'

print 'Getting standard deviation...'
stdev = getstdev(counts, mean)
print 'Done.\n'

plt.hist(counts, bins=100)
x = np.linspace(3350, 3600, 600)                    #
plt.plot(x, mlab.normpdf(x,mean,stdev) * 50000000)  # superimposition of normal fit

ax = plt.gca()
ax.get_yaxis().set_tick_params(direction='out')  # ticks outwards
ax.get_xaxis().set_tick_params(direction='out')

plt.xlabel('Counts')
plt.ylabel('Number of Pixels')

plt.show()

outdata = {'mean':mean, 'stdev':stdev}

print 'Outputting...'
with open('mosaic_bgdata.pkl', 'wb') as outfile:  # saves data to pickle file for later
    cPickle.dump(outdata, outfile)
print 'Done.\n'
