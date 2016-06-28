# Astronomical Image Processing
# David Christopher Ragusa + Aidan Boxall

from astropy.io import fits  # http://www.astropy.org/
import numpy as np
import cPickle, math, os


# Variable Setup
apradius = 3  # minimum radius in pixels
catalog = []

with fits.open("mosaic.fits") as infile:  # with, for automatic closing
    header = infile[0].header
    origdata = np.copy(infile[0].data)

zeropoint = header['MAGZPT']
zeropointerror = header['MAGZRR']
zperrsquared = zeropointerror**2
workdata = np.copy(origdata)
shape = origdata.shape
mask = np.ones(shape)
shape = [i-1 for i in shape]
background = 3423
galaxymin = 3553


# Artifact Mask Coords - (y0, y1, x0, x1)
artifactlist = [[3708,3801,2110,2159],  [3381,3440,2447,2491],  [2286,2336,2108,2153],  [1401,1452,2069,2111],  [424,451,1027,1044],
                [2222,2357,878,934],    [2700,2836,944,1010],   [3203,3418,746,808],    [4077,4115,540,582],    [4015,4050,1442,1475],
                [0,10,969,1719],        [10,26,1286,1526],      [0,47,1634,1646],       [24,117,1411,1454],     [115,140,1289,1540],
                [139,313,1390,1479],    [313,333,1016,1705],    [332,355,1640,1650],    [333,426,1390,1481],    [426,452,1103,1654],
                [424,452,1027,1043],    [452,480,1410,1464],    [480,2562,1427,1450],   [2562,2963,1425,1458],  [2963,2993,1402,1425],
                [2962,3163,1417,1463],  [3073,3355,1244,1583],  [3355,3930,1414,1456],  [3930,4610,1429,1446],  [2697,2807,2552,2569]]


def newfitsfile(name, origdata):
    'Write a new fits file given the location and data.'
    hdu = fits.PrimaryHDU(origdata)
    try:
        hdu.writeto(name)
    except IOError:
        os.remove(name)  # clears if already exists
        hdu.writeto(name)


def fixpixels():
    'Finds dead pixels and replaces them with the average of their neighbours.'

    print 'Fixing dead pixels...'
    deadpoints = np.where(origdata==0)
    deadpoints = zip(deadpoints[0],deadpoints[1])  # list of dead pixel coords

    for i,j in deadpoints:
        toavg = []
        if origdata[i-1, j] != 0:             # if left pixel is not zero,
            toavg.append(origdata[i-1,j])     #  add to toavg
        if origdata[i, j-1] != 0:             # bottom
            toavg.append(origdata[i,j-1])
        if origdata[(i+1)%shape[0], j] != 0:  # right
            toavg.append(origdata[(i+1)%shape[0], j])
        if origdata[i, (j+1)%shape[1]] != 0:  # top
            toavg.append(origdata[i, (j+1)%shape[1]])
        origdata[i,j] = np.mean(toavg)        # fix original data and working data
        workdata[i,j] = np.mean(toavg)
    print 'Done.\n'


def maskedge():
    'Applies an apradius-sized mask around the edges of the image.'

    print 'Masking image edges...'
    for j in xrange(apradius):                        # masks left side
        for i in xrange(shape[0]):
            mask[i,j] = 0
            workdata[i,j] = 0
    for i in xrange(apradius):                        # bottom
        for j in xrange(shape[1]):
            mask[i,j] = 0
            workdata[i,j] = 0
    for j in xrange(shape[1]-apradius+1,shape[1]+1):  # right
        for i in xrange(shape[0]):
            mask[i,j] = 0
            workdata[i,j] = 0
    for i in xrange(shape[0]-apradius+1,shape[0]+1):  # top
        for j in xrange(shape[1]):
            mask[i,j] = 0
            workdata[i,j] = 0
    print 'Done.\n'


def maskartifacts():
    'Applies the mask to artifacts listed in artifactlist.'

    print 'Masking stars/artifacts...'
    for coords in artifactlist:
        for i in xrange(coords[0],coords[1]+1):
            for j in xrange(coords[2],coords[3]+1):  # masks a rectangle
                mask[i,j] = 0
                workdata[i,j] = 0
    print 'Done.\n'


def getdisklist(centre, minradius, maxradius):
    'Returns coords of the disk lying between minradius and maxradius of centre.'

    disk = []
    circlelist = range(-maxradius,maxradius+1)
    for i in circlelist:
        for j in circlelist:
            if minradius**2 <= i**2 + j**2 <= maxradius**2:  # squared distances
                disk.append((centre[0]+i,centre[1]+j))

    return disk


def getapsize(centre):
    'Returns aperture size of object at centre.'
    radius = 3   # start with minimum radius
    mingrad = 26        # tolerance of variation - 1 sigma

    while True:
        disk1 = getdisklist(centre,radius,radius+1)     # 1 pixel from radius
        disk2 = getdisklist(centre,radius+1,radius+2)   # 2 pixels from radius
        for coords in disk1:
            if not mask[coords]:    # if inner disk has hit mask,
                return radius       #  return current radius

        disk1mean = np.mean([origdata[coords] for coords in disk1])
        disk2mean = np.mean([origdata[coords] for coords in disk2])

        if disk2mean >= disk1mean:  # if mean of outer disk is larger or equal
            break
        if disk1mean - disk2mean <= mingrad:  # if difference is smaller than the tolerance
            break
        radius += 1

    return radius


def getmaxlist():
    "Returns a list of current maxima in the 'workdata' image."
    maximum = np.max(workdata)
    maxpoints = np.where(workdata==maximum)
    maxcoords = zip(maxpoints[0],maxpoints[1])  # format into a list of xy coords
    maxcoordscopy = list(maxcoords)             # copy of list

    for coords in maxcoords:
        if not mask[coords]:                # if coords are in mask,
            workdata[coords] = 0            #  remove them from workdata
            maxcoordscopy.remove(coords)    #  remove them from list to return

    return maxcoordscopy


def getmag(count):
    "Returns magnitude from counts."
    return zeropoint - 2.5*math.log10(count)


def geterrormag(count):
    "Returns error on the magnitude."
    return np.sqrt(zperrsquared + (2.5*np.log10(np.exp(1))/count)**2 * count)


def getcirclecount(centre, curapradius):
    "Returns the count of an object centred at 'centre' with radius 'curapradius'."
    hitmask = False
    circle = []
    bgcountlist = []
    apcirclelist = range(-curapradius,curapradius+1)

    for i in apcirclelist:
        for j in apcirclelist:                  # returns coords of disk with radius 'curapradius'
            if i**2 + j**2 <= curapradius**2:
                circle.append((centre[0]+i,centre[1]+j))

    countsum = 0
    bgcountsumlist = []

    for coords in circle:
        try:
            if not mask[coords]:            # break if a coord lies in a masked area
                hitmask = True
            countsum += origdata[coords]    # add counts to list
            workdata[coords] = 0            # remove from workdata
            mask[coords] = 0                # mask area
        except:
            hitmask = True

    if hitmask:
        return False

    for coords in circle:
        maskgrid[coords] = 59000  # white out area on maskgrid

    countsum -= background * len(circle)  # subtract background count level

    return countsum


# MAIN PROGRAM START

fixpixels()

maskedge()

maskartifacts()

print 'Outputting masked mosaic...'
maskpoints = np.where(mask==0)
maskpoints = zip(maskpoints[0],maskpoints[1])  # format into a list of xy coords
maskgrid = np.copy(origdata)
for coords in maskpoints:
    maskgrid[coords] = 0

newfitsfile("maskon.fits", maskgrid)
print 'Done.\n'

print 'Cataloguing stars...'
while True:
    maxcoordslist = getmaxlist()
    while len(maxcoordslist) == 0:
        maxcoordslist = getmaxlist()  # if all maxima lie in the mask, get a new list
    centrecount = origdata[maxcoordslist[0][0], maxcoordslist[0][1]]
    if centrecount <= galaxymin:  # if the count lies below 5 sigma, stop cataloguing
        break
    passes = 0
    for coords in maxcoordslist:
        if mask[coords]:
            curapradius = getapsize(coords)  # if the coord is not masked, get aperture size
            try:
                count = getcirclecount(coords, curapradius)
                if count > 0:  # count can be negative from background subtraction
                    catalog.append({'centre':coords,'peakcount':origdata[coords],'totcount':count,'mag':getmag(count),
                                    'errormag':geterrormag(count),'radius':curapradius})
                    passes += 1
            except:
                pass
    print "%s points found for %s: %s passes" % (len(maxcoordslist), centrecount, passes)

with open('mosaiccatalog.pkl', 'wb') as catalogfile:
    cPickle.dump(catalog, catalogfile)
print 'Done.\n'

print "Catalog length: %s galaxies.\n" % len(catalog)

print 'Outputting mask and sources...'
newfitsfile("galactic_location.fits", maskgrid)  # with mask and sources
print 'Done.\n'
