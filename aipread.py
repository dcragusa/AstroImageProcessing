# Astronomical Image Processing
# David Christopher Ragusa + Aidan Boxall

import cPickle

with open('mosaiccatalog.pkl', 'rb') as catalogfile:
    catalog = cPickle.load(catalogfile)

print ' Location\tRadius\tTot Count\tMax Count\tMagnitude\tMag. Error'
for index,i in enumerate(catalog):
    print " %s\t%s\t%s%s%s\t\t%s\t%s" % (i['centre'],i['radius'],i['totcount'],('\t\t' if index < 4 else '\t'),i['peakcount'],i['mag'],i['errormag'])
