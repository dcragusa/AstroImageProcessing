# Astronomical Image Processing
# David Christopher Ragusa + Aidan Boxall

import scipy.optimize as optimisation
import matplotlib.pyplot as plt
import numpy as np
import cPickle

with open('mosaiccatalog.pkl', 'rb') as catalogfile:  # with, for automatic closing
    catalog = cPickle.load(catalogfile)

mags = [item['mag'] for item in catalog]  # list of magnitudes

values, base = np.histogram(mags, bins=80)
print values
print base
cumulative = np.cumsum(values)  #evaluate the cumulative
cumulative_error = np.sqrt(cumulative)

# Split data into 3 sections: the lower curving section, the middle
# straight line section, and the higher roll-off section.
# This is so we can perform a fit directly in Python.

startbase = []
startcum = []
starterr = []
starterrlist = []
straightbase = []
straightcum = []
straighterrlist = []
straighterr = []
endbase = []
endcum = []
enderrlist = []
enderr = []

for index in range(len(base[:-1])):
    if index in range(0,9):  # plot lower section
        startbase.append(base[index])
        startcum.append(cumulative[index])
        starterrlist.append(cumulative_error[index])
    elif index in range(9,45):  # plot middle straight line section
        straightbase.append(base[index])
        straightcum.append(cumulative[index])
        straighterrlist.append(cumulative_error[index])
    else:  # plot higher section
        endbase.append(base[index])
        endcum.append(cumulative[index])
        enderrlist.append(cumulative_error[index])

# convert all the lists to numpy arrays
startbase = np.array(startbase)
for index, value in enumerate(starterrlist):
    starterr.append(1.5*(1.0/np.log(10))*(value/startcum[index]))
starterr = np.array(starterr)
startcum = np.array(np.log10(startcum))

straightbase = np.array(straightbase)
print straightbase
for index, value in enumerate(straighterrlist):
    straighterr.append(1.5*(1.0/np.log(10))*(value/straightcum[index]))
straighterr = np.array(straighterr)
straightcum = np.array(np.log10(straightcum))

endbase = np.array(endbase)
for index, value in enumerate(enderrlist):
    enderr.append(1.5*(1.0/np.log(10))*(value/endcum[index]))
enderr = np.array(enderr)
endcum = np.array(np.log10(endcum))


def loggraph():

    def func(x,m,c):  # straight line function
        return x**m + c

    params, matrix = optimisation.curve_fit(func, straightbase, straightcum, [0.0,0.0], straighterr)  # the actual optimisation
    m, c = params
    graderr = np.sqrt(matrix[0,0])
    print "Gradient: %s" % m
    print "Intercept: %s" % c
    print "Gradient error: %s" % graderr

    fit = [func(x,m,c) for x in straightbase]

    # plot the three sections, centre in blue
    plt.errorbar(straightbase, straightcum, straighterr, c='blue', fmt=' ')
    plt.errorbar(startbase, startcum, starterr, c='green', fmt=' ')
    plt.errorbar(endbase[:-1], endcum[:-1], enderr[:-1], c='green', fmt=' ')

    plt.plot(straightbase, fit, c='red')  # fit line

    x1,x2,y1,y2 = plt.axis()
    plt.axis((x1-0.5,22,y1,y2))

    plt.ylabel('log(N)')
    plt.xlabel('Apparent Magnitude')

    plt.show()


def normloggraph():

    correctedstraight = []
    for index, value in enumerate(straightcum):
        correctedstraight.append(value-0.6*straightbase[index])

    plt.ylabel('log(N) - 0.6m')
    plt.xlabel('Apparent Magnitude')

    plt.errorbar(straightbase, correctedstraight, straighterr)
    plt.show()

loggraph()
# normloggraph()
