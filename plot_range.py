# =============================================================================
# Script for plotting the data generated from the tempRangeAnalyze.py script.
# This takes reduced file and plots magnetization, energy, and specific heat
# vs. temperature.
#
# Author:           Max Graves
# Last Revision:    26-MAR-2013
# =============================================================================

import subprocess, os, sys, argparse, glob
import DerApproximator as derr
import pylab as pl
import numpy as np
from scipy import interpolate

# =============================================================================
def parseCMD():
    desc= ('Ising 2D script for plotting phase transition.  Takes as its \
            argument the directory holding all of the data files that you ,\
            want to plot.')

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('direcN', help='name of reduced file of interest')
    return parser.parse_args()

def derivatives(x,y):
    '''
    calculate dy by 4-point center differencing using array slices

    frac{y[i-2] - 8y[i-1] + 8[i+1] - y[i+2]}{12h}

    y[0] and y[1] must be defined by lower order methods
    and y[-1] and y[-2] must be defined by lower order methods
    '''

    dy = np.zeros(y.shape,np.float) #we know it will be this size
    h = x[1]-x[0] #this assumes the points are evenely spaced!
    dy[2:-2] = (y[0:-4] - 8*y[1:-3] + 8*y[3:-1] - y[4:])/(12.*h)

    dy[0] = (y[1]-y[0])/(x[1]-x[0])
    dy[1] = (y[2]-y[1])/(x[2]-x[1])
    dy[-2] = (y[-2] - y[-3])/(x[-2] - x[-3])
    dy[-1] = (y[-1] - y[-2])/(x[-1] - x[-2])

    return dy
# =============================================================================
def main():
   
    # load in reduced data from file
    args = parseCMD()
    fileName = args.direcN
    temps, Es, Ms = pl.loadtxt(fileName, unpack=True)

    combine = pl.sort(pl.column_stack((temps,Es)), axis=0)

    EsSorted = pl.array([])
    TsSorted = pl.array([])

    for i in combine:
        TsSorted = pl.append(TsSorted, i[0])
        EsSorted = pl.append(EsSorted, i[1])

    tempsMore = pl.arange(TsSorted[0],TsSorted[-1], 0.01)

    tck = interpolate.splrep(TsSorted, EsSorted, s=23)

    Esmore = interpolate.splev(tempsMore,tck,der=0)

    dy = derivatives(TsSorted,Esmore) 

    # plot energy vs. temp
    fig1 = pl.figure(1)
    p1 = fig1.add_subplot(111)
    pl.scatter(temps, Es)
    pl.grid(True)
    pl.xlabel('Temperature '+r'$[K]$', size=20)
    pl.ylabel('Energy', size=20)
   
    # plot Cv vs. temp
    fig2 = pl.figure(2)
    p2 = fig2.add_subplot(111)
    pl.scatter(tempsMore, dy)
    pl.grid(True)
    pl.xlabel('Temperature '+r'$[K]$', size=20)
    pl.ylabel('Specific Heat', size=20)
    
    # plot magnetization vs. temp
    fig3 = pl.figure(3)
    p3 = fig3.add_subplot(111)
    pl.scatter(temps, pl.absolute(Ms))
    pl.grid(True)
    pl.xlabel('Temperature '+r'$[K]$', size=20)
    pl.ylabel('abs(Magnetization)', size=20)
    
    pl.show()
# =============================================================================
if __name__=='__main__':
    main()
