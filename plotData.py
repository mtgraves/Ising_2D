# =============================================================================
# Plotting script to accompany Ising model code
#
# Author:           Max Graves
# Last Revised:     21-MAR-2013
# =============================================================================

import pylab as pl
import argparse, os, glob, sys

# =============================================================================
# parse the cmd line
# =============================================================================
def parseCMD():
    helpString = ('This script plots the Density vs Yield from the 2D Ising model.\
            It takes a file name as a command line argument.')

    parser = argparse.ArgumentParser(description=helpString)
    parser.add_argument('fileN', help='file to Plot')
    return parser.parse_args()

# =============================================================================
def getD(f):
    """
    CHANGE THIS TO STRIP OUT L,J,H, etc...
    strips the D value out of the data file name.
    """
    k = f[:-4]
    i, t = -1, True
    newt = ''
    while t==True:
        if k[i].isdigit():
            newt += k[i]
        else:
            t = False
        i -= 1
    return newt[::-1]

# =============================================================================
# begin main
# =============================================================================
def main():
    # read in the data
    args = parseCMD()
    fileName = args.fileN
    mcSteps, Es, Ms, E2 = pl.loadtxt(fileName, unpack=True)

    # get parameters from header of data file
    estFile = open(fileName,'r')
    estLines = estFile.readlines();
    T = float(estLines[0].split()[-1])
    L = float(estLines[1].split()[-1])
    H = float(estLines[2].split()[-1])

    # Compute specific heat
    Cv = (pl.average(E2) - (pl.average(Es))**2)/(L*L*T*T)
    print 'C_v: ',Cv

    # plot energy
    fig1 = pl.figure(1)
    ax = fig1.add_subplot(111)
    pl.ylabel('Energy', fontsize=20)
    pl.xlabel('MC steps', fontsize=20)
    ax.plot(mcSteps,Es, marker='o', linewidth=0,
            markerfacecolor='None', markeredgecolor='Orange')

    # plot macroscopic magnetization
    fig2 = pl.figure(2)
    bx = fig2.add_subplot(111)
    pl.ylabel('Magnetization', fontsize=20)
    pl.xlabel('MC steps', fontsize=20)
    bx.plot(mcSteps,Ms, marker='o', linewidth=0,
            markerfacecolor='None', markeredgecolor='Lime')

    pl.show()

# =============================================================================
if __name__=='__main__':
    main()
