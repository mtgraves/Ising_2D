# =============================================================================
# Script for plotting the data generated from the tempRangeSubmit.py script.
# This will plot the ferromagnetic phase transition.
#
# Author:           Max Graves
# Last Revision:    22-MAR-2013
# =============================================================================

import subprocess, os, sys, argparse, glob
import pylab as pl

# =============================================================================
def parseCMD():
    desc= ('Ising 2D script for plotting phase transition.  Takes as its \
            argument the directory holding all of the data files that you ,\
            want to plot.')

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('direcN', help='name of data directory of interest')
    return parser.parse_args()

# =============================================================================
def main():

    args = parseCMD()
    os.chdir(args.direcN)

    # holds all files created by phaseTransSubmit.py
    files = glob.glob("*.dat*")

    # arrays to hold data
    temps, Es, Ms = pl.array([]), pl.array([]), pl.array([])

    # fill arrays from data files
    for f in files:
        estFile = open(f,'r')
        estLines = estFile.readlines();
        tempT = estLines[0].split()[-1]
        print tempT

        temps = pl.append(temps, float(tempT))
        mcSteps, En, Mag = pl.loadtxt(f, unpack=True)

        # take last half of data, no matter size of array
        binAfter = int(0.8*En.size)

        Es = pl.append(Es, pl.average(En[-binAfter:]))
        Ms = pl.append(Ms, pl.average(Mag[-binAfter:]))

    # plot energy vs. temp
    fig1 = pl.figure(1)
    p1 = fig1.add_subplot(111)
    pl.scatter(temps, Es)
    pl.grid(True)
    pl.xlabel('Temperature '+r'$[K]$', size=20)
    pl.ylabel('Energy', size=20)
    
    # plot magnetization vs. temp
    fig2 = pl.figure(2)
    p2 = fig2.add_subplot(111)
    pl.scatter(temps, pl.absolute(Ms))
    pl.grid(True)
    pl.xlabel('Temperature '+r'$[K]$', size=20)
    pl.ylabel('abs(Magnetization)', size=20)
    
    pl.show()
    
# =============================================================================
if __name__=='__main__':
    main()
