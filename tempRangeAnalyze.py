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
    Cvs = pl.array([])

    # fill arrays from data files
    for f in files:
        estFile = open(f,'r')
        estLines = estFile.readlines();
        tempT = float(estLines[0].split()[-1])
        tempL = float(estLines[1].split()[-1])
        tempH = float(estLines[2].split()[-1])
        print 'temp: ',tempT

        temps = pl.append(temps, float(tempT))
        mcSteps, En, Mag, E2 = pl.loadtxt(f, unpack=True)

        # take last half of data, no matter size of array
        binAfter = int(0.95*En.size)

        Cv = (pl.average(E2)-(pl.average(En))**2)/(tempL**2*tempT**2)

        Es = pl.append(Es, pl.average(En[-binAfter:]))
        Ms = pl.append(Ms, pl.average(Mag[-binAfter:]))
        Cvs = pl.append(Cvs, Cv)

    # write temps, Es, Ms to file
    filename = 'ising2D_reduced.txt'
    fid = open(filename, 'w')
    #fid.write('# temp:  %s\n'%T)
    fid.write('# %15s\t%15s\t%15s\n'%('temps','Energies','Magnetism'))
    zipped = zip(temps, Es, Ms)
    pl.savetxt(fid, zipped, fmt='%5.9f\t%5.9f\t%5.9f')
    fid.close()
    print 'Data has been saved to: ',filename

    # plot energy vs. temp
    fig1 = pl.figure(1)
    p1 = fig1.add_subplot(111)
    pl.scatter(temps, Es)
    pl.grid(True)
    pl.xlabel('Temperature '+r'$[K]$', size=20)
    pl.ylabel('Energy', size=20)

    # plot specific heat vs. temp
    fig2 = pl.figure(2)
    p2 = fig2.add_subplot(111)
    pl.scatter(temps, Cvs)
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
