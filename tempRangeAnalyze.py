# =============================================================================
# Script for plotting the data generated from the tempRangeSubmit.py script.
# This will plot the ferromagnetic phase transition.
#
# Note that this has been modified with 'try' statements in order to be
# backwards compatible with old data files that did not contain a column
# for E^2.
#
# Author:           Max Graves
# Last Revision:    24-APR-2013
# =============================================================================

import subprocess, os, sys, argparse, glob
import pylab as pl
import darkPlots as dp

# =============================================================================
def parseCMD():
    desc= ('Ising 2D script for plotting phase transition.  Takes as its \
            argument the directory holding all of the data files that you ,\
            want to plot.  Also, has capabilities of "ghost" plotting')
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('direcN', help='name of data directory of interest')
    parser.add_argument('--include', '-i', type=float,
            default=0.95,
            help='fraction of data to include [0.01-0.99]')
    parser.add_argument('--ghost', '-g', action="store_true", dest="ghost",
            default=False,
            help='Give this flag for ghost plots')

    return parser.parse_args()

# =============================================================================
def main():

    args = parseCMD()
    os.chdir(args.direcN)

    try:
        temps, Es, Ms, Cvs = pl.loadtxt('ising2D_reduced.txt',
                unpack=True)
        print 'Found a reduced file!'
    except:
        print 'Didnt find reduced file.  Creating one'
        
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
            try:
                tempL = float(estLines[1].split()[-1])
            except:
                tempL = 10.
            try:
                tempH = float(estLines[2].split()[-1])
            except:
                tempH = 0.
            print 'temp: ',tempT

            temps = pl.append(temps, float(tempT))
            try:
                mcSteps, En, Mag, E2 = pl.loadtxt(f, unpack=True)
            except:
                mcSteps, En, Mag = pl.loadtxt(f, unpack=True)
                E2 = En**2

            # take last half of data, no matter size of array
            binAfter = int(args.include*En.size)

            Cv = (pl.average(E2)-(pl.average(En))**2)/(tempL**2*tempT**2)
            Es = pl.append(Es, pl.average(En[-binAfter:]))
            Ms = pl.append(Ms, pl.average(Mag[-binAfter:]))
            Cvs = pl.append(Cvs, Cv)

        # write temps, Es, Ms to file
        filename = 'ising2D_reduced.txt'
        fid = open(filename, 'w')
        fid.write('# size:  %s\n'%tempL)
        fid.write('# field:  %s\n'%tempH)
        fid.write('# %15s\t%15s\t%15s\t%15s\n'%('temps','Energies',
            'Magnetism','Specific Heat'))
        zipped = zip(temps, Es, Ms, Cvs)
        pl.savetxt(fid, zipped, fmt='%5.9f\t%5.9f\t%5.9f\t%5.9f')
        fid.close()
        print 'Data has been saved to: ',filename

    # plot energy vs. temp
    fig1 = pl.figure(1)
    p1 = fig1.add_subplot(111)
    if args.ghost:
        pl.scatter(temps, Es, color='white')
        dp.darkPlots(p1)
        pl.grid(True, color='white')
    else:
        pl.scatter(temps, Es, color='Lime')
        pl.grid(True)
    pl.xlabel('Temperature '+r'$[K]$', size=20)
    pl.ylabel('Energy', size=20)

    if args.ghost:
        pl.savefig("2D_EvsT_ghost.png", transparent=True)
    else:
        pl.savefig("2D_EvsT.png")

    # plot specific heat vs. temp
    fig2 = pl.figure(2)
    p2 = fig2.add_subplot(111)
    if args.ghost:
        pl.scatter(temps, Cvs, color='white')
        dp.darkPlots(p2)
        pl.grid(True, color='white')
    else:
        pl.scatter(temps, Cvs, color='Lime')
        pl.grid(True)
    pl.xlabel('Temperature '+r'$[K]$', size=20)
    pl.ylabel('Specific Heat', size=20)

    if args.ghost:
        pl.savefig("2D_CvsT_ghost.png", transparent=True)
    else:
        pl.savefig("2D_CvsT.png")

    # plot magnetization vs. temp
    fig3 = pl.figure(3)
    p3 = fig3.add_subplot(111)
    if args.ghost:
        pl.scatter(temps, pl.absolute(Ms), color='white')
        dp.darkPlots(p3)
        pl.grid(True, color='white')
    else:
        pl.scatter(temps, pl.absolute(Ms), color='Lime')
        pl.grid(True)
    pl.xlabel('Temperature '+r'$[K]$', size=20)
    pl.ylabel('abs(Magnetization)', size=20)

    if args.ghost:
        pl.savefig("2D_MvsT_ghost.png", transparent=True)
    else:
        pl.savefig("2D_MvsT.png")
    
    pl.show()

# =============================================================================
if __name__=='__main__':
    main()
