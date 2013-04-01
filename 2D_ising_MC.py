#==============================================================================
# This classical monte carlo code will simulate 2-D Ising model
#
# Author:               Max Graves
# Last Modified:        21-MAR-2013
#==============================================================================

import random, math, os
import numpy as np
import argparse

#==============================================================================
def parseCMD():
    """
    Parse the command line.
    """
    parser = argparse.ArgumentParser(description = 'simulate 2-D Ising model')
    parser.add_argument('--temp', '-T', type=float,
            help='enter temperature in kelvin')
    parser.add_argument('--field', '-H', type=float, default=0.0,
            help='enter magnetic field strength')
    parser.add_argument('--length', '-L', type=float,
            help='enter the length of the 2-d lattice')
    parser.add_argument('--exchange', '-J', type=float, default=1,
            help='enter the value of the exchange constant')
    parser.add_argument('--sweeps', '-s', type=int,
            help='enter the number of MC sweeps')
    return parser.parse_args()

#==============================================================================
def totalEnergy(aa,L,J):
    """
    Sum products of neighboring spin sites.  Return
    the energy as well as the sum of all of the spins individually.
    """
    w, h, totalSpin, spinSum = 0, 0, 0.0, 0.0
    while (h <= L-1):
        totalSpin += aa[h,w]
        if ((h == 0) and (w==0)):
            #print 'left top row'
            spinSum += (aa[0,1] + aa[0,-1] + aa[1,0] + aa[-1,0])*aa[0,0]
            w += 1
        elif ((h==0) and (w<L-1)):
            #print 'middle top row'
            spinSum += (aa[1,w] + aa[-1,w] + aa[0,w-1] + aa[0,w+1])*aa[0,w]
            w += 1
        elif ((h==0) and (w==L-1)):
            #print 'right top row'
            spinSum += (aa[1,w] + aa[-1,w] + aa[0,0] + aa[0,w-1])*aa[0,w]
            w = 0
            h += 1
        elif ((h<L-1) and (w==0)):
            #print 'left middle'
            spinSum += (aa[h+1,w] + aa[h-1,w] + aa[h,1] + aa[h,-1])*aa[h,0]
            w += 1
        elif ((h<L-1) and (w<L-1)):
            #print 'middle'
            spinSum += (aa[h+1,w] + aa[h-1,w] + aa[h,w+1] + aa[h,w-1])*aa[h,w]
            w += 1
        elif ((h<L-1) and (w==L-1)):
            #print 'right middle'
            spinSum += (aa[h+1,w] + aa[h-1,w] + aa[h,0] + aa[h,w-1])*aa[h,w]
            w = 0
            h += 1
        elif ((h==L-1) and (w==0)):
            #print 'left bottom'
            spinSum += (aa[0,w] + aa[h-1,w] + aa[h,1] + aa[h,-1])*aa[h,w]
            w += 1
        elif ((h==L-1) and (w<L-1)):
            #print 'middle bottom row'
            spinSum += (aa[0,w] + aa[h-1,w] + aa[h,w+1] + aa[h,w-1])*aa[h,w]
            w += 1
        elif ((h==L-1) and (w==L-1)):
            #print 'right bottom'
            spinSum += (aa[0,w] + aa[h-1,w] + aa[h,0] + aa[h,w-1])*aa[h,w]
            break
     
    E = - 0.5 * J * spinSum     # divide by two because of double counting
    return E, totalSpin

# =============================================================================
def energyChange(aa,L,J,w,h):
    """
    Calculate energy change for a single electron flipping.
    """
    spinSum = 0.0
    if ((h == 0) and (w==0)):
        spinSum += (aa[0,1] + aa[0,-1] + aa[1,0] + aa[-1,0])*aa[0,0]
    elif ((h==0) and (w<L-1)):
        spinSum += (aa[1,w] + aa[-1,w] + aa[0,w-1] + aa[0,w+1])*aa[0,w]
    elif ((h==0) and (w==L-1)):
        spinSum += (aa[1,w] + aa[-1,w] + aa[0,0] + aa[0,w-1])*aa[0,w]
    elif ((h<L-1) and (w==0)):
        spinSum += (aa[h+1,w] + aa[h-1,w] + aa[h,1] + aa[h,-1])*aa[h,0]
    elif ((h<L-1) and (w<L-1)):
        spinSum += (aa[h+1,w] + aa[h-1,w] + aa[h,w+1] + aa[h,w-1])*aa[h,w]
    elif ((h<L-1) and (w==L-1)):
        spinSum += (aa[h+1,w] + aa[h-1,w] + aa[h,0] + aa[h,w-1])*aa[h,w]
    elif ((h==L-1) and (w==0)):
        spinSum += (aa[0,w] + aa[h-1,w] + aa[h,1] + aa[h,-1])*aa[h,w]
    elif ((h==L-1) and (w<L-1)):
        spinSum += (aa[0,w] + aa[h-1,w] + aa[h,w+1] + aa[h,w-1])*aa[h,w]
    elif ((h==L-1) and (w==L-1)):
        spinSum += (aa[0,w] + aa[h-1,w] + aa[h,0] + aa[h,w-1])*aa[h,w]
    
    # simply flipping the sign of the energy gives energy change
    delE =  2.0 * J * spinSum     # divide by two because of double counting
    return delE

# =============================================================================
# begin main function
# =============================================================================
def main():

    # assign variables, define constants
    args = parseCMD()
    T,H,L = float(args.temp), float(args.field), int(args.length)
    J, s = float(args.exchange), int(args.sweeps)
    k_B = 1     # = 1.3806503 * pow(10,-23)

    # Create random numpy array of spins as either +/- 1
    latt = (2.0000*np.random.randint(2, size=pow(L,2))-1.0000).reshape(L,L)

    # calculate initial energy and macroscopic magnetism
    E, totalSpin = totalEnergy(latt, L, J)
    #M = np.fabs(1.0*totalSpin/(1.0*pow(L,2)))
    M = 1.0*totalSpin/(1.0*pow(L,2))
    
    # define arrays for mcSteps, Energies, Magnetism
    mcSteps, Es, Ms = np.arange(s), np.array([E]), np.array([M])
    E2 = np.array([E*E])

    # keep track of acceptance
    a = 0       # accepted moves
    r = 0       # rejected moves

    for step in mcSteps:
        # change the spin of one electron randomly
        over, down = random.randint(0,L-1), random.randint(0,L-1)

        # check change in energy
        delE = energyChange(latt, L, J, over, down)

        # calculate Boltzmann factor
        Boltz = np.exp(-1.0*delE/(k_B*T))

        if (delE <= 0):
            reject = False
            E += delE
            if (latt[down,over]==1):
                latt[down,over] = -1
                totalSpin -= 2
            else:
                latt[down,over] = 1
                totalSpin += 2
        else:
            n = random.random()
            if (n <= Boltz):
                E += delE
                reject = False
                if (latt[down,over]==1):
                    latt[down,over] = -1
                    totalSpin -= 2
                else:
                    latt[down,over] = 1
                    totalSpin += 2
            else:
                reject = True

        if reject==True:
            r += 1
        else:
            a += 1

        # calculate magnetism (not absolute value)
        M = 1.0*totalSpin/(1.0*pow(L,2))

        # store magnetism, energy in array
        Es = np.append(Es, E)
        Ms = np.append(Ms, M)
        E2 = np.append(E2, E*E)
    
    print 'acceptance ratio: ', 1.0*a/(r+a)

    if os.path.exists('./data/'):
        os.chdir('./data/')
    else:
        os.mkdir('./data/')
        os.chdir('./data/')

    filename = 'ising2D_L%s_s%s_Temp%s.dat'%(int(L), s, T)
    fid = open(filename, 'w')
    fid.write('# temp:  %s\n'%T)
    fid.write('# size:  %s\n'%L)
    fid.write('# field:  %s\n'%H)
    fid.write('# %15s\t%15s\t%15s\t%15s\n'%('mcSteps','Energies','Magnetism','Energy^2'))
    zipped = zip(mcSteps, Es, Ms, E2)
    np.savetxt(fid, zipped, fmt='%5.9f\t%5.9f\t%5.9f\t%5.9f')
    fid.close()
    print 'Data has been saved to: ',filename

# =============================================================================
if __name__=="__main__":
    main()
