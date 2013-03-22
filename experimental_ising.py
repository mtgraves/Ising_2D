#====================================================================
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#====================================================================
# This classical monte carlo code will simulate 2-D Ising model
#           ** THIS VERSION USES FUNCTIONS **
#                  Author: Max Graves
#                      July 2012
#====================================================================
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#====================================================================

import random
import math
import numpy as np
import argparse
from numpy import *
import pylab

#====================================================================
# Calculate energy of system, summing the products of spins
# of nearest neighbor lattice sites (above, below, R, L)
#====================================================================

def energy(L, J, aa):
    # sum products of vertical neighboring spin sites of lattice
    w = 0
    h = 0
    spinSum = 0.0000
    while (h <= L-1):
        if (w < L-1) and (h < L-1):
            spinSum = spinSum + aa[h+1,w]*aa[h,w]
            print 'w: ', w, 'h: ', h, 'spinSum: ', spinSum
            w += 1
        elif (w == L-1) and (h != L-1):
            spinSum = spinSum + aa[h+1,w]*aa[h,w]
            print 'w: ', w, 'h: ', h, 'spinSum: ', spinSum
            h += 1
            w = 0
        elif (w < L-1) and (h == L-1):
            spinSum = spinSum + aa[h,w]*aa[0,w]
            print 'w: ', w, 'h: ', h, 'spinSum: ', spinSum
            w += 1
        elif (w == L-1) and (h == L-1):
            spinSum = spinSum + aa[h,w]*aa[0,w]
            print 'w: ', w, 'h: ', h, 'spinSum: ', spinSum
            break

    # sum products of horizontal neighboring spin sites of lattice
    w=0
    h=0
    while (w <= L-1):
        if (h < L-1) and (w < L-1):
            spinSum = spinSum + aa[h,w+1]*aa[h,w]
            print 'w: ', w, 'h: ', h, 'spinSum: ', spinSum
            h += 1
        elif (h == L-1) and (w != L-1):
            spinSum = spinSum + aa[h,w+1]*aa[h,w]
            print 'w: ', w, 'h: ', h, 'spinSum: ', spinSum
            w += 1
            h = 0
        elif (h < L-1) and (w == L-1):
            spinSum = spinSum + aa[h,w]*aa[h,0]
            print 'w: ', w, 'h: ', h, 'spinSum: ', spinSum
            h += 1
        elif (h == L-1) and (w == L-1):
            spinSum = spinSum + aa[h,w]*aa[h,0]
            print 'w: ', w, 'h: ', h, 'spinSum: ', spinSum
            break
    return spinSum

#====================================================================
#      begin main function
#====================================================================

def main():

    # set up command line parser options
    parser = argparse.ArgumentParser(description = 'simulate 2-D Ising model')
    parser.add_argument('--temp', '-T', type=float,
            help='enter temperature in kelvin')
    parser.add_argument('--field', '-H', type=float,
            help='enter magnetic field strength')
    parser.add_argument('--length', '-L', type=float,
            help='enter the length of the 2-d lattice')
    parser.add_argument('--exchange', '-J', type=float,
            help='enter the value of the exchange constant')
    parser.add_argument('--sweeps', '-s', type=int,
            help='enter the number of MC sweeps')
    args = parser.parse_args()

    # assign variables from cmd line input
    T = args.temp
    H = args.field
    L = args.length
    J = args.exchange
    s = args.sweeps

    # define Boltzmann constant
    k_B = 1.3806503 * pow(10,-23)

    #====================================================================
    # Create random numpy array of spins as either +/- 1
    #====================================================================

    aa = (2.0000*np.random.randint(2, size=pow(L,2))-1.0000).reshape(L,L)  
    
    SUM = energy(L,J,aa)
    spinSum = SUM

    # calculate energy of the system
    E = - J * spinSum
    print 'initial energy: ', E

    #====================================================================
    # Perform Monte Carlo Sweeps and use Metropolis algorithm
    # to either accept or reject spin flips.
    #====================================================================

    a = 0

    # create array to fill MC sweep number and total energy of system
    arr_a = np.array([])
    arr_E = np.array([])
    arr_a = np.append(arr_a, a)
    arr_E = np.append(arr_E, E)
    print 'array for a: ', arr_a, 'array for E: ', arr_E
    print 'first generated: '
    print aa

    while (a < s):
        # change the spin of one electron randomly chosen
        r = random.randint(0,L-1)
        c = random.randint(0,L-1)
        print 'r:', r, 'c:', c, 'element value: ', aa[c,r]
        if (aa[c,r] == 1):
            aa[c,r] = -1
        else:
            aa[c,r] = 1
        print 'before sampling: '
        print aa
        E_f = -J*aa[c,r]*(aa[c-1,r] + aa[c+1,r] + aa[c,r-1] + aa[c,r+1])
        Boltz = math.exp(-E_f / (k_B * T))
        if (E_f <= 0):
            aa = aa
            Reject = 0
            print 'accepted because negative'
            print 'after sampling: '
            print aa
        else:
            #generate number between 0 and 1, compare to Boltzmann factor
            n = random.random()
            print 'random number: ', n
            if (n <= Boltz):
                aa = aa
                Reject = 0
                print 'accepted by Boltz: ', Boltz
                print 'after sampling: '
                print aa
            else:
                if (aa[c,r] == 1):
                    aa[c,r] = -1
                else:
                    aa[c,r] = 1
                Reject = 1
                print 'rejected by Boltz: ', Boltz
                print 'after sampling: '
                print aa
        print 'neighbors: ', 'bottom: ', aa[c-1,r], 'top: ', aa[c+1,r]
        print 'neighbors: ', 'right: ', aa[c, r+1], 'left: ', aa[c,r-1]
        print 'flip E: ', E_f

        # recalculate energy of system  
        if (Reject == 0):
            aa = energy(L,J,aa)
            E = -J * spinSum
            print 'ENERGY HAS BEEN RECALCULATED TO: ', E
            print '=========================================='
        else:
            print 'ENERGY HAS BEEN LEFT THE SAME: ', E
            print '=========================================='
        a += 1
        arr_a = np.append(arr_a , a)
        arr_E = np.append(arr_E , E)


       #arr = arr.reshape(a+1,2)
    print 'array for a final: '
    print arr_a
    print 'array for E final: '
    print arr_E
    pylab.plot(arr_a , arr_E)
    pylab.show()
    return

#====================================================================
# Standard python boilerplate
#====================================================================
if __name__ == "__main__":
    main()


