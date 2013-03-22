#==============================================================================
# This classical monte carlo code will simulate 2-D Ising model
#
# Author:               Max Graves
# Last Modified:        21-MAR-2013
#==============================================================================

import random, math
import numpy as np
import argparse

#==============================================================================
# parse cmd line 
#==============================================================================
def parseCMD():
    parser = argparse.ArgumentParser(description = 'simulate 2-D Ising model')
    parser.add_argument('--temp', '-T', type=float,
            help='enter temperature in kelvin')
    parser.add_argument('--field', '-H', type=float,
            help='enter magnetic field strength')
    parser.add_argument('--length', '-L', type=float,
            help='enter the length of the 2-d lattice')
    parser.add_argument('--exchange', '-J', type=float, default=1,
            help='enter the value of the exchange constant')
    parser.add_argument('--sweeps', '-s', type=int,
            help='enter the number of MC sweeps')
    return parser.parse_args()

#==============================================================================
# Calculate initial energy and magnetization of system, summing 
# the products of spins of nearest neighbor lattice sites 
#==============================================================================
def calcEnergy(aa,L,J):
    """
    Sum products of neighboring spin sites.  Return
    the energy as well as the macroscopic magnetism.
    """
    # vertical neighbors
    w, h, totalSpin, spinSum = 0, 0, 0.0, 0.0
    while (h <= L-1):
        totalSpin += aa[h,w]
        if ((h == 0) and (w==0)):
            print 'left top row'
            spinSum += (aa[0,1] + aa[0,-1] + aa[1,0] + aa[-1,0])*aa[0,0]
            print 'w: ', w, 'h: ', h, 'spinSum: ', spinSum
            w += 1
        elif ((h==0) and (w<L-1)):
            print 'middle top row'
            spinSum += (aa[1,w] + aa[-1,w] + aa[0,w-1] + aa[0,w+1])*aa[0,w]
            print 'w: ', w, 'h: ', h, 'spinSum: ', spinSum
            w += 1
        elif ((h==0) and (w==L-1)):
            print 'right top row'
            spinSum += (aa[1,w] + aa[-1,w] + aa[0,0] + aa[0,w-1])*aa[0,w]
            print 'w: ', w, 'h: ', h, 'spinSum: ', spinSum
            w = 0
            h += 1
        elif ((h<L-1) and (w==0)):
            print 'left middle'
            spinSum += (aa[h+1,w] + aa[h-1,w] + aa[h,1] + aa[h,-1])*aa[h,0]
            print 'w: ', w, 'h: ', h, 'spinSum: ', spinSum
            w += 1
        elif ((h<L-1) and (w<L-1)):
            print 'middle'
            spinSum += (aa[h+1,w] + aa[h-1,w] + aa[h,w+1] + aa[h,w-1])*aa[h,w]
            print 'w: ', w, 'h: ', h, 'spinSum: ', spinSum
            w += 1
        elif ((h<L-1) and (w==L-1)):
            print 'right middle'
            spinSum += (aa[h+1,w] + aa[h-1,w] + aa[h,0] + aa[h,w-1])*aa[h,w]
            print 'w: ', w, 'h: ', h, 'spinSum: ', spinSum
            w = 0
            h += 1
        elif ((h==L-1) and (w==0)):
            print 'left bottom'
            spinSum += (aa[0,w] + aa[h-1,w] + aa[h,1] + aa[h,-1])*aa[h,w]
            print 'w: ', w, 'h: ', h, 'spinSum: ', spinSum
            w += 1
        elif ((h==L-1) and (w<L-1)):
            print 'middle bottom row'
            spinSum += (aa[0,w] + aa[h-1,w] + aa[h,w+1] + aa[h,w-1])*aa[h,w]
            print 'w: ', w, 'h: ', h, 'spinSum: ', spinSum
            w += 1
        elif ((h==L-1) and (w==L-1)):
            print 'right bottom'
            spinSum += (aa[0,w] + aa[h-1,w] + aa[h,0] + aa[h,w-1])*aa[h,w]
            print 'w: ', w, 'h: ', h, 'spinSum: ', spinSum
            break
     
    # calculate energy of the system
    E = - 0.5 * J * spinSum     # divide by two because of double counting
    M =  math.fabs(totalSpin / pow(L,2))
    return E,M

# =============================================================================
# begin main function
# =============================================================================
def main():

    # assign variables from cmd line input
    args = parseCMD()
    T,H,L,J,s = args.temp, args.field, args.length, args.exchange, args.sweeps

    k_B = 1     # = 1.3806503 * pow(10,-23)

    # Create random numpy array of spins as either +/- 1
    latt = (2.0000*np.random.randint(2, size=pow(L,2))-1.0000).reshape(L,L)

    print latt

    E,M = calcEnergy(latt, L, J)

    print 'E: ',E
    print 'M: ',M
    
    # Perform Monte Carlo Sweeps and use Metropolis algorithm

    # define MC step number to run over (THIS NEEDS TO BE SHIFTED)
    #a = 0
    # create array to fill MC sweep number and total energy of system
    #arr_a = np.array([a])
    #arr_E = np.array([E])
    #arr_M = np.array([M])
    #arr_a = np.append(arr_a, a)
    #arr_E = np.append(arr_E, E)
    #arr_M = np.append(arr_M, M)
    #print 'array for a: ', arr_a, 'array for E: ', arr_E
    #print 'first generated: '
    #print aa

    #while (a < s):
    #    # change the spin of one electron randomly chosen
    #    # 'r' tells you how many columns to count over
    #    r = random.randint(0,L)
    #    # 'c' tells you how many rows to count down 
    #    c = random.randint(0,L)
    #    #print 'r:', r, 'c:', c
    #    if (aa[c,r] == 1):
    #        aa[c,r] = -1
    #    else:
    #        aa[c,r] = 1
    #    #print 'before sampling: '
    #    #print aa 
    #    #E_f = -J*aa[c,r]*(aa[c-1,r] + aa[c+1,r] + aa[c,r-1] + aa[c,r+1])
    #    Boltz = math.exp(-E_f / (k_B * T))
    #    if (E_f < 0):
    #        aa = aa
    #        Reject = 0
    #        #print 'accepted because negative'
    #        #print 'after sampling: '
    #        #print aa
    #    elif (E_f == 0):
    #        aa = aa
    #        Reject = 0
    #        #print 'accepted because zero, no change in energy'
    #        #print 'after sampling: '
    #        #print aa
    #    else:
    #        #generate number between 0 and 1, compare to Boltzmann factor
    #        n = random.random()
    #        #print 'random number: ', n
    #        if (n <= Boltz):
    #            aa = aa
    #            Reject = 0
    #            #print 'accepted by Boltz: ', Boltz
    #            #print 'after sampling: '
    #            #print aa
    #        else:
    #            if (aa[c,r] == 1):
    #                aa[c,r] = -1
    #            else:
    #                aa[c,r] = 1
    #            Reject = 1
    #            #print 'rejected by Boltz: ', Boltz
    #            #print 'after sampling: '
    #            #print aa
    #
    #
    #print 'neighbors: ', 'bottom: ', aa[c-1,r], 'top: ', aa[c+1,r]
    #    #print 'neighbors: ', 'right: ', aa[c, r+1], 'left: ', aa[c,r-1]
    #    #print 'flip E: ', E_f
    #
    #    # recalculate energy and magnetization of spin system  
    #    # sum products of vertical neighboring spin sites of lattice
    #    if (Reject == 0):
    
    #    else:
    #        continue
    #        #print 'ENERGY HAS BEEN LEFT THE SAME: ', E
    #        #print 'MAGNETIZATION WAS LEFT SAME:   ', M
    #    a += 1
    #    arr_a = np.append(arr_a , a)
    #    arr_E = np.append(arr_E , E)
    #    arr_M = np.append(arr_M , M)
    #
    # output results to a file
    #filename = "isingOUT.dat"
    #File = open(filename, "w")
    #File.writelines(arr_a)
    #File.close()

# =============================================================================
if __name__=="__main__":
    main()

