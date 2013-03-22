# =============================================================================
# Script for submitting multiple jobs over range of barrier heights to plot
# energy vs. beta.  Run --help from command line for full description.
#
# Author:           Max Graves
# Last Revision:    9-DEC-2012
# =============================================================================

import subprocess
import os
import sys
import argparse
import pylab as pl
import glob

# =============================================================================
desc= ('PIGS double well energy determination script.  Submits jobs to '+
        'generate data to plot the energy versus beta=2V_max/b^2.  '+
        'This script will submit jobs over a range of barrier heights '+
        'and then move them to a new directory in ../output/dwEnergy/ '+
        'The companion script to this one will create the desired plot '+
        'from the data generated.')

parser = argparse.ArgumentParser(description=desc)
parser.add_argument("-i", "--HMin", type=float,
                default=0.5,
                help="Minimum value of barrier height")
parser.add_argument("-f", "--HMax", type=float,
                default=20.0,
                help="Maximum value of barrier height")
parser.add_argument("-s", "--step", type=float,
                default=0.5,
                help="increment to increase barrier each time")

args = parser.parse_args()

# =============================================================================
def findExe():
    """ checks for PIGS executable """
    
    if os.path.isfile('../PIGS.e'):
        print "PIGS.e executable was found"
        pass
    else:
        print 'Could not find PIGS.e\n'
        print 'Go into parent directory and compile code using Makefile'
        print '...then try again.'
        sys.exit("quitting\n")

# =============================================================================
def main():
    # make sure PIGS.e exists where it should
    findExe()

    # change to parent directory
    os.chdir("..")

    # create array of number of imaginary time step
    cR = pl.arange(args.HMin, (args.HMax+args.step), args.step)

    # submit the executable the correct amount of times
    for i in cR:
        command = ("./PIGS.e", "-M",str(15), "-s",str(3000000),
                "-E",str(1200000),"-c",str(1.0),"-H",str(i), "-p",str(2))
        subprocess.check_call(command)
        print "finished for height: ",i

    print "Now moving files to new directory..."

    os.chdir("./output")
    
    # create /output/dwEnergy subdirectory if it doesn't exist
    if not os.path.exists("./dwEnergy"):
        os.makedirs("./dwEnergy")

    # move all files created into this new directory we just created
    files = glob.glob("*avgBins*")
    for f in files:
        os.rename(f, ('./dwEnergy/'+str(f)))

    outputStr = ("\nThe dwEnergy_analyze.py script existing in this "+
            "directory will crunch all of this data into one numpy " + 
            "array and and plot.")
    
    print outputStr

# =============================================================================
if __name__=='__main__':
    main()
