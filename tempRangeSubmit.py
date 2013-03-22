# =============================================================================
# Script for submitting multiple jobs over range of temperatures for Ising 
# model in order to simulate phase transition.
#
# Author:           Max Graves
# Last Revision:    21-MAR-2013
# =============================================================================

import subprocess, os, sys, argparse, glob
import pylab as pl

# =============================================================================
def parseCMD():
    """
    parse the command line.
    """
    desc= ('Ising 2D phase transition script.  Submits jobs to '+
            'generate data to plot macroscopic magnetization vs temp.  '+
            'This script will submit jobs over a range of temperatures '+
            'and then move them to a new directory in ../data/"uniqueID"/ '+
            'The companion script to this one will create the desired plot '+
            'from the data generated.')

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-i", "--TMin", type=float,
            default=0.5,
            help="Minimum value of temperature")
    parser.add_argument("-f", "--TMax", type=float,
            default=3.0,
            help="Maximum value of temperature")
    parser.add_argument("-e", "--step", type=float,
            default=0.1,
            help="increment to increase temp. each time")
    parser.add_argument("-L", "--Length", type=int,
            default=15,
            help="Number of Lattice Sites per Row")
    parser.add_argument("-J", "--Coupling", type=float,
            default=1.0,
            help="Coupling Constant, J")
    parser.add_argument("-s", "--mcSweeps", type=int,
            default=20000,
            help="Number of MC sweeps to perform")
       
    return parser.parse_args()

# =============================================================================
def findExe():
    """ 
    checks for Ising model code
    """
    if os.path.isfile('./2D_ising_MC.py'):
        print "Ising Driver was found"
        pass
    else:
        print 'Could not find Ising Driver\n'
        print 'Maybe you are in the wrong directory?'
        print '...find the file and try again.'
        sys.exit("quitting\n")

def checkDataDir():
    """
    checks whether the ./data/ directory has .dat files.
    Exits if it does.
    """
    inDir = glob.glob('./data/*dat*')
    if inDir != []:
        print '           ===== ERROR ====='
        print ' The ./data/ directory needs to have no .dat files.\n'
        print ' It may contain subdirectories with .dat files,'
        print ' but not within its top level.  Move into ./data/'
        print ' and move any .dat files into another directory.\n'
        sys.exit()

# =============================================================================
def main():

    # make sure PIGS.e exists where it should
    findExe()

    # check if the data directory has other data files
    checkDataDir()

    # parse cmd line
    args = parseCMD()

    # create array of number of imaginary time step
    temps = pl.arange(args.TMin, (args.TMax+args.step), args.step)

    # submit the executable the correct amount of times
    for temp in temps:
        command = ("python 2D_ising_MC.py"+" -J "+str(args.Coupling)+ " -s "+
                str(args.mcSweeps)+" -L "+str(args.Length)+ " -T "+
                str(temp))
        print command
        subprocess.check_call(command, shell=True)
        print "finished for temp: ",temp

    print "Now moving files to new directory..."

    # create /output/dwEnergy subdirectory if it doesn't exist
    dirName = ('./data/T%s-%s_L%s_s%s/'%(args.TMin,args.TMax,
        args.Length, args.mcSweeps))
    if not os.path.exists(dirName):
        os.mkdir(dirName)

    # move all files created into this new directory we just created
    files = glob.glob("./data/*dat*")
    curD = os.getcwd()
    print curD
    for f in files:
        fN = f[7:]
        #print dirName+str(f)
        os.rename(f, (dirName+str(fN)))

    outputStr = ("\nThe tempRangeAnalyze.py script existing in this "+
            "directory will crunch all of this data into one numpy " + 
            "array and and plot.")
    
    print outputStr

# =============================================================================
if __name__=='__main__':
    main()
