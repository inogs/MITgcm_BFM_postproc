import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates Netcdf from MIT .data hourly files
    ''')


    parser.add_argument(   '--outputdir',"-o",
                                type = str,
                                required = True,
                                help = '/some/path/')

    parser.add_argument(   '--inputdir', '-i',
                                type = str,
                                default = None,
                                required = True,
                                help ='''Directory containg outputs of MITgcm
                                '''
                                )
    parser.add_argument(   '--timelist',
                                type = str,
                                default = None,
                                required = True,
                                help='Filename with times in %Y%m%d-%H:%M:%S format'
                                )

    parser.add_argument(   '--timestep',"-t",
                                type = str,
                                default = None,
                                required = True,
                                help = "timestep of MIT simulation in seconds."
                                )

    parser.add_argument(   '--maskfile',"-m",
                                type = str,
                                default = None,
                                required = True
                                )
    parser.add_argument(   '--varlist',"-v",
                                type = str,
                                default = None,
                                required = True
                                )
    parser.add_argument(   '--freq',"-f",
                                type = str,
                                choices = ['daily','hourly'],
                                required = True
                                )

    return parser.parse_args()

args = argument()

from commons import genUserDateList as DL 
from commons import netcdf4
from commons.mask import Mask
import numpy as np
from datetime import datetime
from commons.utils import addsep, file2stringlist



try:
    from mpi4py import MPI
    comm  = MPI.COMM_WORLD
    rank  = comm.Get_rank()
    nranks =comm.size
    isParallel = True
except:
    rank   = 0
    nranks = 1
    isParallel = False

def readFrame_from_file(filename,Frame,shape):
    jpk,jpj,jpi=shape
    domain_size=jpi*jpj*jpk
    fid=open(filename,'rb')
    fid.seek(domain_size*Frame*4)
    A=np.fromfile(fid,dtype=np.float32,count=domain_size)
    fid.close()
    return A.reshape(jpk,jpj,jpi)


INPUTDIR=addsep(args.inputdir)
OUTDIR=addsep(args.outputdir)
dateformat="%Y%m%d-%H:%M:%S"


TheMask = Mask(args.maskfile)
VARLIST= file2stringlist(args.varlist)
timelist=file2stringlist(args.timelist)

timestep = int(args.timestep)

TimeSteps_in_h = 3600/timestep
ALL_INDEXES = np.arange((len(timelist)))
local_INDEXES = ALL_INDEXES[rank::nranks]


for var in VARLIST:    
    for it in local_INDEXES:
        t = timelist[it]
        if args.freq =='hourly': inputfile = "%s%s.%010d.data" %(INPUTDIR,var, (it+1)*TimeSteps_in_h)
        if args.freq =='daily' : inputfile = "%s%s.%010d.data" %(INPUTDIR,var, (it+1)*TimeSteps_in_h*24)
        outfile   = "%save.%s.%s.nc"  %(OUTDIR,t,var)
        print(outfile)
        M3d = readFrame_from_file(inputfile, 0, TheMask.shape)
        M3d[~TheMask.mask] = 1.e+20
        netcdf4.write_3d_file(M3d, var, outfile, TheMask, compression=True)

