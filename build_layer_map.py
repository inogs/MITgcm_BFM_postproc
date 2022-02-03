import argparse



def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates adriatic medeaf maps
    ''')


    parser.add_argument(   '--outputdir',"-o",
                                type = str,
                                required = True,
                                help = '/some/path/')

    parser.add_argument(   '--inputdir', '-i',
                                type = str,
                                default = None,
                                required = True,
                                help ='''Directory containg merged outputs of MITgcm
                                '''
                                )

    parser.add_argument(   '--rundate',"-d",
                                type = str,
                                default = None,
                                required = True
                                )

    parser.add_argument(   '--maskfile',"-m",
                                type = str,
                                default = None,
                                required = True
                                )
    parser.add_argument(   '--plotlistfile',"-p",
                                type = str,
                                default = None,
                                required = True
                                )
    return parser.parse_args()

args = argument()


import matplotlib
matplotlib.use('Agg')
from commons.Timelist import TimeList
from commons import genUserDateList as DL
import numpy as np
from commons.mask import Mask

from commons.dataextractor import DataExtractor
from layer_integral.mapbuilder import MapBuilder
import pylab as pl
from commons.layer import Layer
from map_gen import map_plotter_basemap_hourly
from datetime import datetime
from commons.utils import addsep

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


RUNDATE=args.rundate

INPUTDIR=addsep(args.inputdir)
OUTPUTDIR=addsep(args.outputdir)

rundate_dt = datetime.strptime(RUNDATE,"%Y%m%d")
datestart_plot = rundate_dt.strftime("%Y%m%d-%H:%M:%S")
dateend        = (rundate_dt+DL.relativedelta(hours=71)).strftime("%Y%m%d-%H:%M:%S")

dateformat="%Y%m%d-%H:%M:%S"


plot_timelist=DL.getTimeList(datestart_plot, dateend, hours=1)
PTL= TimeList(plot_timelist)

TheMask = Mask(args.maskfile)
jpk,jpj,jpi=TheMask.shape
plotlistfile=args.plotlistfile

mb = MapBuilder(plotlistfile,PTL,TheMask,OUTPUTDIR)

for p in mb._MapBuilder__plotlist:
    var=p.varname
    LAYERLIST = p.layerlist

    ALL_INDEXES = np.arange((len(plot_timelist)))
    local_INDEXES = ALL_INDEXES[rank::nranks]

    for iFrame_plot in local_INDEXES:
        dt=plot_timelist[iFrame_plot]
        inputfile="%save.%s.%s.nc"  %(INPUTDIR,dt.strftime(dateformat),var)
        datestr=dt.strftime("%d %h %Y - %H:%M UTC")
        print("rank %d works on %s" %(rank, inputfile))
        De = DataExtractor(TheMask,inputfile,var)

        for il, layer in enumerate(LAYERLIST):
            layer = Layer(layer.top, layer.top)
            map2d = MapBuilder.get_layer_average(De, layer)
            
            clim = p.climlist[il]
            mapdict={'varname':p.varname, 'longname':p.longvarname(), 'clim':clim, 'layer':layer, 'data':map2d, 'date':datestr,'units':p.units()}
            fig,ax=map_plotter_basemap_hourly(mapdict, TheMask)
            outfile = "%save.%02d.%s.%s" % (OUTPUTDIR, ALL_INDEXES[iFrame_plot], p.varname, layer.string())
            print(outfile)
            fig.savefig(outfile + ".png",dpi=86)
            pl.close(fig)
            #import sys
            #sys.exit()

