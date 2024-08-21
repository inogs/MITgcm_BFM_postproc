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
    parser.add_argument(   '--pattern', '-g',
                                type = str,
                                required = False,
                                default = 'ave*nc',
                                help = 'glob search pattern')

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
    parser.add_argument(   '--mapconfig',"-f",
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
import numpy.ma as ma
from commons.mask import Mask
from utilities.mpi_serial_interface import get_mpi_communicator

from commons.dataextractor import DataExtractor
from layer_integral.mapbuilder import MapBuilder
import pylab as pl
from commons.layer import Layer

from map_gen import Map_object

from datetime import datetime
from commons.utils import addsep

try:
    import mpi4py
except:
    pass

comm = get_mpi_communicator()
rank  = comm.Get_rank()
nranks =comm.size





map_obj=Map_object.create_from_file(args.mapconfig)

INPUTDIR=addsep(args.inputdir)
OUTPUTDIR=addsep(args.outputdir)


TL = TimeList.fromfilenames(None, INPUTDIR, args.pattern )

dateformat="%Y%m%d-%H:%M:%S"



TheMask = Mask(args.maskfile)
jpk,jpj,jpi=TheMask.shape
plotlistfile=args.plotlistfile

mb = MapBuilder(plotlistfile,TL,TheMask,OUTPUTDIR)

for p in mb._MapBuilder__plotlist:
    var=p.varname
    LAYERLIST = p.layerlist

    ALL_INDEXES = np.arange(TL.nTimes)
    local_INDEXES = ALL_INDEXES[rank::nranks]

    if var == 'Speed':
        for iFrame_plot in local_INDEXES:
                dt=TL.Timelist[iFrame_plot]
                inputfileU = "%save.%s.%s.nc"  %(INPUTDIR,dt.strftime(dateformat),'U')
                inputfileV = "%save.%s.%s.nc"  %(INPUTDIR,dt.strftime(dateformat),'V')
                datestr=dt.strftime("%d %h %Y - %H:%M UTC")
                print("rank %d works on %s" %(rank, inputfileU), flush=True)
                DeU = DataExtractor(TheMask,inputfileU,'U')
                DeV = DataExtractor(TheMask,inputfileV,'V')

                for il, layer in enumerate(LAYERLIST):
                        layer = Layer(layer.top, layer.top)
                        U2d = MapBuilder.get_layer_average(DeU, layer)
                        V2d = MapBuilder.get_layer_average(DeV, layer)

                        clim = p.climlist[il] #(0.0, 0.4)
                        Udict={'varname':'U', 'longname':'Zonal velocity',      'clim':clim, 'layer':layer, 'data':U2d, 'date':datestr,'units':'m s⁻¹'}
                        Vdict={'varname':'V', 'longname':'Meridional velocity', 'clim':clim, 'layer':layer, 'data':V2d, 'date':datestr,'units':'m s⁻¹'}
                        fig,ax=map_obj.quiver_plotter_basemap_hourly(Udict, Vdict, TheMask)
                        outfile = "%save.%02d.%s.%s" % (OUTPUTDIR, ALL_INDEXES[iFrame_plot], p.varname, layer.string())
                        print(outfile,flush=True)
                        fig.savefig(outfile + ".png",dpi=86)
                        pl.close(fig)
    else:

        for iFrame_plot in local_INDEXES:
            dt=TL.Timelist[iFrame_plot]
            inputfile="%save.%s.%s.nc"  %(INPUTDIR,dt.strftime(dateformat),var)
            datestr=dt.strftime("%d %h %Y - %H:%M UTC")
            print("rank %d works on %s" %(rank, inputfile), flush=True)
            De = DataExtractor(TheMask,inputfile,var)

            for il, layer in enumerate(LAYERLIST):
                layer = Layer(layer.top, layer.top)
                map2d = MapBuilder.get_layer_average(De, layer)

                clim = p.climlist[il]
                mapdict={'varname':p.varname, 'longname':p.longvarname(), 'clim':clim, 'layer':layer, 'data':map2d, 'date':datestr,'units':p.units()}
                fig,ax=map_obj.map_plotter_basemap_hourly(mapdict, TheMask)
                outfile = "%save.%02d.%s.%s" % (OUTPUTDIR, ALL_INDEXES[iFrame_plot], p.varname, layer.string())
                print(outfile,flush=True)
                fig.savefig(outfile + ".png",dpi=86)
                pl.close(fig)
                #import sys
                #sys.exit()

