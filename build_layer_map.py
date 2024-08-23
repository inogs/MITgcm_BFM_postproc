import argparse
from utilities.argparse_types import date_from_str, existing_dir_path, \
    existing_file_path


def argument():
    parser = argparse.ArgumentParser(
        description='Generates adriatic medeaf maps'
    )

    parser.add_argument(   '--outputdir',"-o",
                                type = existing_dir_path,
                                required = True,
                                help = '/some/path/')

    parser.add_argument(   '--inputdir', '-i',
                                type = existing_dir_path,
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

    parser.add_argument(   '--maskfile', "-m",
                                type = existing_file_path,
                                default = None,
                                required = True
                                )
    parser.add_argument(   '--plotlistfile', "-p",
                                type = existing_file_path,
                                default = None,
                                required = True
                                )
    parser.add_argument(   '--mapconfig', "-f",
                                type = existing_file_path,
                                default = None,
                                required = True
                                )
    parser.add_argument(   '--ignore-before', "-b",
                                type = date_from_str,
                                default = None,
                                required = False,
                                help = 'Ignore all the frames before this date'
                                )

    return parser.parse_args()

args = argument()


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from commons.Timelist import TimeList
from commons.mask import Mask
from utilities.mpi_serial_interface import get_mpi_communicator

from commons.dataextractor import DataExtractor
from layer_integral.mapbuilder import MapBuilder
from commons.layer import Layer

from map_gen import Map_object

try:
    import mpi4py
except:
    pass

comm = get_mpi_communicator()
rank = comm.Get_rank()
nranks = comm.size

map_obj = Map_object.create_from_file(args.mapconfig)

TL = TimeList.fromfilenames(None, args.inputdir, args.pattern)

if args.ignore_before is not None:
    valid_frames = [f for f in TL.Timelist if f >= args.ignore_before]
else:
    valid_frames = [f for f in TL.Timelist]

dateformat = "%Y%m%d-%H:%M:%S"

TheMask = Mask(args.maskfile)
jpk, jpj, jpi = TheMask.shape
plotlistfile = args.plotlistfile

mb = MapBuilder(plotlistfile, TL, TheMask, args.outputdir)

for p in mb._MapBuilder__plotlist:
    LAYERLIST = p.layerlist

    # The indices and the dt values of the frames of this process
    local_frames = tuple(enumerate(valid_frames))[rank::nranks]

    if p.varname == 'Speed':
        for frame_index, dt in local_frames:
            input_file_name_U = 'ave.{}.{}.nc'.format(dt.strftime(dateformat), 'U')
            input_file_name_V = 'ave.{}.{}.nc'.format(dt.strftime(dateformat), 'V')
            inputfileU = args.inputdir / input_file_name_U
            inputfileV = args.inputdir / input_file_name_V
            datestr = dt.strftime("%d %h %Y - %H:%M UTC")
            print("rank %d works on %s" %(rank, inputfileU), flush=True)
            DeU = DataExtractor(TheMask,inputfileU, 'U')
            DeV = DataExtractor(TheMask,inputfileV, 'V')

            for il, layer in enumerate(LAYERLIST):
                layer = Layer(layer.top, layer.top)
                U2d = MapBuilder.get_layer_average(DeU, layer)
                V2d = MapBuilder.get_layer_average(DeV, layer)

                clim = p.climlist[il] #(0.0, 0.4)
                Udict={'varname':'U', 'longname':'Zonal velocity',      'clim':clim, 'layer':layer, 'data':U2d, 'date':datestr,'units':'m s⁻¹'}
                Vdict={'varname':'V', 'longname':'Meridional velocity', 'clim':clim, 'layer':layer, 'data':V2d, 'date':datestr,'units':'m s⁻¹'}
                fig, ax=map_obj.quiver_plotter_basemap_hourly(Udict, Vdict, TheMask)
                outfile_name = "ave.{:0>2}.{}.{}.png".format(frame_index, p.varname, layer.string())
                outfile = args.outputdir / outfile_name
                print(outfile, flush=True)
                fig.savefig(outfile, dpi=86)
                plt.close()
    else:
        for frame_index, dt in local_frames:
            if args.ignore_before is not None and dt < args.ignore_before:
                continue
            inputfile = "ave.{}.{}.nc".format(dt.strftime(dateformat), p.varname)
            inputfile_path = args.inputdir / inputfile
            datestr = dt.strftime("%d %h %Y - %H:%M UTC")
            print("rank %d works on %s" %(rank, inputfile_path), flush=True)
            De = DataExtractor(TheMask, inputfile_path, p.varname)

            for il, layer in enumerate(LAYERLIST):
                layer = Layer(layer.top, layer.top)
                map2d = MapBuilder.get_layer_average(De, layer)

                clim = p.climlist[il]
                mapdict={'varname':p.varname, 'longname':p.longvarname(), 'clim':clim, 'layer':layer, 'data':map2d, 'date':datestr,'units':p.units()}
                fig,ax=map_obj.map_plotter_basemap_hourly(mapdict, TheMask)
                outfile_name = "ave.{:0>2}.{}.{}.png".format(frame_index, p.varname, layer.string())
                outfile = args.outputdir / outfile_name
                print(outfile, flush=True)
                fig.savefig(outfile, dpi=86)
                plt.close()

