import argparse
from utilities.argparse_types import path_inside_an_existing_dir, \
    existing_dir_path

def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates adriatic complete meshmask.nc
    and another file with _noRivers suffix in the name
    having tmask without rivers. It uses the river mask 
    included in the mask produced in the A1 phase
    ''')


    parser.add_argument(   '--outfile',"-o",
                                type = path_inside_an_existing_dir,
                                required = True,
                                help = '/some/path/meshmask.nc')

    parser.add_argument(   '--inputdir', '-i',
                                type = existing_dir_path,
                                default = None,
                                required = True,
                                help ='''Directory containg
                                          XC.data
                                          YC.data
                                          DXC.data
                                          DYC.data
                                          hFacC.data
                                '''
                                )

    parser.add_argument(   '--maskfile',"-m",
                                type = str,
                                default = None,
                                required = True,
                                help = "Input maskfile, generated in A1 phase"
                                )
    return parser.parse_args()

args = argument()



import numpy as np
import netCDF4 as NC
from commons import netcdf4

DIR = args.inputdir
tmask=netcdf4.readfile(args.maskfile, 'tmask')

RIVER_MASK_AVAILABLE = True
try:
    rivermask = netcdf4.readfile(args.maskfile, 'rivermask')
except KeyError:
    RIVER_MASK_AVAILABLE = False


CellBottoms = netcdf4.readfile(args.maskfile,'CellBottoms')
nav_lev = netcdf4.readfile(args.maskfile,'depth')
jpk, jpj, jpi = tmask.shape

delZ = np.zeros((jpk,) , np.float32)
delZ[0] = CellBottoms[0]
delZ[1:] = np.diff(CellBottoms)


xC=np.fromfile(DIR / "XC.data" ,dtype=np.float32,count=jpi*jpj).reshape(jpj,jpi)
yC=np.fromfile(DIR / "YC.data" ,dtype=np.float32,count=jpi*jpj).reshape(jpj,jpi)

e1t=np.fromfile(DIR / "DXC.data" ,dtype=np.float32,count=jpi*jpj).reshape(jpj,jpi)
e2t=np.fromfile(DIR / "DYC.data" ,dtype=np.float32,count=jpi*jpj).reshape(jpj,jpi)
e3t_fact=np.fromfile(DIR / "hFacC.data" ,dtype=np.float32,count=jpi*jpj*jpk).reshape(jpk,jpj,jpi)

e3t=np.zeros((1,jpk,jpj,jpi),np.float32)
for jk in range(jpk):
    e3t[0,jk,:,:] = delZ[jk] * e3t_fact[jk,:,:]

tmask=e3t>0



ncOUT=NC.Dataset(args.outfile,"w")

ncOUT.createDimension('x',jpi)
ncOUT.createDimension('y',jpj)
ncOUT.createDimension('z',jpk)
ncOUT.createDimension('time',1)

ncOUT.createDimension('x_a',1)
ncOUT.createDimension('y_a',1)
ncOUT.createDimension('z_a',1)

glamt  =xC
gphit  =yC
nav_lon=xC
nav_lat=yC

ncvar    = ncOUT.createVariable('e1t'   ,'d',('time','z_a', 'y', 'x')  ) ; ncvar[:] = e1t
ncvar    = ncOUT.createVariable('e2t'   ,'d',('time','z_a', 'y', 'x')  ) ; ncvar[:] = e2t
ncvar    = ncOUT.createVariable('e3t'   ,'d',('time',  'z', 'y', 'x'))  ;  ncvar[:] = e3t

ncvar    = ncOUT.createVariable('glamt'   ,'d',('time','z_a', 'y', 'x')) ; ncvar[:] = glamt
ncvar    = ncOUT.createVariable('gphit'   ,'d',('time','z_a', 'y', 'x')) ; ncvar[:] = gphit
ncvar    = ncOUT.createVariable('nav_lat','f',('y','x'))                 ; ncvar[:] = nav_lat
ncvar    = ncOUT.createVariable('nav_lev' ,'f',('z',))                   ; ncvar[:] = nav_lev
ncvar    = ncOUT.createVariable('nav_lon','f',('y','x'))                 ; ncvar[:] = nav_lon
#    float time(time) ;
#    short time_steps(time) ;
ncvar    = ncOUT.createVariable('tmask' ,'d',('time','z', 'y', 'x') )    ; ncvar[:] = tmask 
ncOUT.close()


def write_noRiver_mask(rivermask):
    output_dir = args.outfile.parent
    output_file_name = args.outfile.stem + '_noRivers.nc'
    output_file = output_dir / output_file_name
    with NC.Dataset(output_file,"w") as ncOUT:
        # create meshmask file without rivers in tmask
        tmask_noRivers = tmask - rivermask

        ncOUT.createDimension('x',jpi)
        ncOUT.createDimension('y',jpj)
        ncOUT.createDimension('z',jpk)
        ncOUT.createDimension('time',1)

        ncOUT.createDimension('x_a',1)
        ncOUT.createDimension('y_a',1)
        ncOUT.createDimension('z_a',1)

        glamt  = xC
        gphit  = yC
        nav_lon= xC
        nav_lat= yC

        ncvar    = ncOUT.createVariable('e1t'   ,'d',('time','z_a', 'y', 'x')  ) ; ncvar[:] = e1t
        ncvar    = ncOUT.createVariable('e2t'   ,'d',('time','z_a', 'y', 'x')  ) ; ncvar[:] = e2t
        ncvar    = ncOUT.createVariable('e3t'   ,'d',('time',  'z', 'y', 'x'))  ;  ncvar[:] = e3t

        ncvar    = ncOUT.createVariable('glamt'   ,'d',('time','z_a', 'y', 'x')) ; ncvar[:] = glamt
        ncvar    = ncOUT.createVariable('gphit'   ,'d',('time','z_a', 'y', 'x')) ; ncvar[:] = gphit
        ncvar    = ncOUT.createVariable('nav_lat','f',('y','x'))                 ; ncvar[:] = nav_lat
        ncvar    = ncOUT.createVariable('nav_lev' ,'f',('z',))                   ; ncvar[:] = nav_lev
        ncvar    = ncOUT.createVariable('nav_lon','f',('y','x'))                 ; ncvar[:] = nav_lon
        #    float time(time) ;
        #    short time_steps(time) ;
        ncvar    = ncOUT.createVariable('tmask' ,'d',('time','z', 'y', 'x') )    ; ncvar[:] = tmask_noRivers


if RIVER_MASK_AVAILABLE:
    write_noRiver_mask(rivermask)

