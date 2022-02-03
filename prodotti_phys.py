import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
   Creates products files from MIT chain, similar to the
  MEDSEA_ANALYSISFORECAST_PHY_006_013 ones.
   Standard names are choose from
   http://cfconventions.org/Data/cf-standard-names/30/build/cf-standard-name-table.html.

   Files have been checked from http://puma.nerc.ac.uk/cgi-bin/cf-checker.pl.

   Parallel executable, can be called by mpirun.
   ''',formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(   '--inputdir', '-i',
                                type = str,
                                required = True,
                                help ='The directory wrkdir/MODEL/AVE_FREQ_1/ where chain has run.'
                                )

    parser.add_argument(   '--outputdir',"-o",
                                type = str,
                                required = True,
                                help = 'Path of existing dir')

    parser.add_argument(    '--time',"-t", 
                                type = str,
                                required = True,
                                help = '''Path of input text file with the yyyymmdd list''' )
    parser.add_argument(    '--DType',"-d", 
                                type = str,
                                required = True,
                                help = '''Analysis, simulation , or forecast''',
                                choices = ["an","sm","fc"])
     
    parser.add_argument(    '--bulltime',"-b",
                                type = str,
                                required = True,
                                help = '''The bulletin time a string time in the format yyyymmdd ''')
    parser.add_argument(    '--maskfile', "-m",
                                type = str,
                                required = True,
                                help = '''Path for the maskfile ''')

    return parser.parse_args()


args = argument()

import netCDF4
import numpy as np
import datetime,os
from commons.utils import addsep, file2stringlist
from commons.mask import Mask
from commons.dataextractor import DataExtractor

try:
    from mpi4py import MPI
    comm  = MPI.COMM_WORLD
    rank  = comm.Get_rank()
    nranks =comm.size
except:
    rank   = 0
    nranks = 1



INPUTDIR  = addsep(args.inputdir)
OUTPUTDIR = addsep(args.outputdir)
TIMELIST  = file2stringlist(args.time)
DType     = args.DType
bulletin_date = args.bulltime
maskfile = args.maskfile



TheMask = Mask(maskfile,ylevelsmatvar="gphit", xlevelsmatvar="glamt")
jpk, jpj, jpi = TheMask.shape
nav_lev = TheMask.zlevels
Lon = TheMask.xlevels[0,:].astype(np.float32)
Lat = TheMask.ylevels[:,0].astype(np.float32)
tmask = TheMask.mask



FGROUPS = ['RFVL', 'PSAL', 'TEMP']

if DType == "an": bulletin_type='analysis'
if DType == "sm": bulletin_type='simulation'
if DType == "fc": bulletin_type='forecast'

bulletin_time = datetime.datetime.strptime(bulletin_date,"%Y%m%d")


def readdata(time, var, ndims=3):
    if ndims==3:
        M = np.zeros((24,jpk,jpj,jpi),np.float32)
    else:
        M = np.zeros((24,jpj,jpi), np.float32)
    for iFrame in range(24):
        inputfile = "%save.%s-%02d:00:00.%s.nc" %(INPUTDIR,time,iFrame,var)
        print(inputfile)
        M[iFrame,:] = DataExtractor(TheMask,inputfile,var,dimvar=ndims).values
    return M

def create_Structure(filename, fgroup):
    ref=  'http://medeaf.inogs.it/adriatic'
    inst  ='National Institute of Oceanography and Applied Geophysics - OGS, Italy'
    ncOUT = netCDF4.Dataset(filename,"w",format="NETCDF4")
    ncOUT.createDimension('longitude', jpi)
    ncOUT.createDimension('latitude' ,jpj)
    if (fgroup != 'CO2F') : ncOUT.createDimension('depth'    ,jpk)
    ncOUT.createDimension('time'     ,  24)
    
    setattr(ncOUT,'Conventions'  ,'CF-1.0' )
    setattr(ncOUT,'references'   , ref     )
    setattr(ncOUT,'institution'  , inst    )
    setattr(ncOUT,'source'       , 'MITgcmBFM')
    setattr(ncOUT,'comment'      , ref)
    setattr(ncOUT,'contact'      ,'http://medeaf.inogs.it/contact')
    setattr(ncOUT,'bulletin_date', bulletin_time.strftime("%Y-%m-%d") )
    setattr(ncOUT,'bulletin_type', bulletin_type)
    setattr(ncOUT,'field_type'   , "hourly_mean_beginning_at_time_field")
    
    basename = os.path.basename(filename)
    timestr = basename[:8]
    D = datetime.datetime.strptime(timestr,'%Y%m%d')
    Dref = datetime.datetime(1970,1,1,0,0,0)
    Diff = D-Dref
    
    ncvar = ncOUT.createVariable('time','d',('time',))
    setattr(ncvar,'units',       'seconds since 1970-01-01 00:00:00')
    setattr(ncvar,'long_name'    ,'time')
    setattr(ncvar,'standard_name','time')
    setattr(ncvar,'axis'         ,'T')
    setattr(ncvar,'calendar'     ,'standard')
    ncvar[:] = Diff.days*86400 + np.arange(24)*3600
    
    if (fgroup != 'CO2F') :
        ncvar = ncOUT.createVariable('depth'   ,'f', ('depth',))
        setattr(ncvar,'units'        ,'m')
        setattr(ncvar,'long_name'    ,'depth')
        setattr(ncvar,'standard_name','depth')
        setattr(ncvar,'positive'     ,'down')
        setattr(ncvar,'axis'         ,'Z')
        setattr(ncvar,'valid_min'    ,nav_lev.min())
        setattr(ncvar,'valid_max'    ,nav_lev.max())
        ncvar[:] = nav_lev
    
    ncvar = ncOUT.createVariable('latitude','f' ,('latitude',))
    setattr(ncvar, 'units'        ,'degrees_north')
    setattr(ncvar,'long_name'    ,'latitude')
    setattr(ncvar,'standard_name','latitude')
    setattr(ncvar, 'axis'         ,'Y')
    setattr(ncvar,'valid_min'    , Lat.min())
    setattr(ncvar, 'valid_max'    ,Lat.max())
    ncvar[:]=Lat

    ncvar = ncOUT.createVariable('longitude','f',('longitude',))
    setattr(ncvar, 'units'        ,'degrees_east')
    setattr(ncvar,'long_name'    ,'longitude')
    setattr(ncvar, 'standard_name','longitude')
    setattr(ncvar, 'axis'         ,'X')
    setattr(ncvar, 'valid_min'    , Lon.min())
    setattr(ncvar, 'valid_max'    , Lon.max())
    ncvar[:]=Lon
    
    return ncOUT


def V7_filename(timeobj,FGroup):
    return timeobj.strftime('%Y%m%d_') + "h-OGS--" + FGroup + "-MITgcmBFM-pilot8-b" + bulletin_date +"_" + DType + "-v01.nc"

nGroups = len(FGROUPS)
PROCESSES = np.arange(len(TIMELIST)*nGroups)

for ip in PROCESSES[rank::nranks]:
    (iTime,iFgroup) = divmod(ip,nGroups)
    timestr = TIMELIST[iTime]
    timeobj = datetime.datetime.strptime(timestr,"%Y%m%d")
    FGroup = FGROUPS[iFgroup]
    
    product_file = V7_filename(timeobj, FGroup)
    print("rank =", rank, product_file)
    ncOUT = create_Structure(OUTPUTDIR + product_file,FGroup)
    
    
    if FGroup == 'RFVL':
        setattr(ncOUT,'title','Horizontal Velocity (3D) - Hourly Mean')
        
        ncvar = ncOUT.createVariable('uo', 'f', ('time','depth','latitude','longitude'),zlib=True, fill_value=1.0e+20)
        setattr(ncvar,'missing_value',ncvar._FillValue)
        setattr(ncvar,'units'        ,'"m s-1')
        setattr(ncvar,'long_name'    ,'eastward ocean current velocity')
        setattr(ncvar,'standard_name','eastward_sea_water_velocity')
        setattr(ncvar,'coordinates'  ,'time depth latitude longitude')
        ncvar[:] = readdata(timestr, "U")
        
        
        
        ncvar = ncOUT.createVariable('vo', 'f', ('time','depth','latitude','longitude'),zlib=True, fill_value=1.0e+20)
        setattr(ncvar,'missing_value',ncvar._FillValue)
        setattr(ncvar,'units'        ,'"m s-1')
        setattr(ncvar,'long_name'    ,'northward ocean current velocity')
        setattr(ncvar,'standard_name','northward_sea_water_velocity')
        setattr(ncvar,'coordinates'  ,'time depth latitude longitude')
        ncvar[:] = readdata(timestr, "V")


    if FGroup == 'PSAL':
        setattr(ncOUT,'title','Salinity (3D) - Daily Mean')

        ncvar = ncOUT.createVariable('so', 'f', ('time','depth','latitude','longitude'),zlib=True, fill_value=1.0e+20)
        setattr(ncvar,'missing_value',ncvar._FillValue)
        setattr(ncvar,'units'        ,'0.001')
        setattr(ncvar,'long_name'    ,'salinity')
        setattr(ncvar,'standard_name','sea_water_salinity')
        setattr(ncvar,'coordinates'  ,'time depth latitude longitude')

        ncvar[:] = readdata(timestr, 'S')



    if FGroup == 'TEMP':
        setattr(ncOUT, 'title', "Potential Temperature (3D) - Daily Mean")
        
        ncvar = ncOUT.createVariable('thetao', 'f', ('time','depth','latitude','longitude'),zlib=True, fill_value=1.0e+20)
        setattr(ncvar,'missing_value',ncvar._FillValue)
        setattr(ncvar,'units'        ,'degrees_C')
        setattr(ncvar,'long_name'    ,'potential temperature')
        setattr(ncvar,'standard_name','sea_water_potential_temperature')
        setattr(ncvar,'coordinates'  ,'time depth latitude longitude')
        ncvar[:] = readdata(timestr,"T")


    ncOUT.close()
        
