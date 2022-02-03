import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
   Creates products files from MIT chain, similar to the
  MEDSEA_ANALYSISFORECAST_BGC_006_014 ones.
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



FGROUPS = ['NUTR', 'PFTC', 'BIOL', 'CARB','CO2F']

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
    
    
    if FGroup == 'NUTR':
        setattr(ncOUT,'title','Nitrate, Phosphate, Ammonium and Silicate (3D) - Hourly Mean')
        
        ncvar = ncOUT.createVariable('no3', 'f', ('time','depth','latitude','longitude'),zlib=True, fill_value=1.0e+20)
        setattr(ncvar,'missing_value',ncvar._FillValue)
        setattr(ncvar,'units'        ,'mmol m-3')
        setattr(ncvar,'long_name'    ,'Mole concentration of Nitrate in sea water')
        setattr(ncvar,'standard_name','mole_concentration_of_nitrate_in_sea_water')
        setattr(ncvar,'coordinates'  ,'time depth latitude longitude')
        ncvar[:] = readdata(timestr, "N3n")
        
        
        
        ncvar = ncOUT.createVariable('po4', 'f', ('time','depth','latitude','longitude'),zlib=True, fill_value=1.0e+20)
        setattr(ncvar,'missing_value',ncvar._FillValue)
        setattr(ncvar,'units'        ,'mmol m-3')
        setattr(ncvar,'long_name'    ,'Mole concentration of Phosphate in sea water')
        setattr(ncvar,'standard_name','mole_concentration_of_phosphate_in_sea_water')
        setattr(ncvar,'coordinates'  ,'time depth latitude longitude')
        ncvar[:] = readdata(timestr, "N1p")


        ncvar = ncOUT.createVariable('nh4', 'f', ('time','depth','latitude','longitude'),zlib=True, fill_value=1.0e+20)
        setattr(ncvar,'missing_value',ncvar._FillValue)
        setattr(ncvar,'units'        ,'mmol m-3')
        setattr(ncvar,'long_name'    ,'Mole concentration of Ammonium in sea water')
        setattr(ncvar,'standard_name','mole_concentration_of_ammonium_in_sea_water')
        setattr(ncvar,'coordinates'  ,'time depth latitude longitude')
        ncvar[:] = readdata(timestr, "N4n")

        ncvar = ncOUT.createVariable('si', 'f', ('time','depth','latitude','longitude'),zlib=True, fill_value=1.0e+20)
        setattr(ncvar,'missing_value',ncvar._FillValue)
        setattr(ncvar,'units'        ,'mmol m-3')
        setattr(ncvar,'long_name'    ,'Mole concentration of Silicate in sea water')
        setattr(ncvar,'standard_name','mole_concentration_of_silicate_in_sea_water')
        setattr(ncvar,'coordinates'  ,'time depth latitude longitude')

        ncvar[:] = readdata(timestr, "N5s")

    if FGroup == 'PFTC':
        setattr(ncOUT,'title','Phytoplankton Carbon Biomass, Zooplankton Carbon Biomass and Chlorophyll (3D) - Hourly Mean')

        ncvar = ncOUT.createVariable('phyc', 'f', ('time','depth','latitude','longitude'),zlib=True, fill_value=1.0e+20)
        setattr(ncvar,'missing_value',ncvar._FillValue)
        setattr(ncvar,'units'        ,'mmol m-3')
        setattr(ncvar,'long_name'    ,'Concentration of Phytoplankton Biomass in sea water')
        setattr(ncvar,'standard_name','mole_concentration_of_phytoplankton_expressed_as_carbon_in_sea_water')
        setattr(ncvar,'coordinates'  ,'time depth latitude longitude')
        #CONVERSION from "mgC m-3" to "mmolC m-3"
        # conversion factor: 1/12
        ncvar[:] = readdata(timestr, 'P_c') * (1./12.)

        
        ncvar = ncOUT.createVariable('chl', 'f', ('time','depth','latitude','longitude'),zlib=True, fill_value=1.0e+20)
        setattr(ncvar,'missing_value',ncvar._FillValue)
        setattr(ncvar,'units'        ,'mg m-3')
        setattr(ncvar,'long_name'    ,'Concentration of Chlorophyll in sea water')
        setattr(ncvar,'standard_name','mass_concentration_of_chlorophyll_a_in_sea_water')
        setattr(ncvar,'coordinates'  ,'time depth latitude longitude')        
        ncvar[:] = readdata(timestr, "P_l")

        ncvar = ncOUT.createVariable('zooc', 'f', ('time','depth','latitude','longitude'),zlib=True, fill_value=1.0e+20)
        setattr(ncvar,'missing_value',ncvar._FillValue)
        setattr(ncvar,'units'        ,'mmol m-3')
        setattr(ncvar,'long_name'    ,'Concentration of Zooplankton Biomass in sea water')
        setattr(ncvar,'standard_name','mole_concentration_of_zooplankton_expressed_as_carbon_in_sea_water')
        setattr(ncvar,'coordinates'  ,'time depth latitude longitude')

        
        ncvar[:] = readdata(timestr, "Z_c")* (1./12.)


    if FGroup == 'BIOL':
        setattr(ncOUT, 'title', "Primary Production and Oxygen (3D) - Hourly Mean")
        
        ncvar = ncOUT.createVariable('o2', 'f', ('time','depth','latitude','longitude'),zlib=True, fill_value=1.0e+20)
        setattr(ncvar,'missing_value',ncvar._FillValue)
        setattr(ncvar,'units'        ,'mmol m-3')
        setattr(ncvar,'long_name'    ,'Mole concentration of Dissolved Molecular Oxygen in sea water')
        setattr(ncvar,'standard_name','mole_concentration_of_dissolved_molecular_oxygen_in_sea_water')
        setattr(ncvar,'coordinates'  ,'time depth latitude longitude')
        ncvar[:] = readdata(timestr,"O2o")
        
        ncvar = ncOUT.createVariable('nppv', 'f', ('time','depth','latitude','longitude'),zlib=True, fill_value=1.0e+20)
        setattr(ncvar,'missing_value',ncvar._FillValue)
        setattr(ncvar,'units'        ,'mg m-3 day-1')
        setattr(ncvar,'long_name'    ,'Net Primary Production in sea water')
        setattr(ncvar,'standard_name','net_primary_production_of_biomass_expressed_as_carbon_per_unit_volume_in_sea_water')
        setattr(ncvar,'coordinates'  ,'time depth latitude longitude')
        ncvar[:] = readdata(timestr,"ppn")

        
    if FGroup == 'CARB':
        setattr(ncOUT, 'title',"Dissolved Inorganic Carbon, pH and Alkalinity (3D) - Hourly Mean")

        ncvar = ncOUT.createVariable('ph', 'f', ('time','depth','latitude','longitude'),zlib=True, fill_value=1.0e+20)
        setattr(ncvar,'missing_value',ncvar._FillValue)
        setattr(ncvar,'units'        ,'1')
        setattr(ncvar,'long_name'    ,'PH')
        setattr(ncvar,'standard_name','sea_water_ph_reported_on_total_scale')
        setattr(ncvar,'info'         , 'pH reported on total scale at in situ Temp and Press conditions')
        setattr(ncvar,'coordinates'  ,'time depth latitude longitude')
        ncvar[:] = readdata(timestr, "pH")


        ncvar = ncOUT.createVariable('dissic', 'f', ('time','depth','latitude','longitude'),zlib=True, fill_value=1.0e+20)
        setattr(ncvar,'missing_value',ncvar._FillValue)
        setattr(ncvar,'units'        ,'mol m-3')
        setattr(ncvar,'long_name'    ,"Mole concentration of dissolved inorganic carbon in sea water")
        setattr(ncvar,'standard_name','mole_concentration_of_dissolved_inorganic_carbon_in_sea_water')
        setattr(ncvar,'coordinates'  ,'time depth latitude longitude')
        setattr(ncvar,'info'         , 'In order to calculate DIC in [micro mol / kg of seawater], dissic has to be multiplied by (1.e+6 / seawater density [kg/m3])')
        ncvar[:] = readdata(timestr, "O3c")/(12*1000) # conversion mg/mol


        ncvar = ncOUT.createVariable('talk', 'f', ('time','depth','latitude','longitude'),zlib=True, fill_value=1.0e+20)
        setattr(ncvar,'missing_value',ncvar._FillValue)
        setattr(ncvar,'units'        ,'mol m-3')
        setattr(ncvar,'long_name'    ,"Alkalinity")
        setattr(ncvar,'standard_name','sea_water_alkalinity_expressed_as_mole_equivalent')
        setattr(ncvar,'coordinates'  ,'time depth latitude longitude')
        setattr(ncvar,'info'         , 'In order to calculate ALK in [micro mol / kg of seawater], talk has to be multiplied by (1.e+6 / seawater density [kg/m3])')
        ncvar[:] = readdata(timestr, "O3h")/1000 # conversion mg/mol
        


    if FGroup == 'CO2F':
        setattr(ncOUT, 'title',"Surface partial pressure of CO2 and Surface CO2 flux (2D) - Hourly Mean")
        ncvar = ncOUT.createVariable('fpco2', 'f', ('time','latitude','longitude'),zlib=True, fill_value=1.0e+20)
        setattr(ncvar,'missing_value',ncvar._FillValue)
        setattr(ncvar,'units'        ,'kg m-2 s-1')
        setattr(ncvar,'long_name'    ,"surface downward flux at air-sea interface of carbon dioxide expressed as kg of carbon per square meter per second")
        setattr(ncvar,'standard_name','surface_downward_mass_flux_of_carbon_dioxide_expressed_as_carbon')
        setattr(ncvar,'coordinates'  ,'time latitude longitude')
        ncvar[:] = readdata(timestr, "CO2airflux", ndims=2) *12 * 1.e-6 /86400 # conversion from mmol m-2 day-1 to kg/m2/s
        

        ncvar = ncOUT.createVariable('spco2', 'f', ('time','latitude','longitude'),zlib=True, fill_value=1.0e+20)
        setattr(ncvar,'missing_value',ncvar._FillValue)
        setattr(ncvar,'units'        ,'Pa')
        setattr(ncvar,'long_name'    ,'Surface partial pressure of carbon dioxide in sea water')
        setattr(ncvar,'standard_name','surface_partial_pressure_of_carbon_dioxide_in_sea_water')
        setattr(ncvar,'coordinates'  ,'time latitude longitude')
        ncvar[:] = readdata(timestr, "pCO2",ndims=2) *0.101325 #conversion microatm --> Pascal  1 ppm = 1 microatm = 1.e-6 * 101325 Pa

    ncOUT.close()
        
