import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Calculates Chl statistics on matchups of Model and Sat.
    Main result are
    - BGC_CLASS4_CHL_RMS_SURF_BASIN
    - BGC_CLASS4_CHL_BIAS_SURF_BASIN
    of deliverable CMEMS-Med-biogeochemistry-ScCP-1.0.pdf

    Other similar results are also calculated and dumped in
    the a pickle file provided as argument.
    ''',
    formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(   '--satdir', '-s',
                                type = str,
                                required =False,
                                default = "/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/CCI/MONTHLY_V4/",
                                help = ''' Input satellite dir'''
                                )
    parser.add_argument(   '--inputmodeldir', '-i',
                                type = str,
                                required =True,
                                help = ''' Input model dir, where P_l files are, usually ../wrkdir/POSTPROC/output/AVE_FREQ_2/TMP/'''
                                )
    parser.add_argument(   '--outfile', '-o',
                                type = str,
                                required = True,
                                default = 'export_data_ScMYValidation_plan.pkl',
                                help = 'Output pickle file')
    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                required = True,
                                help = 'Path of the mask file')
    parser.add_argument(   '--coastness', '-c',
                                type = str,
                                required = True,
                                choices = ['coast','open_sea','everywhere'],
                                help = 'definition of mask to apply to the statistics')
 
    parser.add_argument(   '--Timestart', '-t_start',
                                type = str,
                                required = True,
                                help = 'Start time of the timeseriesi, eg. "20220224"')

    parser.add_argument(   '--Time__end', '-t_end',
                                type = str,
                                required = True,
                                help = 'End time of the timeseries, eg. "20230623"')


    return parser.parse_args()


args = argument()
from commons.Timelist import TimeList
from commons.time_interval import TimeInterval
import numpy as np
import os
import Sat.SatManager as Sat
import matchup.matchup as matchup
from commons.dataextractor import DataExtractor
from layer_integral.mapbuilder import MapBuilder
from commons.mask import Mask
from commons.submask import SubMask
from basins import V2 as OGS
from commons.layer import Layer
from commons.utils import addsep
import pickle
#from profiler import DATESTART, DATE__END
import netCDF4

def weighted_mean(Conc, Weight):

    Weight_sum      = np.nansum(Weight) #.sum()
    Mass            = np.nansum(Conc * Weight) #.sum()
    Weighted_Mean   = Mass/Weight_sum
#    Weighted_Std    = np.sqrt(((Conc - Weighted_Mean)**2*Weight).sum()/Weight_sum)
    Weighted_Std    = np.sqrt(np.nansum((Conc - Weighted_Mean)**2*Weight)/Weight_sum)
    return Weighted_Mean, Weighted_Std

#MASKFILE="/g100_work/OGS_prodC/MIT/V1M-dev/V1/devel/wrkdir/POSTPROC/meshmask.nc"
#TheMask=Mask(MASKFILE)
TheMask=Mask(args.maskfile)

Sup_mask = TheMask.cut_at_level(0)
MODEL_DIR= addsep(args.inputmodeldir)
REF_DIR  = addsep(args.satdir)
outfile  = args.outfile
#outfile = "export_data_ScMYValidation_plan_everywhere.pkl"

Timestart=args.Timestart #"20220224"
Time__end=args.Time__end #"20230623"

MASKFILE="/g100_work/OGS_prodC/MIT/V1/prod/wrkdir/POSTPROC/meshmask.nc"
TheMask=Mask(MASKFILE)

TI    = TimeInterval(Timestart,Time__end,"%Y%m%d")
print (TI) 
dateformat ="%Y%m%d"
#dateformat ="%Y%m_d"
sat_TL   = TimeList.fromfilenames(TI, REF_DIR  ,"*.nc", prefix="", dateformat=dateformat)
#model_TL = TimeList.fromfilenames(TI, MODEL_DIR,"*P_l.nc")
#model_TL = TimeList.fromfilenames(TI, MODEL_DIR,"*PFTC*.nc",prefix="",dateformat="%Y%m%d",forceFrequency=None)
model_TL = TimeList.fromfilenames(TI, MODEL_DIR,"20*",prefix="",dateformat="%Y%m%d",forceFrequency=None)
print (sat_TL.Timelist)  
print (" " )
print (model_TL.Timelist)
suffix = os.path.basename(sat_TL.filelist[0])[8:]
suffix_model = "_h-OGS--PFTC-MITgcmBFM-pilot8-b" 


nFrames = model_TL.nTimes
#nFrames = sat_TL.nTimes
#nSUB = len(OGS.P.basin_list)

jpk,jpj,jpi =TheMask.shape
#dtype = [(sub.name, np.bool) for sub in OGS.P]
#SUB = np.zeros((jpj,jpi),dtype=dtype)
#for sub in OGS.P:
#    print sub.name
#    sbmask         = SubMask(sub,maskobject=Sup_mask).mask
#    SUB[sub.name]  = sbmask[0,:,:]

mask20_2D = TheMask.mask_at_level(20.0)
mask0_2D = TheMask.mask_at_level(0.0)
coastmask = mask0_2D

if args.coastness == 'coast':
    coastmask=mask0_2D & (~mask20_2D)
if args.coastness == "open_sea"  : coastmask = mask20_2D
if args.coastness == "everywhere": coastmask = mask0_2D

nSUB=1
BGC_CLASS4_CHL_RMS_SURF_BASIN      = np.zeros((nFrames,nSUB),np.float32)
BGC_CLASS4_CHL_BIAS_SURF_BASIN     = np.zeros((nFrames,nSUB),np.float32)
MODEL_MEAN                         = np.zeros((nFrames,nSUB),np.float32)
SAT___MEAN                         = np.zeros((nFrames,nSUB),np.float32)
MODEL__STD                         = np.zeros((nFrames,nSUB),np.float32)
SAT____STD                         = np.zeros((nFrames,nSUB),np.float32)
BGC_CLASS4_CHL_CORR_SURF_BASIN     = np.zeros((nFrames,nSUB),np.float32)
BGC_CLASS4_CHL_RMS_SURF_BASIN_LOG  = np.zeros((nFrames,nSUB),np.float32)
BGC_CLASS4_CHL_BIAS_SURF_BASIN_LOG = np.zeros((nFrames,nSUB),np.float32)
BGC_CLASS4_CHL_POINTS_SURF_BASIN   = np.zeros((nFrames,nSUB),np.float32)

# This is the surface layer choosen to match satellite chl data
surf_layer = Layer(0,10)

for itime, modeltime in enumerate(model_TL.Timelist):
#for itime, modeltime in enumerate(sat_TL.Timelist):
    print (modeltime)
    CoupledList = sat_TL.couple_with([modeltime])
#    CoupledList = model_TL.couple_with([modeltime])
    try:
        sattime = CoupledList[0][0]
        print (sattime)
    except: 
        print ("except")
        continue
    satfile = REF_DIR + sattime.strftime(dateformat) + suffix
#    modfile = model_TL.filelist[itime]
#    modfile = MODEL_DIR + sattime.strftime(dateformat) + suffix_model + sattime.strftime(dateformat) + "_fc-v01.nc" 
    modfile = MODEL_DIR + sattime.strftime(dateformat) + "/" + sattime.strftime(dateformat) + suffix_model + sattime.strftime(dateformat) + "_fc-v01.nc"
    

#    De         = DataExtractor(TheMask,filename=modfile, varname='P_l')
#    Model      = MapBuilder.get_layer_average(De, surf_layer)
    dset = netCDF4.Dataset(modfile)
    dset.variables['chl'].shape
    mod12h=np.array(dset.variables['chl'])[12,:,:,:]
    Model=mod12h[0,:,:]

#Model      = MapBuilder.get_layer_average(De, surf_layer) 

    #ncIN = NC.netcdf_file(modfile,'r')
    #Model = ncIN.variables['P_i'].data[0,0,:,:].copy()#.astype(np.float64)
    #Model = ncIN.variables['lchlm'].data.copy()
    #ncIN.close()

    Sat16 = Sat.readfromfile(satfile,var='CHL') #.astype(np.float64)


    cloudsLand = (np.isnan(Sat16)) | (Sat16 > 1.e19) | (Sat16<0)
    modelLand  = np.isnan(Model) #lands are nan
#    modelLand  = np.isnan(mod12h)
    nodata     = cloudsLand | modelLand
    selection = ~nodata & coastmask
    M = matchup.matchup(Model[selection], Sat16[selection])

#    for isub, sub in enumerate(OGS.P):
    isub=0
    if (1==1):
#        selection = SUB[sub.name] & (~nodata) & coastmask
        M = matchup.matchup(Model[selection], Sat16[selection])
        
        BGC_CLASS4_CHL_POINTS_SURF_BASIN[itime,isub]  = M.number()
        BGC_CLASS4_CHL_RMS_SURF_BASIN[itime,isub]  = M.RMSE()
        BGC_CLASS4_CHL_BIAS_SURF_BASIN[itime,isub] = M.bias()
        BGC_CLASS4_CHL_CORR_SURF_BASIN[itime,isub] = M.correlation()
        weight = TheMask.area[selection]
        MODEL_MEAN[itime,isub] , MODEL__STD[itime,isub] = weighted_mean( M.Model,weight)
        SAT___MEAN[itime,isub] , SAT____STD[itime,isub] = weighted_mean( M.Ref,  weight)

        Mlog = matchup.matchup(np.log10(Model[selection]), np.log10(Sat16[selection])) #add matchup based on logarithm
        BGC_CLASS4_CHL_RMS_SURF_BASIN_LOG[itime,isub]  = Mlog.RMSE()
        BGC_CLASS4_CHL_BIAS_SURF_BASIN_LOG[itime,isub] = Mlog.bias()

BGC_CLASS4_CHL_EAN_RMS_SURF_BASIN  = BGC_CLASS4_CHL_RMS_SURF_BASIN.mean(axis=0)
BGC_CLASS4_CHL_EAN_BIAS_SURF_BASIN = BGC_CLASS4_CHL_BIAS_SURF_BASIN.mean(axis=0)


LIST   =[i for i in range(11)]
LIST[0]=model_TL.Timelist
#LIST[0]=sat_TL.Timelist
LIST[1]=BGC_CLASS4_CHL_RMS_SURF_BASIN
LIST[2]=BGC_CLASS4_CHL_BIAS_SURF_BASIN
LIST[3]=MODEL_MEAN
LIST[4]=SAT___MEAN
LIST[5]=BGC_CLASS4_CHL_RMS_SURF_BASIN_LOG
LIST[6]=BGC_CLASS4_CHL_BIAS_SURF_BASIN_LOG
LIST[7]=MODEL__STD
LIST[8]=SAT____STD
LIST[9]=BGC_CLASS4_CHL_CORR_SURF_BASIN
LIST[10]=BGC_CLASS4_CHL_POINTS_SURF_BASIN

fid = open(outfile,'wb')
pickle.dump(LIST, fid)
fid.close()
