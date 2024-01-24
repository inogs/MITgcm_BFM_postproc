# OUTPUTS
# images for QUID
# CMEMS-Med-QUID-006-008-V2-V1.0.docx
# Figure IV.3 and table IV.1


import pickle
import matplotlib.pyplot as pl
import matplotlib.dates as mdates
import sys
import numpy as np
import argparse
import matplotlib.ticker as ticker

def argument():
    parser = argparse.ArgumentParser(description = '''
    plot somethings
    ''',
    formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(   '--outdir', '-o',
                            type = str,
                            required =True,
                            default = "./",
                            help = ''' Output image dir'''
                            )

    parser.add_argument(   '--inputfile', '-i',
                            type = str,
                            required = True,
                            default = 'export_data_ScMYValidation_plan.pkl',
                            help = 'Input pickle file')


#    parser.add_argument(   '--var', '-v',
#                                type = str,
#                                required = False,
#                                choices = ['chl','kd'],
#                                help = ''' model var name'''
#                                )


    return parser.parse_args()

args = argument()



fid = open(args.inputfile,'rb')
print (args.inputfile)
print (fid)
LIST = pickle.load(fid)
print (LIST)
#LIST = pickle.load(fid)
fid.close()

TIMES                          = LIST[0]
BGC_CLASS4_CHL_RMS_SURF_BASIN  = LIST[1]
BGC_CLASS4_CHL_BIAS_SURF_BASIN = LIST[2]
BGC_CLASS4_CHL_RMS_SURF_BASIN_LOG = LIST[5]
BGC_CLASS4_CHL_BIAS_SURF_BASIN_LOG= LIST[6]

MODEL_MEAN=LIST[3]
SAT___MEAN=LIST[4]
MODEL_STD = LIST[7]
SAT___STD = LIST[8]
BGC_CLASS4_CHL_CORR_SURF_BASIN= LIST[9]
BGC_CLASS4_CHL_POINTS_SURF_BASIN= LIST[10]

######### FILTER OUT NAN DAYS FROM SAT:

NAN_ARRAY=MODEL_MEAN
nNAN=np.sum(np.isnan(SAT___MEAN))

if (np.sum(np.isnan(SAT___MEAN)) > 0):
    MODEL_MEAN=MODEL_MEAN[~np.isnan(NAN_ARRAY)]
    SAT___MEAN=SAT___MEAN[~np.isnan(NAN_ARRAY)]
    MODEL_STD=MODEL_STD[~np.isnan(NAN_ARRAY)]
    SAT___STD=SAT___STD[~np.isnan(NAN_ARRAY)]
#    NUMBERS=NUMBERS[~np.isnan(NAN_ARRAY)]
 
    BGC_CLASS4_CHL_RMS_SURF_BASIN  = BGC_CLASS4_CHL_RMS_SURF_BASIN[~np.isnan(NAN_ARRAY)]
    BGC_CLASS4_CHL_BIAS_SURF_BASIN = BGC_CLASS4_CHL_BIAS_SURF_BASIN[~np.isnan(NAN_ARRAY)]
    BGC_CLASS4_CHL_RMS_SURF_BASIN_LOG = BGC_CLASS4_CHL_RMS_SURF_BASIN_LOG[~np.isnan(NAN_ARRAY)]
    BGC_CLASS4_CHL_BIAS_SURF_BASIN_LOG= BGC_CLASS4_CHL_BIAS_SURF_BASIN_LOG[~np.isnan(NAN_ARRAY)]
    BGC_CLASS4_CHL_CORR_SURF_BASIN= BGC_CLASS4_CHL_CORR_SURF_BASIN[~np.isnan(NAN_ARRAY)]
    BGC_CLASS4_CHL_POINTS_SURF_BASIN= BGC_CLASS4_CHL_POINTS_SURF_BASIN[~np.isnan(NAN_ARRAY)]


for ii in range(nNAN):
    TIMES.pop(np.argwhere(np.isnan(NAN_ARRAY))[0,0])


#########

from basins import V2 as OGS

nSUB = len(OGS.P.basin_list)
nSUB = 1

for isub,sub in enumerate(OGS.P):
  isub = 0
  if (sub.name=="adr1"):
#  if (isub != 17):
#  if (isub >20):
    print (sub.name)
    fig, ax = pl.subplots()
    ax.plot(TIMES,BGC_CLASS4_CHL_RMS_SURF_BASIN,'-k',label='RMS')
    ax.plot(TIMES,BGC_CLASS4_CHL_BIAS_SURF_BASIN,'-b',label='Bias')
    ax.set_ylabel(sub.name.upper() + ' - CHL [mg/m$^3$]').set_fontsize(14)
    ax.legend(loc="best",labelspacing=0, handletextpad=0,borderpad=0.1)
    leg = pl.gca().get_legend()
    ltext  = leg.get_texts()
    pl.setp(ltext,fontsize=14)
    pl.rc('xtick', labelsize=14)
    pl.rc('ytick', labelsize=14)
#    pl.ylim(-1.3, 1.3)
#    if (args.inputfile == "export_data_ScMYValidation_plan_coast.pkl"):
    if (args.inputfile == "export_data_ScMYValidation_plan_coast.pkl") | (args.inputfile == "export_data_ScMYValidation_plan_coast_V2.pkl" ):
         pl.ylim(-2.00, 5.00)
    else: 
         pl.ylim(-1.00, 2.50)
#    ax.xaxis.set_major_locator(mdates.MonthLocator())
#    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%Y"))
#    ax.xaxis.set_major_locator(mdates.DayLocator(interval=10))
    ax.xaxis.set_major_locator(mdates.MonthLocator())
#    ax.xaxis.set_major_locator(mdates.DayLocator())
#    ax.xaxis.set_major_locator(mdates.WeekdayLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y%m%d"))
    ax.grid(True)
    pl.axhline(y=0.0, color='k', linestyle='-',linewidth=1.5)
    xlabels = ax.get_xticklabels()
    pl.setp(xlabels, rotation=30,fontsize=11)
    ylabels = ax.get_yticklabels()
    pl.setp(ylabels, fontsize=14)
#    #ax.tick_params(direction='left', pad=2)
    #fig.show()
    outfilename_svg=args.outdir+"/"+'chl-RMS-BIAS_' + sub.name + ".svg"
    outfilename_png=args.outdir+"/"+'chl-RMS-BIAS_' + sub.name + ".png"
    print (outfilename_png)
    pl.savefig(outfilename_svg,bbox_inches="tight")
    pl.savefig(outfilename_png,bbox_inches="tight")
    pl.close(fig)

# Histogram for number of obs:
w = 0.9 # bar width
for isub,sub in enumerate(OGS.P):
  isub=0
  if (sub.name=="adr1"):
    print (sub.name)
    fig, ax = pl.subplots()
    ax.bar(TIMES,BGC_CLASS4_CHL_POINTS_SURF_BASIN[:],width=1.0, bottom=0, align="center", color="grey")
    ax.set_ylabel('# of points',fontsize=14)
    if (args.inputfile == "export_data_ScMYValidation_plan_coast.pkl"):
        pl.ylim(0, 10000)
    else:
        pl.ylim(0.0, 70000)
#    ax.xaxis.set_major_locator(mdates.DayLocator(interval=10)) #(byweekday=0, interval=2))
    ax.xaxis.set_major_locator(mdates.MonthLocator()) #(byweekday=0, interval=2))
   # ax.xaxis.set_major_locator(mdates.DayLocator()) 
#    ax.xaxis.set_major_locator(mdates.WeekdayLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y%m%d"))
    xlabels = ax.get_xticklabels()
    pl.setp(xlabels, rotation=30,fontsize=11)
    ax.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
    pl.ticklabel_format(axis='y', style='sci', scilimits=(3,3))
    ylabels = ax.get_yticklabels()
    pl.setp(ylabels, fontsize=14)
    outfile_svg= args.outdir + "/ValidPoints_SATchlsup_validation_" + sub.name + ".svg"
    outfile_png= args.outdir + "/ValidPoints_SATchlsup_validation_" + sub.name + ".png"
    print (outfile_png)
    fig.savefig(outfile_svg,bbox_inches="tight")
    fig.savefig(outfile_png,bbox_inches="tight")
    pl.close(fig)


import sys
sys.exit()

from commons.season import season
S=season()
S.setseasons(["0101", "0501", "0601", "1001"], ["winter","spring","summer","fall"])
from commons import timerequestors
from commons.Timelist import TimeInterval, TimeList
TL=TimeList(TIMES)
from commons.utils import writetable

iSeas=0 # JAN-APR
CLIM_REQ=timerequestors.Clim_season(iSeas,S)
ii,w=TL.select(CLIM_REQ)
RMS__win = np.nanmean(BGC_CLASS4_CHL_RMS_SURF_BASIN[     ii,:],axis=0)
BIAS_win = np.nanmean(BGC_CLASS4_CHL_BIAS_SURF_BASIN[    ii,:],axis=0)
RMSL_win = np.nanmean(BGC_CLASS4_CHL_RMS_SURF_BASIN_LOG[ ii,:],axis=0)
BIASLwin = np.nanmean(BGC_CLASS4_CHL_BIAS_SURF_BASIN_LOG[ii,:],axis=0)

MEAN_MOD_win = np.nanmean(MODEL_MEAN[ii,:],axis=0)
MEAN_REF_win = np.nanmean(SAT___MEAN[ii,:],axis=0)

STD_MOD_win = np.nanmean(MODEL_STD[ii,:],axis=0)
STD_REF_win = np.nanmean(SAT___STD[ii,:],axis=0)
CORR_win    = np.nanmean(BGC_CLASS4_CHL_CORR_SURF_BASIN[ii,:],axis=0)

iSeas=2 # JUN-SEP
CLIM_REQ=timerequestors.Clim_season(iSeas,S)
ii,w=TL.select(CLIM_REQ)
RMS__sum = np.nanmean(BGC_CLASS4_CHL_RMS_SURF_BASIN[     ii,:],axis=0)
BIAS_sum = np.nanmean(BGC_CLASS4_CHL_BIAS_SURF_BASIN[    ii,:],axis=0)
RMSL_sum = np.nanmean(BGC_CLASS4_CHL_RMS_SURF_BASIN_LOG[ ii,:],axis=0)
BIASLsum = np.nanmean(BGC_CLASS4_CHL_BIAS_SURF_BASIN_LOG[ii,:],axis=0)

MEAN_MOD_sum = np.nanmean(MODEL_MEAN[ii,:],axis=0)
MEAN_REF_sum = np.nanmean(SAT___MEAN[ii,:],axis=0)

STD_MOD_sum = np.nanmean(MODEL_STD[ii,:],axis=0)
STD_REF_sum = np.nanmean(SAT___STD[ii,:],axis=0)
CORR_sum    = np.nanmean(BGC_CLASS4_CHL_CORR_SURF_BASIN[ii,:],axis=0)

#mat = np.zeros((nSUB,8),np.float32)
#mat = np.zeros((nSUB,14),np.float32)
mat = np.zeros((nSUB,18),np.float32)

mat[:,0] = RMS__win
mat[:,1] = RMS__sum
mat[:,2] = BIAS_win
mat[:,3] = BIAS_sum
mat[:,4] = RMSL_win
mat[:,5] = RMSL_sum
mat[:,6] = BIASLwin
mat[:,7] = BIASLsum
mat[:,8] = STD_MOD_win
mat[:,9] = STD_REF_win
mat[:,10] = STD_MOD_sum
mat[:,11] = STD_REF_sum
mat[:,12] = CORR_win
mat[:,13] = CORR_sum
#----
mat[:,14] = MEAN_MOD_win
mat[:,15] = MEAN_REF_win
mat[:,16] = MEAN_MOD_sum
mat[:,17] = MEAN_REF_sum

#outfiletable = args.outdir+"/"+"table4.1.dat"
#rows_names=[sub.name for sub in OGS.P.basin_list]
##column_names=['RMSwin','RMSsum','BIASwin','BIASsum', 'RMSLwin','RMSLsum','BIASLwin','BIASLsum']
##column_names=['RMSwin','RMSsum','BIASwin','BIASsum', 'RMSLwin','RMSLsum','BIASLwin','BIASLsum','STD_MODwin','STD_SATwin','STD_MODsum','STD_SATsum','CORRwin','CORRsum']
#column_names=['RMSwin','RMSsum','BIASwin','BIASsum', 'RMSLwin','RMSLsum','BIASLwin','BIASLsum','STD_MODwin','STD_SATwin','STD_MODsum','STD_SATsum','CORRwin','CORRsum','MEAN_MODwin','MEAN_SATwin','MEAN_MODsum','MEAN_SATsum']
#writetable(outfiletable, mat, rows_names, column_names, fmt='%5.3f\t')

