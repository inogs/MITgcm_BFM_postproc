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
# MODEL MEAN and SAT MEAN are calculated in the previous step by "ScMYvalidation_plan_STD_CORR_valid_NAdr.py" and they have NAN at the same timesteps for definition

NAN_ARRAY=MODEL_MEAN.copy()
nNAN=np.sum(np.isnan(SAT___MEAN))

ii= ~np.isnan(NAN_ARRAY)
if (np.sum(~ii) > 0):
    MODEL_MEAN=MODEL_MEAN[ii]
    SAT___MEAN=SAT___MEAN[ii]
    MODEL_STD=MODEL_STD[ii]
    SAT___STD=SAT___STD[ii]
 
    BGC_CLASS4_CHL_RMS_SURF_BASIN  = BGC_CLASS4_CHL_RMS_SURF_BASIN[ii]
    BGC_CLASS4_CHL_BIAS_SURF_BASIN = BGC_CLASS4_CHL_BIAS_SURF_BASIN[ii]
    BGC_CLASS4_CHL_RMS_SURF_BASIN_LOG = BGC_CLASS4_CHL_RMS_SURF_BASIN_LOG[ii]
    BGC_CLASS4_CHL_BIAS_SURF_BASIN_LOG= BGC_CLASS4_CHL_BIAS_SURF_BASIN_LOG[ii]
    BGC_CLASS4_CHL_CORR_SURF_BASIN= BGC_CLASS4_CHL_CORR_SURF_BASIN[ii]
    BGC_CLASS4_CHL_POINTS_SURF_BASIN= BGC_CLASS4_CHL_POINTS_SURF_BASIN[ii]


# Consider defined timeseps only: (discard Nan)
TIMES = [ TIMES[k] for k in range(len(TIMES))  if ii[k] ]

#########

from basins import V2 as OGS

nSUB = len(OGS.P.basin_list)
nSUB = 1

for isub,sub in enumerate(OGS.P):
  isub = 0
  if (sub.name=="adr1"):
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
    if (args.inputfile == "export_data_ScMYValidation_plan_coast.pkl") | (args.inputfile == "export_data_ScMYValidation_plan_coast_V2.pkl" ):
         pl.ylim(-2.00, 5.00)
    else: 
         pl.ylim(-1.00, 2.50)
#    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%Y"))
#    ax.xaxis.set_major_locator(mdates.DayLocator(interval=10))
    ax.xaxis.set_major_locator(mdates.MonthLocator())
#    ax.xaxis.set_major_locator(mdates.WeekdayLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y%m%d"))
    ax.grid(True)
    pl.axhline(y=0.0, color='k', linestyle='-',linewidth=1.5)
    xlabels = ax.get_xticklabels()
    pl.setp(xlabels, rotation=30,fontsize=11)
    ylabels = ax.get_yticklabels()
    pl.setp(ylabels, fontsize=14)
#    #ax.tick_params(direction='left', pad=2)
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
    ax.xaxis.set_major_locator(mdates.MonthLocator()) #(byweekday=0, interval=2))
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

