# OUTPUTS
# images for QUID
# CMEMS-Med-QUID-006-008-V2-V1.0.docx
# Figure IV.2

import pickle
import matplotlib.pyplot as pl
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import sys
import numpy as np
import argparse
from matplotlib.ticker import FormatStrFormatter

def argument():
    parser = argparse.ArgumentParser(description = '''
    plot somethings
    ''',
    formatter_class=argparse.RawTextHelpFormatter
    )

#    parser.add_argument(   '--outdir', '-O',
    parser.add_argument(   '--outdir', '-o',
                            type = str,
                            required =True,
                            default = "./",
                            help = ''' Output image dir'''
                            )
#    parser.add_argument(   '--open_sea_file', '-o',
    parser.add_argument(   '--open_sea_file', '-i',
#    parser.add_argument(   '--inputfile', '-i',
                            type = str,
                            required = True,
                            help = 'Input pickle file for open sea')
    parser.add_argument(   '--coast_file', '-c',
                            type = str,
                            required = False,
                            help = 'Input pickle file for coast')
    parser.add_argument(   '--var', '-v',
                                type = str,
                                required = True,
                                choices = ['chl','kd'],
                                help = ''' model var name'''
                                )



    return parser.parse_args()

args = argument()
#
fid = open(args.open_sea_file,'rb')
LIST = pickle.load(fid)
fid.close()
TIMES,_,_,MODEL_MEAN,SAT___MEAN,_,_,MODEL__STD,SAT____STD,CORR,NUMBERS = LIST

print (TIMES)

# MODEL MEAN and SAT MEAN are calculated in the previous step by "ScMYvalidation_plan_STD_CORR_valid_NAdr.py" and they have NAN at the same timesteps for definition

NAN_ARRAY=MODEL_MEAN.copy()
#NAN_ARRAY=SAT___MEAN.copy()


nNAN=np.sum(np.isnan(SAT___MEAN))

var_label = "CHL [mg/m$^3$]"

ii=~np.isnan(NAN_ARRAY)
#if (np.sum(np.isnan(SAT___MEAN)) > 0):
if (np.sum(~ii) > 0):
    MODEL_MEAN=MODEL_MEAN[ii]
    SAT___MEAN=SAT___MEAN[ii]
    MODEL__STD=MODEL__STD[ii]
    SAT____STD=SAT____STD[ii]
    CORR=CORR[ii]
    NUMBERS=NUMBERS[ii]

# Consider defined timeseps only: (discard Nan)
TIMES = [ TIMES[k] for k in range(len(TIMES))  if ii[k] ]


from basins import V2 as OGS
for isub,sub in enumerate(OGS.P):
  isub=0
#  print sub.name
  if (sub.name=="adr1"):
    print (sub.name)
    fig, ax = pl.subplots()
    ax.plot(TIMES,SAT___MEAN,'og',label=' SAT')
    ax.fill_between(TIMES,SAT___MEAN-SAT____STD,SAT___MEAN+SAT____STD,color='palegreen')
    ax.plot(TIMES,MODEL_MEAN,'-k',label=' RAN')
    ax.plot(TIMES,MODEL_MEAN-MODEL__STD,':k') #,label=' RAN')
    ax.plot(TIMES,MODEL_MEAN+MODEL__STD,':k') #,label=' RAN')

##    ax.fill_between(TIMES,SAT___MEAN[:,isub]-SAT____STD[:,isub],SAT___MEAN[:,isub]+SAT____STD[:,isub],color='palegreen')
##    ax.plot(TIMES,MODEL_MEAN[:,isub],'-k',label='RAN')
##    ax.plot(TIMES,MODEL_MEAN[:,isub]-MODEL__STD[:,isub],':k')
##    ax.plot(TIMES,MODEL_MEAN[:,isub]+MODEL__STD[:,isub],':k')




#    ax.plot(TIMES,model_coast[:,isub],':k',label=' RAN_coast')
#
#   ax.fill_between(TIMES,SAT___MEAN[:,isub]-SAT____STD[:,isub],SAT___MEAN[:,isub]+SAT____STD[:,isub],color='palegreen')
#    ax.plot(TIMES,MODEL_MEAN[:,isub],'-k',label=' RAN')
#    ax.plot(TIMES,MODEL_MEAN[:,isub]-MODEL__STD[:,isub],':k') #,label=' RAN')
#    ax.plot(TIMES,MODEL_MEAN[:,isub]+MODEL__STD[:,isub],':k') 
#
#    ax.set_ylabel(sub.name.upper() + ' - CHL [mg/m$^3$]').set_fontsize(14)
    ax.set_ylabel("%s - %s" %(sub.name.upper(), var_label  ) ).set_fontsize(14)

    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.legend(loc="best",labelspacing=0, handletextpad=0,borderpad=0.1)
    leg = pl.gca().get_legend()
    ltext  = leg.get_texts()
    pl.setp(ltext,fontsize=14)
    pl.rc('xtick', labelsize=14)
    pl.rc('ytick', labelsize=14)
    ax.tick_params(axis='both', labelsize=14)
    ylabels=ax.get_yticklabels()
    pl.setp(ylabels, fontsize=14)
##    pl.ylim(0.0, np.max(MODEL_MEAN[:,isub]+MODEL__STD[:,isub]) * 1.2 )
#    pl.ylim(0.0, max(np.max(SAT___MEAN+SAT____STD)  , np.max(MODEL_MEAN + MODEL__STD)) * 1.5 )
    if (args.open_sea_file == "export_data_ScMYValidation_plan_coast.pkl") | (args.open_sea_file == "export_data_ScMYValidation_plan_coast_V2.pkl" ):
        pl.ylim(0.0, 4.0)
        pl.ylim(0.0, 5.0)
    else:
        pl.ylim(0.0, 2.0)
        pl.ylim(0.0, 2.5)
#    pl.ylim(0.0, np.max(MODEL_MEAN[:,isub]+MODEL__STD[:,isub]) * 1.5 )
    ax.xaxis.set_major_locator(mdates.MonthLocator())
#    ax.xaxis.set_major_locator(mdates.WeekdayLocator())
#    ax.xaxis.set_major_locator(mdates.DayLocator(interval=10))
#    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%Y"))

#    loc = matplotlib.ticker.FixedLocator(matplotlib.dates.date2num(x) )

    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y%m%d"))
    ax.grid(True)
    xlabels = ax.get_xticklabels()
    pl.setp(xlabels, rotation=30, fontsize=11)
    outfilename_svg=args.outdir+"/"+'chl' + sub.name + "_STD.svg"
    outfilename_png=args.outdir+"/"+'chl' + sub.name + "_STD.png"
    pl.savefig(outfilename_svg,bbox_inches="tight")
    pl.savefig(outfilename_png,bbox_inches="tight")
