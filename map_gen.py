import pylab as pl
from mpl_toolkits.basemap import Basemap
import numpy as np

xlim=[12.1,16.1]
ylim=[43.5,45.85]
xC=(xlim[0]+xlim[1])/2
yC=(ylim[0]+ylim[1])/2
map_obj = Basemap(projection='merc',lat_0=xC,lon_0=yC,\
                                 llcrnrlon = xlim[0], \
                                 llcrnrlat = ylim[0], \
                                 urcrnrlon = xlim[1], \
                                 urcrnrlat = ylim[1], \
                                 resolution='h')

def map_plotter_basemap(mapdict,maskobj):
    '''
    Integrates coastlines
    '''
    background_color=(.8, .8, .8)
    fig,ax = pl.subplots()
    ax.set_facecolor(background_color)
    vmin,vmax=mapdict['clim']
    cmap=pl.get_cmap('jet')
    #cmap=pl.get_cmap('gist_rainbow_r')
    #cmap.set_bad(color='w',alpha=1.)
    
    # draw coastlines, country boundaries, fill continents.
    map_obj.drawcoastlines(linewidth=0.5)
    #map.drawcountries(linewidth=0.25)
    map2d=mapdict['data']
    Zm = np.ma.masked_invalid(map2d)
    map_obj.pcolormesh(maskobj.xlevels, maskobj.ylevels, Zm,cmap=cmap,latlon='true',vmin=vmin,vmax=vmax)

    ax.annotate(mapdict['date'][5:] ,xy=(0.9,0.9), xycoords='axes fraction' , fontsize=28);
    parallels = np.arange(43.,46.,.5)
    map_obj.drawparallels(parallels,labels=[1,0,0,0],fontsize=16, dashes=[6,900])
    # draw meridians
    meridians = np.arange(12.,20,1.)
    map_obj.drawmeridians(meridians,labels=[0,0,0,1],fontsize=16,dashes=[6,900])
    #cbar = map_obj.colorbar(cs,location='right',size="5%",pad="2%")    
    #cbar.ax.tick_params(labelsize=5)

     
    fig.set_size_inches(10,10*map_obj.ymax/map_obj.xmax)
    ax.set_position([.10, .08, .88, .88])
    fig.patch.set_alpha(0.0)
    return fig,ax

def map_plotter_basemap_hourly(mapdict,maskobj):
    '''
    Arguments:
    * mapdict  * a dictionary defined as 
                 mapdict={'clim':[vmin, vmax], 'data': 2D numpy ndarray, 'date':string}

    * maskobj  * a Mask object

    Integrates coastlines
    Returns :
    *  fig  * figure handle
    *   ax  * axis handle    
    '''
    ncolors=24
    background_color=(.8, .8, .8)
    fig,ax = pl.subplots()
    ax.set_facecolor(background_color)
    vmin,vmax=mapdict['clim']
    cmap=pl.get_cmap('jet',ncolors)
    #cmap=pl.get_cmap('gist_rainbow_r')
    #cmap.set_bad(color='w',alpha=1.)
    
    # draw coastlines, country boundaries, fill continents.
    map_obj.drawcoastlines(linewidth=0.5)
    #map.drawcountries(linewidth=0.25)
    map2d=mapdict['data']
    Zm = np.ma.masked_invalid(map2d)
    cs=map_obj.pcolormesh(maskobj.xlevels, maskobj.ylevels, Zm,cmap=cmap,latlon='true',vmin=vmin,vmax=vmax)

    layerstr=mapdict['layer'].__repr__().lower()
    ax.annotate(mapdict['date'] ,xy=(0.60,0.93), xycoords='axes fraction' , fontsize=16)
    ax.annotate(layerstr       , xy=(0.77,0.85), xycoords='axes fraction' , ha='center', fontsize=16)
    
    parallels = np.arange(43.,46.,.5)
    map_obj.drawparallels(parallels,labels=[1,0,0,0],fontsize=13, dashes=[6,900])
    # draw meridians
    meridians = np.arange(12.,20,1.)
    map_obj.drawmeridians(meridians,labels=[0,0,0,1],fontsize=13,dashes=[6,900])
    
    nticks = 6
    cbar_ticks_list = np.linspace(vmin, vmax, nticks).tolist()
    cbar_ticks_labels = ["%g" % (t,)  for t in cbar_ticks_list ]

    
    cbar = map_obj.colorbar(cs,location='right',size="5%",pad="2%", ticks=cbar_ticks_list)
    cbar.ax.set_yticklabels(cbar_ticks_labels)  
    cbar.ax.tick_params(labelsize=13)

    title = "%s [ %s ] "  % (mapdict['longname'],mapdict['units'])
    ax.set_title(title, fontsize=16)
    ax.title.set_position((0.5,1.03))

    

    
    #fig.set_size_inches(10,10*map_obj.ymax/map_obj.xmax)
    #ax.set_position([.08, .08, .83, .83])
    fig.set_size_inches(10,10*(map_obj.ymax/map_obj.xmax)* (.83/.93))
    ax.set_position([.08, .03, .83, .93])
    
    fig.patch.set_alpha(0.0)
    return fig,ax

def quiver_basemap_hourly(U,V,maskobj):
    # andrebbe riportata la stessa struttura della map_plotter_basemap_hourly
    
    fig,ax=pl.subplots()
    ax.quiver(maskobj.xlevels, maskobj.ylevels,U,V)
    return fig,ax