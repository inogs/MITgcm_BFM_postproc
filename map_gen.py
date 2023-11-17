import pylab as pl
from mpl_toolkits.basemap import Basemap
import numpy as np
import cmocean as cmo

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
    cs=map_obj.pcolormesh(maskobj.xlevels, maskobj.ylevels, Zm,cmap=cmap,latlon='true',shading='nearest',vmin=vmin,vmax=vmax)

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

def quiver_plotter_basemap_hourly(mapdictU, mapdictV, maskobj, zoomGoT = False):
    '''
    Arguments:
    * mapdict{U/V}  * dictionaries for zonal/meridional velocities defined as
                 mapdict{*}={'clim':[vmin, vmax], 'data': 2D numpy ndarray, 'date':string}

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
    vmin,vmax=mapdictU['clim']
    cmap=cmo.cm.speed

    if zoomGoT == True:
        spaceSteps = int(1)
        '''
        xlim=[13.348,13.840]
        ylim=[45.467,45.850]
        xC=(xlim[0]+xlim[1])/2
        yC=(ylim[0]+ylim[1])/2
        map_obj = Basemap(projection='merc',lat_0=xC,lon_0=yC,\
                                     llcrnrlon = xlim[0], \
                                     llcrnrlat = ylim[0], \
                                     urcrnrlon = xlim[1], \
                                     urcrnrlat = ylim[1], \
                                     resolution='h')
        idxLats = [np.argmin(np.abs(yl - maskobj.ylevels[:,0])) for yl in ylim]
        idxLons = [np.argmin(np.abs(xl - maskobj.xlevels[0,:])) for xl in xlim]

        map_obj.drawcoastlines(linewidth=1)
        #map.drawcountries(linewidth=0.25)
        U2d=mapdictU['data']
        V2d=mapdictV['data']
        Um = np.ma.masked_invalid(U2d)[idxLats[0]:idxLats[1], idxLons[0]:idxLons[1]]
        Vm = np.ma.masked_invalid(V2d)[idxLats[0]:idxLats[1], idxLons[0]:idxLons[1]]
        Zm = np.ma.sqrt(Um**2 + Vm**2)
        lon = maskobj.xlevels[idxLats[0]:idxLats[1], idxLons[0]:idxLons[1]]
        lat = maskobj.ylevels[idxLats[0]:idxLats[1], idxLons[0]:idxLons[1]]

        #cs=map_obj.quiver(maskobj.xlevels[::6, ::6], maskobj.ylevels[::6, ::6], Um[::6, ::6], Vm[::6, ::6], Zm[::6, ::6], cmap=cmap, latlon='true', clim = mapdictU['clim'])
        #cs = map_obj.pcolormesh(lon, lat, Zm, cmap=cmap,latlon='true',shading='nearest', clim = mapdictU['clim'])
        cs = map_obj.pcolormesh(lon, lat, Zm, cmap = cmap, latlon='true',shading='nearest', clim = mapdictU['clim'])
        ColorMatrix = 1. - 1./(np.exp((Zm - 22/24*vmax)/0.001) + 1)
        #qp = map_obj.quiver(lon[::spaceSteps, ::spaceSteps], lat[::spaceSteps, ::spaceSteps], Um[::spaceSteps, ::spaceSteps], Vm[::spaceSteps, ::spaceSteps], latlon='true')
        qp = map_obj.quiver(lon[::spaceSteps, ::spaceSteps], lat[::spaceSteps, ::spaceSteps], Um[::spaceSteps, ::spaceSteps], Vm[::spaceSteps, ::spaceSteps], ColorMatrix[::spaceSteps, ::spaceSteps], latlon='true', cmap = 'Greys_r', clim = mapdictU['clim'])

        layerstr=mapdictU['layer'].__repr__().lower()
        ax.annotate(mapdictU['date'] ,xy=(0.60,0.93), xycoords='axes fraction' , fontsize=16)
        ax.annotate(layerstr       , xy=(0.77,0.85), xycoords='axes fraction' , ha='center', fontsize=16)

        parallels = np.arange(45.5, 45.9, .1)
        map_obj.drawparallels(parallels,labels=[1,0,0,0],fontsize=13, dashes=[6,900])
        # draw meridians
        meridians = np.arange(13.4, 14.0, 0.2)
        map_obj.drawmeridians(meridians,labels=[0,0,0,1],fontsize=13,dashes=[6,900])'''
    else:
        spaceSteps = int(6)

        # draw coastlines, country boundaries, fill continents.
        map_obj.drawcoastlines(linewidth=0.5)
        #map.drawcountries(linewidth=0.25)
        U2d=mapdictU['data']
        V2d=mapdictV['data']
        Um = np.ma.masked_invalid(U2d)
        Vm = np.ma.masked_invalid(V2d)
        Zm = np.ma.sqrt(Um**2 + Vm**2)
        #cs=map_obj.quiver(maskobj.xlevels[::6, ::6], maskobj.ylevels[::6, ::6], Um[::6, ::6], Vm[::6, ::6], Zm[::6, ::6], cmap=cmap, latlon='true', clim = mapdictU['clim'])
        cs=map_obj.pcolormesh(maskobj.xlevels, maskobj.ylevels, Zm, cmap=cmap,latlon='true',shading='nearest', clim = mapdictU['clim'])
        #ColorMatrix = 1. - 1./(np.exp((Zm - 22/24*vmax)/0.001) + 1)
        qp = map_obj.quiver(maskobj.xlevels[::spaceSteps, ::spaceSteps], maskobj.ylevels[::spaceSteps, ::spaceSteps], Um[::spaceSteps, ::spaceSteps], Vm[::spaceSteps, ::spaceSteps], latlon='true') #, ColorMatrix[::6, ::6], cmap = 'Greys_r', clim = mapdictU['clim']

        layerstr=mapdictU['layer'].__repr__().lower()
        ax.annotate(mapdictU['date'] ,xy=(0.60,0.93), xycoords='axes fraction' , fontsize=16)
        ax.annotate(layerstr       , xy=(0.77,0.85), xycoords='axes fraction' , ha='center', fontsize=16)

        parallels = np.arange(43.,46.,.5)
        map_obj.drawparallels(parallels,labels=[1,0,0,0],fontsize=13, dashes=[6,900])
        # draw meridians
        meridians = np.arange(12.,20,1.)
        map_obj.drawmeridians(meridians,labels=[0,0,0,1],fontsize=13,dashes=[6,900])


    nticks = 6
    cbar_ticks_list = np.linspace(vmin, vmax, nticks).tolist()
    cbar_ticks_labels = ["%g" % (t,)  for t in cbar_ticks_list ]


    cbar = map_obj.colorbar(cs, location='right',size="5%",pad="2%", ticks=cbar_ticks_list)
    cbar.ax.set_yticklabels(cbar_ticks_labels)
    cbar.ax.tick_params(labelsize=13)

    varTit = 'Current speed'
    title = "%s [ %s ] "  % (varTit, mapdictU['units'])
    ax.set_title(title, fontsize=16)
    ax.title.set_position((0.5,1.03))




    #fig.set_size_inches(10,10*map_obj.ymax/map_obj.xmax)
    #ax.set_position([.08, .08, .83, .83])
    fig.set_size_inches(10,10*(map_obj.ymax/map_obj.xmax)* (.83/.93))
    ax.set_position([.08, .03, .83, .93])

    fig.patch.set_alpha(0.0)
    return fig,ax

