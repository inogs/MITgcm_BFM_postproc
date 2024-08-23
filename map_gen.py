import pylab as pl
from mpl_toolkits.basemap import Basemap
import numpy as np
import cmocean as cmo
from xml.dom import minidom

from dataclasses import dataclass

@dataclass
class text_options():
    fontsize:int
    xpos:float
    ypos:float
    halignment:str
    @staticmethod
    def build_from_node(Node):
        fontsize=int(Node.attributes['fontsize'].value)
        xpos=float(Node.attributes['xpos'].value)
        ypos=float(Node.attributes['ypos'].value)
        halignment=str(Node.attributes['halignment'].value)
        return text_options(fontsize,xpos,ypos,halignment)
    def annotate(self, ax, text, ):
        ax.annotate(text ,xy=(self.xpos,self.ypos), xycoords='axes fraction' , fontsize=self.fontsize, ha=self.halignment)

@dataclass
class ticks_position():
    start:float
    stop:float
    step:float
    @staticmethod
    def build_from_node(Node):
        start=float(Node.attributes['start'].value)
        stop=float(Node.attributes['stop'].value)
        step=float(Node.attributes['step'].value)
        return ticks_position(start,stop,step)
    def to_array(self):
        return np.arange(self.start, self.stop, self.step)


class Map_object(Basemap):
    def __init__(self, xlim,ylim, datestr, layerstr, parallels, meridians):
        xC=(xlim[0]+xlim[1])/2
        yC=(ylim[0]+ylim[1])/2
        self._parallels=np.array(parallels)
        self._meridians = np.array(meridians)
        self.datestr=datestr
        self.layerstr=layerstr
        super().__init__(projection='merc',lat_0=xC,lon_0=yC,\
                    llcrnrlon = xlim[0], \
                    llcrnrlat = ylim[0], \
                    urcrnrlon = xlim[1], \
                    urcrnrlat = ylim[1], \
                    resolution='h')


    def map_plotter_basemap(self, mapdict,maskobj):
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
        self.drawcoastlines(linewidth=0.5)
        #map.drawcountries(linewidth=0.25)
        map2d=mapdict['data']
        Zm = np.ma.masked_invalid(map2d)
        self.pcolormesh(maskobj.xlevels, maskobj.ylevels, Zm,cmap=cmap,latlon='true',vmin=vmin,vmax=vmax)

        ax.annotate(mapdict['date'][5:] ,xy=(0.9,0.9), xycoords='axes fraction' , fontsize=28);

        self.drawparallels(self._parallels,labels=[1,0,0,0],fontsize=16, dashes=[6,900])
        self.drawmeridians(self._meridians,labels=[0,0,0,1],fontsize=16,dashes=[6,900])
        #cbar = map_obj.colorbar(cs,location='right',size="5%",pad="2%")
        #cbar.ax.tick_params(labelsize=5)


        fig.set_size_inches(10,10*self.ymax/self.xmax)
        ax.set_position([.10, .08, .88, .88])
        fig.patch.set_alpha(0.0)
        return fig,ax

    def map_plotter_basemap_hourly(self,mapdict,maskobj):
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
        self.drawcoastlines(linewidth=0.5)
        #map.drawcountries(linewidth=0.25)
        map2d=mapdict['data']
        Zm = np.ma.masked_invalid(map2d)
        cs=self.pcolormesh(maskobj.xlevels, maskobj.ylevels, Zm,cmap=cmap,latlon='true',shading='nearest',vmin=vmin,vmax=vmax)

        layerstr=mapdict['layer'].__repr__().lower()
        self.datestr.annotate(ax, mapdict['date'])
        self.layerstr.annotate(ax,layerstr)

        self.drawparallels(self._parallels,labels=[1,0,0,0],fontsize=13, dashes=[6,900])
        self.drawmeridians(self._meridians,labels=[0,0,0,1],fontsize=13,dashes=[6,900])

        nticks = 6
        cbar_ticks_list = np.linspace(vmin, vmax, nticks).tolist()
        cbar_ticks_labels = ["%g" % (t,)  for t in cbar_ticks_list ]


        cbar = self.colorbar(cs,location='right',size="5%",pad="2%", ticks=cbar_ticks_list)
        cbar.ax.set_yticklabels(cbar_ticks_labels)
        cbar.ax.tick_params(labelsize=13)

        title = "%s [ %s ] "  % (mapdict['longname'],mapdict['units'])
        ax.set_title(title, fontsize=16)
        ax.title.set_position((0.5,1.03))

        #fig.set_size_inches(10,10*self.ymax/self.xmax)
        #ax.set_position([.08, .08, .83, .83])
        fig.set_size_inches(10,10*(self.ymax/self.xmax)* (.83/.93))
        ax.set_position([.08, .03, .83, .93])

        fig.patch.set_alpha(0.0)
        return fig,ax



    def quiver_plotter_basemap_hourly(self, mapdictU, mapdictV, maskobj):
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


        spaceSteps = int(6)

        # draw coastlines, country boundaries, fill continents.
        self.drawcoastlines(linewidth=0.5)
        #map.drawcountries(linewidth=0.25)
        U2d=mapdictU['data']
        V2d=mapdictV['data']
        Um = np.ma.masked_invalid(U2d)
        Vm = np.ma.masked_invalid(V2d)
        Zm = np.ma.sqrt(Um**2 + Vm**2)
        #cs=self.quiver(maskobj.xlevels[::6, ::6], maskobj.ylevels[::6, ::6], Um[::6, ::6], Vm[::6, ::6], Zm[::6, ::6], cmap=cmap, latlon='true', clim = mapdictU['clim'])
        cs=self.pcolormesh(maskobj.xlevels, maskobj.ylevels, Zm, cmap=cmap,latlon='true',shading='nearest', clim = mapdictU['clim'])
        #ColorMatrix = 1. - 1./(np.exp((Zm - 22/24*vmax)/0.001) + 1)
        qp = self.quiver(maskobj.xlevels[::spaceSteps, ::spaceSteps], maskobj.ylevels[::spaceSteps, ::spaceSteps], Um[::spaceSteps, ::spaceSteps], Vm[::spaceSteps, ::spaceSteps], latlon='true') #, ColorMatrix[::6, ::6], cmap = 'Greys_r', clim = mapdictU['clim']

        layerstr=mapdictU['layer'].__repr__().lower()
        self.datestr.annotate(ax, mapdictU['date'])
        self.layerstr.annotate(ax,layerstr)

        self.drawparallels(self._parallels,labels=[1,0,0,0],fontsize=13, dashes=[6,900])
        self.drawmeridians(self._meridians,labels=[0,0,0,1],fontsize=13,dashes=[6,900])


        nticks = 6
        cbar_ticks_list = np.linspace(vmin, vmax, nticks).tolist()
        cbar_ticks_labels = ["%g" % (t,)  for t in cbar_ticks_list ]


        cbar = self.colorbar(cs, location='right',size="5%",pad="2%", ticks=cbar_ticks_list)
        cbar.ax.set_yticklabels(cbar_ticks_labels)
        cbar.ax.tick_params(labelsize=13)

        varTit = 'Current speed'
        title = "%s [ %s ] "  % (varTit, mapdictU['units'])
        ax.set_title(title, fontsize=16)
        ax.title.set_position((0.5,1.03))


        #fig.set_size_inches(10,10*self.ymax/self.xmax)
        #ax.set_position([.08, .08, .83, .83])
        fig.set_size_inches(10,10*(self.ymax/self.xmax)* (.83/.93))
        ax.set_position([.08, .03, .83, .93])

        fig.patch.set_alpha(0.0)
        return fig,ax


    @staticmethod
    def create_from_file(configfile):
        with open(configfile, 'r') as f:
            xml=minidom.parse(f)
        xmin=float(xml.getElementsByTagName('xlim')[0].attributes['min'].value)
        xmax=float(xml.getElementsByTagName('xlim')[0].attributes['max'].value)
        ymin=float(xml.getElementsByTagName('ylim')[0].attributes['min'].value)
        ymax=float(xml.getElementsByTagName('ylim')[0].attributes['max'].value)
        xlim = [xmin, xmax]
        ylim = [ymin, ymax]


        datestr   = text_options.build_from_node(xml.getElementsByTagName('date')[0])
        layerstr  = text_options.build_from_node(xml.getElementsByTagName('layerstr')[0])
        parallels = ticks_position.build_from_node(xml.getElementsByTagName('parallels')[0])
        meridians = ticks_position.build_from_node(xml.getElementsByTagName('meridians')[0])
        parallel_array=parallels.to_array()
        meridian_array=meridians.to_array()

        return Map_object(xlim, ylim, datestr, layerstr, parallel_array, meridian_array)



