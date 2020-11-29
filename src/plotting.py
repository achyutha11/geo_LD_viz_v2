#!python 

"""
    Functions to plot LD metrics on a global map using Cartopy
"""

import cartopy
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as cf

def cartopy_plot(values, coordinates, map_color=False, title, marker='o', marker_color='red'):
    """
    Plot LD statistics on a global map (projection: equirectangular)
    
    Inputs:
        values: LD-statistic to be compared
        coordinates: List of tuples of latitude and longitude
        map_color: (boolean) True if map with color is preferred
        title: (str)
        marker: (str) Defines type of marker used on map
        marker_color: (str)
        
    Output:
        plot: Matplotlib plot with specified values at specified points on global map
    """
    assert len(values) == len(coordinates)
    
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    if map_color:
        ax.stock_img()
    plt.title(title)
    
    for index, lat_long in enumerate(coordinates):
        plt.plot(lat_long[0], lat_long[1], markersize=(10*values[index]), marker=marker, color=marker_color)
        
        

        
    
    
    
