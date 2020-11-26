#!python 

"""
    Functions to plot LD metrics on a global map using Cartopy
"""

import cartopy
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as cf

def cartopy_plot(values, coordinates):
    """
    Plot LD statistics on a global map
    
    Inputs:
        values: LD-statistic to be compared
        coordinates: List of tuples of latitude and longitude
        
    Output:
        plot: Plot with specified values at required points on global map
    """
    
    
