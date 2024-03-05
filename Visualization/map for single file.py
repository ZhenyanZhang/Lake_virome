import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib as mpl
import matplotlib.patches as mpatches
import os

data=pd.read_excel("D:/OneDrive/Lake_virome/sample_mapping.xlsx")

lat = np.array(data['Latitude'])
lon = np.array(data['Longitude'])
number = np.array(data['number'])

plt.style.use('ggplot')
plt.figure(figsize=(10, 6))

#initial map
map1 = Basemap(projection='robin', lat_0=90, lon_0=0,
               resolution='l', area_thresh=1000.0)

map1.drawcoastlines(linewidth=0.2)
map1.drawcountries(linewidth=0.2)
map1.drawmapboundary(fill_color='white') 
map1.fillcontinents(color='lightgrey', alpha=0.8)

map1.drawmeridians(np.arange(0, 360, 60))
map1.drawparallels(np.arange(-90, 90, 30))

#VALUE
map1.scatter(lon, lat, latlon=True,
             alpha=1, s=number*0.01, c="blue", linewidths=0, marker='o', zorder=10)
plt.colorbar()


plt.show()

plt.savefig(outfigure, dpi=300)
plt.close()